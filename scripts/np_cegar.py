#!/usr/bin/env python3
"""
np_cegar.py
===========
B_n (n <= 8) 의 부분군 H 가 어떤 부울 함수의 NP 안정화군으로 나타나는가를 "증명서가 붙은" 방식으로 정한다.
B_n 전체의 원소(n!*2^n 개)를 표로 만들지 않고, 함수 하나의 안정화군을 n! 번의 사전 조회로 정확히 구한다.

  F(pi x + c) = F(x) (모든 x)  <=>  F(. + c) = F o rho  (rho = pi^-1)
  그래서 평행이동한 함수 F(. + c) 들을 사전에 넣고, 모든 순열 rho 에 대해 F o rho 를 조회하면 안정화군 전체가 나온다.

  1) 닫힘: F = (H 의 꼭짓점 궤도 번호) 의 안정화군 = H^ 가 H 와 같은가.
  2) 반례-유도 (CEGAR): 궤도의 2-색칠 c 에 대한 SAT 를 풀고, 나온 색칠의 함수 f_c 의 안정화군을 정확히 구한다.
     그 위수가 |H| 이면 H 는 나타난다 (f_c 가 증거). 아니면 안정화군의 H 밖 원소 g 로 절 "E_g 의 어떤 쌍은 색이 다르다" 를
     더하고 되풀이한다. SAT 가 모순이면 H 는 나타나지 않는다 (더한 절들이 증거).
  나타남 판정은 f_c 의 안정화군 정확 계산이고, 안 나타남 판정은 실제 원소들이 준 절들의 모순이므로 둘 다 틀릴 수 없다.

사용: from np_cegar import Cube, appears
"""
import itertools
import os
import sys
import numpy as np


class Cube:
    def __init__(self, n):
        self.n = n
        self.V = 1 << n
        self.perms = np.array(list(itertools.permutations(range(n))), dtype=np.int64)   # Bn 과 같은 순서
        self.P = len(self.perms)
        X = np.arange(self.V, dtype=np.int64)
        self.X = X
        dt = np.uint8 if self.V <= 256 else np.uint16
        act = np.zeros((self.P, self.V), dtype=np.int64)
        for i in range(n):
            act += ((X >> i) & 1)[None, :] << self.perms[:, i][:, None]
        self.act = act.astype(dt)
        self.dt = dt

    # ---- 군 원소는 꼭짓점 순열 (길이 V 배열) 로 다룬다 ----
    def gen_vp(self, pidx, c):
        """x -> act[pidx][x] ^ c (Bn 의 부호 a*V+c 와 같은 규약)"""
        return (self.act[pidx].astype(np.int64) ^ c).astype(self.dt)

    def closure(self, gens_vp, limit=200000):
        ident = self.X.astype(self.dt)
        seen = {ident.tobytes(): ident}
        frontier = [ident]
        while frontier:
            new = []
            for a in frontier:
                for g in gens_vp:
                    b = g[a]
                    kb = b.tobytes()
                    if kb not in seen:
                        seen[kb] = b
                        new.append(b)
                        if len(seen) > limit:
                            raise RuntimeError("군이 너무 크다")
            frontier = new
        return seen

    def orbit_labels(self, gens_vp):
        from scipy.sparse import coo_matrix
        from scipy.sparse.csgraph import connected_components
        V = self.V
        if len(gens_vp) == 0:
            return self.X.copy(), V
        rows = np.concatenate([self.X for _ in gens_vp])
        cols = np.concatenate([g.astype(np.int64) for g in gens_vp])
        M = coo_matrix((np.ones(len(rows), dtype=np.int8), (rows, cols)), shape=(V, V))
        k, lab = connected_components(M, directed=False)
        return lab.astype(np.int64), int(k)

    def stab_matches(self, F, limit=None):
        """F 를 보존하는 모든 (rho 번호, c). 원소 h: y -> rho(y ^ c) (곧 g 의 역원)"""
        F = np.asarray(F)
        X = self.X
        T = F[X[:, None] ^ X[None, :]]                        # T[c][y] = F[y ^ c]
        Tc = np.ascontiguousarray(T)
        d = {}
        for c in range(self.V):
            d.setdefault(Tc[c].tobytes(), []).append(c)
        K = np.ascontiguousarray(F[self.act])                  # K[a][y] = F[rho_a(y)]
        out = []
        for a in range(self.P):
            lst = d.get(K[a].tobytes())
            if lst is not None:
                for c in lst:
                    out.append((a, c))
                    if limit is not None and len(out) >= limit:
                        return out
        return out

    def count_stab(self, F):
        F = np.asarray(F)
        X = self.X
        T = np.ascontiguousarray(F[X[:, None] ^ X[None, :]])
        d = {}
        for c in range(self.V):
            kb = T[c].tobytes()
            d[kb] = d.get(kb, 0) + 1
        K = np.ascontiguousarray(F[self.act])
        tot = 0
        for a in range(self.P):
            tot += d.get(K[a].tobytes(), 0)
        return tot

    def h_vp(self, a, c):
        return self.act[a][(self.X ^ c)]


def edge_codes(lab, vp, k):
    L2 = lab[vp.astype(np.int64)]
    mv = L2 != lab
    if not mv.any():
        return np.zeros(0, dtype=np.int64)
    lo = np.minimum(L2[mv], lab[mv])
    hi = np.maximum(L2[mv], lab[mv])
    return np.unique(lo * k + hi)


def appears(cube, gens_vp, order=None, max_iter=20000, per_iter=40, rng=None, verbose=False):
    """반환 dict(status in {'present','absent','notclosed','unknown'}, k, order, iters, f(있다면), clauses)
    order 를 주면 H 의 원소를 나열하지 않는다 (닫힌 군에서는 H = H^ 이므로 'H 의 원소인가' 는 '모든 궤도를 보존하는가' 와 같다)."""
    from pysat.solvers import Cadical153
    rng = rng or np.random.default_rng(1)
    if order is None:
        order = len(cube.closure(gens_vp))
    lab, k = cube.orbit_labels(gens_vp)
    hat = cube.count_stab(lab.astype(np.int64) if cube.V > 256 else lab.astype(np.uint16))
    res = dict(k=k, order=order, hat=hat)
    if hat != order:
        res["status"] = "notclosed"
        return res
    edge_var = {}
    solver = Cadical153()
    nv = [k]

    def evar(e):
        if e not in edge_var:
            nv[0] += 1
            d = nv[0]
            edge_var[e] = d
            a, b = divmod(e, k)
            solver.add_clause([-d, a + 1, b + 1])
            solver.add_clause([-d, -(a + 1), -(b + 1)])
        return edge_var[e]
    nclauses = 0
    seen_clauses = set()
    for it in range(max_iter):
        if not solver.solve():
            solver.delete()
            res.update(status="absent", iters=it, clauses=nclauses)
            return res
        m = solver.get_model()
        val = {abs(x): x > 0 for x in m}
        col = np.array([1 if val.get(i + 1, False) else 0 for i in range(k)], dtype=np.uint8)
        f = col[lab]
        matches = cube.stab_matches(f)
        if len(matches) == order:
            solver.delete()
            res.update(status="present", iters=it, f=f, col=col, clauses=nclauses)
            return res
        idx = np.arange(len(matches))
        if len(matches) > per_iter * 8:
            idx = rng.choice(len(matches), per_iter * 8, replace=False)
        cands = []
        for i in idx:
            a, c = matches[int(i)]
            e = edge_codes(lab, cube.h_vp(a, c), k)
            if len(e):                      # 빈 쌍 집합이면 H^ = H 의 원소
                cands.append(e)
        if not cands and len(matches) > len(idx):
            for (a, c) in matches:
                e = edge_codes(lab, cube.h_vp(a, c), k)
                if len(e):
                    cands.append(e)
                    break
        if not cands:
            raise RuntimeError("안정화군이 H 보다 큰데 H 밖 원소를 못 찾았다")
        cands.sort(key=len)
        added = 0
        for e in cands:
            key = e.tobytes()
            if key in seen_clauses:
                continue
            seen_clauses.add(key)
            solver.add_clause([evar(int(x)) for x in e])
            nclauses += 1
            added += 1
            if added >= per_iter:
                break
        if added == 0:
            raise RuntimeError("새 절이 없다")
        if verbose and it % 50 == 0:
            print("   iter", it, "clauses", nclauses, "stab", len(matches), file=sys.stderr)
    solver.delete()
    res.update(status="unknown", iters=max_iter, clauses=nclauses)
    return res


def vp_from_code(cube, code):
    a, c = divmod(int(code), cube.V)
    return cube.gen_vp(a, c)
