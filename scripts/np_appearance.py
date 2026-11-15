#!/usr/bin/env python3
"""
np_appearance.py
================
B_n 의 부분군 H 가 어떤 부울 함수의 안정화군(NP 안정화군)으로 "나타나는가"를 표지 표 없이 직접 정한다.

  H 의 꼭짓점 궤도들의 집합 Omega (k 개) 위의 2-색칠 c 가 함수 f_c 를 정한다.
  g 가 f_c 를 보존하는 것은 c 가 "g 에 맞다", 곧 g 로 이어진 궤도 쌍 (O(v), O(g v)) 에서 색이 같은 것과 같다.
  H 가 나타난다 <=> 모든 g in B_n \\ H 에 대해 g 로 이어진 쌍 가운데 색이 다른 것이 있는 색칠 c 가 있다.
  g 마다 쌍들의 집합 E_g (궤도가 다른 것만)를 적고, 그 집합이 부분집합 관계로 극소인 것만 남긴 뒤
  "E_g 의 어떤 쌍은 색이 다르다" 는 절들의 충족 가능성을 SAT 으로 푼다.
  E_g 가 공집합인 g 가 있으면 (H 가 닫히지 않았다) 절이 비어 있어 나타나지 않는다.

  python3 -I scripts/np_appearance.py <n> <생성원.tsv> [<켤레류 분포.tsv(aO 대조용)>] [--mus] [--list]
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bn_subgroup_counts as BS  # noqa: E402


class Analyzer:
    def __init__(self, n):
        self.n = n
        self.B = BS.Bn(n)
        B = self.B
        V = B.V
        a_idx = np.arange(B.N) // V
        c_idx = np.arange(B.N) % V
        dtype = np.uint8 if V <= 256 else np.int16
        self.VA = np.empty((B.N, V), dtype=dtype)
        step = 1 << 16
        for s in range(0, B.N, step):
            self.VA[s:s + step] = (B.act[a_idx[s:s + step]] ^ c_idx[s:s + step, None]).astype(dtype)

    def orbit_labels(self, gens):
        """꼭짓점 궤도 번호(0..k-1)와 궤도 수"""
        V = self.B.V
        parent = list(range(V))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        for g in gens:
            vp = self.VA[int(g)]
            for x in range(V):
                rx, ry = find(x), find(int(vp[x]))
                if rx != ry:
                    parent[rx] = ry
        roots = {}
        lab = np.zeros(V, dtype=np.int64)
        for x in range(V):
            r = find(x)
            if r not in roots:
                roots[r] = len(roots)
            lab[x] = roots[r]
        return lab, len(roots)

    def analyze(self, gens, want_mus=False):
        B = self.B
        seen = B.closure(gens)
        order = int(seen.sum())
        lab, k = self.orbit_labels(gens)
        LG = lab[self.VA.astype(np.int64)] if self.VA.dtype != np.int64 else lab[self.VA]   # (N, V) 상의 궤도 번호
        moved = LG != lab[None, :]
        nmoved = moved.sum(axis=1)
        hat = (nmoved == 0)
        hat_order = int(hat.sum())
        res = dict(order=order, k=k, hat_order=hat_order, closed=(hat_order == order))
        if not res["closed"]:
            res["appears"] = False
            return res
        # 쌍 코드 (작은 번호 * k + 큰 번호), 이동하지 않은 꼭짓점은 -1
        lo = np.minimum(LG, lab[None, :])
        hi = np.maximum(LG, lab[None, :])
        codes = np.where(moved, lo * k + hi, -1)
        codes[seen] = -1                               # H 의 원소는 제외 (이동이 없다)
        codes.sort(axis=1)
        dup = np.zeros_like(codes, dtype=bool)
        dup[:, 1:] = codes[:, 1:] == codes[:, :-1]
        codes[dup] = -1
        codes.sort(axis=1)
        rows = np.unique(codes[~seen], axis=0)
        sets = []
        for r in rows:
            e = r[r >= 0]
            if e.size:
                sets.append(frozenset(int(x) for x in e))
        sets.sort(key=len)
        minimal = []
        for s in sets:
            if not any(m <= s for m in minimal):
                minimal.append(s)
        res["n_clauses"] = len(sets)
        res["minimal"] = minimal
        sat, model, mus = solve_sat(k, minimal, want_mus)
        res["appears"] = sat
        res["model"] = model
        res["mus"] = mus
        res["lab"] = lab
        return res


def solve_sat(k, clauses, want_mus=False):
    """색 변수 1..k, 쌍 변수. 반환 (충족 여부, 색칠 또는 None, 극소 모순 부분집합(절 번호 목록) 또는 None)"""
    from pysat.solvers import Cadical153
    edges = sorted(set(e for c in clauses for e in c))
    evar = {e: k + 1 + i for i, e in enumerate(edges)}
    nvars = k + len(edges)
    sel0 = nvars + 1
    s = Cadical153()
    for e, d in evar.items():
        a, b = divmod(e, k)
        # d -> (c_a != c_b)
        s.add_clause([-d, a + 1, b + 1])
        s.add_clause([-d, -(a + 1), -(b + 1)])
    for i, c in enumerate(clauses):
        s.add_clause([-(sel0 + i)] + [evar[e] for e in c])
    assum = [sel0 + i for i in range(len(clauses))]
    ok = s.solve(assumptions=assum)
    if ok:
        m = s.get_model()
        val = {abs(x): x > 0 for x in m}
        col = [1 if val.get(i + 1, False) else 0 for i in range(k)]
        s.delete()
        return True, col, None
    mus = None
    if want_mus:
        core = [a - sel0 for a in s.get_core()]
        core = sorted(core)
        i = 0
        while i < len(core):
            trial = core[:i] + core[i + 1:]
            if not s.solve(assumptions=[sel0 + t for t in trial]):
                core = trial
            else:
                i += 1
        mus = core
    s.delete()
    return False, None, mus


def load_gens(path, B, n):
    out = {}
    for line in open(path, encoding="utf-8"):
        if not line.startswith("G\t"):
            continue
        p = line.rstrip("\n").split("\t")
        j, order, orb, ng = int(p[1]), int(p[2]), int(p[3]), int(p[4])
        vals = [int(t) for t in p[5:]]
        gens = [B.from_images(vals[t * (n + 1):(t + 1) * (n + 1)]) for t in range(ng)]
        out[j] = (order, orb, gens)
    return out


def load_aO(path):
    aO = {}
    for line in open(path, encoding="utf-8"):
        if line.startswith("L\t"):
            p = line.split("\t")
            aO[int(p[1])] = int(p[4])
    return aO


def main():
    pos = [a for a in sys.argv[1:] if not a.startswith("--")]
    n = int(pos[0])
    gens_path = pos[1]
    aO = load_aO(pos[2]) if len(pos) > 2 else None
    want_mus = "--mus" in sys.argv
    A = Analyzer(n)
    G = load_gens(gens_path, A.B, n)
    stats = dict(total=0, closed=0, appear=0, closed_absent=0, mismatch=0)
    absent = []
    for j in sorted(G):
        order, orb, gens = G[j]
        r = A.analyze(gens, want_mus=want_mus)
        if r["order"] != order or r["k"] != orb:
            raise SystemExit("위수나 궤도 수가 다르다 j=%d" % j)
        stats["total"] += 1
        stats["closed"] += r["closed"]
        stats["appear"] += r["appears"]
        if r["closed"] and not r["appears"]:
            stats["closed_absent"] += 1
            absent.append((j, r))
        if aO is not None and r["appears"] != (aO[j] > 0):
            stats["mismatch"] += 1
            print("   표지 표의 aO 와 다르다: j=%d, 이 방법 %s, aO %d" % (j, r["appears"], aO[j]))
    print("n = %d: 부분군 류 %d, 닫힘 %d, 나타남 %d, 닫혔지만 안 나타남 %d, 표지 표(aO > 0)와 다른 것 %s"
          % (n, stats["total"], stats["closed"], stats["appear"], stats["closed_absent"],
             stats["mismatch"] if aO is not None else "(대조 안 함)"))
    if "--list" in sys.argv:
        for j, r in sorted(absent, key=lambda t: -t[1]["k"]):
            print("  j=%d 위수 %d 꼭짓점 궤도 %d, 극소 절 %d개 (서로 다른 절 %d개)%s"
                  % (j, r["order"], r["k"], len(r["minimal"]), r["n_clauses"],
                     "" if r["mus"] is None else ", 극소 모순 부분집합 절 %d개" % len(r["mus"])))


if __name__ == "__main__":
    main()
