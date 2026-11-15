#!/usr/bin/env python3
"""
bn_subgroup_counts.py
=====================
표지 표(GAP)를 쓰지 않고 B_n 의 부분군을 직접 세는 도구. 원소는 부호 있는 순열 (sigma, c): x -> sigma(x) xor c.
  - B_n 의 모든 원소를 부호화 (a * 2^n + c, a = 변수 순열의 번호)해 곱셈표 없이 벡터로 곱한다.
  - 부분군의 생성원에서 원소 전체를 구한다 (너비 우선).
  - B_n 의 위수 2 원소(대합)를 모두 찾고, 위수 4 의 클라인 4군 V4 와 위수 8 의 C2^3 (V8) 를 모두 찾아
    B_n 안의 켤레류로 나눈다 (켤레 작용의 연결 성분).  류 크기에서 |N_B(U)| = |B_n| / (류 크기).
  - 부분군 L 안에 있는 U 와 켤레인 부분군의 수 s(U, L) 를 센다.
scripts/verify_npn_class_identities.py 가 이 수로 "B_n x 1 의 부분군 U 에 대한 고정점 항등식"을 세운다:
  2^{orb(U)} = sum_{칸 (L, N)} W * m_U,   m_U = 2 phi_L(U) (phi = 0), phi_N(U) (phi != 0),
  phi_X(U) = |N_B(U)| s(U, X) / |X|   (X 위의 B_n-집합에서 U 의 고정점 수).
"""
import itertools
import numpy as np


class Bn:
    def __init__(self, n):
        self.n = n
        self.V = 1 << n
        self.perms = list(itertools.permutations(range(n)))   # p[i] = 변수 i 가 가는 곳
        self.P = len(self.perms)
        self.N = self.P * self.V
        pidx = {p: i for i, p in enumerate(self.perms)}
        self.pidx = pidx
        P, V = self.P, self.V
        act = np.zeros((P, V), dtype=np.int64)
        for a, p in enumerate(self.perms):
            for c in range(V):
                r = 0
                for i in range(n):
                    if c >> i & 1:
                        r |= 1 << p[i]
                act[a, c] = r
        self.act = act
        comp = np.zeros((P, P), dtype=np.int64)       # comp[a][b] = p_a o p_b  (p_b 를 먼저)
        for a, pa in enumerate(self.perms):
            for b, pb in enumerate(self.perms):
                comp[a, b] = pidx[tuple(pa[pb[i]] for i in range(n))]
        self.comp = comp
        pinv = np.zeros(P, dtype=np.int64)
        for a, pa in enumerate(self.perms):
            q = [0] * n
            for i in range(n):
                q[pa[i]] = i
            pinv[a] = pidx[tuple(q)]
        self.pinv = pinv
        self.identity = pidx[tuple(range(n))] * V

    def mul(self, g, h):
        """g o h (h 를 먼저 적용). 넘파이 배열(방송 가능)이나 정수."""
        g = np.asarray(g)
        h = np.asarray(h)
        ga, gc = g // self.V, g % self.V
        ha, hc = h // self.V, h % self.V
        return self.comp[ga, ha] * self.V + (gc ^ self.act[ga, hc])

    def inv(self, g):
        g = np.asarray(g)
        a, c = g // self.V, g % self.V
        ia = self.pinv[a]
        return ia * self.V + self.act[ia, c]

    def from_images(self, imgs):
        """GAP 의 내보내기 (꼭짓점 0 의 상, e_1..e_n 의 상) -> 원소 부호"""
        n = self.n
        c = imgs[0]
        p = []
        for i in range(n):
            y = imgs[1 + i] ^ c
            if y == 0 or y & (y - 1):
                raise ValueError("부호 있는 순열이 아니다: %s" % (imgs,))
            p.append(y.bit_length() - 1)
        return self.pidx[tuple(p)] * self.V + c

    def vertex_perm(self, g):
        a, c = int(g) // self.V, int(g) % self.V
        return [int(self.act[a, x]) ^ c for x in range(self.V)]

    def orbits_count(self, gens):
        parent = list(range(self.V))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        for g in gens:
            vp = self.vertex_perm(g)
            for x in range(self.V):
                rx, ry = find(x), find(vp[x])
                if rx != ry:
                    parent[rx] = ry
        return len(set(find(x) for x in range(self.V)))

    def closure(self, gens):
        """생성원으로 만든 부분군의 원소를 참/거짓 배열(길이 |B_n|)로"""
        seen = np.zeros(self.N, dtype=bool)
        seen[self.identity] = True
        frontier = np.array([self.identity], dtype=np.int64)
        gens = [int(g) for g in gens]
        while frontier.size:
            new = np.unique(np.concatenate([self.mul(g, frontier) for g in gens])) if gens else frontier[:0]
            new = new[~seen[new]]
            seen[new] = True
            frontier = new
        return seen


class ElementaryAbelian:
    """B_n 의 위수 2 원소와 C2^2, C2^3 부분군 전부, 그리고 B_n-켤레류 (정확한 부호화, 해시 없음)"""

    def __init__(self, B, ranks=(2, 3), verbose=True):
        self.B = B
        allc = np.arange(B.N, dtype=np.int64)
        sq = B.mul(allc, allc)
        self.invs = allc[(sq == B.identity) & (allc != B.identity)]
        I = len(self.invs)
        self.I = I
        if I >= 2048:
            raise ValueError("대합이 너무 많다")
        self.invlookup = -np.ones(B.N, dtype=np.int64)
        self.invlookup[self.invs] = np.arange(I)
        if verbose:
            print("   B_n 의 위수 2 원소 %d개 (항등원 제외)" % I)
        M = B.mul(self.invs[:, None], self.invs[None, :])
        self.comm = (M == M.T)
        self.prod = M
        np.fill_diagonal(self.comm, False)
        # B_n 의 생성원: 부호 바꾸기 n개와 이웃 호환 n-1개
        n = B.n
        gens = []
        for i in range(n):
            gens.append(B.pidx[tuple(range(n))] * B.V + (1 << i))
        for i in range(n - 1):
            p = list(range(n))
            p[i], p[i + 1] = p[i + 1], p[i]
            gens.append(B.pidx[tuple(p)] * B.V)
        self.bgens = gens
        # 켤레 작용 x -> h x h^-1 이 대합의 번호에 일으키는 순열
        self.conj = []
        for h in gens:
            cx = B.mul(B.mul(h, self.invs), B.inv(h))
            k = self.invlookup[cx]
            assert (k >= 0).all()
            self.conj.append(k)
        self.sub = {}
        if 2 in ranks:
            self._find_v4(verbose)
        if 3 in ranks:
            self._find_v8(verbose)

    @staticmethod
    def _key(S):
        """정렬한 행(각 성분 < 2048)을 정확한 키(바이트열)로"""
        return [bytes(r) for r in np.ascontiguousarray(S.astype(np.int16))]

    def _classes(self, T, verbose, name):
        """T: (m, r) 정렬된 행, 부분군마다 대합 번호들. 켤레류 번호(연결 성분)를 돌려준다"""
        from scipy.sparse import coo_matrix
        from scipy.sparse.csgraph import connected_components
        m = T.shape[0]
        keys = self._key(T)
        index = {k: i for i, k in enumerate(keys)}
        if len(index) != m:
            raise RuntimeError("%s: 같은 부분군이 겹쳐 있다" % name)
        rows, cols = [], []
        for k in self.conj:
            Ti = np.sort(k[T], axis=1)
            kk = self._key(Ti)
            tgt = np.fromiter((index[x] for x in kk), dtype=np.int64, count=m)
            rows.append(np.arange(m))
            cols.append(tgt)
        rows = np.concatenate(rows)
        cols = np.concatenate(cols)
        g = coo_matrix((np.ones(len(rows), dtype=np.int8), (rows, cols)), shape=(m, m))
        ncls, lab = connected_components(g, directed=False)
        sizes = np.bincount(lab, minlength=ncls)
        if verbose:
            print("   %s: 부분군 %d개, B_n-켤레류 %d개" % (name, m, ncls))
        return dict(T=T, lab=lab, sizes=sizes, index=index)

    def _find_v4(self, verbose):
        xs, ys = np.nonzero(np.triu(self.comm, 1))
        zs = self.invlookup[self.prod[xs, ys]]
        ok = zs >= 0
        xs, ys, zs = xs[ok], ys[ok], zs[ok]
        T = np.sort(np.stack([xs, ys, zs], axis=1), axis=1)
        T = np.unique(T, axis=0)
        self.sub[2] = self._classes(T, verbose, "C2^2")

    def _find_v8(self, verbose):
        comm = self.comm
        xs, ys = np.nonzero(np.triu(comm, 1))
        zs = self.invlookup[self.prod[xs, ys]]
        ok = zs >= 0
        xs, ys, zs = xs[ok], ys[ok], zs[ok]
        cc = comm
        chunk = 1500
        uniq = []
        for s in range(0, len(xs), chunk):
            x, y, z = xs[s:s + chunk], ys[s:s + chunk], zs[s:s + chunk]
            pi, w = np.nonzero(cc[x] & cc[y])
            keep = w != z[pi]
            pi, w = pi[keep], w[keep]
            if pi.size == 0:
                continue
            X, Y, Z = x[pi], y[pi], z[pi]
            xw = self.invlookup[self.prod[X, w]]
            yw = self.invlookup[self.prod[Y, w]]
            zw = self.invlookup[self.prod[Z, w]]
            S = np.sort(np.stack([X, Y, Z, w, xw, yw, zw], axis=1), axis=1)
            if (S[:, 1:] == S[:, :-1]).any():
                raise RuntimeError("C2^3: 7 개의 대합이 서로 다르지 않다")
            uniq.append(np.unique(S, axis=0))
        T = np.unique(np.concatenate(uniq), axis=0)
        self.sub[3] = self._classes(T, verbose, "C2^3")

    def class_of(self, r, invs_in_group):
        """부분군의 (항등원이 아닌 원소 중) 대합 번호들이 주어졌을 때 그 부분군의 류 번호"""
        d = self.sub[r]
        t = np.sort(np.array(invs_in_group, dtype=np.int64))
        return int(d["lab"][d["index"][bytes(t.astype(np.int16))]])

    def counts_in(self, r, inL_mask):
        """대합 번호 위의 참/거짓 배열(그 부분군 L 에 들어 있는 대합)로, L 안의 C2^r 부분군을 류별로 센다"""
        d = self.sub[r]
        inside = inL_mask[d["T"]].all(axis=1)
        return np.bincount(d["lab"][inside], minlength=len(d["sizes"]))
