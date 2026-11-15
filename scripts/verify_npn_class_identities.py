#!/usr/bin/env python3
"""
verify_npn_class_identities.py
==============================
G_n = B_n x Z_2 (NPN 군) 의 켤레류마다 "고정점 수 = 안정화군 종류별 NPN 류 수와 치환지표의 합" 항등식을 정확히 세워 검산한다.

  Fix(t) = sum_{K} d_K * pi_K(t),   pi_K(t) = |C_G(t)| |K cap t^G| / |K|        (t 는 G_n 의 켤레류 대표원)
  Fix((g, 0)) = 2^{orb(g)}                                   (orb(g) = g 가 꼭짓점 2^n 개에 일으키는 순환의 수)
  Fix((g, 1)) = 2^{orb(g)}  (g 의 모든 순환이 짝수 길이일 때), 아니면 0
  d_K = 안정화군이 K 인 NPN 류 수.  K = {(g, phi(g)) : g in L}  (phi : L -> Z_2),  N = ker phi.
  pi_K((g,0)) = 2 |C_B(g)| |g^B cap L| / |L|                 (phi = 0, N = L)
              =   |C_B(g)| |g^B cap N| / |N|                 (phi != 0)
  pi_K((g,1)) = 0                                            (phi = 0)
              =   |C_B(g)| (|g^B cap L| - |g^B cap N|) / |N| (phi != 0)
  (B = B_n; g^B 는 g 의 B_n-켤레류.)  그래서 치환지표는 L 과 N 의 B_n-켤레류 분포만 있으면 된다.

두 가지를 한다.
  1. 연립방정식(scripts/solve_npn_pair_resolution.py)이 낸 모든 칸의 값이 G_n 의 켤레류 항등식을 모두 만족하는가.
  2. "직접 센 칸"만 믿고 나머지 칸을 미지수로 두었을 때 이 항등식이 그 칸들을 (유일하게) 정하는가, 그 값이 풀이기의 값과 같은가.
     연립방정식 풀이는 B_n 의 표지 표(꼬인 작용 포함)를 쓰지만 이 방법은 G_n 의 켤레류 항등식만 쓰므로, 표지 표와 따로 얻는 값이다.

  python3 -I scripts/verify_npn_class_identities.py <n> <자료.tsv> <켤레류 분포.tsv> [--fix=고정값.txt]
         [--unknown-orb-min=T] [--gens=생성원.tsv]
    <켤레류 분포.tsv> 는 scripts/export_bn_class_data.g 가, <생성원.tsv> 는 scripts/export_bn_generators.g 가 쓴다.
    --gens=파일      : 켤레류 항등식에 더해, B_n x 1 의 부분군 U = C2^2, C2^3 에 대한 항등식도 쓴다.
                       U 와 켤레인 부분군의 수는 표지 표가 아니라 scripts/bn_subgroup_counts.py 가 직접 센다.
                       (칸 (L, N) 의 U 에 대한 고정점 수: 2 phi_L(U) (phi = 0), phi_N(U) (phi != 0), phi_X(U) = |N_B(U)| s(U,X) / |X|.)
    --fix=파일       : 직접 센 값(칸 = (j, 서명)). 이 파일에 j 가 없는 칸이 미지수가 된다 (n = 6 의 25개 군).
    --unknown-orb-min=T : Kbar 의 꼭짓점 궤도 수가 T 이상인 칸을 미지수로 두고 나머지는 (풀이기의) 값으로 쓴다. --fix 보다 먼저 적용한다.
                          (n <= 5 의 시험용이고, n = 6 에서는 직접 센 칸을 몇 칸이나 다시 얻을 수 있는지 본다.)
"""
import os
import sys
from fractions import Fraction
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import solve_npn_pair_resolution as S  # noqa: E402


def load_classes(path):
    """B 줄, L 줄 -> (n, bsize, B[b] = dict, Lrec[j] = dict)"""
    B, Lrec = {}, {}
    n = bsize = None
    for line in open(path, encoding="utf-8"):
        if line.startswith("# n ="):
            parts = line.replace(",", " ").split()
            n = int(parts[3])
            bsize = int(parts[parts.index("|B_n|") + 2])
            continue
        if line.startswith("#") or not line.strip():
            continue
        p = line.rstrip("\n").split("\t")
        if p[0] == "B":
            B[int(p[1])] = dict(size=int(p[2]), ncyc=int(p[3]), alleven=bool(int(p[4])), pos=p[5], neg=p[6])
        elif p[0] == "L":
            cnt = {}
            for tok in p[5].split():
                b, c = tok.split(":")
                cnt[int(b)] = int(c)
            Lrec[int(p[1])] = dict(order=int(p[2]), orb=int(p[3]), aO=int(p[4]), cnt=cnt)
    return n, bsize, B, Lrec


class Identities:
    """칸 (j, kc) 하나의 치환지표 값을 (b, eps) 마다 돌려준다."""

    def __init__(self, n, bsize, B, Lrec):
        self.n, self.bsize, self.B, self.Lrec = n, bsize, B, Lrec
        self.keys = [(b, e) for b in sorted(B) for e in (0, 1)]
        self.cent = {b: Fraction(bsize, B[b]["size"]) for b in B}
        self.fix = {}
        for b in B:
            self.fix[(b, 0)] = 2 ** B[b]["ncyc"]
            self.fix[(b, 1)] = 2 ** B[b]["ncyc"] if B[b]["alleven"] else 0

    def pi(self, j, kc):
        """{(b, eps): pi} (0 인 것은 뺀다)"""
        Lj = self.Lrec[j]
        out = {}
        if kc == 0:
            for b, c in Lj["cnt"].items():
                out[(b, 0)] = 2 * self.cent[b] * c / Lj["order"]
        else:
            Nk = self.Lrec[kc]
            if 2 * Nk["order"] != Lj["order"]:
                raise ValueError("N 의 위수가 |L|/2 가 아니다: j=%d kc=%d" % (j, kc))
            for b in set(Lj["cnt"]) | set(Nk["cnt"]):
                cN = Nk["cnt"].get(b, 0)
                cL = Lj["cnt"].get(b, 0)
                if cN:
                    out[(b, 0)] = self.cent[b] * cN / Nk["order"]
                if cL - cN:
                    out[(b, 1)] = self.cent[b] * (cL - cN) / Nk["order"]
                if cN > cL:
                    raise ValueError("N 이 L 보다 원소가 많다")
        return out


def rank_solve(M, rhs, nunk):
    """정확한 유리수 가우스-조던. 반환: (모순 행 수, 피벗 열 -> (행, 상수), 자유 열 목록)"""
    rows = [(dict(r), Fraction(b)) for r, b in zip(M, rhs)]
    piv = {}
    inconsistent = 0
    for r, b in rows:
        for c in list(r):
            if c in piv and c in r:
                pr, pb = piv[c]
                f = r[c]
                for cc, a in pr.items():
                    v = r.get(cc, Fraction(0)) - f * a
                    if v == 0:
                        r.pop(cc, None)
                    else:
                        r[cc] = v
                b -= f * pb
        if not r:
            if b != 0:
                inconsistent += 1
            continue
        c0 = min(r)
        a = r[c0]
        r = {c: v / a for c, v in r.items()}
        b = b / a
        for c, (pr, pb) in list(piv.items()):
            if c0 in pr:
                f = pr[c0]
                for cc, aa in r.items():
                    v = pr.get(cc, Fraction(0)) - f * aa
                    if v == 0:
                        pr.pop(cc, None)
                    else:
                        pr[cc] = v
                piv[c] = (pr, pb - f * b)
        piv[c0] = (r, b)
    free = [c for c in range(nunk) if c not in piv]
    return inconsistent, piv, free


def ea_equations(n, bsize, gens_path, C, W, needed):
    """B_n x 1 의 부분군 U = C2^2, C2^3 (B_n-켤레류마다)에 대한 항등식 행.
    반환: [(이름, {칸 (j, 서명): 계수}, 우변 2^{orb(U)})]"""
    import bn_subgroup_counts as BS
    import numpy as np
    B = BS.Bn(n)
    E = BS.ElementaryAbelian(B, verbose=True)
    masks = {}
    nread = 0
    for line in open(gens_path, encoding="utf-8"):
        if not line.startswith("G\t"):
            continue
        p = line.rstrip("\n").split("\t")
        j, order, orb, ng = int(p[1]), int(p[2]), int(p[3]), int(p[4])
        vals = [int(t) for t in p[5:]]
        gens = [B.from_images(vals[k * (n + 1):(k + 1) * (n + 1)]) for k in range(ng)]
        seen = B.closure(gens)
        if int(seen.sum()) != order:
            raise SystemExit("생성원이 만든 군의 위수가 다르다: j=%d, %d != %d" % (j, int(seen.sum()), order))
        if B.orbits_count(gens) != orb:
            raise SystemExit("꼭짓점 궤도 수가 다르다: j=%d" % j)
        masks[j] = seen[E.invs]
        nread += 1
    miss = [j for j in needed if j not in masks]
    if miss:
        raise SystemExit("생성원이 없는 류가 있다: %s" % miss[:10])
    print("   생성원 %d개 류를 읽어 위수와 꼭짓점 궤도 수를 대조했다 (모두 같다)" % nread)
    rows = []
    for r, label in ((2, "C2^2"), (3, "C2^3")):
        d = E.sub[r]
        ncls = len(d["sizes"])
        # 각 류의 대표 부분군의 꼭짓점 궤도 수와 |N_B(U)|
        orbU, normU = [], []
        first = {}
        for i, u in enumerate(d["lab"]):
            first.setdefault(int(u), i)
        for u in range(ncls):
            invs = [int(E.invs[k]) for k in d["T"][first[u]]]
            orbU.append(B.orbits_count(invs))
            if bsize % int(d["sizes"][u]):
                raise SystemExit("류 크기가 |B_n| 의 약수가 아니다")
            normU.append(bsize // int(d["sizes"][u]))
        cnt = {}
        for j in needed:
            cnt[j] = E.counts_in(r, masks[j])
        for u in range(ncls):
            coef = {}
            for (j, sig) in W:
                kc = sig[0]
                if kc == 0:
                    s_ = int(cnt[j][u])
                    if s_:
                        coef[(j, sig)] = Fraction(2 * normU[u] * s_, C[j]["order"])
                else:
                    s_ = int(cnt[kc][u])
                    if s_:
                        coef[(j, sig)] = Fraction(normU[u] * s_, C[kc]["order"])
            rows.append(("%s #%d" % (label, u), coef, 2 ** orbU[u]))
    return rows


def main():
    pos = [a for a in sys.argv[1:] if not a.startswith("--")]
    n, dpath, cpath = int(pos[0]), pos[1], pos[2]
    fixed = None
    orb_min = None
    gens_path = None
    for a in sys.argv[1:]:
        if a.startswith("--gens="):
            gens_path = a[len("--gens="):]
        if a.startswith("--fix="):
            fixed = {}
            for fp in a[len("--fix="):].split(","):
                fixed.update(S.load_fixed(fp))
        if a.startswith("--unknown-orb-min="):
            orb_min = int(a[len("--unknown-orb-min="):])

    C, P, A = S.load(dpath)
    n2, bsize, B, Lrec = load_classes(cpath)
    if n2 != n:
        raise SystemExit("n 이 다르다: %d %d" % (n, n2))
    print("n = %d, |B_n| = %d, B_n 의 켤레류 %d개, 분포를 내보낸 부분군 류 %d개" % (n, bsize, len(B), len(Lrec)))

    # 번호 대조: 같은 표지 표인지 (위수, 꼭짓점 궤도 수, aO)
    bad = [j for j in Lrec if j not in C or (Lrec[j]["order"], Lrec[j]["orb"], Lrec[j]["aO"]) !=
           (C[j]["order"], C[j]["orb"], C[j]["aO"])]
    print("번호 대조 (위수, 꼭짓점 궤도 수, aO) 가 자료와 다른 류: %d개 / %d개" % (len(bad), len(Lrec)))
    if bad:
        raise SystemExit("부분군 류 번호가 자료와 맞지 않는다: %s" % bad[:10])
    if sum(B[b]["size"] for b in B) != bsize:
        raise SystemExit("B_n 의 류 크기 합이 |B_n| 이 아니다")

    # 풀이기의 값
    res = S.solve(C, P, A, False, False, verbose=False, fixed=fixed)
    vars_, det = res["vars"], res["determined"]
    if len(det) != len(vars_):
        raise SystemExit("풀이기가 모든 칸을 정하지 못했다")
    W = {vars_[i]: int(det[i]) for i in det}
    tot_ncls = sum(W.values())
    print("칸 %d개, NPN 류 수 합 %d, 풀이기의 값은 모두 정수" % (len(W), tot_ncls))

    I = Identities(n, bsize, B, Lrec)
    needed = set()
    for (j, sig) in W:
        needed.add(j)
        if sig[0]:
            needed.add(sig[0])
    miss = [j for j in needed if j not in Lrec]
    if miss:
        raise SystemExit("분포가 없는 류가 있다: %s" % miss[:10])

    # 1. 풀이기의 값이 항등식을 모두 만족하는가
    pis = {}
    lhs = defaultdict(Fraction)
    for (j, sig), w in W.items():
        if w == 0:
            continue
        key = (j, sig[0])
        if key not in pis:
            pis[key] = I.pi(j, sig[0])
        for k, v in pis[key].items():
            lhs[k] += w * v
    nbad = 0
    nonint = 0
    for k in I.keys:
        if lhs[k] != I.fix[k]:
            nbad += 1
            if nbad <= 10:
                print("   어긋남 (b=%d, eps=%d): 합 %s, Fix %s" % (k[0], k[1], lhs[k], I.fix[k]))
        if lhs[k].denominator != 1:
            nonint += 1
    print("[1] G_n 의 켤레류 %d개 항등식: 어긋난 것 %d개 (합이 정수가 아닌 것 %d개)" % (len(I.keys), nbad, nonint))
    # 이 항등식들의 가중합이 NPN 류 수 (번사이드)
    burn = sum(Fraction(B[b]["size"] * I.fix[(b, e)]) for b in B for e in (0, 1)) / (2 * bsize)
    print("    번사이드: sum |class| Fix / |G| = %s (표의 NPN 류 수 %d)" % (burn, tot_ncls))

    # 1b. B_n x 1 의 부분군 U = C2^2, C2^3 에 대한 항등식 (표지 표 없이 직접 센 부분군 수)
    ea_rows = []        # (이름, {칸: 계수}, 우변)
    if gens_path:
        ea_rows = ea_equations(n, bsize, gens_path, C, W, needed)
        lhs_bad = 0
        for name, coef, rhs_v in ea_rows:
            tot = sum(W[c] * v for c, v in coef.items())
            if tot != rhs_v:
                lhs_bad += 1
                if lhs_bad <= 10:
                    print("   어긋남 (%s): 합 %s, Fix %s" % (name, tot, rhs_v))
        print("[1b] C2^2, C2^3 부분군 항등식 %d개: 어긋난 것 %d개" % (len(ea_rows), lhs_bad))
        nbad += lhs_bad

    # 2. 미지수 칸을 이 항등식만으로 정한다
    if orb_min is not None:
        unknown = [v for v in vars_ if C[v[0]]["orb"] >= orb_min]
        what = "꼭짓점 궤도 수 >= %d 인 칸" % orb_min
    elif fixed is not None:
        known_j = set(j for (j, _) in fixed)
        unknown = [v for v in vars_ if v[0] not in known_j]
        what = "직접 센 칸이 아닌 칸"
    else:
        return
    ukey = {v: i for i, v in enumerate(unknown)}
    kn_lhs = defaultdict(Fraction)
    for (j, sig), w in W.items():
        if (j, sig) in ukey or w == 0:
            continue
        for k, v in pis[(j, sig[0])].items():
            kn_lhs[k] += w * v
    M, rhs = [], []
    for k in I.keys:
        row = {}
        for v, i in ukey.items():
            val = pis.get((v[0], v[1][0]))
            if val is None:
                val = pis[(v[0], v[1][0])] = I.pi(v[0], v[1][0])
            if k in val:
                row[i] = row.get(i, Fraction(0)) + val[k]
        M.append(row)
        rhs.append(I.fix[k] - kn_lhs[k])
    n_cls_rows = len(M)
    for name, coef, rhs_v in ea_rows:
        row = {}
        kn = Fraction(0)
        for c, v in coef.items():
            if c in ukey:
                row[ukey[c]] = row.get(ukey[c], Fraction(0)) + v
            else:
                kn += W[c] * v
        M.append(row)
        rhs.append(rhs_v - kn)
    incons, piv, free = rank_solve(M, rhs, len(unknown))
    nz = len([r for r in M if r])
    print("[2] %s: %d칸 (Kbar 류 %d개), 방정식 %d개 (켤레류 항등식 %d개 + 부분군 항등식 %d개; 미지수가 들어 있는 것 %d개)"
          % (what, len(unknown), len(set(v[0] for v in unknown)), len(M), n_cls_rows, len(M) - n_cls_rows, nz))
    print("    계수 행렬의 계수(rank) %d / 미지수 %d, 모순 행 %d개" % (len(piv), len(unknown), incons))
    determined = {}
    for c, (r, b) in piv.items():
        if all(cc == c or cc in piv for cc in r):
            determined[c] = b
    nint = sum(1 for x in determined.values() if x.denominator != 1)
    neg = sum(1 for x in determined.values() if x < 0)
    dif = [(unknown[i], determined[i], W[unknown[i]]) for i in determined if determined[i] != W[unknown[i]]]
    print("    항등식만으로 정해진 칸 %d개 / %d개 (정수 아님 %d, 음수 %d); 풀이기의 값과 다른 칸 %d개"
          % (len(determined), len(unknown), nint, neg, len(dif)))
    for d in dif[:10]:
        print("      다름:", d)
    und = [i for i in range(len(unknown)) if i not in determined]
    if und:
        print("    항등식만으로는 따로 정해지지 않는 칸 %d개:" % len(und))
        for i in und[:20]:
            j, sig = unknown[i]
            print("      j=%d (위수 %d, 궤도 %d) 서명 %s, 풀이기 %d" % (j, C[j]["order"], C[j]["orb"], sig, W[unknown[i]]))
        # 정해지지 않은 칸들의 합 조합이 풀이기의 값과 맞는지 (풀이기의 값이 방정식을 만족하므로 항상 참이지만 점검)
    sol_ok = all(sum(a * W[unknown[c]] for c, a in M[r].items()) == rhs[r] for r in range(len(M)))
    print("    풀이기의 값이 이 방정식 %d개를 모두 만족하는가: %s" % (len(M), sol_ok))
    ok = (nbad == 0 and incons == 0 and not dif and sol_ok)
    print("결과:", "일치" if ok else "불일치")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
