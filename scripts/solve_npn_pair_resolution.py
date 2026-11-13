#!/usr/bin/env python3
"""
solve_npn_pair_resolution.py
============================
scripts/npn_pair_resolution.g 가 내보낸 자료로 NPN 안정화군의 세분(어느 Kbar 가 어느 phi 를 갖는 류를 몇 개 갖는가)을 푼다.

NPN 류 하나의 안정화군은 (L, phi): L 은 B_n 의 부분군 (Kbar), phi 는 L 에서 Z_2 로 가는 준동형이다 (반대칭이 없으면 phi = 0).
phi 마다 "서명" sig = (kerClass, t1, t2, t3) 를 적었다. 미지수 W[j, sig] = Kbar 가 켤레류 j 이고 서명이 sig 인 NPN 류의 수. 방정식:
  (D)  sum_sig W[j, sig] = decY[j]                               (Kbar 가 j 인 NPN 류의 수)
  (E0) sum W[j, sig] c0(j, sig, x) = aO[x]                       (c0 = 2 if kerClass = 0 and j = x; 1 if kerClass = x)
  (Ek) sum W[j, sig] ck(j, sig, x) = aOpsi_k[x]   (k = 1, 2, 3)  (ck = 2 if t_k = -1 and j = x; 1 if t_k = x)
  (SYM) W[j, sig] = W[j, sigma(sig)]                              (패리티 대칭 f -> f xor (x_1 xor ... xor x_n))
E0 ~ E3 은 B_n 의 꼬인 작용(psi_0 = 0, psi_1, psi_2, psi_3)의 궤도 분해이다 (npn_pair_resolution.g 머리말 참고).
--basic 이면 (D), (E0) 만 쓰고 phi 를 N = ker phi 의 켤레류로만 구분한다 (서명의 첫 성분만).

  python3 -I scripts/solve_npn_pair_resolution.py <자료.tsv> [--truth] [--basic] [--lp]
    --truth : 자료에 A 줄(B_n x Z_2 표지 표의 참값)이 있으면 참값이 방정식을 만족하는지, 정해진 미지수가 참값과 같은지 본다.
    --lp    : 정해지지 않은 미지수의 범위를 선형계획(scipy)으로 구한다 (실수 근사이므로 참고용).
    --fix=파일 : combine_enum_fix.py 가 쓴 고정값(j s0 s1 s2 s3 W, 한 줄에 하나)을 방정식으로 더한다
                (작은 2-군 몇 개의 W 를 대표 함수의 정확한 대칭군을 직접 세어 얻은 값이다). --basic 과 함께 쓰지 않는다.

정확한 유리수 가우스-조던 소거(fractions)와 "음이 아닌 미지수의 합이 0 이면 모두 0" 전파만 쓴다.
"""
import sys
from fractions import Fraction
from collections import defaultdict


def load_Q(path):
    """Q 줄: j -> 유효한 phi != 0 의 N_B(L)-궤도 수"""
    Q = {}
    for line in open(path, encoding="utf-8"):
        if line.startswith("Q\t"):
            p = line.rstrip("\n").split("\t")
            Q[int(p[1])] = int(p[2])
    return Q


def check_signatures_separate(C, P, Q):
    """서명이 K 의 켤레류(= (Kbar, N_B(L)-궤도 of phi))를 구분하는가: j 마다 서명의 수 = 1 + 궤도 수"""
    bad = []
    for j in sorted(C):
        if C[j]["decY"] <= 0:
            continue
        nsig = len(set(P[j]))
        if nsig != 1 + Q.get(j, 0):
            bad.append((j, nsig, 1 + Q.get(j, 0)))
    return bad


def load(path):
    C, P, A = {}, defaultdict(list), []
    for line in open(path, encoding="utf-8"):
        if line.startswith("#") or not line.strip():
            continue
        p = line.rstrip("\n").split("\t")
        if p[0] == "C":
            j = int(p[1])
            valid = [] if p[7] == "-" else [int(t) for t in p[7].split(",")]
            d = dict(order=int(p[2]), orb=int(p[3]), aO=int(p[4]), decY=int(p[5]), ng=int(p[6]), valid=valid)
            d["aOpsi"] = [int(p[8]), int(p[9]), int(p[10])] if len(p) >= 11 else None
            C[j] = d
        elif p[0] == "P":
            P[int(p[1])].append(tuple(int(t) for t in p[2:6]))
        elif p[0] == "A":
            A.append((int(p[1]), int(p[2]), int(p[3]), int(p[4]), tuple(int(t) for t in p[5:9]) if len(p) >= 9 else None))
    return C, P, A


def build(C, P, basic):
    """변수는 (j, sig). basic 이면 sig = (kerClass,)"""
    occ = [j for j in sorted(C) if C[j]["decY"] > 0]
    vars_ = []
    for j in occ:
        sigs = sorted(set((s[:1] if basic else s) for s in P[j]))
        for s in sigs:
            vars_.append((j, s))
    idx = {v: i for i, v in enumerate(vars_)}
    eqs = []
    for j in occ:
        eqs.append(({idx[v]: 1 for v in vars_ if v[0] == j}, C[j]["decY"], ("D", j)))
    if not basic:
        # 대칭: f -> f xor (x_1 xor ... xor x_n) 는 (L, phi) 인 NPN 류를 (L, phi + psi2|L) 인 NPN 류로 옮기는 전단사이다
        # (패리티 함수 c 는 c o g = c xor psi2(g) 를 모든 g in B_n 에 대해 만족한다). 서명에서 phi + psi2|L 의 서명은
        # (t2 == -1 이면 0, 아니면 t2),  t3,  (kerClass == 0 이면 -1, 아니면 kerClass),  t1  이다.
        for (j, s), i in idx.items():
            kc, t1, t2, t3 = s
            s2 = (0 if t2 == -1 else t2, t3, -1 if kc == 0 else kc, t1)
            if (j, s2) not in idx:
                raise ValueError("대칭으로 옮긴 서명이 변수에 없다: %s %s -> %s" % (j, s, s2))
            i2 = idx[(j, s2)]
            if i2 > i:
                eqs.append(({i: 1, i2: -1}, 0, ("SYM", j)))
    nfam = 1 if basic else 4
    for k in range(nfam):
        rows = defaultdict(dict)
        for (j, s), i in idx.items():
            t = s[0] if k == 0 else s[k]
            if k == 0:
                if t == 0:
                    rows[j][i] = rows[j].get(i, 0) + 2
                else:
                    rows[t][i] = rows[t].get(i, 0) + 1
            else:
                if t == -1:
                    rows[j][i] = rows[j].get(i, 0) + 2
                else:
                    rows[t][i] = rows[t].get(i, 0) + 1
        xs = set(rows) | set(x for x in C if (C[x]["aO"] if k == 0 else C[x]["aOpsi"][k - 1]) != 0)
        for x in sorted(xs):
            rhs = C[x]["aO"] if k == 0 else C[x]["aOpsi"][k - 1]
            eqs.append((rows.get(x, {}), rhs, ("E%d" % k, x)))
    return vars_, idx, eqs


def propagate(nv, eqs):
    val = [None] * nv
    rows = [(dict(r), Fraction(b), nm) for r, b, nm in eqs]
    changed = True
    problems = []
    while changed:
        changed = False
        newrows = []
        for r, b, nm in rows:
            for c in list(r):
                if val[c] is not None:
                    b -= r[c] * val[c]
                    del r[c]
            if not r:
                if b != 0:
                    problems.append((nm, b))
                continue
            positive = all(a > 0 for a in r.values())   # 계수가 모두 양수인 행에서만 "음이 아니다" 로 추론한다
            if positive and b < 0:
                problems.append((nm, b))
            if positive and b == 0:
                for c in r:
                    if val[c] is None:
                        val[c] = Fraction(0)
                        changed = True
                continue
            if len(r) == 1:
                (c, a), = r.items()
                if val[c] is None:
                    val[c] = b / a
                    if val[c] < 0:
                        problems.append((nm, b))
                    changed = True
                continue
            newrows.append((r, b, nm))
        rows = newrows
    return val, rows, problems


def rref(rows):
    piv = {}
    work = [({c: Fraction(v) for c, v in r.items()}, Fraction(b)) for r, b, _ in rows]
    work.sort(key=lambda t: len(t[0]))
    for r, b in work:
        for c in list(r):
            if c in piv and c in r:
                pr, pb = piv[c]
                f = r[c]
                for cc, a in pr.items():
                    nvv = r.get(cc, Fraction(0)) - f * a
                    if nvv == 0:
                        r.pop(cc, None)
                    else:
                        r[cc] = nvv
                b -= f * pb
        if not r:
            if b != 0:
                raise ValueError("모순된 방정식")
            continue
        c0 = min(r)
        a = r[c0]
        r = {c: v / a for c, v in r.items()}
        b = b / a
        for c, (pr, pb) in list(piv.items()):
            if c0 in pr:
                f = pr[c0]
                for cc, aa in r.items():
                    nvv = pr.get(cc, Fraction(0)) - f * aa
                    if nvv == 0:
                        pr.pop(cc, None)
                    else:
                        pr[cc] = nvv
                piv[c] = (pr, pb - f * b)
        piv[c0] = (r, b)
    return piv


def load_fixed(path):
    fixed = {}
    for line in open(path, encoding="utf-8"):
        p = line.split()
        if len(p) == 6:
            fixed[(int(p[0]), tuple(int(t) for t in p[1:5]))] = int(p[5])
    return fixed


def solve(C, P, A, basic, use_truth, verbose=True, fixed=None):
    vars_, idx, eqs = build(C, P, basic)
    if fixed:
        if basic:
            raise ValueError("--fix 는 --basic 과 함께 쓸 수 없다")
        for key, w in sorted(fixed.items()):
            if key not in idx:
                raise ValueError("고정값의 변수가 모형에 없다: %s" % (key,))
            eqs.append(({idx[key]: 1}, w, ("FIX", key[0])))
    nv = len(vars_)
    totNPN = sum(C[j]["decY"] for j in C)
    totNP = sum(C[j]["aO"] for j in C)
    if verbose:
        print("켤레류 %d개, Kbar 로 나타나는 것 %d개, 미지수 %d개, 방정식 %d개 (%s)"
              % (len(C), sum(1 for j in C if C[j]['decY'] > 0), nv, len(eqs), "기본 모형 D+E0" if basic else "꼬인 작용 포함 D+E0..E3+SYM"))
        print("NPN 류 수 %d, B_n-궤도 수(NP) %d, 반대칭 류 수 2 NPN - NP = %d" % (totNPN, totNP, 2 * totNPN - totNP))
    truth = None
    if use_truth and A:
        truth = [Fraction(0)] * nv
        miss = 0
        for ik, ih, an, cnt, sg in A:
            key = (ik, sg[:1] if basic else sg)
            if key not in idx:
                miss += 1
                print("참값의 변수가 모형에 없다:", key)
                continue
            truth[idx[key]] += cnt
        nbad = 0
        for r, b, nm in eqs:
            s = sum(a * truth[c] for c, a in r.items())
            if s != b:
                nbad += 1
                if nbad <= 5:
                    print("   참값이 만족하지 않는 방정식", nm, s, b)
        if verbose:
            print("참값(방법 A) 대조: 모형에 없는 변수 %d개, 참값이 만족하지 않는 방정식 %d개 / %d개" % (miss, nbad, len(eqs)))
    val, rows, problems = propagate(nv, eqs)
    if verbose:
        print("전파: 정해진 미지수 %d개, 남은 방정식 %d개, 모순 %d개" % (sum(1 for x in val if x is not None), len(rows), len(problems)))
        for pr in problems[:10]:
            print("   모순:", pr)
    free_vars = [i for i in range(nv) if val[i] is None]
    piv = rref(rows)
    pivcols = set(piv)
    freecols = [c for c in free_vars if c not in pivcols]
    determined = {i: val[i] for i in range(nv) if val[i] is not None}
    for c, (r, b) in piv.items():
        if all(cc == c or cc in pivcols for cc in r):
            determined[c] = b
    if verbose:
        print("남은 미지수 %d개, 피벗 %d개, 자유도 %d" % (len(free_vars), len(pivcols), len(freecols)))
        print("선형대수와 전파로 정해진 미지수 %d개 / %d개" % (len(determined), nv))
        nonint = [i for i, x in determined.items() if x.denominator != 1]
        neg = [i for i, x in determined.items() if x < 0]
        print("   정해진 값 가운데 정수가 아닌 것 %d개, 음수인 것 %d개" % (len(nonint), len(neg)))
    # 반대칭 류의 몇 퍼센트가 정해졌는가 (phi != 0 인 변수)
    anti_all = 2 * totNPN - totNP
    s_det = sum(x for i, x in determined.items() if vars_[i][1][0] != 0)
    nv_anti = sum(1 for v in vars_ if v[1][0] != 0)
    nd_anti = sum(1 for i in determined if vars_[i][1][0] != 0)
    if verbose:
        print("   정해진 반대칭 류 %s / %d (%.4f%%); 정해진 변수 %d개 / %d개 (phi != 0 인 변수)"
              % (s_det, anti_all, 100.0 * float(s_det) / anti_all if anti_all else 100.0, nd_anti, nv_anti))
    if truth is not None:
        mism = [i for i, x in determined.items() if truth[i] != x]
        if verbose:
            print("정해진 값이 참값과 다른 미지수: %d개" % len(mism))
            for i in mism[:10]:
                print("   ", vars_[i], determined[i], truth[i])
    und = [i for i in range(nv) if i not in determined]
    if verbose:
        print("정해지지 않은 미지수 %d개 (속한 켤레류 j: %d개)" % (len(und), len(set(vars_[i][0] for i in und))))
    return dict(vars=vars_, idx=idx, eqs=eqs, determined=determined, und=und, rows=rows, val=val, truth=truth,
                totNPN=totNPN, totNP=totNP)


def lp_bounds(res):
    import numpy as np
    from scipy.optimize import linprog
    vars_, rows, val, und_all, truth = res["vars"], res["rows"], res["val"], res["und"], res["truth"]
    free = [i for i in range(len(vars_)) if val[i] is None]
    pos = {c: k for k, c in enumerate(free)}
    Aeq = np.zeros((len(rows), len(free)))
    beq = np.zeros(len(rows))
    for r, (row, b, nm) in enumerate(rows):
        for c, a in row.items():
            Aeq[r, pos[c]] = float(a)
        beq[r] = float(b)
    out = {}
    for c in und_all:
        if c not in pos:
            continue
        k = pos[c]
        cost = np.zeros(len(free))
        cost[k] = 1
        lo = linprog(cost, A_eq=Aeq, b_eq=beq, bounds=(0, None), method="highs")
        hi = linprog(-cost, A_eq=Aeq, b_eq=beq, bounds=(0, None), method="highs")
        if lo.status != 0 or hi.status != 0:
            out[c] = None
            continue
        out[c] = (lo.fun, -hi.fun)
    return out


def main():
    path = sys.argv[1]
    use_truth = "--truth" in sys.argv
    basic = "--basic" in sys.argv
    C, P, A = load(path)
    print("자료:", path)
    if not basic and C[next(iter(C))]["aOpsi"] is None:
        print("이 자료에는 꼬인 작용의 분해가 없다. --basic 으로 돌린다.")
        basic = True
    fixed = None
    for a in sys.argv[2:]:
        if a.startswith("--fix="):
            fixed = {}
            for fp in a[len("--fix="):].split(","):
                fixed.update(load_fixed(fp))
            print("고정값 %d개를 방정식으로 더한다" % len(fixed))
    res = solve(C, P, A, basic, use_truth, fixed=fixed)
    Q = load_Q(path)
    if Q:
        bad = check_signatures_separate(C, P, Q)
        print("서명이 K 의 켤레류를 구분하는가 (j 마다 서명의 수 = 1 + N_B(L)-궤도 수): 어긋난 j %d개 / %d개" % (len(bad), len(Q)))
        for b in bad[:10]:
            print("   j = %d: 서명 %d개, 1 + 궤도 수 = %d" % b)
    if "--lp" in sys.argv and res["und"]:
        b = lp_bounds(res)
        tight = [c for c, r in b.items() if r is not None and r[1] - r[0] <= 1e-6 * max(1.0, abs(r[1]))]
        print("선형계획(실수 근사): 범위가 한 점인 미지수 %d개 / %d개" % (len(tight), len(res["und"])))
        if res["truth"] is not None:
            bad = [c for c, r in b.items() if r is not None and not (r[0] - 1e-6 * max(1.0, abs(r[0])) <= float(res["truth"][c]) <= r[1] + 1e-6 * max(1.0, abs(r[1])))]
            print("   참값이 범위 밖인 미지수: %d개" % len(bad))


if __name__ == "__main__":
    main()
