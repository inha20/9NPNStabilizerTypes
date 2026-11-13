#!/usr/bin/env python3
"""
make_npn_k_table.py
===================
scripts/npn_pair_resolution.g 의 자료(C, P, D 줄)를 scripts/solve_npn_pair_resolution.py 로 풀어 NPN 안정화군 K(f) 의 표를 만든다.
B_n x Z_2 의 표지 표를 쓰지 않고, B_n 의 표지 표 하나(원래의 것과 꼬인 작용의 것)와 패리티 대칭만으로 K 의 종류별 NPN 류 수를 구한다.
n <= 5 에서는 npn_stabilizer_types_n<n>.tsv (B_n x Z_2 의 표지 표로 만든 참값)와 행 단위로 견주고, n = 6 에서는 표를 새로 만든다.

  python3 -I scripts/make_npn_k_table.py <자료.tsv> <n> <출력.tsv> [<참값 npn_stabilizer_types_n<n>.tsv>]
        [--fix=고정값1,고정값2,...] [--split=열거결과1,열거결과2,...]

K 의 종류는 (Kbar 의 켤레류 j, 서명)으로 구분한다. 서명이 K 의 켤레류를 구분하는지는 n <= 5 에서 확인했다(K 의 켤레류 수와 같다).
n = 6 에서는 Kbar 의 켤레류 j 하나(이 실행에서는 132)에서만 서명이 K 의 켤레류를 구분하지 못한다(서명 3개, K 의 켤레류 4개).
그 j 의 행은 머리말에 적는다.
--fix : 대표 함수의 정확한 대칭군을 직접 세어 얻은 값(combine_enum_fix.py 가 쓴 파일, 쉼표로 여러 개)을 연립방정식에 더한다.
        n = 6 에서는 이것이 있어야 표가 만들어진다(방정식만으로는 576칸이 남는다). 방정식과 어긋나면 모순으로 멈춘다.
--split : 서명이 K 의 켤레류를 합치는 j 의 행을, 열거 결과(R 줄)의 phi 별 값으로 갈라 K 의 켤레류마다 한 행으로 쓴다
          (그 j 의 nonzero phi 가 모두 서로 다른 N_B(L)-궤도일 때만 가능하다).
"""
import sys
from collections import Counter, defaultdict
from fractions import Fraction

sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
import solve_npn_pair_resolution as S


def load_D(path):
    D = {}
    for line in open(path, encoding="utf-8"):
        if line.startswith("D\t"):
            p = line.rstrip("\n").split("\t")
            D[int(p[1])] = dict(order=int(p[2]), varorb=p[3], d=int(p[4]), vorb=int(p[5]), fixv=int(p[6]), struct=p[7])
    return D


def load_R_per_phi(paths):
    """R 줄 -> {j: [(idx, sig, E, |L|, |N_B(L)|), ...]}"""
    out = defaultdict(list)
    for path in paths:
        for line in open(path, encoding="utf-8"):
            p = line.split()
            if not p or p[0] != "R" or p[-1] == "SKIP":
                continue
            out[int(p[1])].append((int(p[2]), tuple(int(t) for t in p[3:7]), int(p[8]), int(p[10]), int(p[11])))
    return out


def main():
    pos = [a for a in sys.argv[1:] if not a.startswith("--")]
    path, n, outp = pos[0], int(pos[1]), pos[2]
    truthp = pos[3] if len(pos) > 3 else None
    fixed = None
    for a in sys.argv[1:]:
        if a.startswith("--fix="):
            fixed = {}
            for fp in a[len("--fix="):].split(","):
                fixed.update(S.load_fixed(fp))
    split_paths = []
    for a in sys.argv[1:]:
        if a.startswith("--split="):
            split_paths = a[len("--split="):].split(",")
    C, P, A = S.load(path)
    D = load_D(path)
    res = S.solve(C, P, A, False, bool(A), verbose=True, fixed=fixed)
    vars_ = res["vars"]
    det = res["determined"]
    if len(det) != len(vars_):
        print("정해지지 않은 미지수가 %d개 있어 표를 만들지 않는다" % (len(vars_) - len(det)))
        sys.exit(1)
    bsize = None
    # |B_n| = 2^n n!
    import math
    bsize = 2 ** n * math.factorial(n)
    rows = []
    for i, (j, sig) in enumerate(vars_):
        w = det[i]
        if w < 0 or w.denominator != 1:
            print("정수가 아니거나 음수인 값:", vars_[i], w)
            sys.exit(1)
        w = int(w)
        if w == 0:
            continue
        kc = sig[0]
        anti = 0 if kc == 0 else 1
        hj = j if kc == 0 else kc
        ordK = C[j]["order"]
        ordH = C[hj]["order"]
        index = 2 * bsize // ordK
        rows.append(dict(ordK=ordK, index=index, funcs=w * index, ncls=w, anti=anti, ordKb=C[j]["order"], ordH=ordH,
                         j=j, hj=hj, sig=sig))
    Qd = S.load_Q(path)
    merged = [j for j in sorted(Qd) if C[j]["decY"] > 0 and len(set(P[j])) != 1 + Qd[j]]
    nsplit = 0
    if split_paths and merged:
        perphi = load_R_per_phi(split_paths)
        for j in merged:
            if Qd[j] != C[j]["ng"]:
                print("K 를 가를 수 없다: j = %d 의 nonzero phi 가 서로 다른 궤도가 아니다" % j)
                sys.exit(1)
            if j not in perphi:
                print("열거 결과에 j = %d 가 없다" % j)
                sys.exit(1)
            new = []
            for r in [r for r in rows if r["j"] == j]:
                cand = [(idx, sg, E, L, NB) for (idx, sg, E, L, NB) in perphi[j] if tuple(sg) == tuple(r["sig"])]
                tot = sum(Fraction(L * E, 2 * NB) for (_, _, E, L, NB) in cand)
                if tot != r["ncls"]:
                    print("갈랐을 때 합이 맞지 않는다: j = %d, 서명 %s, 합 %s, 표 %d" % (j, r["sig"], tot, r["ncls"]))
                    sys.exit(1)
                if len(cand) > 1:
                    for (idx, sg, E, L, NB) in cand:
                        w = Fraction(L * E, 2 * NB)
                        if w.denominator != 1:
                            print("정수가 아니다", j, sg, w)
                            sys.exit(1)
                        if w > 0:
                            q = dict(r)
                            q["ncls"] = int(w)
                            q["funcs"] = int(w) * r["index"]
                            new.append(q)
                    nsplit += 1
                    rows.remove(r)
            rows.extend(new)
    rows.sort(key=lambda r: (-r["ordK"], r["anti"], -r["ncls"], r["j"], r["hj"]))
    totF = sum(r["funcs"] for r in rows)
    totN = sum(r["ncls"] for r in rows)
    totA = sum(r["ncls"] for r in rows if r["anti"])
    print("K 의 종류 %d개 (반대칭 %d개), NPN 류 수 %d, 함수 수 합 %d (2^(2^n) = %d), 반대칭 류 수 %d"
          % (len(rows), sum(1 for r in rows if r["anti"]), totN, totF, 2 ** (2 ** n), totA))
    if totF != 2 ** (2 ** n):
        print("!!! 함수 수의 합이 맞지 않는다")
    with open(outp, "w", encoding="utf-8") as f:
        f.write("# n = %d, |B_n x Z_2| = %d. NPN 안정화군 K(f) 의 종류(Kbar 의 B_n-켤레류와 phi 의 서명으로 구분)와 종류별 함수 수.\n" % (n, 2 * bsize))
        f.write("# 방법: B_n 의 표지 표(원래의 것과 꼬인 작용의 것)와 패리티 대칭의 연립방정식 (scripts/npn_pair_resolution.g, solve_npn_pair_resolution.py).\n")
        if fixed:
            f.write("# 고정값 %d칸은 대표 함수의 정확한 대칭군을 직접 세어 얻은 값이고(scripts/combine_enum_fix.py), 나머지 칸은 연립방정식이 정했다.\n" % len(fixed))
        f.write("# 열: 1 위수 |K| (= |Kbar|), 2 지수(함수 하나의 NPN 궤도 크기), 3 함수 수, 4 NPN 류 수, 5 반대칭 유무(1 이면 phi != 0, |Kbar| = 2 |H|),\n")
        f.write("#     6 |Kbar|, 7 |H|, 8 Kbar 의 변수 위 궤도 길이, 9 Kbar 의 평행이동 부분의 차원 d, 10 Kbar 의 꼭짓점 궤도 수, 11 Kbar 의 고정 꼭짓점 수, 12 구조(Kbar),\n")
        f.write("#     13 H 의 변수 위 궤도 길이, 14 H 의 평행이동 부분의 차원 d, 15 H 의 꼭짓점 궤도 수, 16 H 의 고정 꼭짓점 수, 17 구조(H),\n")
        f.write("#     18 B_n 켤레류 번호(Kbar 의 것, H 의 것; 이 실행 안에서만 뜻이 있음), 19 서명 (kerClass, t1, t2, t3)\n")
        for j in merged:
            if nsplit:
                f.write("# 주의: Kbar 류 %d 에서는 서명 %d개가 K 의 켤레류 %d개를 나타낸다 (서명이 K 의 켤레류를 다 구분하지 못한다). 열거로 구한 phi 별 값으로 갈라, 서명이 같은 두 행을 따로 적었다.\n"
                        % (j, len(set(P[j])), 1 + Qd[j]))
            else:
                f.write("# 주의: Kbar 류 %d 에서는 서명 %d개가 K 의 켤레류 %d개를 나타낸다 (서명이 K 의 켤레류를 다 구분하지 못한다). 이 류의 행은 그 서명 안의 K 를 합친 값이다.\n"
                        % (j, len(set(P[j])), 1 + Qd[j]))
        for r in rows:
            dk, dh = D.get(r["j"]), D.get(r["hj"])
            f.write("\t".join(str(x) for x in [
                r["ordK"], r["index"], r["funcs"], r["ncls"], r["anti"], r["ordKb"], r["ordH"],
                dk["varorb"], dk["d"], dk["vorb"], dk["fixv"], dk["struct"],
                dh["varorb"], dh["d"], dh["vorb"], dh["fixv"], dh["struct"],
                "%d %d" % (r["j"], r["hj"]), "%d %d %d %d" % tuple(r["sig"])]) + "\n")
    print("표를 썼다:", outp)

    if truthp:
        # 참값(B_n x Z_2 의 표지 표)과 (|K|, 반대칭, |Kbar|, |H|, NPN 류 수) 단위로 견준다
        tr = Counter()
        for line in open(truthp, encoding="utf-8"):
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            tr[(int(p[0]), int(p[4]), int(p[5]), int(p[6]), int(p[3]))] += 1
        mine = Counter((r["ordK"], r["anti"], r["ordKb"], r["ordH"], r["ncls"]) for r in rows)
        print("참값과 견주기 (|K|, 반대칭, |Kbar|, |H|, NPN 류 수) 행 단위: 참값 %d행, 이 표 %d행, 다른 행 %d개"
              % (sum(tr.values()), sum(mine.values()), sum(((tr - mine) + (mine - tr)).values())))
        # 위수별 합계도
        a = Counter(); b = Counter()
        for k, v in tr.items():
            a[(k[0], k[1])] += k[4] * v
        for k, v in mine.items():
            b[(k[0], k[1])] += k[4] * v
        print("   위수·반대칭별 NPN 류 수 합계가 다른 곳: %d개" % sum(1 for k in set(a) | set(b) if a[k] != b[k]))


if __name__ == "__main__":
    main()
