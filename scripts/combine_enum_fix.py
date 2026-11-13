#!/usr/bin/env python3
"""
combine_enum_fix.py
===================
sample_npn_stabilizers.c 의 `enumfix` 모드가 낸 결과(R 줄)를 NPN 류의 수 W[j, 서명] 으로 바꾼다.

R 줄:  R j idx s0 s1 s2 s3 o E total |L| |N_B(L)| nchk
  (j, s0..s3) = (켤레류, 서명), o = L 의 꼭짓점 궤도 수, E = Fix(L, phi) 의 함수 가운데 정확한 안정화군이 K = (L, phi) 인 것의 수,
  total = 2^o (전부 센 함수 수).

변환:  W[j, 서명] = |L| / (2 |N_B(L)|) * (그 서명을 갖는 모든 phi 에 대한 E 의 합)
  (K 의 G_n-켤레가 |B_n| / |N_B(L, phi)| 개이고 궤도의 크기가 2 |B_n| / |L| 이므로 NPN 류 수는 |L| E / (2 |N_B(L, phi)|);
   N_B(L)-궤도 위에서 더하면 phi 전체에 대한 합이 된다.)

  python3 -I scripts/combine_enum_fix.py <자료.tsv> <enum 결과...> [--out=고정값.txt] [--compare]
    --compare : 풀이기(solve_npn_pair_resolution.py)가 이미 정한 값과 비교한다.
    --out     : j s0 s1 s2 s3 W 를 한 줄씩 쓴다 (solve_npn_pair_resolution.py --fix=고정값.txt 에 쓴다).
한 j 의 모든 phi 가 들어 있으면 W 의 합이 decY[j] 와 같은지도 본다.
"""
import os
import sys
from fractions import Fraction
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import solve_npn_pair_resolution as S  # noqa: E402


def read_R(paths):
    rec = defaultdict(lambda: dict(E=0, n=0, L=None, NB=None, o=None, total=None))
    for path in paths:
        for line in open(path, encoding="utf-8"):
            p = line.split()
            if not p or p[0] != "R":
                continue
            if p[-1] == "SKIP":
                continue
            j = int(p[1])
            sig = tuple(int(t) for t in p[3:7])
            o, E, total, L, NB = (int(t) for t in p[7:12])
            r = rec[(j, sig)]
            if r["L"] is not None and (r["L"], r["NB"], r["o"]) != (L, NB, o):
                raise ValueError("같은 j 의 L, N_B(L), o 가 다르다: %s" % (j,))
            r["L"], r["NB"], r["o"], r["total"] = L, NB, o, total
            r["E"] += E
            r["n"] += 1
    return rec


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    opts = [a for a in sys.argv[1:] if a.startswith("--")]
    data, files = args[0], args[1:]
    out = None
    for o in opts:
        if o.startswith("--out="):
            out = o[len("--out="):]
    rec = read_R(files)
    C, P, A = S.load(data)
    Q = S.load_Q(data)
    W = {}
    for (j, sig), r in sorted(rec.items()):
        w = Fraction(r["L"] * r["E"], 2 * r["NB"])
        if w.denominator != 1:
            print("정수가 아니다:", j, sig, w)
            sys.exit(2)
        W[(j, sig)] = int(w)
    js = sorted(set(j for j, _ in W))
    print("열거한 군 %d개, (j, 서명) %d개" % (len(js), len(W)))
    # 1. 한 j 의 W 합 = decY[j]
    nbad = 0
    nfull = 0
    for j in js:
        sigs = set(s for (jj, s) in W if jj == j)
        if sigs != set(P[j]):
            continue   # 이 j 의 모든 서명이 열거에 들어 있지 않다
        nfull += 1
        tot = sum(W[(j, s)] for s in sigs)
        if tot != C[j]["decY"]:
            nbad += 1
            print("  합 불일치 j=%d: 열거 %d, decY %d" % (j, tot, C[j]["decY"]))
    print("서명이 모두 들어 있는 j %d개: W 의 합 != decY 인 j %d개" % (nfull, nbad))
    # 서명 집합이 자료와 같은가 (열거에 없는 서명, 자료에 없는 서명)
    for j in js:
        have = set(s for (jj, s) in W if jj == j)
        extra = have - set(P[j])
        if extra:
            print("  자료에 없는 서명 j=%d: %s" % (j, sorted(extra)))
    if "--compare" in opts:
        basic = False
        res = S.solve(C, P, A, basic, False, verbose=False)
        idx, det = res["idx"], res["determined"]
        ncmp = nmis = nund = 0
        for (j, sig), w in sorted(W.items()):
            if (j, sig) not in idx:
                print("  모형에 없는 변수:", j, sig)
                nmis += 1
                continue
            i = idx[(j, sig)]
            if i in det:
                ncmp += 1
                if det[i] != w:
                    nmis += 1
                    if nmis <= 20:
                        print("  불일치 j=%d sig=%s: 풀이 %s, 열거 %d" % (j, sig, det[i], w))
            else:
                nund += 1
        print("풀이기와 비교: 정해진 값과 대조 %d개, 불일치 %d개, 풀이기가 못 정한 변수 %d개" % (ncmp, nmis, nund))
        if nmis:
            sys.exit(1)
    if out:
        with open(out, "w", encoding="utf-8") as fo:
            for (j, sig), w in sorted(W.items()):
                fo.write("%d %d %d %d %d %d\n" % ((j,) + sig + (w,)))
        print("고정값 %d개를 %s 에 썼다" % (len(W), out))


if __name__ == "__main__":
    main()
