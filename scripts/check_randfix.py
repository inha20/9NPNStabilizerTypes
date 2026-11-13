#!/usr/bin/env python3
"""
check_randfix.py
================
열거하지 못한 군(꼭짓점 궤도가 25개 이상)의 W[j, 서명]을, 균일하게 뽑은 표본(sample_npn_stabilizers.c randfix)으로 견준다.

(L, phi) 마다 L 이 phi 로 고정하는 함수 2^o 개 가운데 안정화군이 정확히 L 인 것의 비율 p = E / 2^o 를 표본으로 어림한다 (Q 줄).
한편 W[j, s] = |L| / (2 |N_B(L)|) * (서명 s 인 phi 에 대한 E 의 합) 이므로 풀이 결과의 W 가 같은 합 E 를 예측한다.
같은 (j, 서명) 의 phi 들에 대해 표본의 sum p̂ 와 예측 sum p 를 z 값으로 견준다 (표준오차 = sqrt(sum p̂ (1 - p̂) / N)).
정확한 값이 아니라 통계적 대조이다. 표본 수가 충분하면 W 가 상대오차 1e-5 쯤까지 맞는지 본다.

  python3 -I scripts/check_randfix.py <자료.tsv> --fix=고정값1,고정값2,... <randfix 출력...> [--zmax=4.5]

고정값 파일들로 풀이기를 돌려 W 를 얻는다 (고정값에 들어 있지 않은 군의 W 는 연립방정식이 정한다).
"""
import os
import sys
import math
from fractions import Fraction
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import solve_npn_pair_resolution as S  # noqa: E402


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    data, qfiles = args[0], args[1:]
    fixed, zmax = {}, 4.5
    for a in sys.argv[1:]:
        if a.startswith("--fix="):
            for fp in a[len("--fix="):].split(","):
                fixed.update(S.load_fixed(fp))
        elif a.startswith("--zmax="):
            zmax = float(a[len("--zmax="):])
    C, P, A = S.load(data)
    res = S.solve(C, P, A, False, False, verbose=False, fixed=fixed)
    idx, det = res["idx"], res["determined"]
    groups = defaultdict(list)
    for path in qfiles:
        for line in open(path, encoding="utf-8"):
            p = line.split()
            if not p or p[0] != "Q":
                continue
            j, ph = int(p[1]), int(p[2])
            sig = tuple(int(t) for t in p[3:7])
            o, exact, nsamp, L, NB = int(p[7]), int(p[8]), int(p[9]), int(p[10]), int(p[11])
            groups[(j, sig)].append((ph, o, exact, nsamp, L, NB))
    rows = []
    bad = 0
    for (j, sig), lst in sorted(groups.items()):
        if (j, sig) not in idx or idx[(j, sig)] not in det:
            print("W 가 정해지지 않은 변수:", j, sig)
            bad += 1
            continue
        W = det[idx[(j, sig)]]
        o, L, NB = lst[0][1], lst[0][4], lst[0][5]
        Eexp = Fraction(W) * 2 * NB / L                     # 서명이 s 인 phi 들의 E 의 합
        pexp = float(Eexp / (Fraction(2) ** o))             # sum_phi p_phi
        phat = sum(e / n for (_, _, e, n, _, _) in lst)
        var = sum((e / n) * (1 - e / n) / n for (_, _, e, n, _, _) in lst)
        sd = math.sqrt(var)
        z = (phat - pexp) / sd if sd > 0 else (0.0 if abs(phat - pexp) < 1e-12 else float("inf"))
        rows.append((j, sig, len(lst), lst[0][3], pexp, phat, sd, z, W))
    print("%5s %-20s %4s %10s %14s %14s %10s %7s %s" % ("j", "서명", "phi수", "표본수", "예측 sum p", "표본 sum p", "표준오차", "z", "W"))
    for (j, sig, k, n, pe, ph, sd, z, W) in rows:
        print("%5d %-20s %4d %10d %14.8f %14.8f %10.2e %7.2f %s" % (j, ",".join(map(str, sig)), k, n, pe, ph, sd, z, W))
    zs = [r[7] for r in rows if math.isfinite(r[7])]
    sz2 = sum(z * z for z in zs)
    print("\n견준 (j, 서명) %d개, sum z^2 = %.1f (자유도 %d), |z| 최대 %.2f, |z| > %.1f 인 것 %d개, 무한대 z %d개, 풀리지 않은 변수 %d개"
          % (len(rows), sz2, len(zs), max((abs(z) for z in zs), default=0.0), zmax, sum(1 for z in zs if abs(z) > zmax),
             sum(1 for r in rows if not math.isfinite(r[7])), bad))
    ok = bad == 0 and all(math.isfinite(r[7]) and abs(r[7]) <= zmax for r in rows)
    print("결과:", "일치 (모든 시험을 통과)" if ok else "불일치가 있다")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
