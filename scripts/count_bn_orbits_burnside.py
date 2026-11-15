#!/usr/bin/env python3
"""
count_bn_orbits_burnside.py
===========================
B_n = Z_2^n x| S_n 가 n변수 부울 함수 2^(2^n)개에 작용할 때의 궤도 수(변수 반전과 순열로 옮겨지는 함수 류의 수)를
Burnside 보조정리로 센다.

  궤도 수 = (1 / |B_n|) * sum_{g in B_n} 2^(g 가 꼭짓점 2^n 개에 만드는 순환의 수)

B_n 의 원소 (pi, c) 를 모두 훑는다(n = 7 이면 645,120개). 순열 pi 마다 평행이동 c 128개를 한 번에 다루고,
순환의 수는 "순환 안의 최솟값" 을 포인터 점프로 구해 최솟값이 자기 자신인 꼭짓점을 센다.
표지 표도 켤레류도 쓰지 않으므로 notes/2026-11-12_stabilizer_types.md 2.1절 표의 궤도 수 열
(표지 표 분해로 구한 값의 합)과 독립인 점검이다.

사용
  python3 -I scripts/count_bn_orbits_burnside.py [N_MAX]      (기본 N_MAX = 7, 약 2초)
출력: n 과 궤도 수, 그리고 n <= 6 에서는 표의 값과 같은지.

2026-11-14. Claude가 쓴 코드이고 저자 검토 전이다.
"""
import itertools
import math
import sys

import numpy as np

# 2.1절 표의 값(표지 표 분해로 구한 것). n = 7 은 표가 없다.
TABLE = {1: 3, 2: 6, 3: 22, 4: 402, 5: 1228158, 6: 400507806843728}


def bn_orbits_on_functions(n):
    nv = 1 << n
    xs = np.arange(nv, dtype=np.int64)
    cs = np.arange(nv, dtype=np.int64)
    rows = np.arange(nv)[:, None]
    total = 0
    # 순환 길이는 g 의 위수 이하이고 B_7 의 원소의 위수는 2^9 보다 훨씬 작다. 점프 9번이면 길이 512까지 덮는다.
    jumps = 9
    for perm in itertools.permutations(range(n)):
        p = np.zeros(nv, dtype=np.int64)          # 변수 순열: 비트 i 가 비트 perm[i] 로 간다
        for i in range(n):
            p |= ((xs >> i) & 1) << perm[i]
        g = p[None, :] ^ cs[:, None]              # (c, x) -> pi(x) xor c
        lab = np.minimum(xs[None, :], g)
        cur = g.copy()                            # cur = g^(2^t)
        for _ in range(jumps):
            lab = np.minimum(lab, lab[rows, cur])
            cur = cur[rows, cur]
        ncyc = (lab == xs[None, :]).sum(axis=1)   # c 마다 순환의 수
        total += sum(1 << int(k) for k in ncyc)
    order = nv * math.factorial(n)
    assert total % order == 0, "Burnside 합이 |B_n| 로 나누어떨어지지 않는다"
    return total // order


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    bad = False
    for n in range(1, nmax + 1):
        v = bn_orbits_on_functions(n)
        note = ""
        if n in TABLE:
            ok = (v == TABLE[n])
            bad |= not ok
            note = "  (표와 같다)" if ok else "  (표와 다르다!)"
        print("n = %d  B_n-궤도 수 = %d%s" % (n, v, note), flush=True)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
