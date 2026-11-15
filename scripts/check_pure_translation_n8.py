#!/usr/bin/env python3
"""
check_pure_translation_n8.py
============================
순수 평행이동군 (H = T = 평행이동만으로 된 닫힌 군) 의 나타남 분류를 n=8 의 실제 B_8 군에서 다시 확인한다.

H = T_S = { x in F_2^8 : sum x_i s_i = 0 } (S = (s_0..s_7) 는 F_2^m 의 서로 다른 0 아닌 벡터, F_2^m 을 생성).
n=8 에서 가능한 m 은 4 (m <= 4 이면 |S| <= 2^m-1 이므로 m >= 4) 부터 8 까지이고, m >= 5 는 5.3절로 모두 나타난다.
m=4 인 S 의 GL(4,2)-류마다
  (a) scripts/pure_translation_absent.py 의 A 부분 (군 연산으로 센 'N(H)/H 에 자유로운 함수의 수')
  (b) np_cegar.appears (B_8 안의 SAT 증거-유도 판정; 나타남은 안정화군의 정확한 계산, 안 나타남은 모순 증거)
의 판정을 비교한다. (a) 에서 자유 함수가 0 개  <=>  (b) 에서 absent 여야 한다.

사용: python3 -u -I scripts/check_pure_translation_n8.py   (약 10초)

2026-11-15. Claude가 쓴 코드이고 저자 검토 전이다.
"""
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

import pure_translation_absent as PA
from np_cegar import Cube, appears
from check_hamming_transvection_family import kernel_basis, vp_translation


def main():
    t0 = time.time()
    res = PA.part_a(maxm=4)
    cases = [(m, n, A, free, S) for (m, n, A, free, S) in res if n == 8 and m == 4]
    print("n=8, m=4 인 S 의 GL(4,2)-류: %d 개" % len(cases), flush=True)
    cube = Cube(8)
    n_absent_pred = 0
    n_absent_obs = 0
    mism = 0
    for (m, n, A, free, S) in cases:
        T = kernel_basis(S, n)
        assert len(T) == n - m
        gens = [vp_translation(n, t).astype(cube.dt) for t in T]
        r = appears(cube, gens, order=1 << (n - m))
        pred_absent = (free == 0)
        obs = r['status']
        n_absent_pred += pred_absent
        n_absent_obs += (obs == 'absent')
        ok = (obs == 'absent') == pred_absent and obs in ('absent', 'present')
        mism += (not ok)
        print("  S=%s |A_S|=%d 자유함수 %d → 예측 %s, appears: %s (k=%s)  %s" % (
            S, A, free, 'absent' if pred_absent else 'present', obs, r.get('k'), '' if ok else '  <<< 불일치'), flush=True)
    print("예측 absent %d 개, np_cegar absent %d 개, 불일치 %d 개 (%.1f초)" % (n_absent_pred, n_absent_obs, mism, time.time() - t0))
    print("결과:", "모두 일치" if mism == 0 else "불일치 있음")


if __name__ == "__main__":
    main()
