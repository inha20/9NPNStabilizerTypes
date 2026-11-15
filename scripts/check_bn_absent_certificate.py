#!/usr/bin/env python3
"""
check_bn_absent_certificate.py
==============================
안 나타나는 군 H 하나를 SAT 와 무관하게 전수로 확인한다. H 로 불변인 함수는 H 의 꼭짓점 궤도마다 값 하나를 고른 것이므로 2^k 개이다.
그 함수 하나하나의 대칭군의 위수를 정확히 세어(np_cegar.Cube.count_stab: 평행이동한 함수와 순열한 함수를 사전에서 맞춰 보는 방식)
위수가 |H| 인 함수가 하나도 없음을 보인다. 없으면 H 는 어떤 함수의 대칭군도 아니다.
np_cegar 의 안 나타남 판정(절들의 모순)과는 다른 방법이므로 독립 점검이다. k = 16 이면 65,536개 함수에 약 1분이 걸린다(2026-11-14 측정 66초).

사용
  python3 -I scripts/check_bn_absent_certificate.py N 표 [--class J | --max-k]
  표: 공개 표(.tsv 나 .tsv.gz). --class J 는 번호가 J 인 류, --max-k(기본)는 안 나타나는 류 가운데 꼭짓점 궤도가 가장 많은 류
  (같으면 번호가 작은 것). 화면에 H 의 구조(평행이동 부분, 궤도 길이)와 위수 분포를 쓴다.

2026-11-14. Claude가 쓴 코드이고 저자 검토 전이다. 2026-11-13의 임시 스크립트 cert16.py 와 explore200.py 를 합쳐 정리한 것이다.
"""
import argparse
import os
import sys
import time
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import np_cegar as C                                   # noqa: E402
from bn_class_io import read_classes                    # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("n", type=int)
    ap.add_argument("table")
    ap.add_argument("--class", dest="cls", type=int, default=None)
    ap.add_argument("--max-k", action="store_true")
    a = ap.parse_args()
    cube = C.Cube(a.n)
    X = cube.X
    cands = list(read_classes(a.table, cube, statuses=("absent",)))
    if a.cls is not None:
        cands = [c for c in cands if c[0] == a.cls]
        if not cands:
            sys.exit("번호 %d 인 안 나타나는 류가 표에 없다" % a.cls)
        pick = cands[0]
    else:
        pick = sorted(cands, key=lambda c: (-c[2], c[0]))[0]
    j, order, k, _, gens, _ = pick
    H = cube.closure(gens)
    assert len(H) == order
    lab, kk = cube.orbit_labels(gens)
    assert kk == k
    trs = sorted(int(np.frombuffer(b, dtype=cube.dt)[0]) for b in H
                 if (np.frombuffer(b, dtype=cube.dt).astype(np.int64) == (X ^ int(np.frombuffer(b, dtype=cube.dt)[0]))).all())
    print("n=%d 번호 %d: |H|=%d, 꼭짓점 궤도 k=%d, 궤도 길이 %s" % (a.n, j, order, k, sorted(Counter(np.bincount(lab).tolist()).items())))
    print("평행이동 부분: %d개, 벡터의 무게 %s" % (len(trs), sorted(bin(c).count("1") for c in trs)))
    t0 = time.time()
    hist = Counter()
    for cb in range(1 << k):
        col = np.array([(cb >> i) & 1 for i in range(k)], dtype=np.uint8)
        hist[cube.count_stab(col[lab])] += 1
        if cb and cb % 16384 == 0:
            print("  %d / %d (%.0f초)" % (cb, 1 << k, time.time() - t0), flush=True)
    bad = hist.get(order, 0)
    print("H 로 불변인 함수 2^%d = %d 개 가운데 대칭군의 위수가 |H| = %d 인 것: %d개, 대칭군의 최소 위수 %d (%.0f초)"
          % (k, 1 << k, order, bad, min(hist), time.time() - t0))
    print("대칭군 위수의 분포(작은 쪽 6개, 위수:함수 수): %s" % ", ".join("%d:%d" % kv for kv in sorted(hist.items())[:6]))
    if bad:
        sys.exit("위수가 |H| 인 함수가 있다: H 는 나타난다 (np_cegar 의 absent 판정과 어긋남)")
    print("결론: 위수가 |H| 인 함수가 없으므로 H 는 어떤 함수의 대칭군도 아니다.")


if __name__ == "__main__":
    main()
