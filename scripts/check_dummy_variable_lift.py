#!/usr/bin/env python3
"""
check_dummy_variable_lift.py
============================
닫힌 군 H <= B_n 가운데 "쓰이지 않는 변수"가 있는 것을 세고, 변수를 하나 덜 쓰는 군과 짝지어지는지 확인한다.

바탕이 되는 사실 (증명은 안정화군 노트 3절 "일반 n")
  H 가 닫힌 군이고 꼭짓점 이름표 lambda_H (꼭짓점 -> H-궤도 번호)가 변수 i 에 의존하지 않으면, 변수 i 의 반전이 H 에 들어 있다.
  쓰이지 않는 변수의 집합을 I, 쓰이는 변수의 집합을 J 라 하면 H = B_I x H_J 이고 H_J 는 B_J 의 닫힌 군이며 J 의 변수가 모두 쓰인다.
  B_I 가 {0,1}^I 위에서 추이적이므로 H 의 꼭짓점 궤도는 (V_I) x (H_J 의 궤도)이고 궤도 수 k 는 같고 궤도 길이는 2^|I| 배이다.
  함수 f 의 안정화군이 B_I x H_J 이면 f 는 I 의 변수에 의존하지 않으므로 f 는 J 의 함수 f_J 이고 Stab(f) = B_I x Stab(f_J) 이다.
  그러므로 B_I x H_J 가 나타남 <=> H_J 가 B_J 에서 나타남. 닫힘도 같다.
  따라서 n 에서 "쓰이지 않는 변수가 정확히 r 개인 닫힌 (나타나는 / 안 나타나는) 류"의 수는 n-r 에서 "모든 변수가 쓰이는 닫힌 (나타나는 / 안 나타나는) 류"의 수와 같다.
  특히 안 나타나는 류의 수 absent(n) 은 n 에 대해 줄지 않고 absent(n) - absent(n-1) 은 모든 변수가 쓰이는 안 나타나는 류의 수이다.

하는 일
  표(공개 표와 같은 모양)의 닫힌 류마다 꼭짓점 궤도를 구하고 쓰이지 않는 변수의 수 r 을 센다.
  (n, r) 별 개수를 (n-r, 0) 의 개수와 견주고, (k, 궤도 길이 / 2^r) 의 모임도 견준다.

사용
  python3 -I scripts/check_dummy_variable_lift.py --pipe 작업폴더 --n7 notes/stabilizer_types_n7_classes.tsv.gz
  작업폴더: run_bn_class_pipeline.sh 의 작업 폴더 (classes_n2.tsv .. classes_n6.tsv 가 있어야 한다). n = 7 의 표는 공개 표.
  시간은 1~2초 (코어 하나).

2026-11-15. Claude가 쓴 코드이고 저자 검토 전이다.
"""
import argparse
import os
import sys
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bn_class_io import read_rows  # noqa: E402


def vertex_perms(n, vals, ng):
    """생성원 ng 개(각 n+1 개의 수)를 꼭짓점 순열(길이 2^n 배열)로 바꾼다."""
    V = 1 << n
    X = np.arange(V, dtype=np.int64)
    out = []
    for t in range(ng):
        im = vals[t * (n + 1):(t + 1) * (n + 1)]
        c = im[0]
        vp = np.full(V, c, dtype=np.int64)
        for i in range(n):
            vp ^= ((X >> i) & 1) * (im[1 + i] ^ c)
        out.append(vp)
    return out


def orbit_labels(n, gens):
    """꼭짓점 궤도마다 가장 작은 꼭짓점 번호를 이름표로 붙인 배열."""
    V = 1 << n
    lab = np.arange(V, dtype=np.int64)
    invs = [np.argsort(g) for g in gens]
    while True:
        old = lab
        for g, gi in zip(gens, invs):
            lab = np.minimum(lab, lab[g])
            lab = np.minimum(lab, lab[gi])
        if np.array_equal(old, lab):
            return lab


def analyse(n, vals, ng):
    """(쓰이지 않는 변수의 수 r, 궤도 수 k, 궤도 길이를 2^r 로 나눠 정렬한 튜플)."""
    V = 1 << n
    X = np.arange(V, dtype=np.int64)
    gens = vertex_perms(n, vals, ng)
    lab = orbit_labels(n, gens)
    r = sum(1 for i in range(n) if np.array_equal(lab[X ^ (1 << i)], lab))
    sizes = sorted(np.bincount(lab)[np.unique(lab)].tolist())
    return r, len(sizes), tuple(s >> r for s in sizes)


def collect(path, n):
    """닫힌 류(나타남, 안 나타남)마다 analyse 의 결과를 모은다. 반환: {상태: [(r, k, 길이들, 위수)]}"""
    res = {"present": [], "absent": []}
    for j, order, k, status, ng, vals, _ in read_rows(path):
        if status not in res:
            continue
        r, kk, lens = analyse(n, vals, ng)
        assert kk == k, (n, j, kk, k)
        res[status].append((r, k, lens, order))
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pipe", required=True)
    ap.add_argument("--n7", required=True)
    a = ap.parse_args()

    data = {}
    for n in range(2, 7):
        data[n] = collect(os.path.join(a.pipe, "classes_n%d.tsv" % n), n)
    data[7] = collect(a.n7, 7)

    ok = True
    for status in ("absent", "present"):
        print("=== %s: n 별, 쓰이지 않는 변수의 수 r 별 류의 수 (r=0 이 모든 변수가 쓰이는 류) ===" % status)
        ess = {}
        for n in range(2, 8):
            cnt = Counter(r for r, _, _, _ in data[n][status])
            ess[n] = cnt[0]
            print("n=%d  합 %4d   " % (n, sum(cnt.values())) + "  ".join("r=%d:%d" % (r, cnt[r]) for r in sorted(cnt)))
        print("모든 변수가 쓰이는 류의 수 (n=2..7): %s" % [ess[n] for n in range(2, 8)])
        tot = {n: len(data[n][status]) for n in range(2, 8)}
        print("%s(n) (n=2..7): %s" % (status, [tot[n] for n in range(2, 8)]))
        diffs = [tot[n] - tot[n - 1] for n in range(3, 8)]
        print("%s(n) - %s(n-1) (n=3..7): %s   (모든 변수가 쓰이는 류의 수 %s 와 같아야 한다)"
              % (status, status, diffs, [ess[n] for n in range(3, 8)]))
        if diffs != [ess[n] for n in range(3, 8)]:
            ok = False
            print("!! 같지 않다")
        # (n, r) 와 (n-r, 0) 견주기
        for n in range(3, 8):
            for r in range(1, n - 1):
                A = Counter((k, lens) for rr, k, lens, _ in data[n][status] if rr == r)
                B = Counter((k, lens) for rr, k, lens, _ in data[n - r][status] if rr == 0)
                if A != B:
                    ok = False
                    print("!! (n=%d, r=%d) 와 (n=%d, r=0) 의 (k, 궤도 길이/2^r) 모임이 다르다: %d 대 %d" % (n, r, n - r, sum(A.values()), sum(B.values())))
        print()
    # 위수 견주기: |H| = |B_r| * |H_J|, |B_r| = 2^r r!
    from math import factorial
    bad = 0
    for status in ("absent", "present"):
        for n in range(3, 8):
            for r in range(1, n - 1):
                A = Counter((k, lens, order // (2 ** r * factorial(r))) for rr, k, lens, order in data[n][status] if rr == r)
                B = Counter((k, lens, order) for rr, k, lens, order in data[n - r][status] if rr == 0)
                if A != B:
                    bad += 1
    print("위수 |H| = 2^r r! |H_J| 까지 넣은 모임 견주기: 어긋난 (상태, n, r) %d개" % bad)
    if bad:
        ok = False
    # 나타남 + 안 나타남 (닫힌 류)
    print("\n=== 닫힌 류 (나타남 + 안 나타남) 의 수, n=2..7 ===")
    closed = [len(data[n]["present"]) + len(data[n]["absent"]) for n in range(2, 8)]
    print(closed)
    print("결과:", "모두 일치" if ok else "불일치가 있다")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
