#!/usr/bin/env python3
"""
analyze_npn_stabilizer_types.py
===============================
npn_stabilizer_types.g 가 쓴 표(notes/npn_stabilizer_types_n*.tsv, npn_projected_types_n*.tsv, *_absent.tsv)와
stabilizer_types.g 의 표(notes/stabilizer_types_n*.tsv)를 읽어 노트에 인용하는 숫자를 모은다.
입력 파일의 열 순서는 각 파일의 머리말 주석을 따른다.

  npn_stabilizer_types_n<n>.tsv      (방법 A, n <= 5)  열 0 |K|, 1 지수, 2 함수 수, 3 NPN 류 수, 4 반대칭 유무, 5 |Kbar|, 6 |H|, ..., 17 H 의 유형
                                     (n = 6 의 표는 방법 C(make_npn_k_table.py)가 쓴 것이고 열 0..6 이 같다. 그 뒤의 열은 다르다)
  npn_projected_types_n<n>.tsv       (방법 B, n <= 6)  열 0 |Kbar|, 1 지수, 2 함수 수, 3 NPN 류 수, 4 H(f) 로 나타나는 B_n-궤도 수, 5 ng, 6 변수 궤도, 7 d, 8 꼭짓점 궤도, 9 고정 꼭짓점, 10 구조
  npn_projected_types_n<n>_absent.tsv                열 0 위수, 1 변수 궤도, 2 d, 3 꼭짓점 궤도, 4 고정 꼭짓점, 5 ng, 6 H(f) 로 나타나는 B_n-궤도 수

출력하는 것
  1. n 별 개수: B_n x Z_2 와 B_n 의 켤레류 수, K 로 나타나는 것, Kbar 로 나타나는 것, H 로 나타나는 것, NPN 류 수(OEIS A000370 과 견줌),
     반대칭이 있는 류 수 (2 NPN - NP 와 견줌), 함수 수의 합.
  2. Kbar 와 H 의 관계: Kbar 로만 나타나는 군(H(f) 로는 나타나지 않는 군), H(f) 로만 나타나는 군.
  3. (n <= 5) 반대칭 유무별 (K 의 켤레류, NPN 류, 함수) 수와 H 의 유형별.
  4. 안정화군이 자명한 함수, Kbar 가 자명한 함수.
  5. 함수가 가장 많은 자명하지 않은 군.

실행: python3 -I scripts/analyze_npn_stabilizer_types.py notes [N_MAX=6] [NA_MAX=5 또는 6(n = 6 의 K 표가 있을 때)]
"""
import csv
import math
import os
import re
import sys
from collections import Counter, defaultdict

NPN_KNOWN = {1: 2, 2: 4, 3: 14, 4: 222, 5: 616126, 6: 200253952527184}              # OEIS A000370
NP_KNOWN = {1: 3, 2: 6, 3: 22, 4: 402, 5: 1228158, 6: 400507806843728}              # B_n-궤도 수 (OEIS A000616)


def rows(path):
    return [r for r in csv.reader(open(path, encoding="utf-8"), delimiter="\t") if r and not r[0].startswith("#")]


def header_class_count(path, label):
    """머리말의 '<label> 의 켤레류 N개' 에서 N 을 읽는다."""
    for line in open(path, encoding="utf-8"):
        m = re.search(re.escape(label) + r" 의 켤레류 (\d+)개", line)
        if m:
            return int(m.group(1))
    return None


def main():
    root = sys.argv[1]
    nmax = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    namax = int(sys.argv[3]) if len(sys.argv) > 3 else 5
    A, B, Babs, R9 = {}, {}, {}, {}
    nclassA, nclassB = {}, {}
    for n in range(1, nmax + 1):
        pb = os.path.join(root, "npn_projected_types_n%d.tsv" % n)
        B[n] = rows(pb)
        Babs[n] = rows(os.path.join(root, "npn_projected_types_n%d_absent.tsv" % n))
        R9[n] = rows(os.path.join(root, "stabilizer_types_n%d.tsv" % n))
        nclassB[n] = header_class_count(pb, "B_n")
        if n <= namax:
            pa = os.path.join(root, "npn_stabilizer_types_n%d.tsv" % n)
            A[n] = rows(pa)
            nclassA[n] = header_class_count(pa, "B_n x Z_2")

    print("[1] n 별 개수")
    print("  n  B_n 켤레류  Kbar 로 나타남  H 로 나타남  Z_2 곱의 켤레류  K 로 나타남   NPN 류 수 (A000370 과)   반대칭 류 (2 NPN - NP 와)   함수 수 합")
    for n in range(1, nmax + 1):
        ncls = sum(int(r[3]) for r in B[n])
        nfun = sum(int(r[2]) for r in B[n])
        anti = "-"
        if n in A:
            anti = "%d (%s)" % (sum(int(r[3]) for r in A[n] if r[4] == "1"), "일치" if sum(int(r[3]) for r in A[n] if r[4] == "1") == 2 * ncls - NP_KNOWN[n] else "불일치")
        print("  %d %9d %12d %14d %16s %12s   %d (%s)   %s   %s" % (
            n, nclassB[n], len(B[n]), len(R9[n]), nclassA.get(n) or "-", len(A[n]) if n in A else "-",
            ncls, "일치" if ncls == NPN_KNOWN[n] else "불일치", anti,
            "2^(2^n) 와 같다" if nfun == 2 ** (2 ** n) else "다르다 %d" % nfun))
        if n in A:
            same = (sum(int(r[3]) for r in A[n]) == ncls and sum(int(r[2]) for r in A[n]) == nfun)
            if not same:
                print("     !!! 방법 A 와 B 의 합이 다르다")

    print("\n[1b] 나타나는 군의 비율과 반대칭 류의 비율")
    print("  n  Kbar 로 나타남/B_n 켤레류   H 로 나타남/B_n 켤레류   K 로 나타남/Z_2 곱의 켤레류   반대칭 류/NPN 류   반대칭이 없는 NPN 류")
    for n in range(1, nmax + 1):
        ncls = sum(int(r[3]) for r in B[n])
        anti = 2 * ncls - NP_KNOWN[n]
        print("  %d  %.4g%%   %.4g%%   %s   %.4g%%   %d" % (
            n, 100.0 * len(B[n]) / nclassB[n], 100.0 * len(R9[n]) / nclassB[n],
            ("%.4g%%" % (100.0 * len(A[n]) / nclassA[n])) if (n in A and nclassA.get(n)) else "-", 100.0 * anti / ncls, ncls - anti))

    print("\n[2] Kbar 와 H(f) 의 관계 (방법 B 표의 5열 = H(f) 로 나타나는 B_n-궤도 수)")
    print("  n  Kbar 로 나타나는 군  그 가운데 H(f) 로는 안 나타나는 것  H(f) 로는 나타나지만 Kbar 로는 안 나타나는 것  Kbar 로 나타나고 ng = 0 (반대칭 없음)인 것")
    for n in range(1, nmax + 1):
        only_k = sum(1 for r in B[n] if int(r[4]) == 0)
        only_h = sum(1 for r in Babs[n] if int(r[6]) > 0)
        ng0 = sum(1 for r in B[n] if int(r[5]) == 0)
        print("  %d %12d %24d %34d %36d" % (n, len(B[n]), only_k, only_h, ng0))
    # 고정 꼭짓점이 있는 Kbar 는 반대칭이 없다 (반대칭 원소는 꼭짓점을 고정하지 못한다)
    bad = sum(1 for n in B for r in B[n] if int(r[9]) > 0 and int(r[5]) > 0)
    print("  점검: 고정 꼭짓점이 있는데 ng > 0 인 Kbar 의 수 = %d (0 이어야 한다)" % bad)
    nonpow = sum(1 for n in B for r in list(B[n]) for _ in [0] if (int(r[5]) + 1) & int(r[5])) + \
        sum(1 for n in Babs for r in Babs[n] if (int(r[5]) + 1) & int(r[5]))
    print("  점검: ng + 1 이 2 의 거듭제곱이 아닌 켤레류의 수 = %d (0 이어야 한다)" % nonpow)
    print("  Kbar 로 나타나는 군의 k = log2(ng + 1) 분포")
    for n in range(1, nmax + 1):
        cnt = Counter(int(r[5]).bit_length() for r in B[n])
        print("  n=%d: %s" % (n, ", ".join("k=%d: %d" % (k, v) for k, v in sorted(cnt.items()))))

    if A:
        print("\n[3] (n <= %d) 반대칭 유무별 K 의 켤레류, NPN 류, 함수" % namax)
        print("  n  반대칭   K 의 켤레류   NPN 류 수          함수 수     함수 비율")
        for n in sorted(A):
            total = 2 ** (2 ** n)
            for flag, name in (("0", "없음"), ("1", "있음")):
                sub = [r for r in A[n] if r[4] == flag]
                nf = sum(int(r[2]) for r in sub)
                print("  %d  %-4s %9d %12d %20d      %.6g%%" % (n, name, len(sub), sum(int(r[3]) for r in sub), nf, 100.0 * nf / total))
        print("  H 의 유형(U 고정 꼭짓점 있음, T0 꼬임, T1 평행이동 있음) 과 반대칭 유무별 NPN 류 수")
        for n in sorted(A):
            c = defaultdict(int)
            for r in A[n]:
                if len(r) > 17 and r[17] in ("U", "T0", "T1"):
                    c[(r[17], r[4])] += int(r[3])
            if c:
                print("  n=%d: %s" % (n, ", ".join("%s/반대칭 %s: %d" % (t, a, v) for (t, a), v in sorted(c.items()))))
            else:
                print("  n=%d: (이 표에는 H 의 유형 열이 없다)" % n)

    print("\n[4] 안정화군이 자명한 함수")
    for n in range(1, nmax + 1):
        k1 = [r for r in B[n] if r[0] == "1"]
        h1 = [r for r in R9[n] if r[0] == "1"]
        print("  n=%d: Kbar 가 자명한 함수 %s개 (NPN 류 %s개); H(f) 가 자명한 함수 %s개 (B_n-궤도 %s개)" % (
            n, k1[0][2] if k1 else 0, k1[0][3] if k1 else 0, h1[0][2] if h1 else 0, h1[0][3] if h1 else 0))
        if k1:
            print("        Kbar 가 자명한 함수의 비율 %.6f%%" % (100.0 * int(k1[0][2]) / 2 ** (2 ** n)))
        if h1:
            # H(f) 가 자명하면 [Kbar : H] <= 2 이므로 Kbar 가 자명하지 않은 함수는 Kbar 의 위수가 2 인 것이다
            d = int(h1[0][2]) - (int(k1[0][2]) if k1 else 0)
            print("        H(f) 는 자명하고 Kbar 는 자명하지 않은 함수 %d개 (H(f) 가 자명한 함수의 %.4g%%; 이때 Kbar 의 위수는 2)" % (
                d, 100.0 * d / int(h1[0][2])))

    print("\n[5] 자명하지 않은 Kbar 가운데 함수가 가장 많은 것 (위수, 함수 수, NPN 류 수, ng, 구조)")
    for n in range(3, nmax + 1):
        sub = sorted((r for r in B[n] if r[0] != "1"), key=lambda r: -int(r[2]))[:5]
        print("  n=%d:" % n)
        for r in sub:
            print("     위수 %5s  함수 %26s  NPN 류 %16s  ng %s  %s" % (r[0], r[2], r[3], r[5], r[10]))

    if A:
        print("\n[6] 위수 |K| 별 K 의 종류 수 (n = 4 이상)")
        for n in sorted(A):
            if n < 4:
                continue
            by = defaultdict(lambda: [0, 0, 0, 0, 0])
            for r in A[n]:
                b = by[int(r[0])]
                b[0] += 1
                b[1] += int(r[3])
                b[2] += int(r[2])
                if r[4] == "1":
                    b[3] += 1
                    b[4] += int(r[3])
            print("  n=%d: 서로 다른 위수 %d개" % (n, len(by)))
            print("     위수   K 의 종류   그 가운데 반대칭   NPN 류 수   (반대칭 류 수)             함수 수")
            for o in sorted(by):
                b = by[o]
                print("   %6d %10d %12d %17d %12d %20d" % (o, b[0], b[3], b[1], b[4], b[2]))


if __name__ == "__main__":
    main()
