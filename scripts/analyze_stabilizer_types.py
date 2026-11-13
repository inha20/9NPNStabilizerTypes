#!/usr/bin/env python3
"""
analyze_stabilizer_types.py
===========================
stabilizer_types.g 가 쓴 표(notes/stabilizer_types_n*.tsv, *_absent.tsv)를 읽어 노트에 인용하는 숫자를 모은다.
입력 파일의 열 순서는 그 파일의 머리말 주석을 따른다.

출력하는 것
  1. n 별 개수: B_n 의 켤레류 수, 안정화군으로 나타나는 것, 닫힌 것, 닫혔지만 나타나지 않는 것.
     함수 수와 B_n-궤도 수의 합이 2^(2^n) 과 알려진 궤도 수에 맞는지.
  2. 유형별(U 고정 꼭짓점 있음, T0 꼬임, T1 평행이동 있음) 켤레류·궤도·함수 수.
  3. 평행이동 부분의 차원 d 별 분포.
  4. 안정화군이 평행이동뿐인 경우(변수 순열 부분이 자명), 안정화군이 자명한 함수.
  5. 고정 꼭짓점이 있는 유형(U)과 S_n 의 부분군 켤레류 수의 비교, 닫혔지만 나타나지 않는 U 유형.
  6. 안정화군이 자명하지 않은 켤레류 가운데 함수가 가장 많은 것.
  7. 안정화군의 변수 순열 부분 Hbar 로 나타나는 S_n 의 켤레류 수 (PLAN.md R5 의 이전 계산과 견주는 값).
  8. 닫혔지만 나타나지 않는 켤레류의 꼭짓점 궤도 수.

실행: python scripts/analyze_stabilizer_types.py notes [N_MAX=6]
"""
import csv
import os
import sys
from collections import Counter, defaultdict

KNOWN_ORBITS = {1: 3, 2: 6, 3: 22, 4: 402, 5: 1228158, 6: 400507806843728}
SN_SUBGROUP_CLASSES = {1: 1, 2: 2, 3: 4, 4: 11, 5: 19, 6: 56}      # S_n 의 부분군 켤레류 수 (GAP 로 센 값과 같다)


def rows(path):
    return [r for r in csv.reader(open(path, encoding="utf-8"), delimiter="\t") if r and not r[0].startswith("#")]


def main():
    root = sys.argv[1]
    nmax = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    occ, ab = {}, {}
    for n in range(1, nmax + 1):
        occ[n] = rows(os.path.join(root, "stabilizer_types_n%d.tsv" % n))
        ab[n] = rows(os.path.join(root, "stabilizer_types_n%d_absent.tsv" % n))

    print("[1] n 별 개수")
    print("  n  |B_n|   켤레류  나타남  닫힘  닫혔지만 안 나타남  닫히지 않아 안 나타남  함수 수 합 == 2^(2^n)  궤도 수 합 == 알려진 값")
    for n in range(1, nmax + 1):
        o, a = occ[n], ab[n]
        nf = sum(int(r[2]) for r in o)
        no = sum(int(r[3]) for r in o)
        ca = sum(1 for r in a if r[6] == "1")
        print("  %d %6d %7d %7d %5d %12d %20d      %s            %s (%d)" % (
            n, 2 ** n * __import__("math").factorial(n), len(o) + len(a), len(o), len(o) + ca, ca, len(a) - ca,
            "예" if nf == 2 ** (2 ** n) else "아니오", "예" if no == KNOWN_ORBITS[n] else "아니오", no))

    print("\n[2] 유형별 (U: 고정 꼭짓점 있음, T0: 고정 꼭짓점도 평행이동도 없음(꼬임), T1: 0 이 아닌 평행이동 있음)")
    print("  n  유형   켤레류        B_n-궤도 수          함수 수      함수 비율")
    for n in range(1, nmax + 1):
        total = 2 ** (2 ** n)
        for t in ("U", "T0", "T1"):
            sub = [r for r in occ[n] if r[10] == t]
            nf = sum(int(r[2]) for r in sub)
            print("  %d  %-3s %8d %17d %24d      %.6g%%" % (n, t, len(sub), sum(int(r[3]) for r in sub), nf, 100.0 * nf / total))

    print("\n[3] 평행이동 부분의 차원 d (안정화군 가운데 입력 반전만으로 된 원소들이 이루는 공간의 차원) 별 (켤레류, 궤도, 함수)")
    for n in range(1, nmax + 1):
        byd = defaultdict(lambda: [0, 0, 0])
        for r in occ[n]:
            v = byd[int(r[7])]
            v[0] += 1
            v[1] += int(r[3])
            v[2] += int(r[2])
        print("  n=%d: %s" % (n, "; ".join("d=%d: %d, %d, %d" % (d, *v) for d, v in sorted(byd.items()))))

    print("\n[4] 안정화군이 평행이동뿐인 경우와 자명한 경우")
    for n in range(1, nmax + 1):
        pt = [r for r in occ[n] if r[4] == "1" and int(r[7]) > 0]
        tr = [r for r in occ[n] if r[0] == "1"]
        print("  n=%d: 평행이동뿐(|Hbar|=1, d>0): 켤레류 %d, 궤도 %d, 함수 %d (d 별: %s); 자명한 안정화군: 함수 %s, 궤도 %s" % (
            n, len(pt), sum(int(r[3]) for r in pt), sum(int(r[2]) for r in pt),
            dict(sorted(Counter(int(r[7]) for r in pt).items())),
            tr[0][2] if tr else 0, tr[0][3] if tr else 0))
    for n in range(1, nmax + 1):
        if occ[n] and any(r[0] == "1" for r in occ[n]):
            tr = [r for r in occ[n] if r[0] == "1"][0]
            print("  n=%d: 자명한 안정화군을 가진 함수의 비율 %.6f%%" % (n, 100.0 * int(tr[2]) / 2 ** (2 ** n)))

    print("\n[5] 고정 꼭짓점이 있는 유형(U)과 S_n 의 부분군 켤레류")
    print("  (U 유형의 B_n-켤레류는 S_n 의 부분군 켤레류와 하나씩 대응한다.)")
    for n in range(1, nmax + 1):
        ou = [r for r in occ[n] if r[10] == "U"]
        au = [r for r in ab[n] if int(r[5]) > 0]
        cu = [r for r in au if r[6] == "1"]
        print("  n=%d: U 유형 나타남 %d + 안 나타남 %d = %d (S_n 의 부분군 켤레류 %d); 안 나타나는 U 유형 가운데 닫힌 것 %d%s" % (
            n, len(ou), len(au), len(ou) + len(au), SN_SUBGROUP_CLASSES[n], len(cu),
            "".join(" [위수 %s, 변수 위 궤도 %s, 꼭짓점 궤도 %s]" % (r[0], r[2], r[4]) for r in cu)))

    print("\n[6] 안정화군이 자명하지 않은 켤레류 가운데 함수가 가장 많은 것 (위수, 함수 수, 궤도 수, 유형, 구조)")
    for n in range(3, nmax + 1):
        sub = sorted((r for r in occ[n] if r[0] != "1"), key=lambda r: -int(r[2]))[:6]
        print("  n=%d:" % n)
        for r in sub:
            print("     위수 %5s  함수 %26s  궤도 %16s  %-2s  %s" % (r[0], r[2], r[3], r[10], r[12]))

    print("\n[7] Hbar(f) 로 나타나는 S_n 의 켤레류 수 (안정화군의 변수 순열 부분의 켤레류)")
    for n in range(1, nmax + 1):
        labels = set(r[6] for r in occ[n])
        print("  n=%d: %d종 (S_n 의 부분군 켤레류 %d종 가운데 %d종은 Hbar 로 나타나지 않는다)" % (
            n, len(labels), SN_SUBGROUP_CLASSES[n], SN_SUBGROUP_CLASSES[n] - len(labels)))

    print("\n[8] 닫혔지만 나타나지 않는 켤레류의 꼭짓점 궤도 수, 그보다 궤도가 많은 닫힌 켤레류가 모두 나타나는지")
    for n in range(2, nmax + 1):
        ca = [r for r in ab[n] if r[6] == "1"]
        mx = max(int(r[4]) for r in ca)
        big_occ = sum(1 for r in occ[n] if int(r[8]) > mx)
        print("  n=%d: 닫혔지만 안 나타나는 %d개의 꼭짓점 궤도 수 분포 %s; 최대 %d. 꼭짓점 궤도가 %d개보다 많은 닫힌 켤레류 %d개는 모두 나타난다." % (
            n, len(ca), dict(sorted(Counter(int(r[4]) for r in ca).items())), mx, mx, big_occ))


if __name__ == "__main__":
    main()
