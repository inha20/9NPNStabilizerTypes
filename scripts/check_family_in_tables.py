#!/usr/bin/env python3
"""
check_family_in_tables.py
=========================
노트 3절 "일반 n (3)"의 전단으로 만든 두 족이 n <= 7 의 표에서 어느 류인지 찾고, 변수를 모두 쓰는 안 나타나는 군 가운데 두 족에 드는 것이 몇 개인지 센다.

두 족 (평행이동 부분 T 의 열 h_1..h_n 은 F_2^m 의 점, K 는 F_2^m 에 하는 작용의 상)
  중심 고정 족: n = 2^m - 1, 열이 F_2^m 의 0 이 아닌 점 전부 (해밍 부호), K = 중심 u0 이 하나로 정해진 전단 전부 (위수 2^(m-1)).
  축 고정 족:   n = 2^(m-1), 열이 F_2^m 의 한 아핀 초평면 (0 을 지나지 않는 것) 의 점 전부 (확장 해밍 부호), K = 축이 하나로 정해진 전단 전부 (위수 2^(m-1)).
둘 다 k = 2^(m-1) + 1 이고 |H| = 2^(n-1) 이다.

하는 일 (scripts/analyze_absent_reduction.py 의 analyse_class 를 쓴다)
  n = 2..7 의 안 나타나는 류 가운데 모든 변수가 쓰이고 (열이 모두 0 이 아님) K != 1 인 것마다
    - 열이 위의 모양인지 (n = 2^m - 1 이고 열이 0 이 아닌 점 전부 / n = 2^(m-1) 이고 열이 한 아핀 초평면 전부),
    - K 의 원소가 모두 선형이고 항등원이 아닌 것은 전단 (g - 1 의 상이 직선) 인지, 그 전단들의 중심이 모두 같은지 (중심 고정) / 고정 초평면이 모두 같은지 (축 고정),
    - |K| = 2^(m-1) 과 k = 2^(m-1) + 1
  을 보고 두 족에 드는 류를 낸다. 그 수와 K != 1 인 류의 수를 견준다.

사용: python3 -I scripts/check_family_in_tables.py --pipe 작업폴더 --n7 notes/stabilizer_types_n7_closed_absent.tsv
      (작업폴더 = scripts/run_bn_class_pipeline.sh 의 작업 폴더). 시간은 1초 안팎.

2026-11-15. Claude가 쓴 코드이고 저자 검토 전이다.
"""
import argparse
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze_absent_reduction import analyse_class  # noqa: E402
from bn_class_io import read_rows  # noqa: E402
from np_cegar import Cube  # noqa: E402


def group_elements(perms, npts):
    ident = tuple(range(npts))
    seen = {ident}
    frontier = [ident]
    gens = [tuple(int(a) for a in g) for g in perms]
    while frontier:
        new = []
        for a in frontier:
            for g in gens:
                b = tuple(g[x] for x in a)
                if b not in seen:
                    seen.add(b)
                    new.append(b)
        frontier = new
    return seen


def is_linear(g):
    N = len(g)
    return g[0] == 0 and all(g[x ^ y] == g[x] ^ g[y] for x in range(N) for y in range(x, N))


def transvection_data(g):
    """g 가 항등원이 아닌 선형 변환일 때 g 가 전단이면 (중심 u, 고정 점들의 집합), 아니면 None."""
    diffs = {g[x] ^ x for x in range(len(g))}
    diffs.discard(0)
    if len(diffs) != 1:
        return None
    u = next(iter(diffs))
    fixed = frozenset(x for x in range(len(g)) if g[x] == x)
    return u, fixed


def column_type(cols, m):
    """열의 모양: 'hamming' (0 이 아닌 점 전부, n = 2^m - 1), 'ext' (한 아핀 초평면 전부, n = 2^(m-1)), 아니면 None."""
    n = len(cols)
    if n == (1 << m) - 1 and sorted(cols) == list(range(1, 1 << m)):
        return "hamming"
    if n == 1 << (m - 1) and len(set(cols)) == n:
        for phi in range(1, 1 << m):
            if all(bin(phi & h).count("1") % 2 == 1 for h in cols):
                return "ext"
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pipe", required=True)
    ap.add_argument("--n7", required=True)
    a = ap.parse_args()
    total = Counter()
    found = []
    for n in range(2, 8):
        path = os.path.join(a.pipe, "classes_n%d.tsv" % n) if n <= 6 else a.n7
        cube = Cube(n)
        for (j, order, k, status, ng, vals, line) in read_rows(path):
            if status != "absent":
                continue
            r = analyse_class(cube, n, order, k, ng, vals)
            if not r["essential"]:
                continue
            if r["K"] == 1:
                total[(n, "K=1")] += 1
                continue
            total[(n, "K>1")] += 1
            m = r["m"]
            ct = column_type(r["cols"], m)
            if ct is None or r["K"] != 1 << (m - 1) or k != (1 << (m - 1)) + 1 or r["pi_order"] != 1:
                continue
            els = group_elements(r["kperms"], 1 << m)
            if not all(is_linear(g) for g in els):
                continue
            ident = tuple(range(1 << m))
            tv = [transvection_data(g) for g in els if g != ident]
            if any(t is None for t in tv):
                continue
            common_centre = len({u for u, f in tv}) == 1
            common_axis = len({f for u, f in tv}) == 1
            fam = None
            if ct == "hamming" and common_centre:
                fam = "중심 고정"
            if ct == "ext" and common_axis:
                fam = "축 고정"
            if fam:
                found.append((n, j, m, fam, order, k))
    print("n 별 (모든 변수가 쓰이는 안 나타나는 류 가운데 K = 1 인 것 / K != 1 인 것)")
    for n in range(2, 8):
        print("  n=%d: K=1 %d개, K!=1 %d개" % (n, total[(n, "K=1")], total[(n, "K>1")]))
    nk1 = sum(v for (n, t), v in total.items() if t == "K=1")
    nk2 = sum(v for (n, t), v in total.items() if t == "K>1")
    print("  합계: K=1 %d개, K!=1 %d개 (모두 %d개)" % (nk1, nk2, nk1 + nk2))
    print()
    print("두 족에 드는 류 (n, 번호, m, 족, |H|, k):")
    for f in found:
        print("  n=%d 번호 %d m=%d %s 족 |H|=%d k=%d" % f)
    print("두 족에 드는 류 %d개 (K != 1 인 %d개 가운데)" % (len(found), nk2))
    print("결과: 끝")


if __name__ == "__main__":
    main()
