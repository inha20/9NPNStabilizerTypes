#!/usr/bin/env python3
"""
circuit_analyze.py
==================
circuit_all_functions.py 가 기록한 함수별 AND 수(도구의 결과)를 B_4-궤도(402개)와 NPN 동치류(222개)로 묶어 본다.

  - 궤도 안에서 AND 수가 일정한 궤도의 수, (최대 - 최소)의 분포
  - 함수마다 직접 합성한 값과 "궤도에서 가장 작은 값"의 차이, 그 합
  - 궤도별 최소 AND 수의 분포

"궤도에서 가장 작은 값"은 그 궤도의 어떤 함수에 대해 도구가 낸 가장 작은 회로를 입력 순열(배선)과
입력 반전(AIG 에서는 비용 0)으로 옮겨 쓰는 경우의 AND 수이다. NPN 쪽은 출력 반전(역시 비용 0)까지 허용한다.

실행: python scripts/circuit_analyze.py notes/circuit_and_counts.txt
"""
import itertools
import sys
from collections import Counter

N = 16
FULL = 1 << N

counts = [None] * FULL
for line in open(sys.argv[1]):
    if line.startswith("#") or not line.strip():
        continue
    chunk, digits = line.split()
    base = int(chunk) * 4096
    for i, ch in enumerate(digits):
        counts[base + i] = int(ch, 16)
assert all(c is not None for c in counts), "빠진 함수가 있다"


def var_perm_on_points(p):
    out = []
    for x in range(N):
        y = 0
        for i in range(4):
            y |= ((x >> p[i]) & 1) << i
        out.append(y)
    return tuple(out)


def act(f, g):
    h = 0
    for x in range(N):
        if (f >> x) & 1:
            h |= 1 << g[x]
    return h


gens = [tuple(x ^ m for x in range(N)) for m in (1, 2, 4, 8)] + [
    var_perm_on_points((1, 0, 2, 3)), var_perm_on_points((0, 2, 1, 3)), var_perm_on_points((0, 1, 3, 2))]


def classes(with_output_negation):
    label = [-1] * FULL
    groups = []
    for f in range(FULL):
        if label[f] != -1:
            continue
        k = len(groups)
        label[f] = k
        members = [f]
        stack = [f]
        while stack:
            u = stack.pop()
            images = [act(u, g) for g in gens]
            if with_output_negation:
                images.append(u ^ (FULL - 1))
            for v in images:
                if label[v] == -1:
                    label[v] = k
                    members.append(v)
                    stack.append(v)
        groups.append(members)
    return groups


def report(name, groups):
    total_direct = sum(counts)
    total_best = 0
    constant = 0
    spread = Counter()
    above = 0
    best_dist_groups = Counter()
    best_dist_funcs = Counter()
    worst = []
    for members in groups:
        vals = [counts[f] for f in members]
        lo, hi = min(vals), max(vals)
        total_best += lo * len(members)
        constant += lo == hi
        spread[hi - lo] += 1
        above += sum(1 for v in vals if v > lo)
        best_dist_groups[lo] += 1
        best_dist_funcs[lo] += len(members)
        worst.append((hi - lo, min(members), lo, hi, len(members)))
    worst.sort(reverse=True)
    print("[%s] %d개" % (name, len(groups)))
    print("  AND 수가 안에서 일정한 것: %d개.  (최대 - 최소)의 분포: %s" % (constant, dict(sorted(spread.items()))))
    print("  직접 합성한 AND 수의 합: %d (함수당 평균 %.3f)" % (total_direct, total_direct / FULL))
    print("  가장 작은 값으로 바꿨을 때의 합: %d (함수당 평균 %.3f), 줄어드는 비율 %.2f%%"
          % (total_best, total_best / FULL, 100.0 * (total_direct - total_best) / total_direct))
    print("  직접 합성한 값이 가장 작은 값보다 큰 함수: %d개 (%.2f%%)" % (above, 100.0 * above / FULL))
    print("  가장 작은 값의 분포 (묶음 수):", dict(sorted(best_dist_groups.items())))
    print("  가장 작은 값의 분포 (함수 수):", dict(sorted(best_dist_funcs.items())))
    print("  차이가 큰 묶음 (차이, 대표 진리표, 최소, 최대, 크기):",
          ["%d: %04x %d~%d (%d개)" % w for w in worst[:5]])


print("함수 65,536개, 도구가 낸 AND 수의 분포:", dict(sorted(Counter(counts).items())))
report("B_4-궤도 (입력 순열과 입력 반전)", classes(False))
report("NPN 동치류 (출력 반전까지)", classes(True))
