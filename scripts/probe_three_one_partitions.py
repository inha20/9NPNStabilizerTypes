#!/usr/bin/env python3
"""
probe_three_one_partitions.py
=============================
4변수 함수 65,536개 전부에서, 변수 하나를 따로 두고 나머지 셋을 묶는 분해(3|1 분할)와
그림 종류 수 t 의 관계를 센다. probe_multiplicity_vs_pictures.py (2|2 분할)의 짝이다.

  - 변수 x_i 를 따로 두면 분해 차트는 2행(x_i = 0, 1) x 8열(나머지 세 변수의 값)이다.
    열 다중도 mu_i = 서로 다른 열의 수 (1..4). mu_i <= 2 이면 f = g(h(나머지 셋), x_i) 꼴로 단순 분해된다.
    (묶음이 변수 하나인 쪽은 다중도가 언제나 2 이하라서 조건이 되지 않는다.)
  - s3 = mu_i <= 2 인 변수 i 의 수,  s2 = 2|2 분할 가운데 단순 분해되는 것의 수 (앞 탐침의 s).
  - t = 그림 종류 수 (세 배치의 그림이 격자 대칭 W 아래 몇 가지인가).

표준 라이브러리만 쓴다.
"""
from collections import Counter

N = 16
FULL = 1 << N


def var_perm_on_points(p):
    out = []
    for x in range(N):
        bits = [(x >> (3 - i)) & 1 for i in range(4)]
        y = 0
        for i in range(4):
            y = (y << 1) | bits[p[i]]
        out.append(y)
    return tuple(out)


def act(f, g):
    h = 0
    for x in range(N):
        if (f >> x) & 1:
            h |= 1 << g[x]
    return h


def orbit_labels(gens):
    label = [-1] * FULL
    k = 0
    for f in range(FULL):
        if label[f] != -1:
            continue
        label[f] = k
        stack = [f]
        while stack:
            u = stack.pop()
            for g in gens:
                v = act(u, g)
                if label[v] == -1:
                    label[v] = k
                    stack.append(v)
        k += 1
    return label, k


neg = [tuple(x ^ m for x in range(N)) for m in (8, 4, 2, 1)]
swapAB = var_perm_on_points((1, 0, 2, 3))
swapCD = var_perm_on_points((0, 1, 3, 2))
swapPairs = var_perm_on_points((2, 3, 0, 1))
swapBC = var_perm_on_points((0, 2, 1, 3))
swapBD = var_perm_on_points((0, 3, 2, 1))
wlab, nw = orbit_labels(neg + [swapAB, swapCD, swapPairs])
assert nw == 805

PARTS = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]


def point(bits):
    return (bits[0] << 3) | (bits[1] << 2) | (bits[2] << 1) | bits[3]


def s2_of(f):
    s = 0
    for P, Q in PARTS:
        rows = []
        for pv in range(4):
            row = []
            for qv in range(4):
                b = [0, 0, 0, 0]
                b[P[0]], b[P[1]] = pv >> 1, pv & 1
                b[Q[0]], b[Q[1]] = qv >> 1, qv & 1
                row.append((f >> point(b)) & 1)
            rows.append(tuple(row))
        cols = set(zip(*rows))
        if min(len(set(rows)), len(cols)) <= 2:
            s += 1
    return s


def mu_three_one(f, i):
    others = [j for j in range(4) if j != i]
    cols = set()
    for v in range(8):
        col = []
        for xi in (0, 1):
            b = [0, 0, 0, 0]
            b[i] = xi
            for k, j in enumerate(others):
                b[j] = (v >> (2 - k)) & 1
            col.append((f >> point(b)) & 1)
        cols.add(tuple(col))
    return len(cols)


joint = Counter()
joint_any = Counter()
violations = 0
mu_dist = Counter()
for f in range(FULL):
    t = len({wlab[f], wlab[act(f, swapBC)], wlab[act(f, swapBD)]})
    mus = [mu_three_one(f, i) for i in range(4)]
    mu_dist[mus[0]] += 1
    s3 = sum(1 for m in mus if m <= 2)
    s2 = s2_of(f)
    joint[(t, s3)] += 1
    joint_any[(t, (s2 > 0) or (s3 > 0))] += 1
    if s3 > 0 and t == 3:
        violations += 1

print("[1] 열 다중도 mu 의 분포 (변수 A 를 따로 둘 때, 함수 수):", dict(sorted(mu_dist.items())))
print("[2] t(그림 종류 수) x s3(3|1 로 단순 분해되는 변수의 수): 함수 수")
print("          s3=0     s3=1     s3=2     s3=3     s3=4")
for t in (1, 2, 3):
    print("    t=%d %8d %8d %8d %8d %8d" % tuple([t] + [joint[(t, s)] for s in range(5)]))
print("[3] 3|1 로 단순 분해되는데 그림이 세 가지인 함수:", violations, "(0 이면 '3|1 단순 분해 => t <= 2')")
print("[4] 2|2 든 3|1 이든 어느 하나로라도 단순 분해되는 함수 / 어느 것으로도 안 되는 함수")
for t in (1, 2, 3):
    print("    t=%d  분해됨 %6d   안 됨 %6d" % (t, joint_any[(t, True)], joint_any[(t, False)]))
