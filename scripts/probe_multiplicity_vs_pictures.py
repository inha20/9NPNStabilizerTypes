#!/usr/bin/env python3
"""
probe_multiplicity_vs_pictures.py
=================================
탐침: 4변수 부울 함수 65,536개 전부에서 '그림 종류 수'(7번 논문 정리 5)와
함수 분해(Ashenhurst-Curtis) 차트의 '다중도'가 어떤 관계인지 전수로 센다.

정의 (입력 x = (A,B,C,D), 진리표는 16비트 정수, 비트 x 가 f(x)):
  - 변수 쌍 분할은 AB|CD, AC|BD, AD|BC 세 가지이다.
  - 분할 P|Q 의 분해 차트는 행을 P의 값, 열을 Q의 값으로 둔 4x4 표이다.
    다중도 쌍 = {서로 다른 행의 수, 서로 다른 열의 수} (순서 없는 쌍).
  - 격자 대칭군 W = AB|CD 분할을 보존하는 변수 순열 8개 x 입력 반전 16개 (위수 128).
  - 그림 종류 수 t(f) = 세 배치의 그림이 W 아래에서 몇 가지인가 (1, 2, 3).
  - d(f) = 세 분할의 다중도 쌍이 몇 가지인가 (1, 2, 3).

확인하는 것:
  (1) 대조: W-궤도 805개, B4-궤도 402개, t = 1/2/3 인 함수 수 (7번 논문의 값과 비교)
  (2) t 를 두 가지 방법으로 계산해 일치하는지 (W-궤도 비교 / 안정화군의 S3 상)
  (3) 같은 그림이면 다중도 쌍이 같다 => d <= t
  (4) t x d 결합 분포 (함수 수, B4-궤도 수)
  (5) 단순 분해(다중도 2 이하)가 되는 분할의 수와 t 의 결합 분포

표준 라이브러리만 쓴다. 실행: python scripts/probe_multiplicity_vs_pictures.py
"""
import itertools
from collections import Counter

N = 16
FULL = 1 << N


def var_perm_on_points(p):
    """변수 순열 p: 새 입력의 i번째 좌표 = 옛 입력의 p[i]번째 좌표. 좌표 0=A(비트 3) ... 3=D(비트 0)."""
    out = []
    for x in range(N):
        bits = [(x >> (3 - i)) & 1 for i in range(4)]
        y = 0
        for i in range(4):
            y = (y << 1) | bits[p[i]]
        out.append(y)
    return tuple(out)


def act(f, g):
    """점 순열 g 로 진리표를 옮긴다: (g.f)(g(x)) = f(x)."""
    h = 0
    for x in range(N):
        if (f >> x) & 1:
            h |= 1 << g[x]
    return h


def compose(g, h):
    return tuple(g[h[x]] for x in range(N))


def orbit_labels(gens):
    """생성원들의 궤도로 65,536개 함수에 번호를 붙인다."""
    label = [-1] * FULL
    sizes = []
    for f in range(FULL):
        if label[f] != -1:
            continue
        k = len(sizes)
        label[f] = k
        stack = [f]
        cnt = 1
        while stack:
            u = stack.pop()
            for g in gens:
                v = act(u, g)
                if label[v] == -1:
                    label[v] = k
                    cnt += 1
                    stack.append(v)
        sizes.append(cnt)
    return label, sizes


def multiplicity_pair(f, P, Q):
    """분할 P|Q (좌표 번호 쌍)의 분해 차트에서 {서로 다른 행 수, 서로 다른 열 수}."""
    rows = []
    for pv in range(4):
        row = []
        for qv in range(4):
            bits = [0, 0, 0, 0]
            bits[P[0]], bits[P[1]] = pv >> 1, pv & 1
            bits[Q[0]], bits[Q[1]] = qv >> 1, qv & 1
            x = (bits[0] << 3) | (bits[1] << 2) | (bits[2] << 1) | bits[3]
            row.append((f >> x) & 1)
        rows.append(tuple(row))
    cols = list(zip(*rows))
    return tuple(sorted((len(set(rows)), len(set(cols)))))


PARTS = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]  # AB|CD, AC|BD, AD|BC
NAMES = ["AB|CD", "AC|BD", "AD|BC"]

# 생성원
neg = [tuple(x ^ m for x in range(N)) for m in (8, 4, 2, 1)]
swapAB = var_perm_on_points((1, 0, 2, 3))
swapCD = var_perm_on_points((0, 1, 3, 2))
swapPairs = var_perm_on_points((2, 3, 0, 1))
swapBC = var_perm_on_points((0, 2, 1, 3))
swapBD = var_perm_on_points((0, 3, 2, 1))
W_gens = neg + [swapAB, swapCD, swapPairs]
B4_gens = W_gens + [swapBC]

print("[1] 궤도 대조")
wlab, wsizes = orbit_labels(W_gens)
blab, bsizes = orbit_labels(B4_gens)
print("    W-궤도 수 :", len(wsizes), "(7번 논문: 805)")
print("    B4-궤도 수:", len(bsizes), "(7번 논문: 402)")

# 방법 1: 세 배치의 그림(f, f.swapBC, f.swapBD)이 속한 W-궤도를 비교
t1 = [0] * FULL
same_pic = [None] * FULL
for f in range(FULL):
    ids = (wlab[f], wlab[act(f, swapBC)], wlab[act(f, swapBD)])
    t1[f] = len(set(ids))
    same_pic[f] = ids
c1 = Counter(t1)
print("    그림 종류 수 1/2/3 인 함수:", c1[1], c1[2], c1[3], "(7번 논문: 3,220 / 36,204 / 26,112)")

# 방법 2 (정리 5): B4-궤도 대표원마다 안정화군이 세 분할에 작용하는 궤도 수
print("[2] 정리 5 방식으로 다시 계산 (B4-궤도 대표원 402개)")
B4 = []
for p in itertools.permutations(range(4)):
    g = var_perm_on_points(p)
    # 이 변수 순열이 세 분할을 어떻게 옮기는지
    img = []
    for (P, Q) in PARTS:
        newP = frozenset(p.index(i) for i in P)
        for j, (P2, Q2) in enumerate(PARTS):
            if newP == frozenset(P2) or newP == frozenset(Q2):
                img.append(j)
    for m in range(N):
        B4.append((tuple(y ^ m for y in g), tuple(img)))
assert len(B4) == 384
reps = {}
for f in range(FULL):
    reps.setdefault(blab[f], f)
t2_by_orbit = {}
for k, f in reps.items():
    parent = [0, 1, 2]
    def find(a):
        while parent[a] != a:
            a = parent[a]
        return a
    for g, img in B4:
        if act(f, g) == f:
            for a in range(3):
                ra, rb = find(a), find(img[a])
                if ra != rb:
                    parent[ra] = rb
    t2_by_orbit[k] = len({find(a) for a in range(3)})
mismatch = sum(1 for f in range(FULL) if t1[f] != t2_by_orbit[blab[f]])
print("    두 방법이 다른 함수 수:", mismatch)

# 다중도
print("[3] 다중도 쌍과 그림")
pairs = [None] * FULL
d = [0] * FULL
viol_same = 0
for f in range(FULL):
    pr = tuple(multiplicity_pair(f, P, Q) for (P, Q) in PARTS)
    pairs[f] = pr
    d[f] = len(set(pr))
    ids = same_pic[f]
    for a in range(3):
        for b in range(a + 1, 3):
            if ids[a] == ids[b] and pr[a] != pr[b]:
                viol_same += 1
print("    '같은 그림인데 다중도 쌍이 다른' 경우:", viol_same, "(0 이면 d <= t)")
print("    d > t 인 함수 수:", sum(1 for f in range(FULL) if d[f] > t1[f]))

print("[4] t(그림 종류 수) x d(다중도 쌍의 종류 수): 함수 수")
joint = Counter((t1[f], d[f]) for f in range(FULL))
print("          d=1      d=2      d=3")
for t in (1, 2, 3):
    print("    t=%d %8d %8d %8d" % (t, joint[(t, 1)], joint[(t, 2)], joint[(t, 3)]))
print("    같은 표를 B4-궤도 수로:")
jorb = Counter((t1[f], d[f]) for f in reps.values())
print("          d=1      d=2      d=3")
for t in (1, 2, 3):
    print("    t=%d %8d %8d %8d" % (t, jorb[(t, 1)], jorb[(t, 2)], jorb[(t, 3)]))

print("[5] 단순 분해(다중도 2 이하)가 되는 분할의 수 s x t: 함수 수")
def simple(pr):
    return sum(1 for p in pr if p[0] <= 2)
js = Counter((t1[f], simple(pairs[f])) for f in range(FULL))
print("          s=0      s=1      s=2      s=3")
for t in (1, 2, 3):
    print("    t=%d %8d %8d %8d %8d" % (t, js[(t, 0)], js[(t, 1)], js[(t, 2)], js[(t, 3)]))

print("[6] 다중도 쌍의 분포 (분할 AB|CD 기준, 함수 수)")
for k, v in sorted(Counter(pairs[f][0] for f in range(FULL)).items()):
    print("    %s: %d" % (k, v))

print("[7] 그림은 다른데(t=3) 다중도 쌍이 모두 같은(d=1) B4-궤도의 대표원 예시 (최대 5개)")
shown = 0
for f in sorted(reps.values()):
    if t1[f] == 3 and d[f] == 1:
        print("    f = 0x%04X, 다중도 쌍 %s" % (f, pairs[f][0],))
        shown += 1
        if shown == 5:
            break
