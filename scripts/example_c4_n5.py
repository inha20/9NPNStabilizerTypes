#!/usr/bin/env python3
"""
example_c4_n5.py
================
입력 반전이 섞인 대칭을 가진 5변수 함수의 구체적인 예 (PLAN.md 2-5).

  f(x1, ..., x5) = (x1 x2 + x1 x3 + x2 x4) + x5 (x1 x2 + x1 x3 + x1 x4 + x2 x3 + x2 x4 + x3 x4)   (덧셈은 mod 2)

앞 괄호는 길이 3인 경로 3-1-2-4 의 이차식, 뒤 괄호는 완전 그래프 K4 의 이차식이다.
f 의 대칭군(입력 반전과 변수 순열)은 위수 4 인 순환군이고, 생성원은
  g : (x1, x2, x3, x4, x5) -> (x3, x4, x2, x1, 1 - x5)
이다. 순열 부분은 4-순환 (1 3 2 4) 이고, 5번 변수를 뒤집어야만 대칭이 된다.

이 스크립트가 하는 일
  1. f 의 대칭군 전체(B_5 의 3,840개 원소)를 직접 구해 위수 4 인 순환군임을 확인한다.
  2. 변수 순열만으로 f 를 보존하는 것(입력 반전 없이)은 (1 2)(3 4) 뿐임을 확인한다.
  3. f 를 4×4 지도 두 장(x5 = 0, 1)으로 그리고, g 가 지도 위에서 하는 일을 확인한다.
  4. f 를 15가지 변수 배치(쌍 분할)로 그려, 그림이 5종류임을 정의대로 센다.
     (check_pictures_by_definition.py 의 Model 을 쓴다.)
  5. 변수 순열만 보는 군 G(f) 와 입력 반전을 허용한 대칭의 순열 부분 H̄(f) 가 쌍 분할 15개에 만드는 궤도 수
     (9 와 5) 를 센다.

표준 라이브러리만 쓴다. 같은 폴더의 check_pictures_by_definition.py 가 필요하다.
"""
import itertools
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_pictures_by_definition import GRAY, Model  # noqa: E402

N = 5
SIZE = 1 << N


def bits(x):
    return tuple((x >> i) & 1 for i in range(N))  # (x1, ..., x5)


def f_value(x1, x2, x3, x4, x5):
    path = (x1 & x2) ^ (x1 & x3) ^ (x2 & x4)
    k4 = (x1 & x2) ^ (x1 & x3) ^ (x1 & x4) ^ (x2 & x3) ^ (x2 & x4) ^ (x3 & x4)
    return path ^ (x5 & k4)


F = [f_value(*bits(x)) for x in range(SIZE)]
TT = sum(F[x] << x for x in range(SIZE))


def act(pi, a, x):
    """(π x)_i = x_{π(i)} 로 두고 π(x) xor a 를 돌려준다. π 는 0부터 센 순열."""
    y = 0
    for i in range(N):
        y |= ((x >> pi[i]) & 1) << i
    return y ^ a


def name(pi):
    """순열 pi(i -> pi[i], 0부터)를 1부터 센 순환 표기로."""
    seen, parts = set(), []
    for i in range(N):
        if i in seen:
            continue
        cyc, j = [], i
        while j not in seen:
            seen.add(j)
            cyc.append(j + 1)
            j = pi[j]
        if len(cyc) > 1:
            parts.append("(" + " ".join(map(str, cyc)) + ")")
    return "".join(parts) or "1"


def show_layers(f_of_cell, row_vars, col_vars, layer_var):
    """지도 두 장을 나란히 문자열로. f_of_cell(c0, c1, l) -> 0/1."""
    head = "행 %s%s, 열 %s%s, 층 %s" % (row_vars[0], row_vars[1], col_vars[0], col_vars[1], layer_var)
    lines = [head,
             "        %s=0              %s=1" % (layer_var, layer_var),
             "      %s%s=%s       %s%s=%s" % (col_vars[0], col_vars[1], "00 01 11 10", col_vars[0], col_vars[1], "00 01 11 10")]
    for c0 in range(4):
        left = " ".join(" %d" % f_of_cell(c0, c1, 0) for c1 in range(4))
        right = " ".join(" %d" % f_of_cell(c0, c1, 1) for c1 in range(4))
        lab = "%d%d" % GRAY[c0]
        lines.append("%s%s=%s   %s     %s%s=%s   %s" % (row_vars[0], row_vars[1], lab, left, row_vars[0], row_vars[1], lab, right))
    return "\n".join(lines)


def main():
    ok = True
    print("f 의 진리표: 십진수 %d, 16진수 0x%08x (비트 x 가 f(x), x 의 i번째 비트가 변수 i). 1 의 개수 %d" % (TT, TT, sum(F)))

    # 1. 대칭군 전체
    stab = []
    for pi in itertools.permutations(range(N)):
        for a in range(SIZE):
            if all(F[act(pi, a, x)] == F[x] for x in range(SIZE)):
                stab.append((pi, a))
    print("\n[1] 대칭군(입력 반전과 변수 순열) 의 위수: %d (B_5 는 3,840)" % len(stab))
    for pi, a in stab:
        flips = "".join("x%d " % (i + 1) for i in range(N) if (a >> i) & 1) or "없음"
        print("    순열 %-8s 반전할 변수: %s" % (name(pi), flips))
    rho = (2, 3, 1, 0, 4)  # 1->3, 2->4, 3->2, 4->1, 5->5
    g = (rho, 1 << 4)
    ok &= len(stab) == 4 and g in stab
    # g 의 거듭제곱
    cur = list(range(SIZE))
    orders = 0
    for k in range(1, 5):
        cur = [act(rho, 1 << 4, y) for y in cur]
        if cur == list(range(SIZE)):
            orders = k
            break
    print("    g = (순열 (1 3 2 4), x5 반전) 의 위수 %d (순환군 C4)" % orders)
    ok &= orders == 4

    # 2. 순열만
    pure = [pi for pi, a in stab if a == 0]
    print("\n[2] 입력 반전 없이 변수 순열만으로 f 를 보존하는 것: %s" % ", ".join(name(pi) for pi in pure))
    ok &= sorted(pure) == sorted([tuple(range(N)), (1, 0, 3, 2, 4)])
    diff = [x for x in range(SIZE) if F[act(rho, 0, x)] != F[x]]
    print("    4-순환 (1 3 2 4) 만 쓰면 f(πx) != f(x) 인 입력이 %d개 (32개 중). 예: x = (%s) 에서 f(x) = %d, f(πx) = %d"
          % (len(diff), ",".join(map(str, bits(diff[0]))), F[diff[0]], F[act(rho, 0, diff[0])]))
    ok &= len(diff) > 0
    print("    5번 변수만 뒤집으면: f(x xor e5) != f(x) 인 입력 %d개" % sum(1 for x in range(SIZE) if F[x ^ 16] != F[x]))

    # 3. 지도
    m = Model(N)

    def point(c0, c1, l):
        g1, g2 = GRAY[c0]
        h1, h2 = GRAY[c1]
        return g1 | (g2 << 1) | (h1 << 2) | (h2 << 3) | (l << 4)

    # Model 의 셀 배치와 같은지 확인
    for idx, cell in enumerate(m.cells):
        if m.base[idx] != point(*cell):
            ok = False
            print("Model 의 셀 배치와 다르다:", cell)
            break
    print("\n[3] 기준 배치의 지도 (행 = (x1 x2), 열 = (x3 x4), 층 = x5; 축의 순서는 Gray 부호 00 01 11 10)")
    print(show_layers(lambda c0, c1, l: F[point(c0, c1, l)], ("x1", "x2"), ("x3", "x4"), "x5"))
    rot_ok = all(F[point(c0, c1, l)] == F[point(c1, (-c0) % 4, 1 - l)] for c0 in range(4) for c1 in range(4) for l in (0, 1))
    rot2_ok = all(F[point(c0, c1, l)] == F[point((-c0) % 4, (-c1) % 4, l)] for c0 in range(4) for c1 in range(4) for l in (0, 1))
    print("    g 가 지도에서 하는 일: 칸 (r, c, 층 l) 의 값 = 칸 (c, -r mod 4, 층 1-l) 의 값 (모든 칸에서): %s" % rot_ok)
    print("    g 를 두 번 하면: 칸 (r, c) 의 값 = 칸 (-r, -c) 의 값, 층은 그대로 (모든 칸에서): %s" % rot2_ok)
    ok &= rot_ok and rot2_ok
    onlyperm_ok = all(F[point(c0, c1, l)] == F[point(c1, (-c0) % 4, l)] for c0 in range(4) for c1 in range(4) for l in (0, 1))
    print("    층을 바꾸지 않고 같은 층 안에서만 90도 돌려 보면 같은가: %s" % onlyperm_ok)
    ok &= not onlyperm_ok

    # 4. 그림 종류
    print("\n[4] 그림 종류 (정의대로: 15가지 변수 배치 각각의 지도를 격자 대칭으로 겹치는 것끼리 묶는다)")
    apply = m._apply

    def canonical(p):
        seen = {p}
        stack = [p]
        while stack:
            q = stack.pop()
            for tabs in m.tables:
                r = apply(tabs, q)
                if r not in seen:
                    seen.add(r)
                    stack.append(r)
        return min(seen)

    groups = {}  # 쌍 분할 -> [(배치 s, 그림)]
    for s, pts in zip(itertools.permutations(range(N)), m.arrangements):
        p = 0
        for c, x in enumerate(pts):
            if (TT >> x) & 1:
                p |= 1 << c
        row = tuple(sorted(i + 1 for i in range(N) if s[i] in (0, 1)))
        col = tuple(sorted(i + 1 for i in range(N) if s[i] in (2, 3)))
        layer = [i + 1 for i in range(N) if s[i] == 4][0]
        key = (tuple(sorted([row, col])), layer)
        groups.setdefault(key, []).append((s, p))
    classes = {}
    for key, items in groups.items():
        canon = {canonical(p) for _, p in items}
        if len(canon) != 1:
            ok = False
            print("같은 쌍 분할인데 그림이 다르다:", key)
        classes.setdefault(canon.pop(), []).append(key)
    T = len(classes)
    print("    쌍 분할 %d개, 그림 종류 T(f) = %d (Model.picture_types 로 센 값: %d)" % (len(groups), T, m.picture_types(TT)))
    ok &= T == 5 and len(groups) == 15 and m.picture_types(TT) == 5
    sizes = sorted(len(v) for v in classes.values())
    print("    종류별 쌍 분할의 수:", sizes, "(g 의 궤도 길이 1, 2, 4, 4, 4 와 같아야 한다)")
    ok &= sizes == [1, 2, 4, 4, 4]
    for k, (canon, keys) in enumerate(sorted(classes.items(), key=lambda kv: (len(kv[1]), sorted(kv[1])))):
        keys = sorted(keys, key=lambda key: (key[1] != 5, key))
        print("\n    그림 %d: 쌍 분할 %d개" % (k + 1, len(keys)))
        for (pair, layer) in keys:
            print("       {%s | %s}, 외톨이 x%d" % (" ".join("x%d" % v for v in pair[0]), " ".join("x%d" % v for v in pair[1]), layer))
        (pair, layer) = keys[0]
        s0, p0 = groups[(pair, layer)][0]
        # 대표 배치에서 지도를 그린다: 변수 i 는 슬롯 s0[i] 에 놓인다
        slot_var = {s0[i]: "x%d" % (i + 1) for i in range(N)}
        print(show_layers(lambda c0, c1, l, p0=p0: (p0 >> m.cells.index((c0, c1, l))) & 1,
                          (slot_var[0], slot_var[1]), (slot_var[2], slot_var[3]), slot_var[4]))
    # 5. 군이 쌍 분할에 만드는 궤도 수 (명제 A2: T(f) = H̄(f) 의 궤도 수)
    def orbit_count(perms):
        keys = list(groups)
        parent = {k: k for k in keys}

        def find(k):
            while parent[k] != k:
                parent[k] = parent[parent[k]]
                k = parent[k]
            return k

        for pi in perms:
            for key in keys:
                pair, layer = key
                moved = (tuple(sorted(tuple(sorted(pi[v - 1] + 1 for v in block)) for block in pair)), pi[layer - 1] + 1)
                a, b = find(key), find(moved)
                if a != b:
                    parent[a] = b
        return len({find(k) for k in keys})

    hbar = sorted({pi for pi, a in stab})
    t_pure = orbit_count(pure)
    t_hbar = orbit_count(hbar)
    print("\n[5] 쌍 분할 15개 위의 궤도 수: 순열만 보는 군 G(f) (위수 %d) %d개, 입력 반전을 허용한 대칭의 순열 부분 H̄(f) (위수 %d) %d개"
          % (len(pure), t_pure, len(hbar), t_hbar))
    ok &= t_pure == 9 and t_hbar == 5 == T

    print("\n전체:", "모두 확인됨" if ok else "확인 실패")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
