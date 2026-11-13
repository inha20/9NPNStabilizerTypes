#!/usr/bin/env python3
"""
check_selfcomplementary_family.py
=================================
notes/2026-11-12_example_c4_n5.md 의 5변수 C4 예와, 그것을 일반화한 구성을 점검한다. 증명이 아니라 점검이다.

구성. m = n-1 개 점 위의 그래프 G 와 순열 ρ 가 "반동형"이다(ij 가 G 의 간선 <=> ρ(i)ρ(j) 가 G 의 간선이 아님).
이때 f(x_1..x_n) = Q_G(x_1..x_m) + x_n * C(w, 2)  (mod 2) 로 두면 (w 는 x_1..x_m 의 1의 개수, Q_G 는 G 의
간선마다 두 변수의 곱을 더한 것) 순열 ρ 와 x_n 의 반전을 함께 한 변환 g 가 f 를 그대로 둔다.
  Q_G(ρx) = Q_G(x) + C(w,2)  이고  C(w,2) 는 ρ 에 대해 불변이므로 f(ρx, ~x_n) = f(x).
자기보수 그래프(자기 여집합과 동형)는 점의 수가 4로 나눈 나머지 0 또는 1일 때만 있으므로
m = 4k 또는 4k+1, 곧 n = 4k+1 또는 4k+2 (n = 5, 6, 9, 10, ...) 에서 쓸 수 있다.

1부  n = 5 의 예(노트의 함수)에 대해, 32개 평행이동 t 마다 f_t(x) = f(x xor t) 의 순열만 쓴 대칭군의
    위수를 구한다. 다상 대칭(순열군을 평행이동으로 공액한 군)이라면 어떤 t 에서 위수가 |H(f)| = 4 가
    되어야 한다. 최댓값이 2 이면 f 의 대칭군은 순열군의 평행이동 공액이 아니다.
2부  m = 4, 5, 8, 9 (n = 5, 6, 9, 10): 점 위의 순환 ρ 를 4-순환들의 곱으로 잡고 ρ 의 쌍 궤도를 따라
    간선을 번갈아 넣어 G 를 만든다(씨앗을 고정해 각 궤도의 시작 위상을 무작위로 정한다).
      (a) 반동형 조건을 모든 점 쌍에서 확인한다.
      (b) 모든 입력 x 에서 f(gx) = f(x) 를 확인한다.
      (c) 평행이동 부분 a = e_n 이 (1+π) 의 상에 들어 있지 않음을 확인한다. 들어 있지 않으면 <g> 는
          순열군의 평행이동 공액이 아니다(고정점 n 에서 (1+π)t 의 n번째 성분은 항상 0).
      (d) n <= 6 에서는 B_n 안의 안정화군 전체를 구하고, 순열 부분 H̄ 가 쌍 분할(M_n)에 만드는 궤도 수 T 를
          센다. 평행이동 공액 시험(1부와 같은 것)도 한다.

필요한 것: numpy. 같은 폴더의 check_general_bounds.py 를 불러 쓴다.
"""
import itertools
import os
import random
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_general_bounds import act, pairings  # noqa: E402


def truth_table(n, edges):
    """f(x) = sum_{ij in edges} x_i x_j + x_{n-1} * C(w,2)  (mod 2), 변수는 0..n-1, 마지막 변수가 x_n."""
    xs = np.arange(1 << n, dtype=np.int64)
    bits = [(xs >> i) & 1 for i in range(n)]
    f = np.zeros(1 << n, dtype=np.int64)
    for i, j in edges:
        f ^= bits[i] & bits[j]
    w = sum(bits[i] for i in range(n - 1))
    f ^= bits[n - 1] & ((w * (w - 1) // 2) & 1)
    return f.astype(np.int8)


def apply_perm(rho, x, n):
    """(ρx)_i = x_{ρ(i)}  (비트 i 는 변수 i)."""
    y = 0
    for i in range(n):
        y |= ((x >> rho[i]) & 1) << i
    return y


def antimorphism_graph(m, rng):
    """m = 4k 또는 4k+1. ρ = (0 1 2 3)(4 5 6 7)... (m = 4k+1 이면 마지막 점은 고정). 간선 집합 E 와 ρ."""
    assert m % 4 in (0, 1)
    rho = list(range(m))
    for b in range(m // 4):
        c = [4 * b + t for t in range(4)]
        for t in range(4):
            rho[c[t]] = c[(t + 1) % 4]
    seen, E = set(), set()
    for i in range(m):
        for j in range(i + 1, m):
            p = frozenset((i, j))
            if p in seen:
                continue
            orbit, q = [], p
            while q not in orbit:
                orbit.append(q)
                q = frozenset(rho[v] for v in q)
            assert len(orbit) % 2 == 0, "쌍 궤도의 길이가 홀수"
            start = rng.randrange(2)  # 이 궤도에서 어느 쪽 위상이 간선인가
            for t, q in enumerate(orbit):
                seen.add(q)
                if (t + start) % 2 == 0:
                    E.add(tuple(sorted(q)))
    return sorted(E), rho


def check_antimorphism(m, E, rho):
    Es = {frozenset(e) for e in E}
    for i in range(m):
        for j in range(i + 1, m):
            inE = frozenset((i, j)) in Es
            img = frozenset((rho[i], rho[j])) in Es
            if inE == img:  # 반동형이면 둘 중 정확히 하나만 간선
                return False
    return True


def stabilizer(f, n):
    """B_n 에서 f(π(x) xor a) = f(x) 인 (π, a) 전부. π(x) 의 i 번 비트는 x 의 π[i] 번 비트."""
    size = 1 << n
    xs = np.arange(size, dtype=np.int32)
    perms = list(itertools.permutations(range(n)))
    P = np.array(perms, dtype=np.int32)
    idx = np.zeros((len(perms), size), dtype=np.int32)
    for i in range(n):
        idx |= ((xs[None, :] >> P[:, i][:, None]) & 1) << i
    out = []
    for a in range(size):
        ok = (f[idx ^ a] == f[None, :]).all(axis=1)
        out.extend((perms[k], a) for k in np.nonzero(ok)[0])
    return out


def translated(f, n, t):
    xs = np.arange(1 << n)
    return f[xs ^ t]


def perm_only_order(f, n):
    """순열만 쓴 대칭군 G(f) 의 위수."""
    size = 1 << n
    xs = np.arange(size, dtype=np.int32)
    perms = list(itertools.permutations(range(n)))
    P = np.array(perms, dtype=np.int32)
    idx = np.zeros((len(perms), size), dtype=np.int32)
    for i in range(n):
        idx |= ((xs[None, :] >> P[:, i][:, None]) & 1) << i
    return int((f[idx] == f[None, :]).all(axis=1).sum())


def orbit_count(perms, n):
    M = pairings(n)
    index = {m: k for k, m in enumerate(M)}
    parent = list(range(len(M)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for p in perms:
        for m in M:
            a, b = find(index[m]), find(index[act(p, m)])
            if a != b:
                parent[a] = b
    return len({find(k) for k in range(len(M))}), len(M)


def conj_test(f, n):
    """모든 평행이동 t 에 대해 f_t 의 순열만 쓴 대칭군의 위수와 그 최댓값."""
    orders = [perm_only_order(translated(f, n, t), n) for t in range(1 << n)]
    return max(orders), orders


def part1():
    print("1부. n = 5 의 노트의 함수 (경로 3-1-2-4, 변수 번호는 0부터: 간선 {2,0},{0,1},{1,3})")
    n = 5
    E = [(0, 2), (0, 1), (1, 3)]
    f = truth_table(n, E)
    tt = sum(int(f[x]) << x for x in range(1 << n))
    print("  진리표 0x%08x, 1의 개수 %d" % (tt, int(f.sum())))
    ok = tt == 0x9ac0e428
    H = stabilizer(f, n)
    mx, orders = conj_test(f, n)
    print("  |H(f)| = %d,  |G(f)| = %d" % (len(H), perm_only_order(f, n)))
    print("  평행이동 32개 각각에서 f(x xor t) 의 순열만 쓴 대칭군 위수의 분포:",
          sorted({o: orders.count(o) for o in set(orders)}.items()))
    print("  최댓값 %d  (|H(f)| = %d 이므로 %s)" % (
        mx, len(H), "평행이동 공액이 아니다" if mx < len(H) else "평행이동 공액일 수 있다"))
    ok = ok and len(H) == 4 and mx == 2
    print("  노트의 진리표와 일치, |H| = 4, 최댓값 2:", "예" if ok else "아니오")
    return ok


def part2():
    print("2부. 일반화 (n = m+1, m = 4k 또는 4k+1)")
    rng = random.Random(20261112)
    allok = True
    for m in (4, 5, 8, 9):
        n = m + 1
        E, rho = antimorphism_graph(m, rng)
        a_ok = check_antimorphism(m, E, rho)
        f = truth_table(n, E)
        # g: (gx)_i = x_{ρ(i)} (i < m), (gx)_{n-1} = 1 - x_{n-1}
        g_ok = True
        for x in range(1 << n):
            y = apply_perm(rho + [m], x, n) ^ (1 << (n - 1))
            if f[y] != f[x]:
                g_ok = False
                break
        # (c) a = e_{n-1} 이 (1+π) 의 상에 있는가
        pi = rho + [m]
        image = {t ^ apply_perm(pi, t, n) for t in range(1 << n)}
        a_not_in_image = (1 << (n - 1)) not in image
        line = "  n=%2d (m=%d)  간선 %2d개  반동형 %s  f(gx)=f(x) %s  e_n 이 (1+π) 의 상 밖에 있음 %s" % (
            n, m, len(E), "예" if a_ok else "아니오", "예" if g_ok else "아니오", "예" if a_not_in_image else "아니오")
        print(line)
        ok = a_ok and g_ok and a_not_in_image
        if n <= 6:
            H = stabilizer(f, n)
            Hbar = sorted({p for p, _ in H})
            T, M = orbit_count(Hbar, n)
            mx, _ = conj_test(f, n)
            print("        |H(f)| = %d,  |H̄(f)| = %d,  T(f) = %d (쌍 분할 %d개),  평행이동 공액 시험의 최댓값 %d" % (
                len(H), len(Hbar), T, M, mx))
            # g 가 정말 H 에 들어 있나
            gp = tuple(rho + [m])
            ok = ok and any(p == gp and a == (1 << (n - 1)) for p, a in H)
            print("        g = (ρ, e_n) 가 H(f) 의 원소: %s" % ("예" if ok else "아니오"))
        allok = allok and ok
    print("2부:", "모두 성립" if allok else "불일치가 있다")
    return allok


def main():
    r1 = part1()
    print()
    r2 = part2()
    print()
    print("전체:", "모두 확인됨" if (r1 and r2) else "확인 실패가 있다")
    return 0 if (r1 and r2) else 1


if __name__ == "__main__":
    sys.exit(main())
