#!/usr/bin/env python3
"""
check_general_n_proofs.py
=========================
notes/2026-11-12_general_n_draft.md 의 G1 증명(고정하는 분할의 수, 최댓값이 위수 2인 순열에서 나온다는 것)과
G1 의 상한을 실제로 주는 함수의 구성을 작은 n 에서 숫자로 점검한다. 증명이 아니라 점검이다.

1부  공식 점검 (N = 6, 8, ..., 60 의 모든 소수 위수 순환 구조, 정수만 쓴다)
     - 위수가 홀수 소수 p 인 g 마다, 같은 점들 위의 위수 2인 순열 τ (맞바꿈 mp 개) 가 g 보다 많은 분할을 고정함
     - 그래서 궤도 수의 상한도 τ 쪽이 크다
     - 맞바꿈 a 개인 대합이 고정하는 분할의 수 φ(a) 가, 고정점이 2개 이상이면 a = 1 에서 최대,
       고정점이 없으면(a = N/2) N >= 10 에서 a = 1 보다 작음
     - I(a) 의 점화식과 I(A+1) <= (2A-1) I(A) (A >= 3)
2부  맞바꿈 하나에서 상한이 나오는 함수: 사슬 모양 함수 f_n 의 영향력(influence)을 직접 계산하고,
     n <= 7 에서는 H̄(f_n) 을 모든 (순열, 반전) 에 대해 직접 구해 {1, (1 2)} 임을 확인한다.
3부  고정점 없는 대합(n = 6, 8)에서 상한이 나오는 함수: 대합 σ = (1 2)(3 4)... 로 불변인 함수를 무작위로
     (씨앗을 고정해) 뽑아, H̄(f) 가 정확히 {1, σ} 인 것을 찾고 그림 종류 수 T(f) 를 센다.

필요한 것: numpy. 같은 폴더의 check_general_bounds.py 를 불러 쓴다.
"""
import itertools
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_general_bounds import (act, double_factorial, fixed_formula_involution,  # noqa: E402
                                  I2, is_prime, pairings, partitions)


def df(m):
    return double_factorial(m)


# ---------------------------------------------------------------- 1부
def fix_involution(N, a):
    """N 개 점 위의 대합(맞바꿈 a 개)이 고정하는 완전 매칭의 수. N 은 짝수."""
    return I2(a) * df(N - 2 * a - 1)


def fix_odd_prime(N, p, c):
    """길이 p(홀수 소수) 순환 c 개와 고정점 N - pc 개인 순열이 고정하는 완전 매칭의 수. N 은 짝수."""
    b = N - p * c
    if c % 2 == 1 or b % 2 == 1:
        return 0
    return df(c - 1) * p ** (c // 2) * df(b - 1)


def part1():
    ok = True
    # I(a) 의 점화식
    for a in range(2, 40):
        if I2(a) != I2(a - 1) + 2 * (a - 1) * I2(a - 2):
            ok = False
            print("점화식 불일치 a =", a)
    for A in range(3, 40):
        if I2(A + 1) > (2 * A - 1) * I2(A):
            ok = False
            print("I(A+1) <= (2A-1) I(A) 불일치 A =", A)
    for A in range(5, 40):
        if I2(A) > df(2 * A - 3):
            ok = False
            print("I(A) <= (2A-3)!! 불일치 A =", A)
    print("점화식, I(A+1) <= (2A-1) I(A) (A >= 3), I(A) <= (2A-3)!! (A >= 5): %s (A < 40)" % ("모두 성립" if ok else "불일치"))

    rows = []
    for N in range(6, 62, 2):
        M = df(N - 1)
        phi = {a: fix_involution(N, a) for a in range(1, N // 2 + 1)}
        # (i) 고정점이 2개 이상인 대합 중에서는 a = 1 이 최대
        for a in range(2, N // 2):
            if phi[a] > phi[1]:
                ok = False
                print("a=%d 가 a=1 보다 큼 (N=%d)" % (a, N))
        # (ii) 고정점이 없는 대합
        fpf_bigger = phi[N // 2] > phi[1]
        if (N >= 10) == fpf_bigger:
            ok = False
            print("고정점 없는 대합 비교가 예상과 다름 N =", N)
        # (iii) 홀수 소수 위수 g 마다 대합 τ 가 더 많이 고정한다
        best_bound = 0
        for p in range(2, N + 1):
            if not is_prime(p):
                continue
            for c in range(1, N // p + 1):
                if p == 2:
                    fg = fix_involution(N, c)
                else:
                    fg = fix_odd_prime(N, p, c)
                bound_g = (M + (p - 1) * fg) // p
                assert (M + (p - 1) * fg) % p == 0
                if p > 2:
                    a = (c // 2) * p if c % 2 == 0 else 1
                    ft = fix_involution(N, a)
                    if ft < fg or 2 * (M + (p - 1) * fg) > p * (M + ft):
                        ok = False
                        print("g 보다 τ 가 적게 고정함: N=%d p=%d c=%d" % (N, p, c))
                best_bound = max(best_bound, bound_g)
        # (iv) 고정점이 1개 이상 있는 순열(홀수 n = N-1 의 경우)만 보면 최댓값은 맞바꿈 하나
        best_with_fixed = 0
        for p in range(2, N + 1):
            if not is_prime(p):
                continue
            for c in range(1, (N - 1) // p + 1):
                fg = fix_involution(N, c) if p == 2 else fix_odd_prime(N, p, c)
                best_with_fixed = max(best_with_fixed, (M + (p - 1) * fg) // p)
        transposition = (M + phi[1]) // 2
        if best_with_fixed != transposition:
            ok = False
            print("홀수 n 의 최댓값이 맞바꿈 하나가 아님: N =", N)
        rows.append((N, best_bound, transposition, (M + phi[N // 2]) // 2))
    print("  N | 항등이 아닌 소수 위수 순열 하나의 최대 궤도 수 | 맞바꿈 하나 | 고정점 없는 대합")
    for N, best_bound, t, f in rows[:6]:
        print("%3d | %d | %d | %d" % (N, best_bound, t, f))
    print("...  (N = 6..60 모두 검사)")
    print("1부:", "모두 성립" if ok else "불일치가 있다")
    return ok


# ---------------------------------------------------------------- 2부, 3부
def chain_function(n):
    """G_0 = x1 and x2, G_j = x_{j+2} or G_{j-1} (j 홀수), x_{j+2} and G_{j-1} (j 짝수). 진리표를 0/1 배열로."""
    size = 1 << n
    xs = np.arange(size)
    bit = lambda i: ((xs >> (i - 1)) & 1).astype(np.uint8)  # 변수 x_i (i = 1..n)
    g = bit(1) & bit(2)
    for j in range(1, n - 1):
        g = (bit(j + 2) | g) if j % 2 == 1 else (bit(j + 2) & g)
    return g


def influences(f, n):
    xs = np.arange(1 << n)
    return [int((f != f[xs ^ (1 << i)]).sum()) for i in range(n)]


def perm_part_of_symmetry(f, n):
    """(순열 π, 반전 a) 가 f(π(x) xor a) = f(x) 를 모든 x 에서 만족하는 π 의 목록. π(x) 의 i 번 비트는 x 의 π[i] 번 비트."""
    size = 1 << n
    xs = np.arange(size, dtype=np.int32)
    perms = list(itertools.permutations(range(n)))
    P = np.array(perms, dtype=np.int32)
    idx = np.zeros((len(perms), size), dtype=np.int32)
    for i in range(n):
        idx |= ((xs[None, :] >> P[:, i][:, None]) & 1) << i
    ok = np.zeros(len(perms), dtype=bool)
    for a in range(size):
        ok |= (f[idx ^ a] == f[None, :]).all(axis=1)
    return [perms[k] for k in np.nonzero(ok)[0]]


def orbit_count_on_pairings(group, n):
    """순열들의 집합(군)이 M_n 에 만드는 궤도 수."""
    M = pairings(n)
    index = {m: k for k, m in enumerate(M)}
    parent = list(range(len(M)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for perm in group:
        for m in M:
            a, b = find(index[m]), find(index[act(perm, m)])
            if a != b:
                parent[a] = b
    return len({find(k) for k in range(len(M))}), len(M)


def part2():
    ok = True
    print("사슬 모양 함수 f_n 의 영향력 (변수 x1..xn 순서):")
    for n in range(3, 13):
        f = chain_function(n)
        inf = influences(f, n)
        distinct = len(set(inf[2:])) == n - 2 and inf[0] == inf[1] and inf[0] not in inf[2:]
        # 공식: Inf_{j+2} = 2^n q_j 2^{-(n-2-j)}, q_1 = 3/4, q_j = 1 - q_{j-1}/2.  Inf_1 = Inf_2 = 2^n 2^{-(n-1)} = 2
        q = {1: (3, 4)}
        for j in range(2, n - 1):
            a, b = q[j - 1]
            q[j] = (2 * b - a, 2 * b)  # 1 - a/(2b) = (2b - a)/(2b)
        formula_ok = inf[0] == 2 and all(
            inf[j + 1] * q[j][1] * 2 ** (n - 2 - j) == 2 ** n * q[j][0] for j in range(1, n - 1))
        if not (distinct and formula_ok):
            ok = False
        print("  n=%2d  %s  %s" % (n, inf, "서로 다름(x1=x2 제외), 공식과 일치" if distinct and formula_ok else "확인 실패"))
    print("H̄(f_n) 직접 계산 (n <= 7):")
    for n in range(4, 8):
        f = chain_function(n)
        group = perm_part_of_symmetry(f, n)
        expected = {tuple(range(n)), (1, 0) + tuple(range(2, n))}
        T, M = orbit_count_on_pairings(group, n)
        size = double_factorial(n - 1) if n % 2 == 0 else double_factorial(n)
        bound = (size + fixed_formula_involution(n, 1)) // 2
        good = set(group) == expected and T == bound
        ok = ok and good
        print("  n=%d  |H̄|=%d  T(f)=%d  공식값 %d  %s" % (n, len(group), T, bound, "일치" if good else "불일치"))
    print("2부:", "모두 성립" if ok else "불일치가 있다")
    return ok


def sigma_invariant_random(n, rng):
    """σ = (1 2)(3 4)...(변수 2i-1 과 2i 를 맞바꿈) 에 대해 f(σ(x)) = f(x) 인 함수를 무작위로 뽑는다."""
    size = 1 << n
    f = np.zeros(size, dtype=np.uint8)
    for x in range(size):
        y = 0
        for i in range(n // 2):
            lo, hi = (x >> (2 * i)) & 1, (x >> (2 * i + 1)) & 1
            y |= (hi << (2 * i)) | (lo << (2 * i + 1))
        if x <= y:
            f[x] = f[y] = rng.integers(2)
    return f


def hex_table(f, n):
    """진리표를 16진수로 (맨 앞 = x1..xn 이 모두 1인 입력, 맨 뒤 = 모두 0인 입력)."""
    digits = []
    for base in range((1 << n) - 4, -4, -4):
        digits.append("%x" % int("".join(str(int(f[x])) for x in range(base + 3, base - 1, -1)), 2))
    return "".join(digits)


def part3():
    ok = True
    for n in (6, 8):
        sigma = tuple(i ^ 1 for i in range(n))
        size = double_factorial(n - 1)
        bound = (size + fixed_formula_involution(n, n // 2)) // 2
        rng = np.random.default_rng(20261112)
        found, tried = None, 0
        while found is None and tried < 50:
            tried += 1
            f = sigma_invariant_random(n, rng)
            group = perm_part_of_symmetry(f, n)
            if set(group) == {tuple(range(n)), sigma}:
                found = (f, group)
        if found is None:
            ok = False
            print("n=%d: 50번 안에 찾지 못함" % n)
            continue
        f, group = found
        T, M = orbit_count_on_pairings(group, n)
        good = T == bound
        ok = ok and good
        print("n=%d  씨앗 20261112 의 %d 번째 함수에서 찾음" % (n, tried))
        print("     진리표(16진수, 맨 앞 = x1..x%d 가 모두 1인 입력): %s" % (n, hex_table(f, n)))
        print("     |H̄(f)| = %d (= <(1 2)(3 4)...>)  T(f) = %d  공식값 %d  %s" % (len(group), T, bound, "일치" if good else "불일치"))
    print("3부:", "모두 성립" if ok else "불일치가 있다")
    return ok


def main():
    results = [part1(), part2(), part3()]
    print("전체:", "모두 성립" if all(results) else "불일치가 있다")
    return 0 if all(results) else 1


if __name__ == "__main__":
    sys.exit(main())
