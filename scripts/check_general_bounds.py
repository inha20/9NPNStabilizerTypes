#!/usr/bin/env python3
"""
check_general_bounds.py
=======================
notes/2026-11-12_general_n_draft.md 의 공식을 작은 n 에서 직접 세어 확인한다.

변수 쌍 분할의 집합 M_n (n 이 홀수이면 외톨이 하나 포함) 위에서, 순환 구조마다
  - 그 순열이 고정하는 분할의 수 (직접 셈 / 공식)
  - 그 순열이 생성하는 순환군의 궤도 수
를 구하고, "항등이 아닌 순열 하나가 만들 수 있는 궤도 수의 최댓값"을 출력한다.
그림 종류 수 T(f) 가 최댓값 |M_n| 이 아니면 이 값 이하이다 (정리 G1).

표준 라이브러리만 쓴다. n = 5..10 은 직접 세고 (n = 10 은 몇 초), n = 11, 12 는 공식만 쓴다.
"""
import itertools
import sys
from math import comb


def double_factorial(m):
    r = 1
    while m > 1:
        r *= m
        m -= 2
    return r


def pairings(n):
    """M_n 의 원소를 frozenset(frozenset) 으로 만든다."""
    def rec(rest, single_used):
        if not rest:
            yield []
            return
        i = rest[0]
        if n % 2 == 1 and not single_used:
            for tail in rec(rest[1:], True):
                yield [frozenset([i])] + tail
        for j in rest[1:]:
            others = [x for x in rest[1:] if x != j]
            for tail in rec(others, single_used):
                yield [frozenset([i, j])] + tail
    return [frozenset(p) for p in rec(list(range(n)), False)]


def partitions(n, maxpart=None):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield []
        return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - k, k):
            yield [k] + rest


def perm_of_type(ctype):
    perm, pos = [], 0
    for length in ctype:
        perm += [pos + (k + 1) % length for k in range(length)]
        pos += length
    return perm


def act(perm, pairing):
    return frozenset(frozenset(perm[i] for i in block) for block in pairing)


def I2(a):
    """a 개의 맞바꿈 쌍을 '제 짝끼리' 또는 '두 쌍을 엇갈려' 맞추는 방법의 수."""
    return sum(comb(a, 2 * j) * double_factorial(2 * j - 1) * 2 ** j for j in range(a // 2 + 1))


def fixed_formula_involution(n, a):
    b = n - 2 * a
    return I2(a) * (double_factorial(b - 1) if n % 2 == 0 else double_factorial(b))


def fixed_formula_odd_prime(n, p, c):
    """길이 p 순환 c 개와 고정점 b = n - pc 개인 순열이 고정하는 분할의 수."""
    b = n - p * c
    if c % 2 == 1:
        return 0
    inner = double_factorial(c - 1) * p ** (c // 2)
    if n % 2 == 0:
        return inner * double_factorial(b - 1) if b % 2 == 0 else 0
    return inner * double_factorial(b) if b % 2 == 1 else 0


def is_prime(p):
    return p > 1 and all(p % q for q in range(2, int(p ** 0.5) + 1))


def main():
    ok = True
    print(" n | |M_n| | 항등이 아닌 순열 하나의 최대 궤도 수 | 그 순열 | 맞바꿈 하나 | 고정점 없는 대합(짝수 n)")
    for n in range(5, 13):
        size = double_factorial(n - 1) if n % 2 == 0 else double_factorial(n)
        best = (0, None)
        direct = n <= 10
        if direct:
            M = pairings(n)
            assert len(M) == size
            index = {m: k for k, m in enumerate(M)}
        for ctype in partitions(n):
            lengths = set(ctype) - {1}
            if len(lengths) != 1:
                continue
            p = lengths.pop()
            if not is_prime(p):
                continue
            c = ctype.count(p)
            formula = fixed_formula_involution(n, c) if p == 2 else fixed_formula_odd_prime(n, p, c)
            if direct:
                perm = perm_of_type(ctype)
                fixed = sum(1 for m in M if act(perm, m) == m)
                if fixed != formula:
                    ok = False
                    print("   공식 불일치: n=%d 순환 구조 %s 직접 %d 공식 %d" % (n, ctype, fixed, formula))
                # 궤도 수도 직접 센다
                seen, orbits = set(), 0
                for m in M:
                    if m in seen:
                        continue
                    orbits += 1
                    x = m
                    while x not in seen:
                        seen.add(x)
                        x = act(perm, x)
                if orbits * p != size + (p - 1) * fixed:
                    ok = False
                    print("   궤도 수 불일치: n=%d %s" % (n, ctype))
            orbits_formula = (size + (p - 1) * formula) // p
            if orbits_formula > best[0]:
                best = (orbits_formula, "길이 %d 순환 %d개" % (p, c))
        transposition = (size + fixed_formula_involution(n, 1)) // 2
        fpf = (size + fixed_formula_involution(n, n // 2)) // 2 if n % 2 == 0 else None
        print("%2d | %5d | %5d | %s | %d | %s%s" % (n, size, best[0], best[1], transposition,
                                                 fpf if fpf is not None else "-", "" if direct else "  (공식만)"))
    print("모두 일치한다." if ok else "일치하지 않는 항목이 있다.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
