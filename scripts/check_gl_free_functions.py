#!/usr/bin/env python3
"""
check_gl_free_functions.py
==========================
GL(m,2) 와 AGL(m,2) 가 부울 함수 위에서 하는 작용에 "자유로운" 함수(안정화군이 자명한 함수)가 있는지를 m = 1..4 에서
군의 모든 원소를 훑어 센다. m = 5 는 scripts/agl52_orbits.c 의 전수 열거(notes/gl52_orbits_output.txt, agl52_orbits_output.txt)를 읽어 적는다.

왜 필요한가 (notes/2026-11-12_stabilizer_types.md 3절의 "일반 n").
  - 해밍-전단 족(축 고정, n = 2^r): AGL(r,2) 가 F_2^r 의 점 2^r 개 위의 함수에 자유로운 함수가 없으면 안 나타난다.
  - 중심 고정 족(n = 2^m - 1): GL(m-1,2) 가 F_2^(m-1) 의 0 이 아닌 점 2^(m-1) - 1 개 위의 함수에 자유로운 함수가 없으면 안 나타난다.
이 스크립트가 두 질문의 답을 작은 m 에서 세어 준다.

방법: 군의 원소를 점 위의 순열로 만들고, 항등원이 아닌 원소 p 마다 "p 의 사이클 위에서 상수인 함수" 전체에 표시를 한다
(p 가 고정하는 함수가 정확히 그것들이다). 어떤 항등원이 아닌 원소가 고정하는 함수는 모두 표시되고, 표시되지 않은 함수가 자유로운 함수이다.

실행: python3 -I scripts/check_gl_free_functions.py   (수 초)

2026-11-15. Claude가 쓴 코드이고 저자 검토 전이다.
"""
import itertools
import sys
import time

import numpy as np


def span(vecs):
    s = {0}
    for v in vecs:
        s |= {x ^ v for x in s}
    return s


def gl_matrices(m):
    """GL(m,2) 의 모든 원소를 열 벡터 m 개의 튜플로 낸다 (열 i = 기본벡터 e_i 의 상)."""
    def rec(cols, sp):
        if len(cols) == m:
            yield tuple(cols)
            return
        for v in range(1, 1 << m):
            if v not in sp:
                yield from rec(cols + [v], sp | {x ^ v for x in sp})
    yield from rec([], {0})


def apply(cols, v):
    out = 0
    i = 0
    while v:
        if v & 1:
            out ^= cols[i]
        v >>= 1
        i += 1
    return out


def cycles_of(perm):
    seen = [False] * len(perm)
    out = []
    for i in range(len(perm)):
        if not seen[i]:
            c = []
            j = i
            while not seen[j]:
                seen[j] = True
                c.append(j)
                j = perm[j]
            out.append(c)
    return out


def count_free(perms, npts):
    """perms: 점 0..npts-1 위의 순열들(항등원 제외). 어떤 순열에도 고정되지 않는 함수의 수와 전체 함수의 수를 낸다."""
    total = 1 << npts
    marked = np.zeros(total, dtype=bool)
    for p in perms:
        arr = np.zeros(1, dtype=np.int64)
        for c in cycles_of(p):
            cm = 0
            for j in c:
                cm |= 1 << j
            arr = np.concatenate([arr, arr | cm])
        marked[arr] = True
    return int(total - marked.sum()), total


def run_gl(m):
    """GL(m,2) 가 0 이 아닌 점 2^m - 1 개 위의 함수에 하는 작용."""
    npts = (1 << m) - 1
    perms = []
    order = 0
    for cols in gl_matrices(m):
        order += 1
        p = [apply(cols, v) - 1 for v in range(1, 1 << m)]
        if p != list(range(npts)):
            perms.append(p)
    free, total = count_free(perms, npts)
    return order, npts, free, total


def run_agl(m):
    """AGL(m,2) 가 F_2^m 의 점 2^m 개 위의 함수에 하는 작용."""
    npts = 1 << m
    perms = []
    order = 0
    ident = list(range(npts))
    for cols in gl_matrices(m):
        lin = [apply(cols, v) for v in range(npts)]
        for t in range(npts):
            order += 1
            p = [lin[v] ^ t for v in range(npts)]
            if p != ident:
                perms.append(p)
    free, total = count_free(perms, npts)
    return order, npts, free, total


def main():
    print("=== A. GL(m,2) 가 0 이 아닌 점 2^m - 1 개 위의 함수에 하는 작용 (중심 고정 족의 몫 GL(m-1,2) 에 쓴다) ===")
    gl_free = {}
    for m in range(1, 5):
        t0 = time.time()
        order, npts, free, total = run_gl(m)
        gl_free[m] = free
        print(f"  m={m}: |GL(m,2)| = {order}, 점 {npts}개, 함수 2^{npts} = {total}개 가운데 자유로운 함수 {free}개 ({time.time() - t0:.1f}초)")
    print("  m=5: scripts/agl52_orbits.c -DLINEAR_ONLY 의 전수 열거 (notes/gl52_orbits_output.txt, 점 0 의 값이 자유로워 함수 2^32 개 = 31개 점 위의 함수 2^31 개 x 2):")
    n5 = None
    with open("notes/gl52_orbits_output.txt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("FREE_ORBITS"):
                n5 = int(line.split()[1])
            if line.startswith("SUMMARY"):
                summ = line.strip()
    print(f"    {summ}")
    print(f"    자유 궤도 {n5}개 (각 궤도의 크기는 |GL(5,2)| = 9999360, 0 이 아닌 점 31개 위의 함수로는 {n5 // 2}개 궤도)")
    assert n5 is not None and n5 > 0
    assert gl_free[1] > 0 and gl_free[2] == 0 and gl_free[3] == 0 and gl_free[4] == 0

    print()
    print("=== B. AGL(m,2) 가 F_2^m 의 점 2^m 개 위의 함수에 하는 작용 (해밍-전단 족의 몫 AGL(r,2) 에 쓴다) ===")
    agl_free = {}
    for m in range(1, 5):
        t0 = time.time()
        order, npts, free, total = run_agl(m)
        agl_free[m] = free
        print(f"  m={m}: |AGL(m,2)| = {order}, 점 {npts}개, 함수 2^{npts} = {total}개 가운데 자유로운 함수 {free}개 "
              f"(|AGL| {'>' if order > total else '<='} 2^(2^m): {order > total}) ({time.time() - t0:.1f}초)")
    n5a = None
    with open("notes/agl52_orbits_output.txt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("FREE_ORBITS"):
                n5a = int(line.split()[1])
            if line.startswith("SUMMARY"):
                summa = line.strip()
    print("  m=5: scripts/agl52_orbits.c 의 전수 열거 (notes/agl52_orbits_output.txt):")
    print(f"    {summa}")
    print(f"    자유 궤도 {n5a}개 (|AGL(5,2)| = 319979520 < 2^32 이므로 세기만으로는 가려지지 않았다)")
    assert n5a == 0
    assert agl_free[1] > 0 and agl_free[2] == 0 and agl_free[3] == 0 and agl_free[4] == 0

    print()
    print("정리: 자유로운 함수가 없는 곳 = AGL(m,2) 는 m = 2, 3, 4, 5, GL(m,2) 는 m = 2, 3, 4. GL(5,2) 에는 있다 (자유 궤도 40개).")
    print("결과: 끝")


if __name__ == "__main__":
    sys.path.insert(0, ".")
    main()
