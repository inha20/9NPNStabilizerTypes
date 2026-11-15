#!/usr/bin/env python3
"""
check_agl_free_witnesses_nauty.py
=================================
scripts/pure_translation_absent.py 가 m = 6, 7 에서 찾은 "AGL(m,2) 에서 안정화군이 자명한 함수"를 nauty 로 따로 확인한다.
(GAP 의 Stabilizer(AGL(m,2), 집합, OnSets) 는 m = 7 부터 너무 느려서 GAP 으로는 m = 6 까지만 했다: scripts/check_agl_free_witnesses.g.
 m = 8 은 평면이 690,880 개라서 nauty 의 밀집 그래프 표현이 메모리를 넘어 이 방법으로 하지 못한다. m = 8 은 파이썬의 정확한 되짚기
 (m <= 7 에서 nauty 와 같은 값을 낸 것으로 시험한 바로 그 함수)로만 확인한 것이다.)

방법: m >= 3 일 때 F_2^m 의 아핀 평면(점 4개 {x, x+u, x+w, x+u+w})의 모임은 점집합 위의 Steiner 4중계(SQS(2^m))이고
그 자기동형군은 AGL(m,2) 이다. 그래서 함수 g 의 AGL(m,2) 안 안정화군의 위수 = 다음 그래프의 색칠을 보존하는 자기동형군의 위수이다.
  꼭짓점 = 점 2^m 개 + 아핀 평면마다 하나, 간선 = 평면과 그 네 점, 색 = (g = 0 인 점 | g = 1 인 점 | 평면).
이 방법을 먼저 m = 3, 4, 5 에서 많은 함수(무작위 + 특별한 함수)로 scripts/pure_translation_absent.py 의 되짚기와 견주어 시험하고
(위수가 모두 같아야 한다), 그 뒤 m = 6, 7, 8 의 함수에 쓴다.

실행: python3 -u -I scripts/check_agl_free_witnesses_nauty.py   (pynauty 필요: pip install pynauty; 약 2분)

2026-11-15. Claude가 쓴 코드이고 저자 검토 전이다.
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pure_translation_absent import stab_elements  # noqa: E402

import pynauty  # noqa: E402


def planes(m):
    N = 1 << m
    out = []
    for x in range(N):
        for y in range(x + 1, N):
            for z in range(y + 1, N):
                t = x ^ y ^ z
                if t > z:
                    out.append((x, y, z, t))
    return out


_PLANES = {}


def nauty_order(gt, m):
    """g 의 AGL(m,2) 안 안정화군의 위수 (그래프의 색칠 보존 자기동형군의 위수)."""
    N = 1 << m
    if m not in _PLANES:
        _PLANES[m] = planes(m)
    P = _PLANES[m]
    adj = {}
    for i, q in enumerate(P):
        v = N + i
        adj[v] = list(q)
    g = pynauty.Graph(N + len(P), directed=False, adjacency_dict=adj,
                      vertex_coloring=[set(int(x) for x in range(N) if gt[x] == 0),
                                       set(int(x) for x in range(N) if gt[x] == 1),
                                       set(range(N, N + len(P)))])
    gens, g1, g2, orbits, norb = pynauty.autgrp(g)
    return int(round(g1 * 10 ** g2))


def special_functions(m):
    N = 1 << m
    out = {}
    out["g=0"] = np.zeros(N, dtype=np.int64)
    pt = np.zeros(N, dtype=np.int64)
    pt[0] = 1
    out["한 점의 지시함수"] = pt
    sub = np.array([1 if (x >> 2) == 0 else 0 for x in range(N)], dtype=np.int64)
    out["부분공간(앞 2 좌표)의 지시함수"] = sub
    out["이차식 x1x2+x2x3+..."] = np.array([sum(((x >> i) & 1) & ((x >> (i + 1)) & 1) for i in range(m - 1)) & 1 for x in range(N)], dtype=np.int64)
    out["선형함수 x1"] = np.array([x & 1 for x in range(N)], dtype=np.int64)
    out["무게 1 인 점들의 지시함수"] = np.array([1 if bin(x).count("1") == 1 else 0 for x in range(N)], dtype=np.int64)
    return out


def main():
    print("=== 시험: m = 3, 4, 5 에서 nauty 의 위수와 pure_translation_absent.stab_elements 의 위수 견주기 ===", flush=True)
    bad = 0
    tot = 0
    for m in (3, 4, 5):
        t0 = time.time()
        N = 1 << m
        # m = 5 에서는 대칭이 큰 함수의 안정화군(수억 개의 원소)을 되짚기로 나열할 수 없으므로 무작위 함수만 쓴다
        funcs = list(special_functions(m).items()) if m <= 4 else []
        for seed in range(40):
            funcs.append(("무작위 씨앗 %d" % seed, np.random.default_rng(7 * 10 ** 6 + seed).integers(0, 2, size=N)))
        # 안정화군이 큰 함수도 섞는다: 무작위 함수를 몇 개의 좌표만 쓰게 접는다
        for seed in range(10 if m <= 4 else 0):
            r = np.random.default_rng(8 * 10 ** 6 + seed)
            small = r.integers(0, 2, size=1 << (m - 2))
            funcs.append(("앞 %d 좌표만 쓰는 무작위 %d" % (m - 2, seed), np.array([small[x & ((1 << (m - 2)) - 1)] for x in range(N)], dtype=np.int64)))
        orders = {}
        for name, gt in funcs:
            a = nauty_order(gt, m)
            b = len(stab_elements(gt, m))
            tot += 1
            orders[b] = orders.get(b, 0) + 1
            if a != b:
                bad += 1
                print("  불일치 m=%d %s: nauty %d, 되짚기 %d" % (m, name, a, b))
        print("  m=%d: 함수 %d개 시험, 위수 분포 %s (%.1f초)" % (m, len(funcs), dict(sorted(orders.items())), time.time() - t0), flush=True)
    print("  불일치 %d개 / %d개" % (bad, tot))
    assert bad == 0

    print()
    print("=== m = 6, 7: pure_translation_absent.py 가 찾은 함수 (씨앗 규칙 3*10^6 + 씨앗) 의 안정화군 위수 ===", flush=True)
    for m in (6, 7):
        t0 = time.time()
        seed = 0
        while True:
            gt = np.random.default_rng(3 * 10 ** 6 + seed).integers(0, 2, size=1 << m)
            if len(stab_elements(gt, m)) == 1:
                break
            seed += 1
        o = nauty_order(gt, m)
        print("  m=%d: 씨앗 %d, 값이 1 인 점 %d개, nauty 의 안정화군 위수 = %d (%.1f초)" % (m, seed, int(gt.sum()), o, time.time() - t0), flush=True)
        assert o == 1
    print("(m = 8 은 평면 690,880 개의 그래프가 nauty 의 밀집 표현으로 메모리를 넘어 이 스크립트로 하지 못한다.)")
    print("결과: 모두 일치")


if __name__ == "__main__":
    main()
