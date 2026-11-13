#!/usr/bin/env python3
"""
circuit_exact_compare.py
========================
정확 합성으로 구한 NPN 동치류별 최솟값(circuit_exact.py)을 도구가 낸 값(circuit_metrics.py)과 견준다.
기준은 반전이 비용 0인 세 가지(and, l2, l3)이다. 이 기준들에서는 최솟값이 NPN 동치류 안에서 같다.

하는 일
  1. 대표의 최적 회로를 옮겨 쓰는 것이 실제로 되는지 확인한다 (탐색 노트 5.4절).
     동치류의 대표마다 정확 합성이 낸 회로를 하나 두고, 입력선을 바꿔 끼우고 입력과 출력을 반전한 회로를
     게이트 단위로 계산해, 65,536개 함수 각각의 진리표가 그대로 나오는지 본다.
     (반전은 이 세 기준에서 비용이 0이므로, 옮긴 회로의 크기는 대표의 최솟값과 같다.)
  2. 최솟값의 분포(동치류 수와 함수 수)와 65,536개 전체의 합.
  3. 도구의 값이 최솟값보다 얼마나 큰가 (흐름 A, B). 도구의 값이 최솟값보다 작은 함수가 있으면 오류로 멈춘다.

실행: python scripts/circuit_exact_compare.py notes/circuit_metrics.txt.gz notes/circuit_exact_and.txt notes/circuit_exact_l2.txt notes/circuit_exact_l3.txt
      정확 합성 파일은 있는 것만 주면 된다. 하한만 구해진 동치류가 있으면 그 기준은 범위로 적는다.
준비: pip install numpy
"""
import itertools
import os
import re
import sys
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import circuit_metrics_analyze as an

N = 16
FULL = 1 << N
ONES = FULL - 1


def load_exact(path):
    rows = {}
    metric = None
    for line in open(path):
        parts = line.split()
        if len(parts) < 6:
            continue
        metric = parts[2]
        rows[int(parts[0], 16)] = (int(parts[3]), parts[4], parts[6] if len(parts) > 6 else "")   # 마지막 줄이 유효
    return metric, rows


def support(f):
    return [i for i in range(4) if any(((f >> x) & 1) != ((f >> (x ^ (1 << i))) & 1) for x in range(N))]


GATE = re.compile(r"^([A-Z])=([01]+)\(([A-Za-z]+)\)$")


def run_circuit(circuit, words):
    """회로 문자열을 16비트 낱말들(입력 a, b, ... 의 값) 위에서 계산한다. 노드 진리표의 비트 번호는 첫째 입력이 가장 높은 자리."""
    gates = [GATE.match(g).groups() for g in circuit.split(";")]
    val = {chr(ord("a") + i): w for i, w in enumerate(words)}
    for name, tt, fanins in reversed(gates):
        k = len(fanins)
        out = 0
        for idx in range(1 << k):
            if tt[len(tt) - 1 - idx] == "1":
                term = ONES
                for pos, fi in enumerate(fanins):
                    v = val[fi]
                    term &= v if (idx >> (k - 1 - pos)) & 1 else (v ^ ONES)
                out |= term
        val[name] = out
    return val[gates[0][0]]


def verify_transport(rows):
    """대표의 회로에 입력 순열·입력 반전·출력 반전을 붙인 회로가 65,536개 함수를 모두 내는지 게이트 단위로 확인한다.
    돌려주는 값: (회로를 계산해 확인한 함수의 수, 회로가 적히지 않아 계산하지 않은 함수의 수)."""
    maps = [an.point_map(p, a) for p in itertools.permutations(range(4)) for a in range(N)]
    done = bytearray(FULL)
    checked = skipped = 0
    for rep, (value, status, circuit) in sorted(rows.items()):
        sup = support(rep)
        has_circuit = bool(circuit) and GATE.match(circuit.split(";")[0]) is not None
        for t in maps:
            inv = [0] * N
            for x in range(N):
                inv[t[x]] = x
            h = 0
            for x in range(N):
                if (rep >> x) & 1:
                    h |= 1 << t[x]
            for comp in (0, ONES):
                target = h ^ comp
                if done[target]:
                    continue
                done[target] = 1
                if not has_circuit:
                    skipped += 1
                    continue
                # 옮긴 회로: 대표 회로의 j번째 입력에, 점 y 에서 (t 의 역으로 옮긴 점의 sup[j]번 비트)를 넣는다
                words = []
                for i in sup:
                    w = 0
                    for y in range(N):
                        w |= ((inv[y] >> i) & 1) << y
                    words.append(w)
                got = run_circuit(circuit, words) ^ comp
                assert got == target, "옮긴 회로의 진리표가 다르다: 대표 %04x -> %04x" % (rep, target)
                checked += 1
    assert checked + skipped == FULL
    return checked, skipped


def main():
    vals = an.load(sys.argv[1])
    canon = an.canon_npn
    reps, labels = np.unique(canon, return_inverse=True)
    sizes = np.bincount(labels)
    for path in sys.argv[2:]:
        metric, rows = load_exact(path)
        print("=" * 96)
        print("기준 %s  (%s)" % (metric, os.path.basename(path)))
        print("=" * 96)
        missing = [r for r in reps if int(r) not in rows]
        if missing:
            print("  아직 값이 없는 동치류 %d개. 건너뛴다." % len(missing))
            continue
        open_rows = {r: v for r, v in rows.items() if v[1] != "exact"}
        checked, skipped = verify_transport(rows)
        print("  옮겨 쓰기 확인: 함수 %d개에서 대표의 회로를 옮긴 회로를 계산해 진리표가 같음을 보았다 (불일치 없음). "
              "나머지 %d개는 회로가 적히지 않은 경우이다 (게이트가 필요 없는 함수, LUT 하나짜리, 하한만 구해진 동치류)." % (checked, skipped))
        lo = np.array([rows[int(r)][0] for r in reps], dtype=np.int64)            # 최솟값 또는 그 하한
        exact_known = np.array([rows[int(r)][1] == "exact" for r in reps])
        e = lo[labels]
        for fl in an.FLOWS:
            c = vals[(fl, metric)]
            assert (c >= e).all(), "도구의 값이 최솟값(또는 하한)보다 작은 함수가 있다: 흐름 %s" % fl
        if open_rows:
            # 하한만 있는 동치류: 상한은 도구가 낸 가장 작은 값
            best_tool = np.full(len(reps), 10 ** 9, dtype=np.int64)
            for fl in an.FLOWS:
                np.minimum.at(best_tool, labels, vals[(fl, metric)])
            print("  최솟값이 정해지지 않은 동치류 %d개 (대표: 증명된 하한 ~ 도구가 낸 가장 작은 값):" % len(open_rows))
            for i, r in enumerate(reps):
                if not exact_known[i]:
                    print("    %04x: %d ~ %d  (함수 %d개)" % (r, lo[i], best_tool[i], sizes[i]))
            hi = np.where(exact_known, lo, best_tool)
            print("  최솟값이 정해진 동치류 %d개의 분포 (동치류 수): %s" % (
                exact_known.sum(), dict(sorted(Counter(lo[exact_known].tolist()).items()))))
            print("  최솟값이 정해진 동치류 %d개의 분포 (함수 수)  : %s" % (
                exact_known.sum(), dict(sorted(Counter(e[exact_known[labels]].tolist()).items()))))
            print("  65,536개 전체의 최솟값 합: %d 이상 %d 이하" % ((lo * sizes).sum(), (hi * sizes).sum()))
        else:
            print("  최솟값의 분포 (동치류 수): %s" % dict(sorted(Counter(lo.tolist()).items())))
            print("  최솟값의 분포 (함수 수)  : %s" % dict(sorted(Counter(e.tolist()).items())))
            print("  65,536개 전체의 최솟값 합: %d (함수당 평균 %.3f)" % (e.sum(), e.sum() / FULL))
        for fl in an.FLOWS:
            c = vals[(fl, metric)]
            diff = c - e
            cls_min = np.full(len(reps), 10 ** 9, dtype=np.int64)
            cls_max = np.zeros(len(reps), dtype=np.int64)
            np.minimum.at(cls_min, labels, c)
            np.maximum.at(cls_max, labels, c)
            known = exact_known
            print("  흐름 %s: 도구의 합 %d, 최솟값(하한) 합보다 %+.2f%%.  (도구 - 최솟값)의 분포(함수 수): %s" % (
                fl, c.sum(), 100.0 * (c.sum() - e.sum()) / e.sum(), dict(sorted(Counter(diff.tolist()).items()))))
            print("          최솟값이 정해진 동치류 %d개 가운데: 모든 함수에서 최솟값 %d개, 일부 함수에서만 %d개, 어느 함수에서도 아님 %d개" % (
                known.sum(), (known & (cls_max == lo)).sum(), (known & (cls_min == lo) & (cls_max > lo)).sum(),
                (known & (cls_min > lo)).sum()))
        print()


if __name__ == "__main__":
    main()
