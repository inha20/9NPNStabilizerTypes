#!/usr/bin/env python3
"""
compare_enumeration.py
======================
enumerate_all_functions.c 의 전수 열거 결과를 GAP 의 표지 표 계산(picture_types_by_symmetry.g 의 출력)과 비교한다.

사용법: python scripts/compare_enumeration.py notes/picture_types_by_symmetry_output.txt part0.txt [part1.txt ...]

군의 종류는 (위수, 변수 위의 궤도 길이, 쌍 분할 위의 궤도 수)로 맞춘다.
n <= 5 에서는 이 세 값이 S_n 의 부분군 켤레류를 구별한다 (스크립트가 확인한다).
"""
import re
import sys
from collections import defaultdict


def read_parts(paths):
    n = None
    counts = {"hbar": defaultdict(int), "g": defaultdict(int)}
    checked = 0
    covered = 0
    for path in paths:
        for line in open(path, encoding="utf-8"):
            if line.startswith("# NV"):
                w = line.split()
                n = int(w[2])
                covered += int(w[9]) - int(w[8])
                checked += int(w[11])
            elif line.startswith("#") or not line.strip():
                continue
            else:
                kind, order, pattern, orbits, count = line.split()
                counts[kind][(int(order), pattern, int(orbits))] += 2 * int(count)
    return n, counts, covered, checked


def read_gap(path, n):
    rows = []
    active = False
    for line in open(path, encoding="utf-8"):
        m = re.match(r"n = (\d+):", line)
        if m:
            active = int(m.group(1)) == n
            continue
        if active and re.match(r"\s+\d+ \|", line):
            c = [x.strip() for x in line.split("|")]
            pattern = "".join(re.findall(r"\d+", c[3]))
            rows.append((int(c[1]), pattern, int(c[4]), int(c[5]), int(c[7]), c[2]))
    return rows


def main():
    gap_path, parts = sys.argv[1], sys.argv[2:]
    n, counts, covered, checked = read_parts(parts)
    half = 1 << ((1 << n) - 1)
    print("n = %d: 훑은 함수 %d개 (절반 %d개 가운데), 거르지 않고 다시 확인한 함수 %d개" % (n, covered, half, checked))
    ok = covered == half
    rows = read_gap(gap_path, n)
    keys = [(o, p, t) for (o, p, t, _, _, _) in rows]
    if len(set(keys)) != len(keys):
        print("경고: (위수, 궤도 길이, 궤도 수)가 같은 켤레류가 있어 종류를 구별하지 못한다")
        ok = False
    for kind, col, label in (("hbar", 3, "Hbar(f) (입력 반전 허용)"), ("g", 4, "G(f) (변수 순열만)")):
        got = counts[kind]
        total = sum(got.values())
        print("\n[%s] 함수 수 합 %d (2^(2^n) = %d)" % (label, total, 1 << (1 << n)))
        ok = ok and total == 1 << (1 << n)
        bad = 0
        expected_keys = set()
        for row in rows:
            key = (row[0], row[1], row[2])
            expected_keys.add(key)
            mark = "일치" if got.get(key, 0) == row[col] else "불일치"
            bad += mark == "불일치"
            print("  위수 %3d  궤도 %-6s  쌍 분할 궤도 %2d  %-10s  전수 열거 %12d  표지 표 %12d  %s"
                  % (row[0], row[1], row[2], row[5], got.get(key, 0), row[col], mark))
        extra = [k for k in got if k not in expected_keys]
        if extra:
            print("  표지 표에 없는 종류:", extra)
        ok = ok and bad == 0 and not extra
        if kind == "hbar":
            dist = defaultdict(int)
            for (o, p, t), c in got.items():
                dist[t] += c
            print("  그림 종류 수 T 의 분포:", dict(sorted(dist.items())))
    print("\n모두 일치한다." if ok else "\n일치하지 않는 항목이 있다.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
