#!/usr/bin/env python3
"""
circuit_native_abc_compare.py
=============================
Yosys 안의 ABC 로 잰 값과, 소스에서 따로 빌드한 ABC 로 같은 스크립트를 돌려 얻은 값이 같은지 열마다 견준다
(notes/2026-11-12_circuit_metrics.md 4절 "도구의 판이 달라도 같은 값이 나오는가").

두 파일은 circuit_metrics.py 의 출력과 같은 모양이다. 줄마다 "묶음 흐름 기준 문자열"이고,
문자열은 함수마다 16진수 두 자리(묶음 b 의 i 번째 문자는 함수 b*4096 + i 의 값)이다.
한쪽 파일에 없는 (묶음, 흐름, 기준)은 견주지 않고 따로 적는다.

직접 빌드한 ABC 쪽 파일은 이렇게 만든다.
    CIRCUIT_ABC=<빌드한 abc 실행 파일의 경로> python scripts/circuit_metrics.py <출력 파일> 1 6

실행: python scripts/circuit_native_abc_compare.py <Yosys 판 파일> <직접 빌드한 판 파일>
      (파일 이름이 .gz 로 끝나면 gzip 으로 푼다)
출력: 열마다 완전히 같은지, 다르면 어긋난 함수의 수와 값 차이의 범위, 합의 차이, 열 합계에 대한 비율.
"""
import gzip
import sys


def load(path):
    opener = gzip.open if path.endswith(".gz") else open
    cols = {}
    with opener(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            batch, flow, crit, s = line.split()
            key = (int(batch), flow, crit)
            if key in cols:
                raise SystemExit("%s: 같은 열이 두 번 있다: %s" % (path, key))
            if len(s) % 2:
                raise SystemExit("%s: 문자열 길이가 홀수이다: %s" % (path, key))
            cols[key] = [int(s[i:i + 2], 16) for i in range(0, len(s), 2)]
    return cols


def main():
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    ref = load(sys.argv[1])
    nat = load(sys.argv[2])
    print("Yosys 판: %s  (열 %d개)" % (sys.argv[1], len(ref)))
    print("직접 빌드한 판: %s  (열 %d개)" % (sys.argv[2], len(nat)))
    print("직접 빌드한 판의 묶음: %s" % sorted({k[0] for k in nat}))
    only_nat = sorted(k for k in nat if k not in ref)
    if only_nat:
        raise SystemExit("Yosys 판에 없는 열이 있다: %s" % only_nat)

    order = {"A": 0, "B": 1}
    same = 0
    diff_keys = []
    n_functions = 0
    for key in sorted(nat, key=lambda k: (k[0], order.get(k[1], 9), k[2])):
        a, b = ref[key], nat[key]
        if len(a) != len(b):
            raise SystemExit("열의 길이가 다르다: %s (%d, %d)" % (key, len(a), len(b)))
        n_functions = len(a)
        d = [y - x for x, y in zip(a, b)]
        n_diff = sum(1 for v in d if v)
        if n_diff == 0:
            same += 1
            continue
        diff_keys.append(key)
        print("다른 열 묶음 %d 흐름 %s 기준 %s: 어긋난 함수 %d개, 값 차이(직접 빌드 - Yosys) %d..%d, "
              "합의 차이 %+d (열 합계 %d의 %.4f%%)"
              % (key[0], key[1], key[2], n_diff, min(d), max(d), sum(d), sum(a), abs(sum(d)) / sum(a) * 100))

    print("열 %d개(묶음마다 함수 %d개) 가운데 완전히 같은 열 %d개, 다른 열 %d개"
          % (len(nat), n_functions, same, len(diff_keys)))
    if diff_keys:
        print("다른 열의 기준: %s" % sorted({k[2] for k in diff_keys}))
        sums = [sum(nat[k]) - sum(ref[k]) for k in diff_keys]
        print("합의 차이의 범위: %+d..%+d" % (min(sums), max(sums)))
        print("어긋난 값의 크기는 모두 %d 이하" % max(abs(v) for k in diff_keys
                                                        for v in (y - x for x, y in zip(ref[k], nat[k]))))


if __name__ == "__main__":
    main()
