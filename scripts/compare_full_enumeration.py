#!/usr/bin/env python3
"""
compare_full_enumeration.py
===========================
전수 열거(enumerate_full_stabilizers.c)가 센 "안정화군의 종류별 함수 수"와
표지 표 계산(stabilizer_types.g)의 표를 견준다.

두 쪽이 군의 종류를 적는 방식이 달라서, 둘 다 같은 여섯 값으로 묶어 합친 뒤 견준다.
    (위수, 순열 부분의 위수, 평행이동 부분의 차원 d, 순열 부분의 변수 위 궤도 길이, 꼭짓점 위 궤도 수, 고정 꼭짓점 수)
켤레류가 다르면서 이 여섯 값이 같은 것들은 한 묶음으로 합쳐진다. 여섯 값은 켤레류에서 정해지므로,
표지 표에 틀린 값이 있으면 합친 수도 어긋나는 것이 보통이다.

실행: python scripts/compare_full_enumeration.py <stabilizer_types_n?.tsv> <열거 출력> [<열거 출력> ...]
      (열거 출력은 ./enumfull? <부분> <부분 수> 의 출력. 여러 부분으로 나눠 돌렸으면 모두 준다.)
"""
import csv
import re
import sys
from collections import defaultdict


def pattern_of(text):
    """'[ 1, 1, 2 ]' -> 112 (작은 것부터 이어 붙인 자릿수)"""
    nums = sorted(int(t) for t in re.findall(r"\d+", text))
    return int("".join(str(t) for t in nums))


def load_table(path):
    n = None
    groups = defaultdict(int)
    orbits = defaultdict(int)
    for line in open(path, encoding="utf-8"):
        if line.startswith("# n ="):
            n = int(re.match(r"# n = (\d+)", line).group(1))
        if not line.strip() or line.startswith("#"):
            continue
        r = line.rstrip("\n").split("\t")
        order, functions, norbits = int(r[0]), int(r[2]), int(r[3])
        hbord, varorb, d, vorbits, fixv = int(r[4]), pattern_of(r[5]), int(r[7]), int(r[8]), int(r[9])
        key = (order, hbord, d, varorb, vorbits, fixv)
        groups[key] += functions
        orbits[key] += norbits
    return n, groups, orbits


def load_enum(paths):
    groups = defaultdict(int)
    nv = None
    for path in paths:
        for line in open(path, encoding="utf-8"):
            if line.startswith("# NV"):
                this = int(line.split()[2])
                if nv not in (None, this):
                    raise SystemExit("변수 수가 다른 파일이 섞였다")
                nv = this
            if not line.startswith("full "):
                continue
            t = line.split()
            order, hbord, d, pat, vorbits, fixv, pairorbits, count = (int(x) for x in t[1:9])
            groups[(order, hbord, d, pat, vorbits, fixv)] += 2 * count   # 최상위 비트가 0 인 절반만 훑었다
    return nv, groups


def main():
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    n, table, table_orbits = load_table(sys.argv[1])
    nv, enum = load_enum(sys.argv[2:])
    if n != nv:
        raise SystemExit("표는 n = %s, 열거는 n = %s 이다" % (n, nv))
    total_t, total_e = sum(table.values()), sum(enum.values())
    print("n = %d: 함수 수의 합: 표지 표 %d, 전수 열거 %d (2^(2^n) = %d)" % (n, total_t, total_e, 2 ** (2 ** n)))
    keys = sorted(set(table) | set(enum))
    bad = 0
    print("%8s %8s %3s %8s %8s %8s | %14s %14s" % ("위수", "순열부분", "d", "변수궤도", "꼭짓점궤도", "고정점", "표지 표", "전수 열거"))
    for k in keys:
        a, b = table.get(k, 0), enum.get(k, 0)
        flag = "" if a == b else "   <-- 다르다"
        if a != b:
            bad += 1
        print("%8d %8d %3d %8d %8d %8d | %14d %14d%s" % (k + (a, b, flag)))
    print("묶음 %d개 가운데 다른 것: %d개" % (len(keys), bad))
    print("모두 일치한다." if bad == 0 and total_t == total_e == 2 ** (2 ** n) else "일치하지 않는다.")
    return 0 if bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
