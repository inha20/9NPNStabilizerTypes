#!/usr/bin/env python3
"""
export_unexplained_for_gap.py
=============================
analyze_unexplained_classes.py 가 다루는 류(주석표에서 설명됨 = 0)마다, H 의 생성원과 "정규화군 몫 판정법으로는 덮이지 않는 색칠"
(자유 색칠, K = H 의 안정화군 K 와 N(H) 의 교집합이 H 인 색칠)의 꼭짓점 집합 S 를 GAP 파일로 쓴다.
identify_unexplained_overgroups.g 가 읽어서 S 를 보존하는 B_n 의 부분군 K 의 구조를 적는다.

꼭짓점 v (0 .. 2^n - 1) 는 GAP 의 점 v + 1 이다. B_n 의 생성원은 x -> x xor 2^i (i = 0..n-1)과 비트 자리의 순열 둘
(호환 (0 1)과 n-순환)이다. 한 생성원은 점 1..2^n 의 상의 목록으로 쓴다.

사용: python3 -I scripts/export_unexplained_for_gap.py N 표.tsv.gz 주석표.tsv 출력.g

2026-11-15. Claude가 쓴 코드이고 저자 검토 전이다. 2026-11-15의 임시 스크립트를 정리한 것이다.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import np_cegar as C                                    # noqa: E402
from bn_class_io import read_classes                    # noqa: E402
import analyze_unexplained_classes as A                 # noqa: E402


def point_list(vp):
    return "[" + ",".join(str(int(v) + 1) for v in vp) + "]"


def main():
    n, table, annotated, outpath = int(sys.argv[1]), sys.argv[2], sys.argv[3], sys.argv[4]
    cube = C.Cube(n)
    V = cube.V
    assert V <= 256
    X = np.arange(V)
    VA = (cube.act[:, None, :].astype(np.uint8) ^ np.arange(V, dtype=np.uint8)[None, :, None]).reshape(-1, V)
    todo = set(A.read_unexplained(annotated))

    def bitperm(perm):
        vp = np.zeros(V, dtype=np.int64)
        for i in range(n):
            vp += ((X >> i) & 1) << perm[i]
        return vp

    gens_bn = [X ^ (1 << i) for i in range(n)]
    gens_bn.append(bitperm([1, 0] + list(range(2, n))))
    gens_bn.append(bitperm(list(range(1, n)) + [0]))
    with open(outpath, "w", encoding="utf-8") as out:
        out.write("# analyze_unexplained_classes.py 의 류마다 H 의 생성원과 자유 색칠의 꼭짓점 집합 S (export_unexplained_for_gap.py 가 쓴다)\n")
        out.write("NBITS := %d;\nBNGENS := [%s];\nDATA := [];\n" % (n, ",\n".join(point_list(g) for g in gens_bn)))
        for j, order, k, status, gens, _ in read_classes(table, cube, statuses=("absent",)):
            if j not in todo:
                continue
            r = A.analyze(cube, VA, j, order, k, gens)
            lab, _ = cube.orbit_labels(gens)
            free = [x for x in r["rows"] if x["sKN"] == order]
            out.write("Add(DATA, rec(j := %d, order := %d, k := %d, NH := %d,\n gens := [%s],\n free := [\n"
                      % (j, order, k, r["NH"], ",\n".join(point_list(g) for g in gens)))
            parts = []
            for x in free:
                S = [v + 1 for v in range(V) if x["col"][lab[v]] == 1]
                parts.append("rec(col := \"%s\", sK := %d, S := [%s])"
                             % ("".join(str(int(c)) for c in x["col"]), x["sK"], ",".join(map(str, S))))
            out.write(",\n".join(parts) + "\n]));\n")
    print("%d류를 %s 에 썼다" % (len(todo), outpath))


if __name__ == "__main__":
    main()
