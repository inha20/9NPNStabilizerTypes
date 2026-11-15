#!/usr/bin/env python3
"""
analyze_bn_absent_classes.py
============================
닫혀 있지만 안 나타나는 B_n 의 켤레류(표에서 상태가 absent 인 줄)마다 "정규화군 몫 판정법"을 적용하고 구조를 적은 표를 만든다.

정규화군 몫 판정법 (안 나타남의 충분조건이고 증명이 한 줄이다)
  H 가 닫힌 군이고 꼭짓점 궤도가 k 개라고 하자. 정규화군 N(H) 는 H 의 궤도들을 서로 바꾸고, H 가 닫혀 있으므로
  몫 P = N(H)/H 는 궤도 집합(k 개의 점) 위에 충실하게 작용한다(궤도를 모두 제자리에 두는 원소는 H^ = H 에 속한다).
  궤도의 2-색칠 c 에 대응하는 함수 f_c (궤도마다 값 하나)는 H 로 불변이다. P 의 단위원이 아닌 원소 p 가 c 를 보존하면
  (c∘p = c) p 의 대표원 g ∈ N(H) \\ H 가 f_c 를 보존하므로 f_c 의 대칭군은 H 보다 크다.
  따라서 P 가 2^k 개의 색칠 가운데 "단위원 밖의 어떤 원소에도 보존되지 않는" 색칠(자유 궤도)을 하나도 못 가지면 H 는 나타나지 않는다.
  거꾸로는 성립하지 않는다(f_c 의 대칭군이 N(H) 밖의 원소로 커질 수 있다). 자유 궤도가 있는지는 SAT 로 정한다.
  표의 "설명됨 = 1" 은 자유 궤도가 없어서 이 판정법만으로 안 나타남이 증명되는 류이다.

  구현: H 가 닫혀 있으면 N(H) 는 H 의 궤도 분할을 보존하는 B_n 의 원소 전체와 같다(그런 g 는 H^ = H 를 켤레로 보존한다).
  그래서 B_n 의 모든 원소(n! 2^n 개)를 한 번에 검사한다. 몫이 궤도 위에 충실하게 작용한다는 것(|N/H| = |P|)은 매번 확인한다.

쓰는 것: 표 파일(번호 = 입력 표의 번호). 열: 번호, 위수, k, 궤도 길이(개수x길이), 평행이동 부분의 위수 |C|, 순열 부분의 위수 |Sigma|,
        |N(H)/H|, 설명됨(1/0). |C| 는 H 속 평행이동(x -> x xor c)의 수이고 |Sigma| = |H| / |C| 이다.
화면에는 개수와 k 분포를 쓴다.

사용
  python3 -I scripts/analyze_bn_absent_classes.py N 표 [출력표.tsv]
  표: merge_bn_class_tables.py 가 낸 공개 표(.tsv 나 .tsv.gz). 안 나타남 줄만 분석하고, 표 전체의 개수도 같이 쓴다.
  n = 7 의 184류에 약 45초(코어 하나, 메모리 최대 약 0.4GB; 2026-11-14 측정).

2026-11-14. Claude가 쓴 코드이고 저자 검토 전이다. 2026-11-13의 임시 스크립트
(normalizer_criterion.py, absent_classes.py, analyze_absent_lib.py)를 하나로 정리하고 번호를 입력 표의 번호로 바꿨다.
"""
import argparse
import os
import sys
import time
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import np_cegar as C                                    # noqa: E402
from bn_class_io import read_classes, read_rows, summarize, describe   # noqa: E402


def free_coloring_exists(k, perms):
    """perms: 궤도 집합 위의 순열들(행). 단위원이 아닌 모든 순열에 보존되지 않는 2-색칠이 있는가?
    절: 각 순열 p 마다 "c[i] != c[p[i]] 인 i 가 있다". 변수 e_ab 는 c[a] xor c[b]."""
    from pysat.solvers import Cadical153
    ident = np.arange(k)
    cl = set()
    for p in perms:
        if (p == ident).all():
            continue
        e = np.nonzero(p != ident)[0]
        cl.add(tuple(sorted(set((min(int(i), int(p[i])), max(int(i), int(p[i]))) for i in e))))
    mins = []
    for c in sorted(cl, key=len):                       # 극소 절만 남긴다
        s = set(c)
        if not any(set(m) <= s for m in mins):
            mins.append(c)
    solver = Cadical153()
    ev, nv = {}, k
    for c in mins:
        lits = []
        for e in c:
            if e not in ev:
                nv += 1
                ev[e] = nv
                a, b = e
                solver.add_clause([-nv, a + 1, b + 1])
                solver.add_clause([-nv, -(a + 1), -(b + 1)])
            lits.append(ev[e])
        solver.add_clause(lits)
    ok = solver.solve()
    solver.delete()
    return ok


def comp(sizes):
    c = Counter(sizes)
    return "+".join("%dx%d" % (v, s) for s, v in sorted(c.items(), reverse=True))


def analyze(cube, VA, gens, order):
    V, X = cube.V, cube.X
    lab, k = cube.orbit_labels(gens)
    labu = lab.astype(np.uint8)
    first = np.array([np.nonzero(lab == o)[0][0] for o in range(k)])
    rep = first[lab]
    LG = labu[VA]                                        # (원소 수, V): 각 원소가 보내는 라벨
    norm = (LG == LG[:, rep]).all(axis=1)                # 궤도를 궤도로 보내는 원소 = N(H) (H 가 닫혔으므로)
    idx = np.nonzero(norm)[0]
    assert len(idx) % order == 0
    perms = np.unique(LG[idx][:, first], axis=0)         # 궤도 집합 위에 유도된 순열
    assert len(perms) == len(idx) // order, "몫이 궤도 위에 충실하게 작용하지 않는다(H 가 닫혀 있지 않다?)"
    ntrans = 0
    for b in cube.closure(gens):
        vp = np.frombuffer(b, dtype=cube.dt).astype(np.int64)
        if (vp == (X ^ vp[0])).all():
            ntrans += 1
    sizes = tuple(sorted(np.bincount(lab).tolist(), reverse=True))
    free = free_coloring_exists(k, perms)
    fixed = np.ones(V, dtype=bool)                       # 모든 생성원이 제자리에 두는 꼭짓점
    for g in gens:
        fixed &= (g.astype(np.int64) == X)
    typ = "U" if fixed.any() else ("T1" if ntrans > 1 else "T0")   # U: 고정 꼭짓점이 있다, T1: 0이 아닌 평행이동이 있다, T0: 둘 다 아니다
    return dict(order=order, k=k, sizes=sizes, ntrans=ntrans, sigma=order // ntrans, NH=len(idx) // order, explained=not free, typ=typ, perms=perms)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("n", type=int)
    ap.add_argument("table")
    ap.add_argument("out", nargs="?", default=None)
    ap.add_argument("--kmin", type=int, default=0, help="꼭짓점 궤도가 kmin 개 이상인 류만 분석한다(빠르게 큰 k 만 볼 때)")
    ap.add_argument("--gap-out", default=None, help="분석한 류마다 몫 P 가 궤도 집합에 일으키는 순열군을 GAP 파일로 쓴다(identify_quotient_groups.g 가 읽는다)")
    a = ap.parse_args()
    cube = C.Cube(a.n)
    V = cube.V
    assert V <= 256
    VA = (cube.act[:, None, :].astype(np.uint8) ^ np.arange(V, dtype=np.uint8)[None, :, None]).reshape(-1, V)
    allrows = [(r[0], r[1], r[2], r[3]) for r in read_rows(a.table)]
    print("n=%d 표 전체: %s" % (a.n, describe(summarize(allrows)).replace("\n", "\n     ")))
    t0 = time.time()
    rows = []
    gap = open(a.gap_out, "w", encoding="utf-8") if a.gap_out else None
    if gap:
        gap.write("# n=%d: 몫 P = N(H)/H 가 H 의 꼭짓점 궤도(점 1..k) 위에 일으키는 순열군. analyze_bn_absent_classes.py --gap-out 이 쓴다.\nQUOTIENTS := [];\n" % a.n)
    for j, order, k, status, gens, _ in read_classes(a.table, cube, statuses=("absent",)):
        if k < a.kmin:
            continue
        r = analyze(cube, VA, gens, order)
        assert r["k"] == k
        r["j"] = j
        if gap:
            gap.write("Add(QUOTIENTS, rec(j := %d, order := %d, k := %d, NH := %d, gens := [\n" % (j, order, k, r["NH"]))
            gap.write(",\n".join("[" + ",".join(str(int(x) + 1) for x in p) + "]" for p in r["perms"]))
            gap.write("\n]));\n")
        r.pop("perms")
        rows.append(r)
    if gap:
        gap.close()
    rows.sort(key=lambda r: (-r["k"], -r["order"], r["NH"], r["j"]))
    expl = [r for r in rows if r["explained"]]
    unexpl = [r for r in rows if not r["explained"]]
    print("안 나타나는 류 %d개 분석 (%.0f초)" % (len(rows), time.time() - t0))
    print("  정규화군 몫 판정법으로 설명되는 류 %d개, 설명 못 하는 류 %d개" % (len(expl), len(unexpl)))
    print("  설명되는 류의 k 분포 %s" % sorted(Counter(r["k"] for r in expl).items()))
    print("  설명 못 하는 류의 k 분포 %s" % sorted(Counter(r["k"] for r in unexpl).items()))
    tc = Counter(r["typ"] for r in rows)
    print("  유형별 개수 U(고정 꼭짓점 있음) %d, T0(꼬임) %d, T1(평행이동 있음) %d" % (tc["U"], tc["T0"], tc["T1"]))
    for r in rows[:6]:
        print("  j=%d 위수 %d k=%d 궤도 %s |C|=%d |Sigma|=%d |N(H)/H|=%d 설명됨=%d"
              % (r["j"], r["order"], r["k"], comp(r["sizes"]), r["ntrans"], r["sigma"], r["NH"], r["explained"]))
    if a.out:
        with open(a.out, "w", encoding="utf-8") as f:
            f.write("# n=%d 닫혔지만 안 나타나는 군의 켤레류 %d개 (번호는 %s 의 번호와 같다)\n" % (a.n, len(rows), os.path.basename(a.table)))
            f.write("# 번호\t위수\t꼭짓점 궤도 수\t궤도 길이(개수x길이)\t평행이동 부분 위수 |C|\t순열 부분 위수 |Sigma|\t|N(H)/H|\t정규화군 몫만으로 설명됨(1/0)\t유형(U/T0/T1)\n")
            f.write("# 설명됨 = 1 : 몫 N(H)/H 가 궤도의 2-색칠 위에서 자유 궤도를 갖지 않아, 모든 H-불변 함수의 대칭군이 H 보다 크다는 것이 증명된다 (scripts/analyze_bn_absent_classes.py 의 머리말)\n")
            f.write("# 유형: U = 모든 원소가 함께 고정하는 꼭짓점이 있다, T1 = 0이 아닌 평행이동이 있다(|C| > 1), T0 = 둘 다 아니다 (notes/2026-11-12_stabilizer_types.md 1절)\n")
            for r in rows:
                f.write("%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%s\n" % (r["j"], r["order"], r["k"], comp(r["sizes"]), r["ntrans"], r["sigma"], r["NH"], 1 if r["explained"] else 0, r["typ"]))


if __name__ == "__main__":
    main()
