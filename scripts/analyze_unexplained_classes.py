#!/usr/bin/env python3
"""
analyze_unexplained_classes.py
==============================
정규화군 몫 판정법이 설명하지 못하는 류(analyze_bn_absent_classes.py 가 낸 표에서 설명됨 = 0)가 왜 안 나타나는지 본다.
또 "안 나타남" 판정을 SAT 와 별개로, 색칠 전수로 다시 확인할 수 있다(--verify-all).

바탕이 되는 사실 (증명이 몇 줄이다)
  H 가 닫힌 군이고 꼭짓점 궤도가 O_1..O_k 라고 하자. 궤도의 2-색칠 c 마다 H-불변 함수 f_c 가 하나 있고, H 의 안정화군이
  H 와 같은 함수가 있는가가 곧 "H 가 나타나는가" 이다.
  궤도 집합의 분할 pi 마다, 같은 블록의 궤도에 같은 이름을 붙인 꼭짓점 이름표의 안정화군을 H_pi 라 하자(H <= H_pi).
  (보조정리) c 가 pi 의 모든 블록에서 상수이면 f_c 는 이름표의 함수이므로 H_pi 로 불변이다. 따라서 H_pi != H 이면 f_c 의
  안정화군은 H 보다 크다.
  (역) f_c 의 안정화군에 H 밖의 원소 g 가 있으면, 그림 Gamma_g(궤도 i 와 g 가 i 에서 보내는 궤도를 잇는다)의 연결성분을 블록으로
  하는 분할 pi_g 에서 c 가 상수이고 g 는 H_{pi_g} 에 들어 있다. 그러므로
      c 가 H 보다 큰 안정화군을 갖는다  <=>  H_pi != H 인 분할 pi 가 있고 c 가 pi 의 블록에서 상수이다.
  H_pi != H 인 분할(실현된 분할)의 모임은 더 거친 쪽으로 닫혀 있고, 덮을 때는 가장 고운(극소) 분할만 보면 된다.
  H 가 나타나지 않음 <=> 모든 2-색칠이 어떤 극소 실현 분할의 블록에서 상수이다.
  - 정규화군 몫 판정법은 실현된 분할 가운데 몫 P = N(H)/H 의 원소 p 의 순환 분할(N형)만 쓴 경우이다.
  - 병합 그래프 판정법: 두 궤도만 합친 분할(쌍 병합)이 실현되는 쌍을 변으로 하는 그래프가 이분 그래프가 아니면
    (홀수 순환이 있으면) 모든 색칠이 어떤 변의 양끝에서 같은 색이므로 H 는 나타나지 않는다.

하는 일 (설명 못 하는 류마다)
  1) B_n 의 모든 원소를 훑어서 모든 2-색칠(보색 제외) 의 안정화군 K_c 의 위수와 K_c 와 N(H) 의 교집합의 위수를 구한다.
     교집합이 H 이면(자유 색칠) 정규화군 몫 판정법으로는 덮이지 않는 색칠이다.
  2) 모든 분할 pi 에 대해 H_pi 의 위수를 (np_cegar.Cube.count_stab 으로) 구해 극소 실현 분할을 찾고, N형인지 가른다.
     2-블록 분할은 1)의 값과 같아야 한다(두 방법의 상호 점검).
  3) 쌍 병합 그래프와 그 이분성을 적는다.
  4) N형 분할(몫 P 의 p != 1 의 순환 분할)과 실현된 쌍 병합만 써서 덮이지 않는 색칠이 몇 개인지 센다(0개이면 이 둘로 충분하다).
--pair-graph-all 은 안 나타나는 모든 류에 병합 그래프 판정법을 적용해 정규화군 몫 판정법(주석표의 설명됨 열)과 견준다.
--verify-all 은 안 나타나는 모든 류에 대해 모든 2-색칠의 안정화군 위수를 count_stab 으로 구하고 H 의 위수보다 큰지 본다
  (SAT 를 쓰지 않는 확인; k 가 --kmax 를 넘는 류는 건너뛴다).

사용
  python3 -I scripts/analyze_unexplained_classes.py N 표.tsv.gz 주석표.tsv                       (n=7: 12류, 약 21초)
  python3 -I scripts/analyze_unexplained_classes.py N 표.tsv.gz 주석표.tsv --pair-graph-all      (n=7: 184류, 약 5초)
  python3 -I scripts/analyze_unexplained_classes.py N 표.tsv.gz 주석표.tsv --verify-all [--kmax 16]   (n=7: 184류, 51,152색칠, 약 66초)
  표: merge_bn_class_tables.py 의 공개 표, 주석표: analyze_bn_absent_classes.py 의 출력(설명됨 열이 0인 줄을 분석한다).
  시간은 코어 하나, 메모리는 최대 약 0.5GB (2026-11-15 측정).

2026-11-15. Claude가 쓴 코드이고 저자 검토 전이다.
"""
import argparse
import os
import sys
import time
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import np_cegar as C                                    # noqa: E402
from bn_class_io import read_classes                    # noqa: E402


def set_partitions(k):
    """{0..k-1} 의 모든 집합 분할을 블록 번호 튜플(처음 나온 순서로 0,1,2,...)로 낸다."""
    def rec(i, cur, m):
        if i == k:
            yield tuple(cur)
            return
        for b in range(m + 1):
            cur.append(b)
            yield from rec(i + 1, cur, m + 1 if b == m else m)
            cur.pop()
    yield from rec(0, [], 0)


def canon(lbl):
    m = {}
    out = []
    for x in lbl:
        if x not in m:
            m[x] = len(m)
        out.append(m[x])
    return tuple(out)


def cycle_partition(p):
    k = len(p)
    lbl = [-1] * k
    nb = 0
    for i in range(k):
        if lbl[i] == -1:
            j = i
            while lbl[j] == -1:
                lbl[j] = nb
                j = int(p[j])
            nb += 1
    return canon(lbl)


def refines(a, b):
    """분할 a 가 분할 b 보다 곱거나 같은가(a 의 각 블록이 b 의 한 블록 안에 든다)."""
    seen = {}
    for x, y in zip(a, b):
        if seen.setdefault(x, y) != y:
            return False
    return True


def const_on_blocks(col, pi):
    seen = {}
    for c, b in zip(col, pi):
        if seen.setdefault(b, c) != c:
            return False
    return True


def blocks_all(pi, sizes):
    bl = {}
    for i, b in enumerate(pi):
        bl.setdefault(b, []).append(i)
    return "".join("(" + "+".join(str(sizes[i]) for i in v) + ")" for v in bl.values())


def merged_only(pi, sizes):
    bl = {}
    for i, b in enumerate(pi):
        bl.setdefault(b, []).append(i)
    out = ["+".join(str(sizes[i]) for i in v) for v in bl.values() if len(v) > 1]
    return " | ".join(out)


def bipartite(k, edges):
    adj = {i: [] for i in range(k)}
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
    color = {}
    for s in range(k):
        if s in color:
            continue
        color[s] = 0
        stack = [s]
        while stack:
            u = stack.pop()
            for v in adj[u]:
                if v not in color:
                    color[v] = 1 - color[u]
                    stack.append(v)
                elif color[v] == color[u]:
                    return False
    return True


def coloring(k, mask):
    return np.array([0] + [(mask >> (i - 1)) & 1 for i in range(1, k)], dtype=np.uint8)


def read_annotated(path):
    """주석표에서 번호 -> 설명됨(1/0) 사전을 읽는다."""
    out = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            out[int(p[0])] = int(p[7])
    return out


def read_unexplained(path):
    return [j for j, e in read_annotated(path).items() if e == 0]


def analyze(cube, VA, j, order, k, gens):
    V = cube.V
    lab, k2 = cube.orbit_labels(gens)
    assert k2 == k
    labu = lab.astype(np.uint8)
    sizes = np.bincount(lab).tolist()
    first = np.array([np.nonzero(lab == o)[0][0] for o in range(k)])
    rep = first[lab]
    LG = labu[VA]
    inH = (LG == labu[None, :]).all(axis=1)              # 모든 궤도를 제자리에 두는 원소 = H^ (닫혔으므로 H)
    assert int(inH.sum()) == order, "H 가 닫혀 있지 않다?"
    norm = (LG == LG[:, rep]).all(axis=1)                # N(H)
    idx = np.nonzero(norm)[0]
    NH = len(idx) // order
    perms = np.unique(LG[idx][:, first], axis=0)
    assert len(perms) == NH
    ident = np.arange(k)
    cyc = {cycle_partition(p) for p in perms if not (p == ident).all()}

    # 1) 모든 2-색칠(보색 제외)의 안정화군
    rows = []
    for mask in range(1 << (k - 1)):
        col = coloring(k, mask)
        f = col[lab]
        stab = (f[VA] == f[None, :]).all(axis=1)
        sK = int(stab.sum())
        sKN = int((stab & norm).sum())
        assert not (inH & ~stab).any()
        if mask == 0:
            blocks = None
        else:
            code = (lab[None, :].astype(np.uint16) * k + LG[stab].astype(np.uint16))
            parent = list(range(k))

            def find(a):
                while parent[a] != a:
                    parent[a] = parent[parent[a]]
                    a = parent[a]
                return a
            for cd in np.unique(code):
                a, b = divmod(int(cd), k)
                ra, rb = find(a), find(b)
                if ra != rb:
                    parent[ra] = rb
            blocks = canon([find(i) for i in range(k)])
        rows.append(dict(mask=mask, col=col, sK=sK, sKN=sKN, blocks=blocks))

    # 2) 모든 분할의 H_pi 위수 (count_stab)
    realized = {}
    for pi in set_partitions(k):
        realized[pi] = cube.count_stab(np.array(pi, dtype=np.uint16)[lab])
    assert realized[tuple(range(k))] == order
    for r in rows:                                       # 두 방법의 상호 점검
        assert realized[canon(r["col"])] == r["sK"], "count_stab 과 전원소 훑기가 다르다"
    R = [pi for pi, v in realized.items() if v > order]
    minimal = [pi for pi in R if not any(q != pi and refines(q, pi) for q in R)]
    covered = all(any(const_on_blocks(r["col"], pi) for pi in minimal) for r in rows)
    assert covered == all(r["sK"] > order for r in rows)

    # 3) 쌍 병합 그래프
    edges, edges_n, pair_parts = [], [], set()
    for a in range(k):
        for b in range(a + 1, k):
            pi = canon([a if m == b else m for m in range(k)])
            if realized[pi] > order:
                edges.append((a, b))
                pair_parts.add(pi)
                if pi in cyc:
                    edges_n.append((a, b))
    # 4) N형 분할(모든 p != 1 의 순환 분할)과 실현된 쌍 병합만 쓰면 모든 색칠이 덮이는가, 안 덮이는 색칠은 몇 개인가
    simple = set(cyc) | pair_parts
    n_uncovered = sum(1 for r in rows if not any(const_on_blocks(r["col"], pi) for pi in simple))
    return dict(j=j, order=order, k=k, sizes=sizes, NH=NH, rows=rows, minimal=minimal, cyc=cyc, edges=edges, edges_n=edges_n,
                bip=bipartite(k, edges), realized=realized, n_uncovered=n_uncovered)


def report(r):
    sizes = r["sizes"]
    k, order = r["k"], r["order"]
    free = [x for x in r["rows"] if x["sKN"] == order]
    print("j=%d |H|=%d k=%d 궤도 %s |N(H)/H|=%d" % (r["j"], order, k, "+".join(str(s) for s in sorted(sizes, reverse=True)), r["NH"]))
    print("  색칠 %d개(보색 제외) 가운데 N(H)로 덮이지 않는(자유) 색칠 %d개; 모든 색칠의 안정화군이 H 보다 큼: %s"
          % (len(r["rows"]), len(free), all(x["sK"] > order for x in r["rows"])))
    for x in free:
        print("   자유 색칠 %s: |K|/|H|=%d, K 의 궤도 %s, K∩N(H)=H 이므로 H 는 K 안에서 자기 정규화군"
              % ("".join(str(int(c)) for c in x["col"]), x["sK"] // order,
                 "전체(상수 함수)" if x["blocks"] is None else blocks_all(x["blocks"], sizes)))
    mins = []
    for pi in sorted(r["minimal"], key=lambda t: (max(t), t)):
        mins.append("%s%s" % (merged_only(pi, sizes), "[N형]" if pi in r["cyc"] else "[비N형]"))
    print("  극소 실현 분할 %d개: %s" % (len(r["minimal"]), "; ".join(mins)))
    print("  쌍 병합 그래프: 변 %d개(N형 %d개), 이분 그래프 %s; 변 %s"
          % (len(r["edges"]), len(r["edges_n"]), "예" if r["bip"] else "아니오 (홀수 순환 있음 -> 병합 그래프 판정법으로 안 나타남)",
             ", ".join("%d-%d" % (sizes[a], sizes[b]) for a, b in r["edges"])))
    print("  N형 분할과 쌍 병합만으로 덮이지 않는 색칠: %d개 (보색 제외 %d개 가운데)" % (r["n_uncovered"], len(r["rows"])))


def verify_all(cube, table, kmax):
    t0 = time.time()
    ncls = ncol = 0
    skipped = []
    bad = []
    ratios = []
    for j, order, k, status, gens, _ in read_classes(table, cube, statuses=("absent",)):
        if k > kmax:
            skipped.append((j, k))
            continue
        lab, k2 = cube.orbit_labels(gens)
        assert k2 == k
        assert cube.count_stab(lab.astype(np.uint16)) == order, "닫혀 있지 않다"
        mn = None
        for mask in range(1 << (k - 1)):
            cnt = cube.count_stab(coloring(k, mask).astype(np.uint16)[lab])
            ncol += 1
            if cnt <= order:
                bad.append((j, mask, cnt))
            mn = cnt if mn is None else min(mn, cnt)
        ratios.append(mn / order)
        ncls += 1
    print("전수 확인: 안 나타나는 류 %d개, 색칠(보색 제외) %d개 (%.0f초)" % (ncls, ncol, time.time() - t0))
    print("  H 와 안정화군이 같은 색칠이 있는 류: %d개 %s" % (len(bad), bad[:5]))
    print("  류마다 (색칠의 안정화군 위수의 최솟값)/|H| 의 최솟값 %.3g, 최댓값 %.3g" % (min(ratios), max(ratios)))
    if skipped:
        print("  건너뛴 류(k > %d): %s" % (kmax, skipped))
    return not bad


def pair_graph_all(cube, table, annotated):
    """안 나타나는 모든 류에 병합 그래프 판정법(쌍 병합 그래프가 이분 그래프가 아님)을 적용하고
    정규화군 몫 판정법(주석표의 설명됨 열)과 견준다."""
    expl = read_annotated(annotated)
    t0 = time.time()
    tab = Counter()
    left = []
    for j, order, k, status, gens, _ in read_classes(table, cube, statuses=("absent",)):
        lab, k2 = cube.orbit_labels(gens)
        assert k2 == k
        edges = []
        for a in range(k):
            for b in range(a + 1, k):
                lp = np.array(canon([a if m == b else m for m in range(k)]), dtype=np.uint16)[lab]
                if cube.count_stab(lp) > order:
                    edges.append((a, b))
        nonbip = not bipartite(k, edges)
        tab[(expl[j], nonbip)] += 1
        if expl[j] == 0 and not nonbip:
            left.append((j, k))
    print("병합 그래프 판정법을 안 나타나는 류 %d개에 적용 (%.0f초)" % (sum(tab.values()), time.time() - t0))
    print("  정규화군 몫 판정법으로 설명됨 %d개: 병합 그래프 판정법으로도 설명 %d개, 아닌 것 %d개"
          % (tab[(1, True)] + tab[(1, False)], tab[(1, True)], tab[(1, False)]))
    print("  정규화군 몫 판정법으로 설명 못 함 %d개: 병합 그래프 판정법으로 설명 %d개, 아닌 것 %d개 %s"
          % (tab[(0, True)] + tab[(0, False)], tab[(0, True)], tab[(0, False)], left))
    print("  병합 그래프 판정법만으로 설명 %d개; 두 판정법 중 하나로 설명 %d개; 둘 다 못 하는 류 %d개"
          % (tab[(1, True)] + tab[(0, True)], tab[(1, True)] + tab[(1, False)] + tab[(0, True)], tab[(0, False)]))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("n", type=int)
    ap.add_argument("table")
    ap.add_argument("annotated")
    ap.add_argument("--verify-all", action="store_true")
    ap.add_argument("--pair-graph-all", action="store_true")
    ap.add_argument("--kmax", type=int, default=12)
    a = ap.parse_args()
    cube = C.Cube(a.n)
    V = cube.V
    assert V <= 256
    if a.verify_all:
        sys.exit(0 if verify_all(cube, a.table, a.kmax) else 1)
    if a.pair_graph_all:
        pair_graph_all(cube, a.table, a.annotated)
        return
    todo = set(read_unexplained(a.annotated))
    VA = (cube.act[:, None, :].astype(np.uint8) ^ np.arange(V, dtype=np.uint8)[None, :, None]).reshape(-1, V)
    t0 = time.time()
    res = []
    for j, order, k, status, gens, _ in read_classes(a.table, cube, statuses=("absent",)):
        if j in todo:
            r = analyze(cube, VA, j, order, k, gens)
            report(r)
            res.append(r)
    print("정규화군 몫 판정법이 설명 못 하는 류 %d개를 분석 (%.0f초)" % (len(res), time.time() - t0))
    nbip = sum(1 for r in res if not r["bip"])
    nonn = sum(1 for r in res if any(pi not in r["cyc"] for pi in r["minimal"]))
    allfree = sum(1 for r in res if any(x["sKN"] == r["order"] for x in r["rows"]))
    print("  자유 색칠이 있는 류 %d개(설명 못 하는 류와 같아야 한다)" % allfree)
    print("  극소 실현 분할에 비N형이 있는 류 %d개" % nonn)
    print("  쌍 병합 그래프가 이분 그래프가 아닌(홀수 순환이 있는) 류 %d개" % nbip)
    nsimple = sum(1 for r in res if r["n_uncovered"] > 0)
    print("  N형 분할과 쌍 병합을 함께 써도 덮이지 않는 색칠이 있는 류 %d개 %s" % (nsimple, [r["j"] for r in res if r["n_uncovered"] > 0]))
    print("  k 분포 %s" % sorted(Counter(r["k"] for r in res).items()))


if __name__ == "__main__":
    main()
