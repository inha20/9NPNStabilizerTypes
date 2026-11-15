#!/usr/bin/env python3
"""
analyze_absent_reduction.py
===========================
닫혔지만 안 나타나는 군 H <= B_n 을 "평행이동 부분 T 로 줄이기"로 읽는다 (노트 3절 '줄이기').

  T = H 의 평행이동 부분 (T = H ∩ Z_2^n),  m = n - dim T,
  M = T 의 검사 행렬 (m x n, 열 h_1..h_n 은 F_2^m 의 벡터; T = { x : sum x_i h_i = 0 }),
  Pi_T = 같은 열을 가진 좌표끼리의 순열군 (열 값 v 의 중복도 mu_v 에 대해 곱 prod mu_v!),
  K = H 가 V/T = F_2^m (= 검사값의 공간) 위에 하는 작용의 상 (아핀 변환군).
H 가 닫혀 있으면 H 는 T x| Pi_T 를 품고, |H| = |T| |Pi_T| |K| 이다. 꼭짓점 궤도는 K 의 F_2^m 위 궤도와 같고 k = (K 의 궤도 수) 이다.
이 스크립트는 모든 변수가 쓰이는 (영 열이 없는) 안 나타나는 군마다
  (1) k = K 의 궤도 수,  (2) |H| = |T| |Pi_T| |K|,  (3) 작용이 잘 정의되는가 (M g(x) 가 M x 로만 정해지는가)
를 점검하고, (n, m, |Pi_T|, |K|, k) 별 개수와 'K = 1 이고 Pi_T = 1' (순수 평행이동군) 의 개수를 낸다.

사용: python3 -I scripts/analyze_absent_reduction.py --pipe 작업폴더 --n7 notes/stabilizer_types_n7_closed_absent.tsv
      (작업폴더 = scripts/run_bn_class_pipeline.sh 의 작업 폴더. classes_n2.tsv .. classes_n6.tsv 가 있다)
      시간은 1초 안팎.

2026-11-15. Claude가 쓴 코드이고 저자 검토 전이다.
"""
import argparse
import os
import sys
from collections import Counter
from math import factorial

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

from bn_class_io import read_rows, gens_to_vp
from np_cegar import Cube


def nullspace_gf2(rows, n):
    """정수 비트마스크 rows (F_2^n 의 벡터들) 에 모두 직교하는 벡터들의 기저."""
    # 행 기사다리꼴
    piv = {}
    for r in rows:
        v = r
        while v:
            lead = v.bit_length() - 1
            if lead in piv:
                v ^= piv[lead]
            else:
                piv[lead] = v
                break
    # 기약 기사다리꼴
    leads = sorted(piv)
    for l in leads:
        for l2 in leads:
            if l2 != l and (piv[l2] >> l) & 1:
                piv[l2] ^= piv[l]
    free = [i for i in range(n) if i not in piv]
    basis = []
    for f in free:
        v = 1 << f
        for l in leads:
            if (piv[l] >> f) & 1:
                v |= 1 << l
        basis.append(v)
    return basis


def group_size(perms, npts):
    """점 npts 개 위의 순열(튜플/배열) 생성원들이 만드는 군의 위수 (원소를 모두 나열)."""
    ident = tuple(range(npts))
    seen = {ident}
    frontier = [ident]
    gens = [tuple(int(a) for a in g) for g in perms]
    while frontier:
        new = []
        for a in frontier:
            for g in gens:
                b = tuple(g[x] for x in a)
                if b not in seen:
                    seen.add(b)
                    new.append(b)
        frontier = new
    return len(seen)


def orbit_count(perms, npts):
    parent = list(range(npts))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for g in perms:
        for x in range(npts):
            a, b = find(x), find(int(g[x]))
            if a != b:
                parent[a] = b
    return len({find(x) for x in range(npts)})


def analyse_class(cube, n, order, k, ng, vals):
    gens = gens_to_vp(cube, vals, ng)
    Hset = cube.closure(gens, limit=2_000_000)
    assert len(Hset) == order, (len(Hset), order)
    X = cube.X.astype(np.int64)
    # 평행이동 부분
    T = []
    for vp in Hset.values():
        v = vp.astype(np.int64)
        c = int(v[0])
        if np.array_equal(v, X ^ c):
            T.append(c)
    dimT = len(T).bit_length() - 1
    assert len(T) == 1 << dimT
    # 검사 행렬: T 에 직교하는 벡터 r_1..r_m ; 열 h_i = (r_1[i], ..., r_m[i])
    Tbasis = []
    piv = {}
    for t in T:
        v = t
        while v:
            lead = v.bit_length() - 1
            if lead in piv:
                v ^= piv[lead]
            else:
                piv[lead] = v
                Tbasis.append(t)
                break
    rows = nullspace_gf2(Tbasis, n)
    m = len(rows)
    assert m == n - dimT
    cols = [sum(((rows[j] >> i) & 1) << j for j in range(m)) for i in range(n)]
    essential = all(h != 0 for h in cols)
    mult = Counter(cols)
    pi_order = 1
    for v, mu in mult.items():
        pi_order *= factorial(mu)
    # 검사값 syn[x] = M x
    syn = np.zeros(cube.V, dtype=np.int64)
    for i in range(n):
        syn ^= ((X >> i) & 1) * cols[i]
    rep = np.full(1 << m, -1, dtype=np.int64)
    for x in range(cube.V - 1, -1, -1):
        rep[syn[x]] = x
    assert (rep >= 0).all()          # M 은 위로의 사상 (행이 독립)
    kperms = []
    for g in gens:
        gv = g.astype(np.int64)
        perm = syn[gv[rep]]           # y -> M g(rep(y))
        # 잘 정의되는가: 모든 x 에서 M g(x) = perm[M x]
        assert np.array_equal(syn[gv], perm[syn]), "작용이 잘 정의되지 않음"
        kperms.append(perm)
    Ksize = group_size(kperms, 1 << m)
    korb = orbit_count(kperms, 1 << m)
    return dict(dimT=dimT, m=m, essential=essential, pi_order=pi_order, K=Ksize, korb=korb,
                ok_k=(korb == k), ok_order=(order == (1 << dimT) * pi_order * Ksize),
                cols=cols, kperms=kperms)        # cols, kperms 는 scripts/check_family_in_tables.py 가 쓴다


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pipe", required=True)
    ap.add_argument("--n7", required=True)
    ap.add_argument("--nmax", type=int, default=7)
    args = ap.parse_args()
    bad = 0
    summary = {}
    for n in range(2, args.nmax + 1):
        path = os.path.join(args.pipe, "classes_n%d.tsv" % n) if n <= 6 else args.n7
        cube = Cube(n)
        cnt = Counter()
        ess = 0
        tot = 0
        pure = 0
        for (j, order, k, status, ng, vals, line) in read_rows(path):
            if status not in (None, "absent") or (status is None and False):
                continue
            if n <= 6 and status != "absent":
                continue
            tot += 1
            r = analyse_class(cube, n, order, k, ng, vals)
            if not r["ok_k"] or not r["ok_order"]:
                bad += 1
                print("  어긋남: n=%d j=%d %s" % (n, j, r))
            if not r["essential"]:
                continue
            ess += 1
            key = (r["m"], r["pi_order"], r["K"], k)
            cnt[key] += 1
            if r["K"] == 1 and r["pi_order"] == 1:
                pure += 1
        summary[n] = (tot, ess, pure)
        print("=== n=%d: 안 나타나는 류 %d개, 모든 변수가 쓰이는 류 %d개 (그 가운데 순수 평행이동군 %d개) ===" % (n, tot, ess, pure))
        print("    (m, |Pi_T|, |K|, k) : 개수")
        for key in sorted(cnt):
            print("    %-24s : %d" % (str(key), cnt[key]))
        print("    m 별:", dict(sorted(Counter({m: sum(v for (mm, p, K, kk), v in cnt.items() if mm == m) for m in {mm for (mm, p, K, kk) in cnt}}).items())))
        print("    |K| 별:", dict(sorted(Counter({K: sum(v for (mm, p, KK, kk), v in cnt.items() if KK == K) for K in {KK for (mm, p, KK, kk) in cnt}}).items())))
        print("    |Pi_T| 별:", dict(sorted(Counter({p: sum(v for (mm, pp, K, kk), v in cnt.items() if pp == p) for p in {pp for (mm, pp, K, kk) in cnt}}).items())), flush=True)
    print("점검 (k = K 의 궤도 수, |H| = |T||Pi_T||K|, 작용이 잘 정의됨): 어긋난 류 %d개" % bad)
    print("결과:", "모두 일치" if bad == 0 else "어긋남 있음")


if __name__ == "__main__":
    main()
