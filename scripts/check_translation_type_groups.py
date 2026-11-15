#!/usr/bin/env python3
"""
check_translation_type_groups.py
================================
K = 1 인 닫힌 군 ("평행이동형 군" H = T x| Pi_T) 의 나타남 판정을 공개 표와 B_8 의 SAT 판정에 견준다.

바탕이 되는 사실 (증명은 안정화군 노트 3절 "일반 n (2)")
  닫힌 군 H 의 평행이동 부분을 T = ker M (M: F_2^n -> F_2^m 위로의 사상, 열 h_1..h_n) 이라 하고, H 가 잉여류의 공간 F_2^m 위에 하는
  작용의 상을 K 라 하자. K = 1 이면 H = T x| Pi_T 이다 (Pi_T = 같은 열을 가진 좌표끼리의 순열군). 이것은 T 의 모든 잉여류를 제자리에
  두는 원소 전부이므로 닫혀 있고, 거꾸로 열이 주어진 T 마다 이 군은 K = 1 인 닫힌 군이다.
  모든 변수가 쓰이면 열이 0 이 아니다. 열의 중복도 mu(v) = #{i : h_i = v} (v != 0, mu 의 지지가 F_2^m 을 생성) 의 GL(m,2)-류가
  H 의 켤레류를 정한다. A_mu = {A in GL(m,2) : mu(Av) = mu(v)},  P_mu = F_2^m x| A_mu,  N(H)/H = P_mu,  꼭짓점 궤도 수는 2^m,
  |H| = 2^(n-m) prod_v mu(v)!.
  H 가 나타남 <=> P_mu 에 자유 함수 (안정화군이 자명한 g: F_2^m -> F_2) 가 있다.

하는 일
  A. m = 1..4, n = 2..NMAX 에서 모든 변수가 쓰이는 mu 의 GL(m,2)-류를 나열하고, P_mu 의 자유 함수를 2^m 개 점 위의 함수 전부로 가려
     안 나타나는 류를 낸다 (군 연산으로 세는 것이고 SAT 를 쓰지 않는다).
  B. 공개 표의 닫힌 류 가운데 K = 1 이고 모든 변수가 쓰이는 것을 scripts/analyze_absent_reduction.py 의 analyse_class 로 골라
     (n, m, |Pi_T|, 상태) 별 개수를 A 의 예측과 견준다. n = 2..6 은 나타남과 안 나타남 모두, n = 7 은 안 나타남만이다.
     m >= 5 인 K = 1 류는 모두 나타나야 한다 (노트 3절: m = 5 는 mu 가 상수일 때만, 곧 n = 31 c 에서만 안 나타나고, m >= 6 은 없다).
  C. n = 8 에서 m <= 4 인 mu 의 류마다 H = T x| Pi_T 를 B_8 안에서 만들어 np_cegar.appears (SAT 증거 유도 판정) 의 상태를 A 의 예측과 견준다.
     appears 는 H 의 위수가 2^(n-m) prod mu(v)! 와 같은지(닫힘)도 점검한다.

사용
  python3 -u -I scripts/check_translation_type_groups.py --pipe 작업폴더 --n7 notes/stabilizer_types_n7_closed_absent.tsv [--nmax 9] [--parts ABC]
  작업폴더: scripts/run_bn_class_pipeline.sh 의 작업 폴더 (classes_n2.tsv .. classes_n6.tsv 가 있다). 시간은 A 와 B 가 합쳐 수 초, C 가 90초쯤 (코어 하나; python-sat 필요).
  (--nmax 9 는 A 의 중복도 류를 n = 9 까지 나열한다. C 는 --nmax >= 8 이 필요하다.)

2026-11-15. Claude가 쓴 코드이고 저자 검토 전이다.
"""
import argparse
import itertools
import os
import sys
import time
from collections import Counter
from math import factorial

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pure_translation_absent as PA  # noqa: E402
from analyze_absent_reduction import analyse_class  # noqa: E402
from bn_class_io import read_rows  # noqa: E402
from check_hamming_transvection_family import kernel_basis, vp_coord_perm, vp_translation  # noqa: E402
from np_cegar import Cube, appears  # noqa: E402


def all_multisets(N, n):
    """N 개의 점 위에서 크기 n 인 다중집합 전부를 중복도 벡터(행)로 낸다. 모양 (C(N+n-1, n), N)."""
    combos = np.array(list(itertools.combinations_with_replacement(range(N), n)), dtype=np.int16)
    R = combos.shape[0]
    M = np.zeros((R, N), dtype=np.int8)
    rows = np.arange(R)
    for j in range(n):
        M[rows, combos[:, j]] += 1
    return M


def spanning_rows(M, m):
    """지지가 F_2^m 을 생성하는 행만 남긴다 (점 v 는 열 v-1)."""
    N = (1 << m) - 1
    rk = np.zeros(1 << N, dtype=np.int8)
    for mask in range(1 << N):
        rk[mask] = PA.rank_gf2([v for v in range(1, N + 1) if (mask >> (v - 1)) & 1])
    supp = ((M > 0).astype(np.int64) * (1 << np.arange(N, dtype=np.int64))).sum(axis=1)
    return M[rk[supp] == m]


def orbit_reps(M, n, m):
    """M (모두 지지가 F_2^m 을 생성하는 중복도 벡터) 의 GL(m,2)-궤도마다 한 행(궤도의 가장 작은 행 번호)을 낸다.
    GL(m,2) 는 기본 전단 v -> v + bit_j(v) e_i (i != j) 들이 만든다 (모두 대합이다)."""
    R, N = M.shape
    base = (n + 1) ** np.arange(N, dtype=np.int64)
    keys = M.astype(np.int64) @ base
    order = np.argsort(keys, kind="stable")
    skeys = keys[order]
    assert (np.diff(skeys) > 0).all()
    pts = np.arange(1, N + 1, dtype=np.int64)
    nbr = []
    for i in range(m):
        for j in range(m):
            if i == j:
                continue
            img = pts ^ (((pts >> j) & 1) << i)
            k2 = M[:, img - 1].astype(np.int64) @ base
            pos = np.searchsorted(skeys, k2)
            assert (skeys[pos] == k2).all()
            nbr.append(order[pos])
    lab = np.arange(R)
    while True:
        new = lab.copy()
        for idx in nbr:
            new = np.minimum(new, lab[idx])
        if (new == lab).all():
            break
        lab = new
    return np.nonzero(lab == np.arange(R))[0]


def part_a(nmax):
    """반환: 목록 [(n, m, mu(길이 2^m-1 튜플), |A_mu|, 자유 함수의 수)] (mu 는 GL(m,2)-류의 대표)."""
    out = []
    for m in range(1, 5):
        t0 = time.time()
        GP = PA.perm_table(PA.gl_imgs(m), m)
        N = (1 << m) - 1
        ts = np.arange(1 << m)
        cnt = 0
        for n in range(max(2, m), nmax + 1):
            M = spanning_rows(all_multisets(N, n), m)
            if M.shape[0] == 0:
                continue
            for r in orbit_reps(M, n, m).tolist():
                mu = M[r].astype(np.int64)
                mu_full = np.zeros(1 << m, dtype=np.int64)
                mu_full[1:] = mu
                ok = (mu_full[GP[:, 1:]] == mu_full[1:][None, :]).all(axis=1)
                A = GP[ok]
                elems = (A[:, None, :] ^ ts[None, :, None]).reshape(-1, 1 << m)
                free = int((~PA.nonfree_functions(elems, m)).sum())
                out.append((n, m, tuple(int(x) for x in mu), int(A.shape[0]), free))
                cnt += 1
        print("  m=%d: GL(m,2) 위수 %d, mu 의 류 %d개 (n = 2..%d, %.1f초)" % (m, GP.shape[0], cnt, nmax, time.time() - t0), flush=True)
    return out


def pi_order_of(mu):
    p = 1
    for x in mu:
        p *= factorial(x)
    return p


def table_k1(n, path, statuses):
    """표에서 K = 1 이고 모든 변수가 쓰이는 닫힌 류: [(번호, m, |Pi_T|, 상태)]."""
    cube = Cube(n)
    out = []
    for (j, order, k, status, ng, vals, line) in read_rows(path):
        if status not in statuses:
            continue
        r = analyse_class(cube, n, order, k, ng, vals)
        assert r["ok_k"] and r["ok_order"], (n, j, r)
        if r["essential"] and r["K"] == 1:
            out.append((j, r["m"], r["pi_order"], status))
    return out


def part_c(pred):
    """n = 8, m <= 4 인 mu 마다 B_8 안에서 H = T x| Pi_T 를 만들어 appears 와 견준다."""
    n = 8
    cube = Cube(n)
    cases = [p for p in pred if p[0] == n]
    mism = 0
    nabs_pred = 0
    nabs_obs = 0
    for (nn, m, mu, aord, free) in cases:
        cols = []
        for v, c in enumerate(mu, start=1):
            cols += [v] * c
        assert len(cols) == n
        T = kernel_basis(cols, n)
        assert len(T) == n - m
        gens = [vp_translation(n, t).astype(cube.dt) for t in T]
        for v in range(1, (1 << m)):
            pos = [i for i in range(n) if cols[i] == v]
            for a, b in zip(pos, pos[1:]):
                pi = list(range(n))
                pi[a], pi[b] = b, a
                gens.append(vp_coord_perm(n, pi).astype(cube.dt))
        order = (1 << (n - m)) * pi_order_of(mu)
        r = appears(cube, gens, order=order)
        obs = r["status"]
        pred_absent = (free == 0)
        nabs_pred += pred_absent
        nabs_obs += (obs == "absent")
        ok = obs in ("absent", "present") and ((obs == "absent") == pred_absent)
        mism += (not ok)
        if not ok or pred_absent:
            print("  m=%d mu=%s |A|=%d 자유함수 %d → 예측 %s, appears: %s (k=%s)%s" % (
                m, {v + 1: c for v, c in enumerate(mu) if c}, aord, free, "absent" if pred_absent else "present",
                obs, r.get("k"), "" if ok else "  <<< 불일치"), flush=True)
    print("n=8, m <= 4 인 mu 의 류 %d개: 예측 absent %d 개, np_cegar absent %d 개, 불일치 %d 개" % (len(cases), nabs_pred, nabs_obs, mism))
    return mism


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pipe", required=True)
    ap.add_argument("--n7", required=True)
    ap.add_argument("--nmax", type=int, default=9)
    ap.add_argument("--parts", default="ABC")
    a = ap.parse_args()
    bad = 0

    print("=== A. m = 1..4: 열의 중복도 mu 의 류와 P_mu 의 자유 함수 ===", flush=True)
    pred = part_a(a.nmax)
    pa = Counter()
    pt = Counter()
    for (n, m, mu, aord, free) in pred:
        pt[(n, m)] += 1
        if free == 0:
            pa[(n, m)] += 1
    ns = list(range(2, a.nmax + 1))
    print("\nn 별 (모든 변수가 쓰이는 K = 1 인 닫힌 류의 수 / 그 가운데 안 나타나는 류의 수), m <= 4")
    print("  n:          " + " ".join("%4d" % n for n in ns))
    for m in range(1, 5):
        print("  m=%d 전체:    " % m + " ".join("%4d" % pt[(n, m)] for n in ns))
        print("  m=%d 안 나타남: " % m + " ".join("%4d" % pa[(n, m)] for n in ns))
    print("  m<=4 전체:    " + " ".join("%4d" % sum(pt[(n, m)] for m in range(1, 5)) for n in ns))
    print("  m<=4 안 나타남: " + " ".join("%4d" % sum(pa[(n, m)] for m in range(1, 5)) for n in ns))
    dist = Counter()
    for (n, m, mu, aord, free) in pred:
        if free == 0 and all(c <= 1 for c in mu):
            dist[n] += 1
    print("  (그 가운데 열이 서로 다른 것, 곧 순수 평행이동군) 안 나타남: " + " ".join("%4d" % dist[n] for n in ns))
    print()

    if "B" in a.parts:
        print("=== B. 공개 표와 견주기 (K = 1, 모든 변수가 쓰임) ===", flush=True)
        for n in range(2, 8):
            t0 = time.time()
            if n <= 6:
                rows = table_k1(n, os.path.join(a.pipe, "classes_n%d.tsv" % n), ("absent", "present"))
            else:
                rows = table_k1(n, a.n7, ("absent",))
            tab = Counter((m, p, st) for (_, m, p, st) in rows if m <= 4)
            tab_hi = Counter(st for (_, m, p, st) in rows if m >= 5)
            prd = Counter()
            for (nn, m, mu, aord, free) in pred:
                if nn != n:
                    continue
                st = "present" if free > 0 else "absent"
                if n == 7 and st == "present":
                    continue
                prd[(m, pi_order_of(mu), st)] += 1
            same = (tab == prd)
            hi_ok = tab_hi.get("absent", 0) == 0
            nabs = sum(v for (m, p, st), v in tab.items() if st == "absent")
            print("  n=%d: 표의 K=1 류 %d개 (m<=4 에서 안 나타남 %d, m>=5 %d개 가운데 안 나타남 %d), 예측과 %s, m>=5 는 모두 나타남: %s (%.1f초)" % (
                n, len(rows), nabs, sum(tab_hi.values()), tab_hi.get("absent", 0), "같다" if same else "다르다", hi_ok, time.time() - t0), flush=True)
            if not same:
                bad += 1
                for key in sorted(set(tab) | set(prd)):
                    if tab[key] != prd[key]:
                        print("      어긋남 (m, |Pi_T|, 상태) = %s: 표 %d, 예측 %d" % (key, tab[key], prd[key]))
            if not hi_ok:
                bad += 1
        print()

    if "C" in a.parts:
        if a.nmax < 8:
            print("C 는 --nmax >= 8 이 필요하다")
        else:
            print("=== C. n = 8 의 실제 B_8 군에서 SAT 판정과 견주기 (m <= 4) ===", flush=True)
            t0 = time.time()
            bad += part_c(pred)
            print("  (%.1f초)" % (time.time() - t0))
            print()

    print("결과:", "모두 일치" if bad == 0 else "어긋남 있음")


if __name__ == "__main__":
    main()
