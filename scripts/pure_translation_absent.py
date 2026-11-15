#!/usr/bin/env python3
"""
pure_translation_absent.py
==========================
평행이동만으로 이루어진 닫힌 군 T <= Z_2^n (순수 평행이동 군)이 언제 안 나타나는지 모든 n 에 대해 가린다.

바탕이 되는 사실 (증명은 안정화군 노트 3절 "일반 n")
  T 를 Z_2^n 의 부분공간(이진 부호)이라 하고 m = n - dim T, 몫 사상 M: Z_2^n -> Z_2^n / T = F_2^m 이라 하자.
  M x = sum_i x_i h_i 로 쓰면 열 h_1..h_n 이 F_2^m 을 생성하고 T = {x : sum x_i h_i = 0} 이다.
  (1) T 가 닫힌 군이고 모든 변수가 쓰임 <=> h_i 가 서로 다른 0 이 아닌 벡터이다. 집합 S = {h_i} 를 GL(m,2) 의 작용으로 나눈 류가 T 의 류이다.
  (2) N(T)/T 는 F_2^m 위에 아핀 변환으로 작용하고 P_S = F_2^m x| Stab_{GL(m,2)}(S) 이다.
  (3) F_2^m 의 함수 g 의 주기군이 0 이면 f = g 합성 M 의 안정화군은 Stab_{P_S}(g) 의 원상이다. 주기가 있으면 안정화군이 T 보다 크다.
      그러므로 T 가 나타남 <=> P_S 의 안정화군이 자명한 g: F_2^m -> F_2 (자유 함수)가 있다.
  (4) 꼭짓점 궤도 수는 k = 2^m 이다.
  (5) 아핀 변환 (L, t) 가 g 의 안정화군에 들어 있고 P_S 에 들어 있는 것은 L 이 S 를 보존하는 것과 같다(평행이동은 모두 P_S 의 원소이므로).
      그래서 g 가 P_S 에서 자유 <=> g 의 주기가 0 이고 g 의 안정화군의 선형 부분이 Stab(S) 와 항등원에서만 만난다.

하는 일
  A. m = 1..4 에 대해 F_2^m 에서 0 을 뺀 집합의 부분집합 S 를 GL(m,2) 로 나눠 모든 류를 구하고, F_2^m 을 생성하는 S 마다 P_S 의 자유 함수가 있는지 가린다.
     (2^m 개 점 위의 함수 전부를 검사한다. 비자명한 원소가 고정하는 함수는 그 원소의 순환 분할에서 상수인 함수이므로
      순환 분할마다 상수인 함수를 표시해 표시되지 않은 함수가 있는지 본다.)
  B. m = 5..8 에 대해 아핀 변환군 AGL(m,2) 에서 안정화군이 자명한 함수가 있는지를 무작위 함수로 본다. 안정화군은 자기상관 함수를
     보존해야 한다는 조건으로 가지를 친 되짚기로 정확히 센다(이 계산을 m = 3, 4 에서 모든 아핀 변환을 훑는 방법과 견주어 시험한다).
     m = 5 는 안정화군의 위수 분포를 적고, GL(5,2) 의 극대 부분군 다섯 류마다 그 안에서 자유인 함수를 찾는다.
     (AGL(5,2) 에 자유 함수가 없다는 것은 이 표본이 아니라 scripts/agl52_orbits.c 의 전수 열거로 증명된다. 표본은 그 결과와 맞는지의 점검이다.)
  C. m >= 9 일 때 |AGL(m,2)| <= 2^(2^(m-2)) 임을 확인한다 (아핀 변환 g != 1 이 고정하는 함수의 비율이 2^(-2^(m-2)) 이하이다).
  D. 표(공개 표)에서 n <= 7 의 순수 평행이동이고 모든 변수가 쓰이는 안 나타나는 류를 뽑아 A 의 예측과 견준다 (n, m, k, |N(H)/H|).

사용
  python3 -u -I scripts/pure_translation_absent.py [--parts ABCD] --pipe 작업폴더 --n7 표.tsv.gz --annot7 주석표.tsv
  (D 를 생략하려면 --parts ABC). 시간은 1분 안팎 (코어 하나).

2026-11-15. Claude가 쓴 코드이고 저자 검토 전이다.
"""
import argparse
import math
import os
import sys
import time
from collections import Counter, defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bn_class_io import read_rows  # noqa: E402
from check_dummy_variable_lift import analyse  # noqa: E402


# ----------------------------------------------------------------------------------------------
# A. m <= 4: 부분집합의 류와 자유 함수
# ----------------------------------------------------------------------------------------------
def gl_imgs(m):
    """GL(m,2) 의 모든 원소를 기저 e_i = 1<<i 의 상(정수 m 개)으로 나열한다. 모양 (|GL|, m)."""
    N = 1 << m
    out = []

    def rec(imgs, span):
        if len(imgs) == m:
            out.append(imgs)
            return
        for v in range(1, N):
            if v not in span:
                rec(imgs + [v], span | {s ^ v for s in span})

    rec([], {0})
    return np.array(out, dtype=np.int64)


def perm_table(imgs, m):
    """각 원소가 F_2^m 의 점 x = 0..2^m-1 에 주는 상. 모양 (|GL|, 2^m)."""
    xs = np.arange(1 << m)
    P = np.zeros((imgs.shape[0], 1 << m), dtype=np.int64)
    for i in range(m):
        P ^= ((xs >> i) & 1)[None, :] * imgs[:, i][:, None]
    return P


def subset_canon(m, GP):
    """F_2^m 에서 0 을 뺀 집합의 부분집합(점 v 가 비트 v-1)마다 GL(m,2) 궤도의 가장 작은 마스크."""
    N = (1 << m) - 1
    masks = np.arange(1 << N, dtype=np.int64)
    canon = masks.copy()
    lo_bits = min(N, 8)
    hi_bits = N - lo_bits
    lo_idx = np.arange(1 << lo_bits)
    hi_idx = np.arange(1 << hi_bits)
    lo_mask = (1 << lo_bits) - 1
    for g in range(GP.shape[0]):
        bitval = 1 << (GP[g, 1:] - 1)
        lo = np.zeros(1 << lo_bits, dtype=np.int64)
        for j in range(lo_bits):
            lo |= ((lo_idx >> j) & 1) * bitval[j]
        hi = np.zeros(1 << hi_bits, dtype=np.int64)
        for j in range(hi_bits):
            hi |= ((hi_idx >> j) & 1) * bitval[lo_bits + j]
        np.minimum(canon, lo[masks & lo_mask] | hi[masks >> lo_bits], out=canon)
    return masks, canon


def rank_gf2(vecs):
    basis = []
    for v in vecs:
        for b in basis:
            v = min(v, v ^ b)
        if v:
            basis.append(v)
    return len(basis)


def nonfree_functions(elems, m):
    """elems: (E, 2^m) 점의 순열들. 항등원이 아닌 원소가 고정하는 함수(길이 2^m 의 0/1 수열을 정수로 본 것)를 True 로 표시한 배열을 낸다."""
    K = 1 << m
    mark = np.zeros(1 << K, dtype=bool)
    ident = np.arange(K, dtype=np.int64)
    parts = set()
    CH = 40000
    for s in range(0, elems.shape[0], CH):
        q = elems[s:s + CH].copy()
        lab = np.broadcast_to(ident, q.shape).copy()
        for _ in range(m):                          # m 번 제곱하면 길이 2^m 이하의 순환을 덮는다
            lab = np.minimum(lab, np.take_along_axis(lab, q, axis=1))
            q = np.take_along_axis(q, q, axis=1)
        keep = ~(lab == ident).all(axis=1)
        for row in np.unique(lab[keep], axis=0):
            parts.add(row.tobytes())
    for key in parts:
        lab = np.frombuffer(key, dtype=np.int64)
        blocks = {}
        for x, l in enumerate(lab.tolist()):
            blocks[l] = blocks.get(l, 0) | (1 << x)
        funcs = np.zeros(1, dtype=np.int64)
        for bm in blocks.values():
            funcs = np.concatenate([funcs, funcs | bm])
        mark[funcs] = True
    return mark


def part_a(maxm=4):
    """반환: 목록 [(m, n, |A_S|, 자유 함수의 수, S 의 점 목록)] (S 는 F_2^m 을 생성하는 류의 대표)."""
    res = []
    for m in range(1, maxm + 1):
        t0 = time.time()
        imgs = gl_imgs(m)
        GP = perm_table(imgs, m)
        N = (1 << m) - 1
        masks, canon = subset_canon(m, GP)
        reps = np.nonzero(canon == masks)[0]
        nspan = 0
        for mask in reps.tolist():
            S = [v for v in range(1, N + 1) if (mask >> (v - 1)) & 1]
            if rank_gf2(S) < m:
                continue
            nspan += 1
            member = np.zeros(N + 1, dtype=bool)
            member[S] = True
            ok = member[GP[:, S]].all(axis=1)
            A = GP[ok]
            ts = np.arange(1 << m)
            elems = (A[:, None, :] ^ ts[None, :, None]).reshape(-1, 1 << m)
            mark = nonfree_functions(elems, m)
            res.append((m, len(S), A.shape[0], int((~mark).sum()), S))
        print("  m=%d: GL(m,2) 위수 %d, F_2^m 을 생성하는 부분집합의 류 %d개 (%.1f초)" % (m, imgs.shape[0], nspan, time.time() - t0), flush=True)
    return res


# ----------------------------------------------------------------------------------------------
# B. m >= 5: AGL(m,2) 에서 안정화군이 자명한 함수
# ----------------------------------------------------------------------------------------------
def wht(a):
    a = np.array(a, dtype=np.int64)
    n = len(a)
    h = 1
    while h < n:
        a = a.reshape(-1, 2, h)
        x, y = a[:, 0, :].copy(), a[:, 1, :].copy()
        a[:, 0, :], a[:, 1, :] = x + y, x - y
        a = a.reshape(n)
        h *= 2
    return a


def autocorr(gt, m):
    """a[v] = sum_x (-1)^(g(x) + g(x+v))."""
    f = 1 - 2 * np.asarray(gt, dtype=np.int64)
    W = wht(f)
    return wht(W * W) // (1 << m)


def stab_elements(gt, m):
    """아핀 변환 x -> L x + t 가운데 g(Lx + t) = g(x) 인 것을 (Lx 배열, t) 의 목록으로 낸다.
    선형 부분 L 은 자기상관 함수를 보존해야 한다는 조건으로 가지를 친 되짚기로 찾는다."""
    N = 1 << m
    gt = np.asarray(gt, dtype=np.int64)
    a = autocorr(gt, m)
    cls = defaultdict(list)
    for v in range(1, N):
        cls[int(a[v])].append(v)
    size = {v: len(cls[int(a[v])]) for v in range(1, N)}
    basis, span = [], {0}
    for v in sorted(range(1, N), key=lambda v: (size[v], v)):
        if v not in span:
            basis.append(v)
            span |= {s ^ v for s in span}
            if len(basis) == m:
                break
    found = []

    def rec(i, img_of):
        if i == m:
            Lx = np.array([img_of[x] for x in range(N)], dtype=np.int64)
            for t in range(N):
                if np.array_equal(gt[Lx ^ t], gt):
                    found.append((Lx, t))
            return
        b = basis[i]
        imgspan = set(img_of.values())
        for w in cls[int(a[b])]:
            if w in imgspan:
                continue
            new = {}
            ok = True
            for s, ims in img_of.items():
                x2, y2 = s ^ b, ims ^ w
                if a[y2] != a[x2]:
                    ok = False
                    break
                new[x2] = y2
            if ok:
                merged = dict(img_of)
                merged.update(new)
                rec(i + 1, merged)

    rec(0, {0: 0})
    return found


def brute_stab_order(gt, m, GP):
    """모든 아핀 변환을 훑어서 센다 (m <= 4)."""
    gt = np.asarray(gt, dtype=np.int64)
    tot = 0
    for t in range(1 << m):
        tot += int((gt[GP ^ t] == gt).all(axis=1).sum())
    return tot


# GL(5,2) 의 극대 부분군 다섯 류 (ATLAS): 점 안정화군, 초평면 안정화군, 직선 안정화군, 평면 안정화군, 싱어 순환군의 정규화군 31:5
def lin_imgs(Lx, m):
    return tuple(int(Lx[1 << i]) for i in range(m))


def compose(L1, L2, m):
    """(L1 합성 L2) 의 기저 상. L 은 기저 상의 튜플."""
    def ap(L, x):
        r = 0
        for i in range(m):
            if (x >> i) & 1:
                r ^= L[i]
        return r
    return tuple(ap(L1, L2[i]) for i in range(m))


def singer_normalizer():
    """x^5 + x^2 + 1 이 정하는 GF(32) 에서 x 곱하기(C)와 제곱(F)이 생성하는 군 (위수 155)."""
    def mulx(a):
        a <<= 1
        return a ^ 0b100101 if a & 32 else a
    C = (2, 4, 8, 16, 5)
    pw = [1]
    for _ in range(9):
        pw.append(mulx(pw[-1]))
    F = tuple(pw[2 * i] for i in range(5))
    ident = (1, 2, 4, 8, 16)
    group = {ident}
    frontier = [ident]
    while frontier:
        nxt = []
        for g in frontier:
            for h in (C, F):
                x = compose(h, g, 5)
                if x not in group:
                    group.add(x)
                    nxt.append(x)
        frontier = nxt
    return group


def maximal_tests(m=5):
    sing = singer_normalizer()
    return [
        ("점 안정화군 2^4:L4(2)", lambda L: L[0] == 1),
        ("초평면 안정화군 2^4:L4(2)", lambda L: all((L[i] & 1) == (1 if i == 0 else 0) for i in range(m))),
        ("직선 안정화군 2^6:(S3 x L3(2))", lambda L: L[0] < 4 and L[1] < 4),
        ("평면 안정화군 2^6:(S3 x L3(2))", lambda L: L[0] < 8 and L[1] < 8 and L[2] < 8),
        ("싱어 정규화군 31:5", lambda L: L in sing),
    ], len(sing)


def part_b():
    print("  [시험] 되짚기가 모든 아핀 변환 훑기와 같은 값을 내는가 (무작위 함수와 대칭이 있는 함수)")
    rng = np.random.default_rng(12345)
    bad = 0
    for m in (3, 4):
        GP = perm_table(gl_imgs(m), m)
        cnt = 0
        for _ in range(60):
            gt = rng.integers(0, 2, size=1 << m)
            if rng.random() < 0.5:                  # 대칭이 있는 함수: 아핀 변환으로 옮긴 함수와 섞는다
                gt = np.where(rng.random(1 << m) < 0.5, gt, gt[GP[rng.integers(GP.shape[0])]])
            e, b = len(stab_elements(gt, m)), brute_stab_order(gt, m, GP)
            cnt += 1
            if e != b:
                bad += 1
        print("    m=%d: %d개 함수 견줌" % (m, cnt), flush=True)
    print("    어긋난 함수 %d개" % bad)
    assert bad == 0

    print("  [m=5] AGL(5,2) 의 안정화군 위수 분포 (무작위 함수)")
    t0 = time.time()
    NS = 6000
    dist = Counter()
    for seed in range(NS):
        gt = np.random.default_rng(10 ** 6 + seed).integers(0, 2, size=32)
        dist[len(stab_elements(gt, 5))] += 1
    print("    표본 %d개 (%.1f초): " % (NS, time.time() - t0) + ", ".join("%d:%d" % kv for kv in sorted(dist.items())))
    print("    가장 작은 위수 %d. 안정화군이 자명한 함수의 궤도가 있다면 그 궤도는 모든 함수의 |AGL(5,2)|/2^32 = %.4f 를 차지하므로 표본 %d개에서 하나도 안 나올 확률은 %.2e 이하이다."
          % (min(dist), agl_order(5) / 2 ** 32, NS, (1 - agl_order(5) / 2 ** 32) ** NS))
    print("    (정확한 답은 scripts/agl52_orbits.c 의 전수 열거 notes/agl52_orbits_output.txt: 궤도 382개, 안정화군 위수의 최솟값 4, 자유 궤도 0개. 표본은 그 결과와 맞는지의 점검이다.)")

    print("  [m=5] GL(5,2) 의 극대 부분군마다, 주기가 0 이고 안정화군의 선형 부분이 그 군과 항등원에서만 만나는 함수")
    tests, singer_order = maximal_tests()
    assert singer_order == 155, singer_order
    for name, inM in tests:
        found = None
        for seed in range(20000):
            gt = np.random.default_rng(2 * 10 ** 6 + seed).integers(0, 2, size=32)
            el = stab_elements(gt, 5)
            if any(t == 0 and False for _, t in el):
                continue
            lin = [lin_imgs(Lx, 5) for Lx, t in el]
            ident = (1, 2, 4, 8, 16)
            if len(set(lin)) != len(lin):           # 같은 L 이 둘 이상 = 주기가 0 이 아니다
                continue
            if all((L == ident) or (not inM(L)) for L in lin):
                found = (seed, len(el))
                break
        print("    %-34s %s" % (name, ("씨앗 %d 의 함수 (안정화군 위수 %d)" % found) if found else "못 찾았다 (!)"))
        assert found, name

    print("  [m=6..8] AGL(m,2) 에서 안정화군이 자명한 함수")
    for m in (6, 7, 8):
        t0 = time.time()
        seed = 0
        while True:
            gt = np.random.default_rng(3 * 10 ** 6 + seed).integers(0, 2, size=1 << m)
            if len(stab_elements(gt, m)) == 1:
                break
            seed += 1
        hexs = "".join("%x" % int("".join(str(int(b)) for b in gt[i:i + 4][::-1]), 2) for i in range(0, 1 << m, 4))
        print("    m=%d: 씨앗 %d 의 무작위 함수의 안정화군이 자명하다 (%.1f초)" % (m, seed, time.time() - t0), flush=True)


def agl_order(m):
    o = 1 << m
    for i in range(m):
        o *= (1 << m) - (1 << i)
    return o


def part_c():
    print("  m   log2|AGL(m,2)|   2^(m-2)   |AGL| <= 2^(2^(m-2)) ?   m^2+m <= 2^(m-2) ?")
    ok = True
    for m in range(5, 21):
        o = agl_order(m)
        ok_exact = o <= 2 ** (2 ** (m - 2))
        ok_simple = (m * m + m) <= 2 ** (m - 2)
        print("  %-3d %-16.2f %-9d %-22s %s" % (m, math.log2(o), 2 ** (m - 2), ok_exact, ok_simple))
        if m >= 9 and not (ok_exact and ok_simple):
            ok = False
    return ok


# ----------------------------------------------------------------------------------------------
# D. 표와 견주기
# ----------------------------------------------------------------------------------------------
def read_annot(path):
    d = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            d[int(p[0])] = dict(order=int(p[1]), k=int(p[2]), C=int(p[4]), Sigma=int(p[5]), NH=int(p[6]), expl=int(p[7]), typ=p[8])
    return d


def table_pure_absent(n, tablepath, annotpath):
    """표에서 순수 평행이동이고 모든 변수가 쓰이는 안 나타나는 류의 (n, m, k, |N/H|) 목록."""
    ann = read_annot(annotpath)
    out = []
    for j, order, k, status, ng, vals, _ in read_rows(tablepath):
        if status != "absent":
            continue
        pure = True
        for t in range(ng):
            im = vals[t * (n + 1):(t + 1) * (n + 1)]
            c = im[0]
            if any((im[1 + i] ^ c) != (1 << i) for i in range(n)):
                pure = False
                break
        if not pure:
            continue
        r, kk, lens = analyse(n, vals, ng)
        if r != 0:
            continue
        d = int(round(math.log2(order)))
        assert (1 << d) == order
        a = ann[j]
        assert a["order"] == order and a["k"] == k
        out.append((n, n - d, k, a["NH"]))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parts", default="ABCD")
    ap.add_argument("--pipe")
    ap.add_argument("--n7")
    ap.add_argument("--annot7")
    a = ap.parse_args()

    res = None
    if "A" in a.parts or "D" in a.parts:
        print("=== A. m = 1..4: 부분집합 S 의 류와 자유 함수 ===", flush=True)
        res = part_a(4)
        print()
        print("S 의 크기 n 별, F_2^m 을 생성하는 류의 수(전체)와 그 가운데 안 나타나는 류의 수")
        allc, absc = Counter(), Counter()
        for m, n, aord, nfree, S in res:
            allc[(n, m)] += 1
            if nfree == 0:
                absc[(n, m)] += 1
        ns = sorted({n for n, _ in allc})
        print("  n:          " + " ".join("%3d" % n for n in ns))
        for m in range(1, 5):
            print("  m=%d 전체:    " % m + " ".join("%3d" % allc[(n, m)] for n in ns))
            print("  m=%d 안 나타남: " % m + " ".join("%3d" % absc[(n, m)] for n in ns))
        print("  합친 안 나타나는 류 (n 별): " + " ".join("%3d" % sum(absc[(n, m)] for m in range(1, 5)) for n in ns))
        print()
        print("안 나타나는 류의 목록 (m, n, |Stab(S)|, |P_S| = 2^m |Stab(S)|, k = 2^m)")
        for m, n, aord, nfree, S in sorted(res, key=lambda r: (r[0], r[1], r[2])):
            if nfree == 0:
                print("  m=%d n=%2d |A|=%5d |P|=%6d k=%2d  S=%s" % (m, n, aord, aord << m, 1 << m, S))
        mmax = max((m for m, n, aord, nfree, S in res if nfree == 0), default=0)
        nmax = max((n for m, n, aord, nfree, S in res if nfree == 0), default=0)
        print("\nm <= 4 에서 안 나타나는 류가 있는 가장 큰 m = %d (k = %d), 가장 큰 n = %d" % (mmax, 1 << mmax, nmax))
        print()

    if "B" in a.parts:
        print("=== B. m >= 5: AGL(m,2) 에서 안정화군이 자명한 함수 ===", flush=True)
        part_b()
        print()

    if "C" in a.parts:
        print("=== C. m >= 9: 세어서 보이기 ===")
        okc = part_c()
        print("m = 9..20 에서 부등식이 성립한다: %s" % okc)
        print("(m >= 21 은 m^2+m <= 2^(m-2) 가 m 이 늘면 계속 성립하므로 같다.)")
        print()

    if "D" in a.parts:
        print("=== D. 표와 견주기 (n <= 7) ===")
        pred = Counter()
        for m, n, aord, nfree, S in res:
            if nfree == 0 and n <= 7:
                pred[(n, m, 1 << m, aord << m)] += 1
        found = Counter()
        for n in range(2, 8):
            if n <= 6:
                tab = os.path.join(a.pipe, "classes_n%d.tsv" % n)
                ann = os.path.join(a.pipe, "annot_n%d.tsv" % n)
            else:
                tab, ann = a.n7, a.annot7
            for key in table_pure_absent(n, tab, ann):
                found[key] += 1
        keys = sorted(set(pred) | set(found))
        print("  (n, m, k, |N(H)/H|): 예측 / 표")
        bad = 0
        for kk in keys:
            flag = "" if pred[kk] == found[kk] else "   <-- 다르다"
            if pred[kk] != found[kk]:
                bad += 1
            print("  %s: %d / %d%s" % (kk, pred[kk], found[kk], flag))
        print("순수 평행이동이고 모든 변수가 쓰이는 안 나타나는 류: 예측 %d개, 표 %d개, 다른 항목 %d개"
              % (sum(pred.values()), sum(found.values()), bad))
        print("결과:", "모두 일치" if bad == 0 else "불일치가 있다")
        sys.exit(0 if bad == 0 else 1)


if __name__ == "__main__":
    main()
