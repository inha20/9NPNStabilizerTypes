#!/usr/bin/env python3
"""
check_npn_bruteforce.py
=======================
n <= 4 에서 모든 함수(2^(2^n)개)의 NPN 안정화군을 정의대로 구해 npn_stabilizer_types.g 의 표(방법 A)와 견준다.
표지 표를 쓰지 않는 독립 검증이다.

함수 f 마다 B_n 의 모든 원소 g = (pi, c)  (x -> pi(x xor c)) 를 시험한다.
    H(f) = { g : f o g = f }          (대칭)
    A(f) = { g : f o g = not f }      (반대칭)
    Kbar(f) = H(f) 합 A(f)             |A| 는 0 이거나 |H| 이다.
그리고 군의 종류를 다음 열두 값으로 적는다 (표지 표 쪽 TSV 에서도 같은 값을 읽어 묶는다).
    반대칭 유무, |H|, H 의 순열 부분의 위수, H 의 평행이동 부분의 차원 d, H 의 순열 부분의 변수 위 궤도 유형,
    H 의 꼭짓점 위 궤도 수, H 의 고정 꼭짓점 수, |Kbar|, Kbar 의 순열 부분의 위수, Kbar 의 평행이동 부분의 차원,
    Kbar 의 순열 부분의 변수 위 궤도 유형, Kbar 의 꼭짓점 위 궤도 수.
묶음마다 함수 수와 NPN 류 수를 견준다. 그 뒤 표의 각 행에 적힌 대표 함수(진리표 16진수)의 안정화군을 정의대로 구해
행의 값(반대칭 유무, |H|, |Kbar| 와 위의 열두 값)과 견준다 (n = 5 도 한다). NPN 류 수는 (묶음 안 함수들의 |Kbar| 의 합) / (2 |B_n|) 이다
(함수 하나의 NPN 궤도 크기가 2 |B_n| / |K| 이고 K 와 Kbar 의 위수가 같으므로).

실행: python3 -I scripts/check_npn_bruteforce.py <n> <npn_stabilizer_types_n?.tsv>                     (n <= 4: 파이썬이 직접 계산)
      python3 -I scripts/check_npn_bruteforce.py <n> <npn_stabilizer_types_n?.tsv> <열거 출력> [...]   (n = 3, 4, 5: C 열거 출력과 견줌)
      (열거 출력은 ./enumnpn? <부분> <부분 수> 의 출력. 여러 부분으로 나눠 돌렸으면 모두 준다.)
"""
import itertools
import math
import re
import sys
from collections import defaultdict

import numpy as np


def build_group(n):
    """B_n 의 원소 (pi, c, 꼭짓점 사상). 꼭짓점 사상은 x -> pi(x xor c), pi(y) 의 i번째 좌표 = y 의 pi[i]번째 좌표."""
    NP = 1 << n
    elems = []
    for pi in itertools.permutations(range(n)):
        for c in range(NP):
            vm = []
            for x in range(NP):
                y = x ^ c
                r = 0
                for i in range(n):
                    r |= ((y >> pi[i]) & 1) << i
                vm.append(r)
            elems.append((pi, c, tuple(vm)))
    return elems


def components(size, links):
    par = list(range(size))

    def find(a):
        while par[a] != a:
            par[a] = par[par[a]]
            a = par[a]
        return a

    for a, b in links:
        ra, rb = find(a), find(b)
        if ra != rb:
            par[ra] = rb
    return [find(a) for a in range(size)]


def pattern_of_sizes(sizes):
    return int("".join(str(s) for s in sorted(sizes)))


def var_pattern(n, perm_set):
    comp = components(n, [(i, pi[i]) for pi in perm_set for i in range(n)])
    cnt = defaultdict(int)
    for r in comp:
        cnt[r] += 1
    return pattern_of_sizes(cnt.values())


def vertex_data(NP, vms):
    comp = components(NP, [(x, vm[x]) for vm in vms for x in range(NP)])
    cnt = defaultdict(int)
    for r in comp:
        cnt[r] += 1
    fixed = sum(1 for x in range(NP) if all(vm[x] == x for vm in vms))
    return len(cnt), fixed


def group_key(n, elems, idx_list):
    """원소 목록 (B_n 의 원소 번호들) 이 이루는 군의 (위수, 순열 부분의 위수, d, 변수 위 궤도 유형, 꼭짓점 궤도 수, 고정 꼭짓점 수)"""
    NP = 1 << n
    ident = tuple(range(n))
    perm_set = {elems[e][0] for e in idx_list}
    ntrans = sum(1 for e in idx_list if elems[e][0] == ident)
    vorb, fixv = vertex_data(NP, [elems[e][2] for e in idx_list])
    return (len(idx_list), len(perm_set), int(math.log2(ntrans)), var_pattern(n, perm_set), vorb, fixv)


def brute_force(n):
    NP = 1 << n
    FULL = (1 << NP) - 1
    elems = build_group(n)
    nfun = 1 << NP
    F = np.arange(nfun, dtype=np.uint64)
    nelem = len(elems)
    stabH = np.zeros((nfun, nelem), dtype=bool)
    stabA = np.zeros((nfun, nelem), dtype=bool)
    NF = np.uint64(FULL) ^ F
    for e, (_, _, vm) in enumerate(elems):
        res = np.zeros(nfun, dtype=np.uint64)
        for x in range(NP):
            res |= ((F >> np.uint64(vm[x])) & np.uint64(1)) << np.uint64(x)
        stabH[:, e] = (res == F)
        stabA[:, e] = (res == NF)
    groups = defaultdict(int)       # 키 -> 함수 수
    ksum = defaultdict(int)         # 키 -> |Kbar| 의 합
    for f in range(nfun):
        H = np.flatnonzero(stabH[f]).tolist()
        A = np.flatnonzero(stabA[f]).tolist()
        if A and len(A) != len(H):
            raise SystemExit("|A| != |H| at f = %d" % f)
        hk = group_key(n, elems, H)
        kk = group_key(n, elems, H + A) if A else hk
        key = (1 if A else 0,) + hk + (kk[0], kk[1], kk[2], kk[3], kk[4])
        groups[key] += 1
        ksum[key] += len(H) + len(A)
    return groups, ksum


def load_rows(path):
    return [line.rstrip("\n").split("\t") for line in open(path, encoding="utf-8") if line.strip() and not line.startswith("#")]


def check_representatives(n, path):
    """표의 각 행의 대표 함수에서 H, A 를 정의대로 구해 행의 값과 견준다."""
    NP = 1 << n
    FULL = (1 << NP) - 1
    elems = build_group(n)
    VM = np.array([vm for (_, _, vm) in elems], dtype=np.int64)
    bad = 0
    nrows = 0
    for r in load_rows(path):
        nrows += 1
        if r[21] == "-":
            print("  대표 함수가 없는 행:", r[:8])
            bad += 1
            continue
        f = int(r[21], 16)
        bits = np.array([(f >> x) & 1 for x in range(NP)], dtype=np.int8)
        img = bits[VM]                      # img[e][x] = f(g_e(x))
        H = np.flatnonzero((img == bits).all(axis=1)).tolist()
        A = np.flatnonzero((img == 1 - bits).all(axis=1)).tolist()
        hk = group_key(n, elems, H)
        kk = group_key(n, elems, H + A) if A else hk
        key = (1 if A else 0,) + hk + (kk[0], kk[1], kk[2], kk[3], kk[4])
        row_key = (int(r[4]), int(r[6]), int(r[7]), int(r[12]), pattern_of_text(r[9]), int(r[14]), int(r[16]),
                   int(r[5]), int(r[8]), int(r[13]), pattern_of_text(r[10]), int(r[15]))
        if key != row_key or len(H) + len(A) != int(r[0]):
            print("  대표 함수가 행과 다르다:", r[21], key, row_key, len(H) + len(A), r[0])
            bad += 1
    print("대표 함수 점검: 행 %d개 가운데 다른 것 %d개" % (nrows, bad))
    return bad


def pattern_of_text(text):
    nums = sorted(int(t) for t in re.findall(r"\d+", text))
    return int("".join(str(t) for t in nums))


def load_table(path):
    groups = defaultdict(int)
    classes = defaultdict(int)
    for line in open(path, encoding="utf-8"):
        if not line.strip() or line.startswith("#"):
            continue
        r = line.rstrip("\n").split("\t")
        antisym, kb, h, hb, kbb = int(r[4]), int(r[5]), int(r[6]), int(r[7]), int(r[8])
        pat_h, pat_k = pattern_of_text(r[9]), pattern_of_text(r[10])
        d_h, d_k, vo_h, vo_k, fixv = int(r[12]), int(r[13]), int(r[14]), int(r[15]), int(r[16])
        key = (antisym, h, hb, d_h, pat_h, vo_h, fixv, kb, kbb, d_k, pat_k, vo_k)
        groups[key] += int(r[2])
        classes[key] += int(r[3])
    return groups, classes


def load_enum(paths, n):
    """enumerate_npn_stabilizers.c 의 출력. 최상위 비트가 0 인 절반만 훑었으므로 2배 한다."""
    groups = defaultdict(int)
    ksum = defaultdict(int)
    nv = None
    for path in paths:
        for line in open(path, encoding="utf-8"):
            if line.startswith("# NV"):
                this = int(line.split()[2])
                if nv not in (None, this):
                    raise SystemExit("변수 수가 다른 파일이 섞였다")
                nv = this
            if not line.startswith("npn "):
                continue
            t = [int(x) for x in line.split()[1:]]
            antisym = t[0]
            hk, kk, count = t[1:7], t[7:13], t[13]
            key = (antisym, hk[0], hk[1], hk[2], hk[3], hk[4], hk[5], kk[0], kk[1], kk[2], kk[3], kk[4])
            groups[key] += 2 * count
            ksum[key] += 2 * count * kk[0]
    if nv != n:
        raise SystemExit("표는 n = %s, 열거는 n = %s 이다" % (n, nv))
    return groups, ksum


def main():
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    n = int(sys.argv[1])
    if len(sys.argv) > 3:
        bf_groups, bf_ksum = load_enum(sys.argv[3:], n)
        how = "C 열거"
    else:
        if n > 4:
            raise SystemExit("n <= 4 만 직접 계산한다 (n = 5 는 enumerate_npn_stabilizers.c 의 출력을 함께 준다)")
        bf_groups, bf_ksum = brute_force(n)
        how = "전수 계산"
    tab_groups, tab_classes = load_table(sys.argv[2])
    gnpn = 2 * math.factorial(n) * 2 ** n
    total_bf, total_tab = sum(bf_groups.values()), sum(tab_groups.values())
    print("n = %d: 함수 수의 합: 표지 표 %d, %s %d (2^(2^n) = %d)" % (n, total_tab, how, total_bf, 2 ** (2 ** n)))
    keys = sorted(set(bf_groups) | set(tab_groups))
    bad = 0
    ncls_bf = 0
    print("묶음 열쇠 = (반대칭, |H|, |Hbar|, d(H), 변수궤도(H), 꼭짓점궤도(H), 고정점(H), |Kbar|, |Kbar의 순열 부분|, d(Kbar), 변수궤도(Kbar), 꼭짓점궤도(Kbar))")
    print("%-52s | %14s %14s | %10s %10s" % ("열쇠", "함수: 표지 표", how, "류: 표지 표", how))
    for k in keys:
        fa, fb = tab_groups.get(k, 0), bf_groups.get(k, 0)
        ca = tab_classes.get(k, 0)
        s = bf_ksum.get(k, 0)
        cb = s // gnpn
        if s % gnpn:
            print("!!! 류 수가 정수가 아니다:", k, s, gnpn)
            bad += 1
        ncls_bf += cb
        flag = "" if (fa == fb and ca == cb) else "   <-- 다르다"
        if flag:
            bad += 1
        print("%-52s | %14d %14d | %10d %10d%s" % (str(k), fa, fb, ca, cb, flag))
    print("NPN 류 수: 표지 표 %d, %s %d" % (sum(tab_classes.values()), how, ncls_bf))
    print("묶음 %d개 가운데 다른 것: %d개" % (len(keys), bad))
    bad_reps = check_representatives(n, sys.argv[2])
    ok = (bad == 0 and bad_reps == 0 and total_bf == total_tab == 2 ** (2 ** n) and sum(tab_classes.values()) == ncls_bf)
    print("모두 일치한다." if ok else "일치하지 않는다.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
