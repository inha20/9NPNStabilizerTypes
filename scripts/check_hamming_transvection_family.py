#!/usr/bin/env python3
"""
check_hamming_transvection_family.py
====================================
전단(transvection)으로 만든 안 나타나는 군의 두 족을 실제 B_n 의 군으로 만들어 확인한다 (노트 3절 "일반 n (3)", "(4)").

[axis] 해밍-전단 족 (축 고정; r >= 2, n = 2^r, m = r+1)
  U = { h_i = i + 2^r : i = 0..n-1 }  (F_2^m 에서 마지막 좌표가 1 인 벡터 전부, 곧 W0 의 여집합, W0 = {마지막 좌표 0})
  T = { x in F_2^n : sum_i x_i h_i = 0 }  (확장 해밍 부호; 평행이동 부분)
  E = { pi_u : i -> i XOR u (u in W0) }   (좌표 위 정칙 작용, 순수 순열; F_2^m 위에서는 축 W0 을 점마다 고정하는 전단들)
  H = T x| E  (위수 2^(n-m) * 2^r = 2^(n-1))
  H 의 꼭짓점 궤도: 길이 |T| 인 궤도 2^r 개 (짝수 무게 T-잉여류, 각각 E 가 고정) + 길이 2^(n-1) 인 궤도 1 개 (홀수 무게 전부) → k = n + 1.
  H 는 닫혀 있고, N(H)/H 는 AGL(r,2) 이며 2^r 개의 짝수 궤도 위에는 자연스럽게, 마지막 궤도 위에는 자명하게 작용한다.
  AGL(r,2) 에 자유로운 함수가 없으면 H 는 나타나지 않는다 (r=2,3,4: |AGL(r,2)| > 2^(2^r); r=5: scripts/agl52_orbits.c 의 전수 열거).

[centre] 중심 고정 족 (m >= 2, n = 2^m - 1)
  S = F_2^m 의 0 이 아닌 점 전부 (열 h_i = i+1), T = 해밍 부호 [n, n-m] 의 평행이동군,
  K = { v -> v + psi(v) u0 : psi(u0) = 0 } (중심 u0 = 1 을 고정하고 축이 달라지는 전단들, 위수 2^(m-1)) 를 좌표 순열로 옮긴 군, H = T x| K~ (위수 2^(n-1)).
  꼭짓점 궤도: K 의 F_2^m 위 궤도 = {0}, {u0}, 쌍 {v, v+u0} 2^(m-1)-1 개 → k = 2^(m-1) + 1.
  H 는 닫혀 있고, N(H)/H = Z_2 x GL(m-1,2) (Z_2 는 u0 만큼의 평행이동으로 {0} 과 {u0} 를 맞바꾼다, GL(m-1,2) 는 F_2^m/<u0> 의 0 이 아닌 점 위의 쌍들을 치환한다).
  GL(m-1,2) 가 그 점들 위의 함수 가운데 자유로운 것을 갖지 않으면 H 는 나타나지 않는다 (m=3: GL(2,2) 가 점 3개, m=4: GL(3,2) 가 점 7개, m=5: GL(4,2) 가 점 15개.
  자유로운 함수가 없다는 것은 scripts/check_gl_free_functions.py 가 센다. GL(5,2) 에는 자유로운 함수가 있어서 m=6 (n=63) 은 이 판정이 말하지 못한다).
  m=2 (n=3) 도 이 판정이 말하지 못하지만 노트 2.6절의 가장 작은 예 <(x2 x3), x -> x xor 111> 이다.

r=2,3 (centre 는 m=2,3) 은 B_n 전체를 표로 만들어 (np_cegar.Cube) 2^k 가지 색칠 전부의 안정화군 위수를 정확히 구한다.
더 큰 매개변수는 B_n 을 표로 만들 수 없으므로 궤도 구조와 위수를 직접 확인하고 (닫힘은 --nauty, 정규화군은 --gap), 안 나타남은 위의 계수 판정으로 닫는다.

사용: python3 -I scripts/check_hamming_transvection_family.py [--family axis|centre] [--r 2 3 4] [--gap /path/to/gap_plain.sh] [--nauty]
      (--nauty 는 pynauty 가 필요하다: pip install pynauty)
      axis: r=2,3 이 약 25초, r=4 가 약 35초(nauty 포함, GAP 의 N(H) 계산은 1초 미만).
      centre: m=2,3 이 약 25초, m=4 (n=15) 가 8~13분 (대부분 nauty 가 Q_15 를 다루는 시간; 다른 계산과 코어를 나누면 13분쯤).
      저장된 출력: notes/check_hamming_transvection_family_axis_output.txt, notes/check_hamming_transvection_family_centre_output.txt
      (명령: --family axis --r 2 3 4 --gap <gap_plain.sh> --nauty 와 --family centre --r 2 3 4 --gap <gap_plain.sh> --nauty).
      axis 의 r=5 (n=32) 는 B_32 를 다룰 수 없어 이 스크립트가 하지 않는다. 그 경우는 노트 3절의 증명과 scripts/agl52_orbits.c 에 기댄다.

2026-11-15. Claude가 쓴 코드이고 저자 검토 전이다.
"""
import argparse
import os
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np


def kernel_basis(cols, n):
    """F_2^m 의 벡터 h_0..h_{n-1} (정수 비트마스크) 에 대해 {x : sum x_i h_i = 0} 의 기저 (n비트 정수 목록)."""
    # 가우스 소거: 각 i 에 대해 (h_i, e_i) 쌍을 줄여 h 가 0 이 되는 조합을 모은다
    basis = {}            # 선행 비트 -> (값, 조합)
    ker = []
    for i in range(n):
        v, comb = cols[i], 1 << i
        while v:
            lead = v.bit_length() - 1
            if lead in basis:
                bv, bc = basis[lead]
                v ^= bv
                comb ^= bc
            else:
                basis[lead] = (v, comb)
                break
        if v == 0:
            ker.append(comb)
    return ker


def vp_translation(n, t):
    X = np.arange(1 << n, dtype=np.int64)
    return (X ^ t).astype(np.int64)


def vp_coord_perm(n, pi):
    """좌표 순열 pi (i -> pi[i]) 가 만드는 꼭짓점 순열: x 의 i번 비트가 pi[i]번 비트로 간다."""
    X = np.arange(1 << n, dtype=np.int64)
    Y = np.zeros_like(X)
    for i in range(n):
        Y |= ((X >> i) & 1) << pi[i]
    return Y


def family_generators(r):
    n = 1 << r
    m = r + 1
    cols = [i + (1 << r) for i in range(n)]
    T = kernel_basis(cols, n)
    assert len(T) == n - m, (len(T), n - m)
    gens = [vp_translation(n, t) for t in T]
    perms = []
    for j in range(r):
        u = 1 << j
        pi = [i ^ u for i in range(n)]
        perms.append(pi)
        gens.append(vp_coord_perm(n, pi))
    return n, m, T, perms, gens


def centre_family_generators(m):
    """중심 고정 족 (m >= 2, n = 2^m - 1): S = F_2^m 의 0 이 아닌 점 전부 (열 h_i = i+1), T = 해밍 부호의 평행이동군,
    K = { v -> v + psi(v) u0 : psi(u0) = 0 } (중심 u0 = 1 을 고정하고 축이 달라지는 전단들, 위수 2^(m-1)) 를 좌표 순열로 옮긴 군.
    H = T x| K~ (위수 2^(n-m) 2^(m-1) = 2^(n-1))."""
    n = (1 << m) - 1
    cols = [i + 1 for i in range(n)]
    T = kernel_basis(cols, n)
    assert len(T) == n - m, (len(T), n - m)
    gens = [vp_translation(n, t) for t in T]
    perms = []
    for j in range(1, m):                       # psi_j(v) = v 의 j번 비트; v -> v xor (psi_j(v) * 1)
        pi = [((i + 1) ^ (((i + 1) >> j) & 1)) - 1 for i in range(n)]
        assert sorted(pi) == list(range(n))
        perms.append(pi)
        gens.append(vp_coord_perm(n, pi))
    return n, m, T, perms, gens


FAMILIES = {
    # 이름: (생성기, k 의 예상값(n, 매개변수), 설명)
    "axis": (family_generators, lambda n, r: n + 1, "해밍-전단 족 (축 고정, n = 2^r, k = n + 1)"),
    "centre": (centre_family_generators, lambda n, m: (1 << (m - 1)) + 1, "중심 고정 족 (n = 2^m - 1, k = 2^(m-1) + 1)"),
}


def components(gens, V):
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    rows = np.concatenate([np.arange(V, dtype=np.int64) for _ in gens])
    cols = np.concatenate([g.astype(np.int64) for g in gens])
    M = coo_matrix((np.ones(len(rows), dtype=np.int8), (rows, cols)), shape=(V, V))
    k, lab = connected_components(M, directed=False)
    return int(k), lab.astype(np.int64)


def nauty_aut_order(n, lab):
    """꼭짓점 색칠 lab (길이 2^n) 을 입힌 정육면체 Q_n 의 자동동형군의 위수 (pynauty). Aut(Q_n) = B_n 이므로
    이것은 B_n 안에서 이 색칠(곧 이름표 함수)을 보존하는 군의 위수이다."""
    from pynauty import Graph, autgrp
    V = 1 << n
    adj = {v: [v ^ (1 << i) for i in range(n) if (v >> i) & 1 == 0] for v in range(V)}   # 변마다 한 번만 적는다
    classes = {}
    for v in range(V):
        classes.setdefault(int(lab[v]), set()).add(v)
    g = Graph(V, directed=False, adjacency_dict=adj, vertex_coloring=[classes[c] for c in sorted(classes)])
    gens, m1, m2, orbits, norb = autgrp(g)
    return int(round(m1 * 10 ** m2))


def gap_normalizer(gap, n, T, perms):
    """B_n = Z_2 wr S_n 을 2n 점 위의 순열군으로 놓고 H 와 N(H) 를 GAP 으로 계산한다 (점 i 와 i+n 이 좌표 i 의 값 0/1)."""
    def cyc_swap(t):
        cs = [f"({i+1},{i+1+n})" for i in range(n) if (t >> i) & 1]
        return "".join(cs) if cs else "()"

    def cyc_perm(pi):
        img = [pi[i] + 1 for i in range(n)] + [pi[i] + 1 + n for i in range(n)]
        return "PermList(" + str(img) + ")"

    lines = [f"n := {n};",
             "cyc := PermList(Concatenation(Concatenation([2..n],[1]), Concatenation([n+2..2*n],[n+1])));",
             "tr := PermList(Concatenation([2,1],[3..n],[n+2,n+1],[n+3..2*n]));",
             "sw := (1,n+1);",
             "B := Group(cyc, tr, sw);",
             "Print(\"|B_n| = \", Size(B), \"\\n\");"]
    gens = [cyc_swap(t) for t in T] + [cyc_perm(pi) for pi in perms]
    lines.append("H := Group(" + ", ".join(gens) + ");")
    lines.append("Print(\"|H| = \", Size(H), \"  (예상 2^(n-1) = \", 2^(n-1), \")\\n\");")
    lines.append("Nn := Normalizer(B, H);")
    lines.append("Print(\"|N(H)| = \", Size(Nn), \"  |N(H)/H| = \", Size(Nn)/Size(H), \"\\n\");")
    lines.append("QUIT;")
    src = "\n".join(lines) + "\n"
    with tempfile.NamedTemporaryFile("w", suffix=".g", delete=False) as f:
        f.write(src)
        path = f.name
    try:
        out = subprocess.run([gap, "-A", "-q", "-o", "3g", path], capture_output=True, text=True, timeout=600)
        return out.stdout + out.stderr
    except subprocess.TimeoutExpired:
        return "GAP 시간 초과 (600 s)"
    finally:
        os.unlink(path)


def run_small(r, nauty=False, fam="axis"):
    """B_n 전체를 표로 만들어 정확히 확인 (n <= 8)."""
    from np_cegar import Cube, appears
    gen_fn, k_fn, desc = FAMILIES[fam]
    n, m, T, perms, gens_vp = gen_fn(r)
    cube = Cube(n)
    gens = [g.astype(cube.dt) for g in gens_vp]
    H = cube.closure(gens)
    order = len(H)
    lab, k = cube.orbit_labels(gens)
    sizes = np.bincount(lab)
    print(f"[{desc}, 매개변수 {r}] n={n}, m={m}, |T|=2^{len(T)}, |H|={order} (= 2^{n-1}: {order == 2**(n-1)}), 꼭짓점 궤도 k={k} (예상 {k_fn(n, r)})")
    print(f"       궤도 길이(길이, 개수): {sorted((int(s), int((sizes == s).sum())) for s in set(sizes.tolist()))}")
    assert k == k_fn(n, r)
    # 닫힘: 궤도 번호 함수의 안정화군 = H ?
    st = cube.count_stab(lab)
    print(f"       닫힘: |Stab(궤도 번호 함수)| = {st}  {'== |H| (닫힘)' if st == order else '!= |H| (닫혀 있지 않음)'}")
    assert st == order
    if nauty:
        an = nauty_aut_order(n, lab)
        print(f"       nauty: |Aut(색칠한 정육면체)| = {an}  {'== |H|' if an == order else '!= |H|'}")
        assert an == order
    # 2^k 가지 색칠 전부
    ratios = {}
    hit = 0
    for mask in range(1 << k):
        col = np.array([(mask >> j) & 1 for j in range(k)], dtype=np.int64)
        f = col[lab]
        s_ = cube.count_stab(f)
        assert s_ % order == 0
        ratios[s_ // order] = ratios.get(s_ // order, 0) + 1
        if s_ == order:
            hit += 1
    print(f"       2^{k} = {1 << k} 가지 색칠 중 안정화군 위수가 |H| 인 것: {hit} 개  → {'나타남' if hit else '나타나지 않음'}")
    print(f"       |Stab(f_c)|/|H| 의 분포 {{비: 색칠 수}}: {dict(sorted(ratios.items()))}")
    assert hit == 0
    # 증거-유도 판정기 (SAT) 로도 확인
    res = appears(cube, gens, order=order)
    print(f"       np_cegar.appears → status={res['status']}, k={res['k']}")
    assert res['status'] == 'absent'
    return n, order, k


def run_large(r, gap=None, nauty=False, fam="axis"):
    gen_fn, k_fn, desc = FAMILIES[fam]
    n, m, T, perms, gens_vp = gen_fn(r)
    V = 1 << n
    k, lab = components(gens_vp, V)
    sizes = np.bincount(lab)
    print(f"[{desc}, 매개변수 {r}] n={n}, m={m}, |T|=2^{len(T)}, 생성원 {len(gens_vp)} 개, |H| = 2^{len(T) + len(perms)} = 2^{n-1}, 꼭짓점 궤도 k={k} (예상 {k_fn(n, r)})")
    print(f"       궤도 길이(길이, 개수): {[(int(s), int((sizes == s).sum())) for s in sorted(set(sizes.tolist()))]}")
    assert k == k_fn(n, r)
    from math import prod
    if fam == "axis":
        q = (1 << r) * prod((1 << r) - (1 << i) for i in range(r))
        print(f"       |AGL({r},2)| = {q},  2^(2^{r}) = {2 ** (1 << r)}  → {'|AGL| > 2^(2^r): 자유로운 함수가 있을 수 없다' if q > 2 ** (1 << r) else '계수로는 닫히지 않음'}")
    else:
        q = 2 * prod((1 << (m - 1)) - (1 << i) for i in range(m - 1))
        pts = (1 << (m - 1)) - 1
        print(f"       |Z_2 x GL({m-1},2)| = {q},  GL({m-1},2) 가 작용하는 점 {pts}개 위의 함수 2^{pts} = {2 ** pts} 개 (자유 함수가 있는지는 scripts/check_gl_free_functions.py 가 센다)")
    if nauty:
        an = nauty_aut_order(n, lab)
        print(f"       닫힘 (nauty): |Aut(색칠한 정육면체 Q_{n})| = {an}  {'== |H| = 2^%d (닫힘)' % (n - 1) if an == 2 ** (n - 1) else '!= |H|'}")
        assert an == 2 ** (n - 1)
    if gap:
        print("       GAP:", gap_normalizer(gap, n, T, perms).strip().replace("\n", "\n            "))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--family", choices=sorted(FAMILIES), default="axis", help="axis: 해밍-전단 족 (n = 2^r), centre: 중심 고정 족 (n = 2^m - 1)")
    ap.add_argument("--r", type=int, nargs="+", default=None, help="매개변수 (axis 는 r, centre 는 m)")
    ap.add_argument("--gap", default=None, help="gap_plain.sh 경로 (주면 N(H)/H 를 GAP 으로 계산)")
    ap.add_argument("--nauty", action="store_true", help="닫힘을 pynauty (색칠한 정육면체의 자동동형군) 로도 확인")
    args = ap.parse_args()
    params = args.r if args.r else ([2, 3, 4] if args.family == "axis" else [2, 3, 4])
    for r in params:
        n = (1 << r) if args.family == "axis" else (1 << r) - 1
        if n <= 8:
            run_small(r, args.nauty, args.family)
            if args.gap:
                n_, m_, T, perms, _ = FAMILIES[args.family][0](r)
                print("       GAP:", gap_normalizer(args.gap, n_, T, perms).strip().replace("\n", "\n            "))
        else:
            run_large(r, args.gap, args.nauty, args.family)
    print("결과: 끝")


if __name__ == "__main__":
    main()
