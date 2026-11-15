#!/usr/bin/env python3
"""
check_n3_minimal_absent_example.py
==================================
n = 3 에서 닫혔지만 안 나타나는 군 4개 가운데, 정규화군 몫 판정법이 설명하지 못하는 하나(켤레류 번호 31, 위수 4, 꼭짓점 궤도 길이 4, 2, 2)를
저장소의 다른 코드를 쓰지 않고 처음부터 확인한다. B_3 의 원소 48개와 부울 함수 256개를 모두 훑는다. 손으로 따라갈 수 있는 가장 작은 예이다.

확인하는 것
  1) B_3 = { x -> pi(x) xor c } 의 원소 48개. H = < (x2 x3), x -> x xor 111 > 의 위수 4와 꼭짓점 궤도 길이.
  2) 부울 함수 256개 모두의 안정화군을 구해, 안정화군이 정확히 H 인 함수가 없음.
  3) H 로 불변인 함수 8개의 안정화군 위수.
  4) 이유. z = (x -> x xor 111) 은 B_3 의 모든 원소와 가환이므로 B_3 는 몫 {0,1}^3 / <z> (4점)에 작용한다.
     이 작용의 상은 S_4 (위수 24), 핵은 {id, z} 이다. H 는 z 를 품으므로 H 의 상은 호환 하나이다.
     z 로 불변인 함수 f 는 4점 위의 2-색칠이고, f 의 안정화군은 그 색칠의 S_4 안 안정화군의 원상이다.
     4점의 2-색칠의 안정화군은 S_a x S_b (위수 24, 6, 4)뿐이므로 위수 2인 군(호환)은 안 나온다. 그래서 H 는 안 나타난다.

사용:  python3 -I scripts/check_n3_minimal_absent_example.py
시간: 1초 미만.

2026-11-15. Claude가 쓴 코드이고 저자 검토 전이다.
"""
from itertools import permutations, product

N = 3
V = 1 << N                      # 꼭짓점 0..7, 꼭짓점의 i 번째 비트(0부터)가 변수 x_{i+1}
FULL = V - 1                    # 111


def permute_bits(x, p):
    """비트 j 를 비트 p[j] 로 옮긴다."""
    y = 0
    for j in range(N):
        if (x >> j) & 1:
            y |= 1 << p[j]
    return y


def compose(a, b):              # a 다음 b 가 아니라 a(b(x))
    return tuple(a[b[x]] for x in range(V))


def generate(gens):
    ident = tuple(range(V))
    seen = {ident}
    todo = [ident]
    while todo:
        g = todo.pop()
        for h in gens:
            k = compose(h, g)
            if k not in seen:
                seen.add(k)
                todo.append(k)
    return frozenset(seen)


def orbits(group):
    seen, out = set(), []
    for x in range(V):
        if x not in seen:
            o = sorted({g[x] for g in group})
            seen.update(o)
            out.append(o)
    return out


def main():
    Bn = frozenset(tuple(permute_bits(x, p) ^ c for x in range(V)) for p in permutations(range(N)) for c in range(V))
    print(f"1) B_3 의 원소 수 {len(Bn)} (48이어야 한다)")
    assert len(Bn) == 48

    swap23 = tuple(permute_bits(x, (0, 2, 1)) for x in range(V))        # x2 와 x3 을 맞바꾼다
    z = tuple(x ^ FULL for x in range(V))                               # 모든 입력 반전
    H = generate([swap23, z])
    ob = sorted(len(o) for o in orbits(H))
    print(f"   H = <(x2 x3), x -> x xor 111> 의 위수 {len(H)}, 꼭짓점 궤도 길이 {ob}")
    assert len(H) == 4 and ob == [2, 2, 4]

    # 2) 함수 256개의 안정화군
    stab = {}
    for f in range(1 << V):
        bits = [(f >> x) & 1 for x in range(V)]
        stab[f] = frozenset(g for g in Bn if all(bits[g[x]] == bits[x] for x in range(V)))
    exact = [f for f in stab if stab[f] == H]
    print(f"2) 부울 함수 {len(stab)}개 가운데 안정화군이 정확히 H 인 함수: {len(exact)}개 (0이어야 한다)")
    assert not exact

    # 3) H 로 불변인 함수
    inv = [f for f in stab if H <= stab[f]]
    ords = sorted(len(stab[f]) for f in inv)
    print(f"3) H 로 불변인 함수 {len(inv)}개의 안정화군 위수 {ords}")
    assert len(inv) == 8 and ords == [8, 8, 12, 12, 12, 12, 48, 48]

    # 4) 몫 {0,1}^3/<z> 와 S_4
    cls = {}
    for x in range(V):
        cls[x] = min(x, x ^ FULL)
    reps = sorted(set(cls.values()))
    assert len(reps) == 4
    idx = {r: i for i, r in enumerate(reps)}

    def act(g):                 # g 가 4점 위에 하는 순열 (z 와 가환이어서 잘 정의된다)
        out = []
        for r in reps:
            out.append(idx[cls[g[r]]])
        for x in range(V):      # 잘 정의됨 확인
            assert idx[cls[g[x]]] == out[idx[cls[x]]]
        return tuple(out)

    img = {g: act(g) for g in Bn}
    image = set(img.values())
    kernel = [g for g in Bn if img[g] == (0, 1, 2, 3)]
    print(f"4) B_3 가 4점 {{0,1}}^3/<z> 에 하는 작용: 상의 위수 {len(image)} (S_4 이면 24), 핵의 위수 {len(kernel)} (<z> 이면 2)")
    assert len(image) == 24 and len(kernel) == 2 and z in kernel

    Hbar = {img[h] for h in H}
    moved = sorted(i for i in range(4) if any(p[i] != i for p in Hbar))
    print(f"   H 의 상: 위수 {len(Hbar)}, 움직이는 점 {moved} (호환이면 점 둘)")
    assert len(Hbar) == 2 and len(moved) == 2

    # 4점의 2-색칠의 S_4 안 안정화군
    s4_orders = {}
    for col in product((0, 1), repeat=4):
        s = [p for p in image if all(col[p[i]] == col[i] for i in range(4))]
        s4_orders.setdefault(len(s), 0)
        s4_orders[len(s)] += 1
    print(f"   4점의 2-색칠 16개의 S_4 안 안정화군 위수별 개수 {dict(sorted(s4_orders.items()))} (위수 24, 6, 4뿐이다)")
    assert set(s4_orders) == {24, 6, 4}

    # z 로 불변인 함수의 안정화군 = 색칠의 안정화군의 원상
    for f in range(1 << V):
        if z in stab[f]:
            col = [None] * 4
            for x in range(V):
                col[idx[cls[x]]] = (f >> x) & 1
            s = [p for p in image if all(col[p[i]] == col[i] for i in range(4))]
            assert len(stab[f]) == 2 * len(s)
    print("   z 로 불변인 함수마다 (안정화군의 위수) = 2 x (색칠의 S_4 안 안정화군의 위수): 모두 맞다")
    print("끝: H 는 S_4 의 호환의 원상이고, 위수 2인 군은 4점의 2-색칠의 안정화군이 못 되므로 H 는 안 나타난다.")


if __name__ == "__main__":
    main()
