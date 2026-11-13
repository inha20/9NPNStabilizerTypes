#!/usr/bin/env python3
"""
check_pictures_by_definition.py
===============================
그림 종류 수 T(f)를 정의대로 센다. 군론(안정화군, 표지 표, 명제 A2)을 쓰지 않는다.

  1. 함수 f 를 n! 가지 변수 배치 각각으로 지도에 그린다 (모델 A: 변수를 둘씩 묶어 4-순환 축에
     Gray 순서 00, 01, 11, 10 으로 놓고, n 이 홀수이면 마지막 변수를 층으로 둔다).
  2. 격자 대칭(축의 한 칸 이동, 축 뒤집기, 축 맞바꿈, 층 뒤집기가 만드는 군)으로 겹치는 그림끼리 묶는다.
  3. 묶음의 수가 T(f) 이다.

쓰임:
  python scripts/check_pictures_by_definition.py notes/class_representatives_n5.txt
      export_class_representatives.g 가 내보낸 함수들에서 T 를 다시 세어 예측값과 대조한다.
  python scripts/check_pictures_by_definition.py --selftest4
      n = 4 의 65,536개 함수 전부에서 T 를 세어 7번 논문의 3,220 / 36,204 / 26,112 와 비교한다.

두 경우 모두 먼저 명제 A1 을 확인한다: 격자 대칭군의 원소가 모두 "기준 쌍 분할을 보존하는 부호 순열"이고
그 개수가 2^n * 2^k * k! (k = n 을 2 로 나눈 몫) 이다.

표준 라이브러리만 쓴다.
"""
import itertools
import sys
from collections import Counter
from math import factorial

GRAY = [(0, 0), (0, 1), (1, 1), (1, 0)]


class Model:
    def __init__(self, n):
        self.n = n
        k, odd = divmod(n, 2)
        self.k, self.odd = k, odd
        cells = list(itertools.product(range(4), repeat=k))
        if odd:
            cells = [c + (l,) for c in cells for l in (0, 1)]
        self.cells = cells
        index = {c: i for i, c in enumerate(cells)}

        def point(c):
            x = 0
            for a in range(k):
                g1, g2 = GRAY[c[a]]
                x |= (g1 << (2 * a)) | (g2 << (2 * a + 1))
            if odd:
                x |= c[k] << (n - 1)
            return x

        self.base = [point(c) for c in cells]  # 기준 배치에서 셀 -> 입력점

        # 격자 대칭의 생성원 (셀의 치환). 다른 축의 이동·뒤집기는 축 맞바꿈과의 켤레로 얻어진다.
        fns = [
            lambda c: ((c[0] + 1) % 4,) + c[1:],
            lambda c: ((-c[0]) % 4,) + c[1:],
        ]
        for a in range(k - 1):
            fns.append(lambda c, a=a: c[:a] + (c[a + 1], c[a]) + c[a + 2:])
        if odd:
            fns.append(lambda c: c[:k] + (1 - c[k],))
        self.gens = [tuple(index[g(c)] for c in cells) for g in fns]
        self.tables = [self._byte_tables(g) for g in self.gens]

        # 배치 s 에서 셀 c 에 놓이는 입력점: 기준 입력점의 좌표를 s 로 섞는다
        self.arrangements = []
        for s in itertools.permutations(range(n)):
            self.arrangements.append(
                [sum(((xb >> s[i]) & 1) << i for i in range(n)) for xb in self.base]
            )

    def _byte_tables(self, g):
        tabs = []
        for b in range(len(self.cells) // 8):
            t = [0] * 256
            for v in range(256):
                r = 0
                for j in range(8):
                    if (v >> j) & 1:
                        r |= 1 << g[8 * b + j]
                t[v] = r
            tabs.append(t)
        return tabs

    @staticmethod
    def _apply(tabs, p):
        r = 0
        shift = 0
        for t in tabs:
            r |= t[(p >> shift) & 255]
            shift += 8
        return r

    def check_A1(self):
        """격자 대칭군을 전부 만들어, 원소가 모두 p0 를 보존하는 부호 순열이고 개수가 맞는지 본다."""
        n, k, odd = self.n, self.k, self.odd
        ident = tuple(range(len(self.cells)))
        group = {ident}
        stack = [ident]
        while stack:
            h = stack.pop()
            for g in self.gens:
                hg = tuple(g[h[c]] for c in range(len(h)))
                if hg not in group:
                    group.add(hg)
                    stack.append(hg)
        pair_of = {}
        for a in range(k):
            pair_of[2 * a] = pair_of[2 * a + 1] = a
        if odd:
            pair_of[n - 1] = "single"
        negations = set()
        for h in group:
            pt = [0] * (1 << n)
            for c, xb in enumerate(self.base):
                pt[xb] = self.base[h[c]]
            c0 = pt[0]
            pi = []
            for i in range(n):
                e = pt[1 << i] ^ c0
                if e == 0 or e & (e - 1):
                    return False, "부호 순열이 아닌 원소가 있다"
                pi.append(e.bit_length() - 1)
            for x in range(1 << n):
                y = c0
                for i in range(n):
                    if (x >> i) & 1:
                        y ^= 1 << pi[i]
                if pt[x] != y:
                    return False, "부호 순열이 아닌 원소가 있다"
            for a in range(k):
                if pair_of[pi[2 * a]] != pair_of[pi[2 * a + 1]] or pair_of[pi[2 * a]] == "single":
                    return False, "쌍 분할 p0 를 보존하지 않는 원소가 있다"
            if pi == list(range(n)):
                negations.add(c0)
        expected = (1 << n) * (1 << k) * factorial(k)
        ok = len(group) == expected and len(negations) == (1 << n)
        return ok, "|W| = %d (기대 %d), 입력 반전 %d개 (기대 %d)" % (len(group), expected, len(negations), 1 << n)

    def picture_types(self, tt):
        """진리표 tt 인 함수의 그림 종류 수를 정의대로 센다."""
        pics = set()
        for pts in self.arrangements:
            p = 0
            for c, x in enumerate(pts):
                if (tt >> x) & 1:
                    p |= 1 << c
            pics.add(p)
        count = 0
        apply = self._apply
        tables = self.tables
        while pics:
            start = pics.pop()
            count += 1
            seen = {start}
            stack = [start]
            while stack:
                p = stack.pop()
                for tabs in tables:
                    q = apply(tabs, p)
                    if q not in seen:
                        seen.add(q)
                        stack.append(q)
            pics -= seen
        return count


def selftest4():
    m = Model(4)
    ok, msg = m.check_A1()
    print("[A1] n = 4:", "통과" if ok else "실패", "-", msg)
    dist = Counter(m.picture_types(tt) for tt in range(1 << 16))
    print("[n=4 전수] T = 1/2/3 인 함수:", dist[1], dist[2], dist[3], "(7번 논문: 3220 36204 26112)")
    return ok and (dist[1], dist[2], dist[3]) == (3220, 36204, 26112)


def check_file(path):
    rows = []
    n = None
    for line in open(path, encoding="utf-8"):
        if line.startswith("#"):
            if line.startswith("# n ="):
                n = int(line.split("=")[1])
            continue
        if line.strip():
            j, size, orbits, t, tt = line.split()
            rows.append((int(j), int(size), int(orbits), int(t), int(tt)))
    m = Model(n)
    ok, msg = m.check_A1()
    print("[A1] n = %d:" % n, "통과" if ok else "실패", "-", msg)
    bad = 0
    seen_t = Counter()
    for j, size, orbits, t, tt in rows:
        got = m.picture_types(tt)
        seen_t[got] += 1
        if got != t:
            bad += 1
            print("    불일치: 켤레류 %d, 예측 %d, 정의대로 센 값 %d" % (j, t, got))
    print("[그림 비교] n = %d: 함수 %d개(안정화군 유형마다 하나) 가운데 불일치 %d개" % (n, len(rows), bad))
    print("    정의대로 센 T 값별 유형 수:", dict(sorted(seen_t.items())))
    return ok and bad == 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    good = selftest4() if sys.argv[1] == "--selftest4" else check_file(sys.argv[1])
    print("모두 일치한다." if good else "일치하지 않는 항목이 있다.")
    sys.exit(0 if good else 1)
