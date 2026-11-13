#!/usr/bin/env python3
"""
circuit_sample_orbits.py
========================
같은 NPN 동치류(입력 순열·입력 반전·출력 반전으로 옮겨지는 함수들)의 함수를 도구에 넣으면 도구의 값이 얼마나 달라지는지를
표본으로 잰다 (PLAN.md 6절 C8의 (나), 5변수로 넓히기). 4변수는 65,536개 전수 측정이 이미 있으므로
(circuit_metrics.py) 같은 스크립트를 4변수에도 돌려 두 방법이 같은 결론을 내는지 확인하는 데 쓴다.

하는 일
  1. 무작위 함수 M개를 뽑는다 (모든 변수에 의존하는 것만, 시드 고정).
  2. 함수마다 같은 동치류의 함수 K개를 만든다: 자기 자신과 무작위 NPN 변환 K-1개.
     (변환은 군 B_n x Z_2 에서 고르게 뽑으므로 각각이 그 궤도의 무작위 원소이다. 대칭이 있는 함수에서는 겹칠 수 있다.)
  3. 그 K개와, 동치류의 표준형(진리표가 가장 작은 함수) 하나를 ABC 로 합성해 비용을 잰다.
       and_A : 흐름 A (strash 뒤 resyn2 세 번)의 2입력 AND 수
       and_B : 흐름 B (strash, dc2, compress2rs 두 번, dc2)의 2입력 AND 수
       l2, l3, l4 : 흐름 B 의 결과를 `dch -f; if -K k -a` 로 옮긴 k입력 LUT 수 (입력이 둘 이상인 노드만 센다)
     circuit_metrics.py 와 같은 명령이다. 흐름 B 의 `and` 는 4변수 전수 측정에서 동치류 안에서 늘 일정했다.
  4. 함수마다 B_n 안의 안정화군 위수와 NPN 군 안의 안정화군 위수도 적는다 (stabilizer_types.g 의 표와 이어 보려고).

실행: python scripts/circuit_sample_orbits.py run <n: 4 또는 5> <M> <K> <시드> <출력 파일> [<나눌 수> <나머지>] [-abc <경로>] [-kind <종류>]
      -kind npn   (기본) 변환 K-1개를 군 B_n x Z_2 에서 고른다.
      -kind perm  변수 순열만 쓴다 (입력·출력 반전 없음).   -kind neg  입력 반전만 쓴다 (변수 순열·출력 반전 없음).
      -kind order 같은 함수를 최소항 순서만 바꿔 쓴다. ABC 는 이 순서에 반응하지 않으므로 (BLIF 를 읽는 순간 같은 AIG 가 된다) 값이 모두 같고,
                  그래서 "동치류 안의 값이 다른 것은 입력 이름과 극성이 다르기 때문"임을 보이는 점검으로만 쓴다.
      python scripts/circuit_sample_orbits.py analyze <출력 파일> [<출력 파일> ...]
      ABC 실행 파일은 환경 변수 CIRCUIT_ABC 로도 줄 수 있다. 마지막 두 인자를 주면 함수 번호가 (번호 % 나눌 수 == 나머지)인 것만 한다.
출력: 한 줄에 함수 하나. "번호 진리표 표준형 |Stab_B| |Stab_NPN|" 뒤에 비용 다섯 가지가 이어진다.
      각 비용은 "변환 K개의 값(쉼표로 이음);표준형의 값" 꼴이다. 맨 앞의 변환은 항등이다 (함수 자신).
준비: pip install numpy. ABC 는 소스에서 빌드한다 (circuit_exact.py 머리말 참고).
"""
import itertools
import os
from math import comb, factorial
import random
import subprocess
import sys
import tempfile

import numpy as np

RESYN2 = ["balance", "rewrite", "refactor", "balance", "rewrite", "rewrite -z", "balance",
          "refactor -z", "rewrite -z", "balance"]
COMPRESS2RS = ["balance -l", "resub -K 6 -l", "rewrite -l", "resub -K 6 -N 2 -l", "refactor -l",
               "resub -K 8 -l", "balance -l", "resub -K 8 -N 2 -l", "rewrite -l", "resub -K 10 -l",
               "rewrite -z -l", "resub -K 10 -N 2 -l", "balance -l", "resub -K 12 -l", "refactor -z -l",
               "resub -K 12 -N 2 -l", "rewrite -z -l", "balance -l"]
FLOW_A = RESYN2 * 3
FLOW_B = ["dc2"] + COMPRESS2RS * 2 + ["dc2"]
METRICS = ["and_A", "and_B", "l2", "l3", "l4"]
BATCH = 100          # ABC 한 번에 맡기는 함수 수


class Space:
    """n변수 함수의 공간과 NPN 변환."""

    def __init__(self, n):
        self.n = n
        self.size = 1 << n                       # 진리표의 길이
        self.full = (1 << self.size) - 1
        maps = []
        for perm in itertools.permutations(range(n)):
            for c in range(self.size):
                row = []
                for x in range(self.size):
                    y = 0
                    for i in range(n):           # y 의 i번째 비트 = (x ^ c) 의 perm[i]번째 비트
                        y |= (((x ^ c) >> perm[i]) & 1) << i
                    row.append(y)
                maps.append(row)
        self.maps = np.array(maps, dtype=np.int64)       # (n! 2^n, 2^n)
        self.pow2 = (np.uint64(1) << np.arange(self.size, dtype=np.uint64))

    def bits(self, f):
        return np.array([(f >> x) & 1 for x in range(self.size)], dtype=np.uint64)

    def images(self, f):
        """B_n 의 모든 원소가 f 를 옮긴 함수(정수)의 배열."""
        return (self.bits(f)[self.maps] * self.pow2).sum(axis=1)

    def act(self, f, g, e):
        h = int(self.images_one(f, g))
        return h ^ self.full if e else h

    def images_one(self, f, g):
        return (self.bits(f)[self.maps[g]] * self.pow2).sum()

    def canonical_and_stabilizers(self, f):
        im = self.images(f)
        comp = im ^ np.uint64(self.full)
        canon = int(min(im.min(), comp.min()))
        stab_b = int((im == np.uint64(f)).sum())
        stab_npn = stab_b + int((comp == np.uint64(f)).sum())
        return canon, stab_b, stab_npn

    def support_size(self, f):
        return sum(1 for i in range(self.n)
                   if any(((f >> x) & 1) != ((f >> (x ^ (1 << i))) & 1) for x in range(self.size)))


def blif_text(space, f, shuffle_seed=None):
    names = " ".join("x%d" % i for i in range(space.n))
    lines = [".model f", ".inputs " + names, ".outputs y", ".names " + names + " y"]
    cubes = ["".join(str((x >> i) & 1) for i in range(space.n)) + " 1" for x in range(space.size) if (f >> x) & 1]
    if shuffle_seed is not None:
        random.Random(shuffle_seed).shuffle(cubes)       # 같은 함수, 최소항을 늘어놓는 순서만 다르다
    lines += cubes
    lines.append(".end")
    return "\n".join(lines) + "\n"


def synthesize(space, funcs, abc):
    """funcs: (진리표(정수), 최소항 순서의 시드 또는 None) 목록. 각각에 대해 METRICS 의 값을 돌려준다 (사전의 목록)."""
    out = []
    for start in range(0, len(funcs), BATCH):
        chunk = funcs[start:start + BATCH]
        with tempfile.TemporaryDirectory() as d:
            os.mkdir(os.path.join(d, "in"))
            os.mkdir(os.path.join(d, "out"))
            cmds = []
            for j, (f, sh) in enumerate(chunk):
                with open(os.path.join(d, "in", "%d.blif" % j), "w") as fh:
                    fh.write(blif_text(space, f, sh))
                for tag, flow in (("A", FLOW_A), ("B", FLOW_B)):
                    cmds += ["read_blif in/%d.blif" % j, "strash"] + flow + ["write_aiger out/%d.%s.aig" % (j, tag)]
                for k in (2, 3, 4):
                    cmds += ["read_aiger out/%d.B.aig" % j, "dch -f", "if -K %d -a" % k,
                             "write_blif out/%d.l%d.blif" % (j, k)]
            with open(os.path.join(d, "s.abc"), "w") as fh:
                fh.write("\n".join(cmds) + "\n")
            subprocess.run([os.path.abspath(abc), "-q", "source s.abc"], cwd=d, stdin=subprocess.DEVNULL,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            for j in range(len(chunk)):
                row = {}
                for tag in ("A", "B"):
                    with open(os.path.join(d, "out", "%d.%s.aig" % (j, tag)), "rb") as fh:
                        row["and_" + tag] = int(fh.readline().split()[5])
                for k in (2, 3, 4):
                    luts = 0
                    for line in open(os.path.join(d, "out", "%d.l%d.blif" % (j, k))):
                        if line.startswith(".names") and len(line.split()) >= 4:
                            luts += 1
                    row["l%d" % k] = luts
                out.append(row)
    return out


def run(argv):
    n, m_funcs, k_var, seed, out_path = int(argv[0]), int(argv[1]), int(argv[2]), int(argv[3]), argv[4]
    rest = argv[5:]
    abc = os.environ.get("CIRCUIT_ABC")
    if "-abc" in rest:
        i = rest.index("-abc")
        abc = rest[i + 1]
        rest = rest[:i] + rest[i + 2:]
    kind = "npn"
    if "-kind" in rest:
        i = rest.index("-kind")
        kind = rest[i + 1]
        rest = rest[:i] + rest[i + 2:]
    if kind not in ("npn", "order", "perm", "neg"):
        sys.exit("-kind 는 npn, order, perm, neg 중 하나")
    parts, part = (int(rest[0]), int(rest[1])) if len(rest) >= 2 else (1, 0)
    if not abc:
        sys.exit("ABC 실행 파일의 경로가 없다 (-abc 또는 환경 변수 CIRCUIT_ABC)")
    space = Space(n)
    rng = random.Random(seed)
    funcs = []
    while len(funcs) < m_funcs:
        f = rng.getrandbits(space.size)
        if space.support_size(f) == n:
            funcs.append(f)
    group_n = len(space.maps)
    with open(out_path, "w") as out:
        out.write("# n=%d M=%d K=%d seed=%d part=%d/%d kind=%s\n" % (n, m_funcs, k_var, seed, part, parts, kind))
        batch, meta = [], []
        for idx, f in enumerate(funcs):
            vrng = random.Random(seed * 1000003 + idx)       # 변환은 함수 번호만으로 정해진다 (나눠 돌려도 같다)
            if idx % parts != part:
                continue
            canon, stab_b, stab_npn = space.canonical_and_stabilizers(f)
            if kind == "npn":
                variants = [(f, None)]
                for _ in range(k_var - 1):
                    variants.append((space.act(f, vrng.randrange(group_n), vrng.getrandbits(1)), None))
            elif kind in ("perm", "neg"):           # 군의 한 부분만 쓴다: 변수 순열만, 또는 입력 반전만 (출력 반전은 쓰지 않는다)
                variants = [(f, None)]
                for _ in range(k_var - 1):
                    g = vrng.randrange(factorial(n)) * space.size if kind == "perm" else vrng.randrange(space.size)
                    variants.append((space.act(f, g, 0), None))
            else:                                   # 같은 함수를 최소항 순서만 바꿔 K번 (ABC 는 이 순서에 반응하지 않는다)
                variants = [(f, None)] + [(f, vrng.getrandbits(60)) for _ in range(k_var - 1)]
            batch.append(variants + [(canon, None)])
            meta.append((idx, f, canon, stab_b, stab_npn))
            if len(batch) * (k_var + 1) >= BATCH:
                flush(space, abc, batch, meta, k_var, out)
                batch, meta = [], []
        if batch:
            flush(space, abc, batch, meta, k_var, out)


def flush(space, abc, batch, meta, k_var, out):
    flat = [t for variants in batch for t in variants]
    res = synthesize(space, flat, abc)
    pos = 0
    for variants, (idx, f, canon, stab_b, stab_npn) in zip(batch, meta):
        rows = res[pos:pos + len(variants)]
        pos += len(variants)
        cells = []
        for metric in METRICS:
            vals = [r[metric] for r in rows]
            cells.append(",".join(str(v) for v in vals[:-1]) + ";" + str(vals[-1]))
        out.write("%d %x %x %d %d %s\n" % (idx, f, canon, stab_b, stab_npn, " ".join(cells)))
    out.flush()


def load(paths):
    data = []
    for path in paths:
        for line in open(path):
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            row = {"idx": int(p[0]), "f": int(p[1], 16), "canon": int(p[2], 16), "stab_b": int(p[3]),
                   "stab_npn": int(p[4])}
            for metric, cell in zip(METRICS, p[5:]):
                vs, c = cell.split(";")
                row[metric] = ([int(v) for v in vs.split(",")], int(c))
            data.append(row)
    data.sort(key=lambda r: r["idx"])
    return data


def exp_min_of_k(costs, k):
    """K개 가운데 무작위로 (중복 없이) k개를 고를 때 최솟값의 기댓값."""
    c = sorted(costs)
    kk = len(c)
    return sum(c[i] * comb(kk - 1 - i, k - 1) for i in range(kk)) / comb(kk, k)


def mean_se(xs):
    m = len(xs)
    mu = sum(xs) / m
    sd = (sum((x - mu) ** 2 for x in xs) / (m - 1)) ** 0.5 if m > 1 else 0.0
    return mu, sd / m ** 0.5


def analyze(paths):
    data = load(paths)
    if not data:
        sys.exit("읽은 줄이 없다")
    m = len(data)
    k = len(data[0]["and_B"][0])
    print("함수 %d개, 함수마다 변환 %d개(자기 자신 포함)와 표준형 하나" % (m, k))
    triv = sum(1 for r in data if r["stab_b"] == 1)
    print("B_n 안의 안정화군이 자명한 함수 %d개 (%.1f%%), NPN 군 안의 안정화군이 자명한 함수 %d개" %
          (triv, 100.0 * triv / m, sum(1 for r in data if r["stab_npn"] == 1)))
    ks = [kk for kk in (1, 2, 4, 8, 16, 32, 64, 128, 256) if kk <= k]
    for metric in METRICS:
        base = [r[metric][0][0] for r in data]                          # 함수 자신의 값
        mean_all = sum(sum(r[metric][0]) / k for r in data) / m           # 궤도의 무작위 원소의 값(평균)
        const = sum(1 for r in data if len(set(r[metric][0] + [r[metric][1]])) == 1)
        const_vars = sum(1 for r in data if len(set(r[metric][0])) == 1)
        canon = [r[metric][1] for r in data]
        print("\n[%s]" % metric)
        print("  함수 자신의 값 평균 %.3f, 변환 %d개 전체의 평균 %.3f, 표준형의 값 평균 %.3f" %
              (sum(base) / m, k, mean_all, sum(canon) / m))
        print("  변환 %d개의 값이 모두 같은 함수 %d개 (%.1f%%); 표준형의 값까지 모두 같은 함수 %d개 (%.1f%%)" %
              (k, const_vars, 100.0 * const_vars / m, const, 100.0 * const / m))
        print("  궤도의 무작위 원소 하나를 합성하는 것에 견준 값 (함수마다 변환 %d개에서 무작위로 k개를 고를 때 최솟값의 기댓값):" % k)
        for kk in ks:
            gains = [(sum(r[metric][0]) / k) - exp_min_of_k(r[metric][0], kk) for r in data]
            best = sum(exp_min_of_k(r[metric][0], kk) for r in data) / m
            mu, se = mean_se(gains)
            print("    k = %2d: 최솟값의 평균 %.3f  (한 번 합성보다 %.3f 작음, 표준오차 %.3f, 상대 %.2f%%)" %
                  (kk, best, mu, se, 100.0 * mu / mean_all))
        print("    표준형만 합성: 평균 %.3f  (한 번 합성보다 %.2f%% 작음; 음수면 큼)" %
              (sum(canon) / m, 100.0 * (mean_all - sum(canon) / m) / mean_all))
        diffs = [(sum(r[metric][0]) / k) - r[metric][1] for r in data]
        mu, se = mean_se(diffs)
        print("    함수마다 (변환 평균 - 표준형) 의 평균 %.3f, 표준오차 %.3f" % (mu, se))
        spread = {}
        for r in data:
            d = max(r[metric][0]) - min(r[metric][0])
            spread[d] = spread.get(d, 0) + 1
        print("    변환 %d개의 (최댓값 - 최솟값) 분포: %s" % (k, dict(sorted(spread.items()))))


if __name__ == "__main__":
    if len(sys.argv) >= 2 and sys.argv[1] == "run":
        run(sys.argv[2:])
    elif len(sys.argv) >= 3 and sys.argv[1] == "analyze":
        analyze(sys.argv[2:])
    else:
        sys.exit(__doc__)
