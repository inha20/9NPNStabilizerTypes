#!/usr/bin/env python3
"""
circuit_multi_output.py
=======================
함수 둘을 한 회로로 함께 구현할 때 얼마나 겹쳐 쓰이는지를 표본으로 잰다 (탐색 노트 4.3절).
"궤도의 대표 회로를 옮겨 쓰는" 방식은 함수마다 회로를 따로 두므로 비용이 C(f) + C(g) 이다.
함께 합성한 회로의 크기 C(f, g) 가 이보다 얼마나 작은지가, 그 방식이 여러 함수를 구현할 때 놓치는 몫이다.
같은 궤도의 두 함수가 다른 궤도의 두 함수보다 더 많이 겹치는지도 본다.

비용은 AND 수(흐름 B). C(f), C(g) 는 circuit_metrics.py 의 전수 측정값을 쓴다.
f 는 네 변수에 모두 의존하는 함수 가운데에서 무작위로 고른다. g 는 무리마다 다르게 고른다.

  G1 같은 궤도, 무작위   : g 는 f 에 무작위 입력 순열·입력 반전을 가한 것 (f 와 다른 것).
  G2 다른 궤도, 같은 크기: g 는 f 와 다른 궤도에 있고 C(g) = C(f) 인 함수 가운데 무작위.
  G3 입력 하나 반전      : g(x) = f(x 의 한 비트를 뒤집은 것).
  G4 입력 둘 맞바꿈      : g 는 f 의 두 변수를 맞바꾼 것 (f 와 다른 것).

실행: python scripts/circuit_multi_output.py notes/circuit_metrics.txt.gz [무리마다의 쌍 수, 기본 1000] [seed, 기본 1]
      환경 변수 CIRCUIT_ABC 에 따로 빌드한 ABC 의 경로를 주면 그것을 쓰고, 없으면 Yosys 안의 ABC 를 쓴다.
준비: pip install numpy yowasp-yosys
"""
import os
import random
import subprocess
import sys
import tempfile
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import circuit_metrics as cm
import circuit_metrics_analyze as an

N = 16
FULL = 1 << N


def run_abc(d, cmds):
    native = os.environ.get("CIRCUIT_ABC")
    if native:
        with open(os.path.join(d, "s.abc"), "w") as fh:
            fh.write("\n".join(cmds) + "\n")
        subprocess.run([os.path.abspath(native), "-q", "source s.abc"], cwd=d, stdin=subprocess.DEVNULL,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    else:
        with open(os.path.join(d, "s.abc"), "w") as fh:
            fh.write("\n".join(["backup"] + cmds + ["restore", "strash"]) + "\n")
        with open(os.path.join(d, "d.v"), "w") as fh:
            fh.write("module dummy(input a, input b, output y); assign y = a & b; endmodule\n")
        with open(os.path.join(d, "s.ys"), "w") as fh:
            fh.write("read_verilog d.v\nproc\nopt\ntechmap\nabc -script s.abc\n")
        subprocess.run(["yowasp-yosys", "-q", "-s", "s.ys"], cwd=d)


def joint_sizes(function_sets):
    """함수 묶음마다, 그 함수들을 출력으로 하는 회로 하나를 흐름 B 로 최적화한 AND 수."""
    out = []
    with tempfile.TemporaryDirectory() as d:
        os.mkdir(os.path.join(d, "in"))
        os.mkdir(os.path.join(d, "out"))
        cmds = []
        for k, fs in enumerate(function_sets):
            with open(os.path.join(d, "in", "%d.blif" % k), "w") as fh:
                fh.write(".model m\n.inputs x0 x1 x2 x3\n.outputs %s\n" % " ".join("y%d" % j for j in range(len(fs))))
                for j, f in enumerate(fs):
                    fh.write(".names x0 x1 x2 x3 y%d\n" % j)
                    for x in range(N):
                        if (f >> x) & 1:
                            fh.write("%d%d%d%d 1\n" % (x & 1, (x >> 1) & 1, (x >> 2) & 1, (x >> 3) & 1))
                fh.write(".end\n")
            cmds += ["read_blif in/%d.blif" % k, "strash"] + cm.FLOWS["B"] + ["write_aiger out/%d.aig" % k]
        run_abc(d, cmds)
        for k in range(len(function_sets)):
            with open(os.path.join(d, "out", "%d.aig" % k), "rb") as fh:
                out.append(int(fh.readline().split()[5]))
    return out


def main():
    vals = an.load(sys.argv[1])
    single = vals[("B", "and")]
    pairs_per_group = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
    rng = random.Random(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
    canon = an.canon_b4
    full = [f for f in range(FULL) if len(cm.support(f)) == 4]
    by_cost = {}
    for f in full:
        by_cost.setdefault(int(single[f]), []).append(f)
    perms = an.PERMS                      # (384, 65536): 변환마다 함수가 어디로 가는지
    weights = an.WEIGHTS
    neg_one = [k for k in range(len(weights)) if weights[k] == 1]
    # 반전 없이 두 변수만 맞바꾸는 변환: 무게 0 이고, 자기 자신과 합성하면 항등이며, 항등이 아닌 것 가운데 호환(6개)
    ident = [k for k in range(len(weights)) if weights[k] == 0 and (perms[k] == an.F).all()][0]
    zero = [k for k in range(len(weights)) if weights[k] == 0 and k != ident]
    transpositions = []
    for k in zero:
        p = perms[k].astype("int64")
        if (p[p] == an.F).all():
            # 호환과 (호환 둘의 곱)을 가른다: 변수 하나짜리 함수 x_i (진리표) 가 움직이는 개수가 2 이면 호환
            moved = sum(1 for lit in (0xAAAA, 0xCCCC, 0xF0F0, 0xFF00) if p[lit] != lit)
            if moved == 2:
                transpositions.append(k)
    assert len(neg_one) == 4 * 24 and len(transpositions) == 6

    groups = {"G1 같은 궤도, 무작위": [], "G2 다른 궤도, 같은 크기": [], "G3 입력 하나 반전": [], "G4 입력 둘 맞바꿈": []}
    while len(groups["G1 같은 궤도, 무작위"]) < pairs_per_group:
        f = rng.choice(full)
        g = int(perms[rng.randrange(len(weights))][f])
        if g != f:
            groups["G1 같은 궤도, 무작위"].append((f, g))
    while len(groups["G2 다른 궤도, 같은 크기"]) < pairs_per_group:
        f = rng.choice(full)
        g = rng.choice(by_cost[int(single[f])])
        if canon[g] != canon[f]:
            groups["G2 다른 궤도, 같은 크기"].append((f, g))
    pure_neg = [k for k in neg_one if all(int(perms[k][lit]) in (lit, lit ^ 0xFFFF) for lit in (0xAAAA, 0xCCCC, 0xF0F0, 0xFF00))]
    assert len(pure_neg) == 4
    while len(groups["G3 입력 하나 반전"]) < pairs_per_group:
        f = rng.choice(full)
        g = int(perms[rng.choice(pure_neg)][f])
        assert g != f
        groups["G3 입력 하나 반전"].append((f, g))
    while len(groups["G4 입력 둘 맞바꿈"]) < pairs_per_group:
        f = rng.choice(full)
        g = int(perms[rng.choice(transpositions)][f])
        if g != f:
            groups["G4 입력 둘 맞바꿈"].append((f, g))

    print("무리마다 %d쌍. 비용은 AND 수(흐름 B). 따로 = C(f) + C(g), 함께 = 두 출력을 한 회로로 합성한 값." % pairs_per_group)
    print("%-22s %9s %9s %16s %9s %12s  %s" % ("무리", "따로 평균", "함께 평균", "줄어든 평균(표준오차)", "줄어든 비율", "안 줄어든 쌍", "줄어든 AND 수의 분포"))
    for name, pairs in groups.items():
        joint = joint_sizes(pairs)
        sep = [int(single[f]) + int(single[g]) for f, g in pairs]
        gain = [s - j for s, j in zip(sep, joint)]
        n = len(pairs)
        mean = sum(gain) / n
        se = (sum((x - mean) ** 2 for x in gain) / (n - 1) / n) ** 0.5
        print("%-22s %9.3f %9.3f %9.3f (%.3f) %8.2f%% %12d  %s" % (
            name, sum(sep) / n, sum(joint) / n, mean, se,
            100.0 * sum(gain) / sum(sep), sum(1 for x in gain if x <= 0), dict(sorted(Counter(gain).items()))))

    # 한 궤도의 함수 전부를 한꺼번에 구현하는 예: 최소항 16개 (4-16 디코더)
    minterms = [1 << x for x in range(N)]
    assert len({int(canon[m]) for m in minterms}) == 1
    j = joint_sizes([minterms])[0]
    print()
    print("예: 최소항 16개(한 궤도 전부, 4-16 디코더). 따로 = %d, 함께 = %d." % (sum(int(single[m]) for m in minterms), j))


if __name__ == "__main__":
    main()
