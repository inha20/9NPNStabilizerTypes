#!/usr/bin/env python3
"""
circuit_stronger_scripts.py
===========================
circuit_analyze.py 가 보인 "궤도 안의 차이"가 약한 합성 스크립트 때문인지 확인한다.

기준 흐름(balance/rewrite/refactor 묶음 세 번)에서 AND 수가 궤도의 최솟값보다 컸던 함수를 뽑아,
더 강한 스크립트 두 가지로 다시 합성한다.
  강한 1: 같은 묶음을 열 번
  강한 2: ABC 의 compress2rs 에 해당하는 순서(resub 포함)를 두 번, 앞뒤로 dc2
그리고 차이가 없던 함수도 일부 뽑아, 강한 스크립트가 궤도의 최솟값보다 더 작은 회로를 찾는지 본다.

실행: python scripts/circuit_stronger_scripts.py notes/circuit_and_counts.txt [차이 있는 함수 표본 수=600] [차이 없는 함수 표본 수=200]
준비: pip install yowasp-yosys
"""
import json
import os
import random
import subprocess
import sys
import tempfile
from collections import Counter

N = 16
FULL = 1 << N
RESYN2 = ["balance", "rewrite", "refactor", "balance", "rewrite", "rewrite -z", "balance",
          "refactor -z", "rewrite -z", "balance"]
COMPRESS2RS = ["balance -l", "resub -K 6 -l", "rewrite -l", "resub -K 6 -N 2 -l", "refactor -l",
               "resub -K 8 -l", "balance -l", "resub -K 8 -N 2 -l", "rewrite -l", "resub -K 10 -l",
               "rewrite -z -l", "resub -K 10 -N 2 -l", "balance -l", "resub -K 12 -l", "refactor -z -l",
               "resub -K 12 -N 2 -l", "rewrite -z -l", "balance -l"]
SCRIPTS = {
    "강한 1 (묶음 열 번)": ["strash"] + RESYN2 * 10,
    "강한 2 (compress2rs 두 번, dc2)": ["strash", "dc2"] + COMPRESS2RS * 2 + ["dc2"],
}

counts = [None] * FULL
for line in open(sys.argv[1]):
    if line.startswith("#") or not line.strip():
        continue
    chunk, digits = line.split()
    for i, ch in enumerate(digits):
        counts[int(chunk) * 4096 + i] = int(ch, 16)
n_gap = int(sys.argv[2]) if len(sys.argv) > 2 else 600
n_flat = int(sys.argv[3]) if len(sys.argv) > 3 else 200


def var_perm_on_points(p):
    return tuple(sum(((x >> p[i]) & 1) << i for i in range(4)) for x in range(N))


def act(f, g):
    h = 0
    for x in range(N):
        if (f >> x) & 1:
            h |= 1 << g[x]
    return h


gens = [tuple(x ^ m for x in range(N)) for m in (1, 2, 4, 8)] + [
    var_perm_on_points((1, 0, 2, 3)), var_perm_on_points((0, 2, 1, 3)), var_perm_on_points((0, 1, 3, 2))]
best = [None] * FULL
for f in range(FULL):
    if best[f] is not None:
        continue
    members, stack, seen = [f], [f], {f}
    while stack:
        u = stack.pop()
        for g in gens:
            v = act(u, g)
            if v not in seen:
                seen.add(v)
                members.append(v)
                stack.append(v)
    lo = min(counts[m] for m in members)
    for m in members:
        best[m] = lo

rng = random.Random(20261112)
gap = [f for f in range(FULL) if counts[f] > best[f]]
flat = [f for f in range(FULL) if counts[f] == best[f] and counts[f] >= 3]
gap_sample = rng.sample(gap, min(n_gap, len(gap)))
flat_sample = rng.sample(flat, min(n_flat, len(flat)))
print("기준 흐름에서 궤도의 최솟값보다 컸던 함수 %d개 가운데 %d개, 같았던 함수 가운데 %d개를 뽑았다."
      % (len(gap), len(gap_sample), len(flat_sample)))


def synthesize(funcs, script):
    with tempfile.TemporaryDirectory() as d:
        with open(os.path.join(d, "f.v"), "w") as fh:
            for f in funcs:
                fh.write("module f_%04x(input [3:0] x, output y); wire [15:0] t = 16'h%04x; assign y = t[x]; endmodule\n" % (f, f))
        with open(os.path.join(d, "s.abc"), "w") as fh:
            fh.write("\n".join(script) + "\n")
        with open(os.path.join(d, "s.ys"), "w") as fh:
            fh.write("read_verilog f.v\nproc\nopt\ntechmap\nopt\nabc -script s.abc\nopt_clean\ntee -q -o stat.json stat -json\n")
        subprocess.run(["yowasp-yosys", "-q", "-s", "s.ys"], cwd=d, check=True)
        stat = json.load(open(os.path.join(d, "stat.json")))
    return {int(name.split("_")[-1], 16): info.get("num_cells_by_type", {}).get("$lut", 0)
            for name, info in stat["modules"].items()}


for name, script in SCRIPTS.items():
    res = synthesize(gap_sample + flat_sample, script)
    still = sum(1 for f in gap_sample if res[f] > best[f])
    equal = sum(1 for f in gap_sample if res[f] == best[f])
    below = sum(1 for f in gap_sample if res[f] < best[f])
    base_gap = sum(counts[f] - best[f] for f in gap_sample)
    new_gap = sum(max(res[f] - best[f], 0) for f in gap_sample)
    worse = sum(1 for f in gap_sample if res[f] > counts[f])
    print("[%s]" % name)
    print("  차이 있던 함수 %d개: 여전히 최솟값보다 큼 %d, 최솟값과 같아짐 %d, 최솟값보다 작아짐 %d (기준보다 나빠진 것 %d)"
          % (len(gap_sample), still, equal, below, worse))
    print("  최솟값과의 차이의 합: 기준 %d -> %d" % (base_gap, new_gap))
    fb = sum(1 for f in flat_sample if res[f] < best[f])
    fa = sum(1 for f in flat_sample if res[f] > best[f])
    print("  차이 없던 함수 %d개: 최솟값보다 작아짐 %d, 커짐 %d" % (len(flat_sample), fb, fa))
