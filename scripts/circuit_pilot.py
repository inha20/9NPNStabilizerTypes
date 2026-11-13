#!/usr/bin/env python3
"""
circuit_pilot.py
================
회로 최적화 탐색의 예비 실험 (PLAN.md 6절 "둘째 갈래"의 C2, C3).

4변수 부울 함수의 B_4-궤도 402개에서 궤도마다 함수를 몇 개씩 뽑아, 기존 합성 도구(Yosys 안의 ABC)로
AND/NOT 회로를 합성하고 게이트 수가 궤도 안에서 어떻게 갈리는지 본다.

  측정 1 (AND 수)   : ABC 의 AIG 최적화(balance/rewrite/refactor 를 세 번 되풀이) 뒤의 2입력 AND 수
  측정 2 (NOT 수)   : 같은 회로를 AND2 와 INV 셀로 옮겼을 때의 인버터 수
  측정 3 (대표+변환): 궤도 대표 함수의 회로에 입력 순열(배선, 비용 0)과 입력 반전(인버터)을 붙여
                      같은 함수를 얻을 때의 게이트 수. 필요한 입력 인버터의 수는 대표를 그 함수로 옮기는
                      변환들 가운데 반전하는 입력이 가장 적은 것으로 센다.

도구가 낸 값은 최적값이 아닐 수 있다. 이 스크립트는 "도구의 결과"를 재는 것이다.

준비: pip install yowasp-yosys   (Yosys 와 ABC 가 들어 있다)
실행: python scripts/circuit_pilot.py [궤도마다 뽑을 함수 수, 기본 8] [결과 CSV 경로]
"""
import itertools
import json
import os
import random
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict

N = 16
FULL = 1 << N
K = int(sys.argv[1]) if len(sys.argv) > 1 else 8
CSV = sys.argv[2] if len(sys.argv) > 2 else None

RESYN2 = ["balance", "rewrite", "refactor", "balance", "rewrite", "rewrite -z", "balance",
          "refactor -z", "rewrite -z", "balance"]
AIG_SCRIPT = ["strash"] + RESYN2 * 3


def var_perm_on_points(p):
    out = []
    for x in range(N):
        y = 0
        for i in range(4):
            y |= ((x >> p[i]) & 1) << i
        out.append(y)
    return tuple(out)


def act(f, g):
    h = 0
    for x in range(N):
        if (f >> x) & 1:
            h |= 1 << g[x]
    return h


# B_4 의 원소 384개: (점의 치환, 반전하는 입력의 수)
B4 = []
for p in itertools.permutations(range(4)):
    base = var_perm_on_points(p)
    for c in range(N):
        B4.append((tuple(b ^ c for b in base), bin(c).count("1")))

gens = [tuple(x ^ m for x in range(N)) for m in (1, 2, 4, 8)] + [
    var_perm_on_points((1, 0, 2, 3)), var_perm_on_points((0, 2, 1, 3)), var_perm_on_points((0, 1, 3, 2))]
label = [-1] * FULL
orbits = []
for f in range(FULL):
    if label[f] != -1:
        continue
    label[f] = len(orbits)
    members = [f]
    stack = [f]
    while stack:
        u = stack.pop()
        for g in gens:
            v = act(u, g)
            if label[v] == -1:
                label[v] = len(orbits)
                members.append(v)
                stack.append(v)
    orbits.append(sorted(members))
assert len(orbits) == 402

rng = random.Random(20261112)
sample = []
for members in orbits:
    rep = members[0]
    rest = members[1:]
    chosen = [rep] + (rest if len(rest) <= K - 1 else rng.sample(rest, K - 1))
    sample.append(chosen)
functions = [f for chosen in sample for f in chosen]
print("궤도 %d개, 뽑은 함수 %d개 (궤도마다 최대 %d개, 대표는 진리표가 가장 작은 함수)" % (len(orbits), len(functions), K))


def synthesize(funcs, use_cells):
    """Yosys+ABC 로 합성해 함수별 셀 수를 돌려준다. use_cells=True 이면 AND2/INV 셀로 옮긴다."""
    with tempfile.TemporaryDirectory() as d:
        with open(os.path.join(d, "f.v"), "w") as fh:
            for f in funcs:
                fh.write("module f_%04x(input [3:0] x, output y); wire [15:0] t = 16'h%04x; assign y = t[x]; endmodule\n" % (f, f))
        with open(os.path.join(d, "s.abc"), "w") as fh:
            fh.write("\n".join(AIG_SCRIPT + (["map -a"] if use_cells else [])) + "\n")
        abc = "abc -g AND -script s.abc" if use_cells else "abc -script s.abc"
        with open(os.path.join(d, "s.ys"), "w") as fh:
            fh.write("read_verilog f.v\nproc\nopt\ntechmap\nopt\n%s\nopt_clean\ntee -q -o stat.json stat -json\n" % abc)
        subprocess.run(["yowasp-yosys", "-q", "-s", "s.ys"], cwd=d, check=True)
        stat = json.load(open(os.path.join(d, "stat.json")))
    out = {}
    for name, info in stat["modules"].items():
        out[int(name.split("_")[-1], 16)] = info.get("num_cells_by_type", {})
    return out


aig = synthesize(functions, False)
cells = synthesize(functions, True)


def support_size(f):
    n = 0
    for i in range(4):
        if any(((f >> x) & 1) != ((f >> (x ^ (1 << i))) & 1) for x in range(N)):
            n += 1
    return n


and_aig = {f: (0 if support_size(f) <= 1 else aig[f].get("$lut", 0)) for f in functions}
and_cell = {f: cells[f].get("$_AND_", 0) for f in functions}
not_cell = {f: cells[f].get("$_NOT_", 0) for f in functions}
print("AND 수가 두 흐름(AIG / AND2·INV 셀)에서 다른 함수: %d개" % sum(1 for f in functions if and_aig[f] != and_cell[f]))

and_const = 0
and_spread = Counter()
not_spread = Counter()
tot_spread = Counter()
transform_cmp = Counter()
dist_and = Counter()
rows = []
for chosen in sample:
    rep = chosen[0]
    a = [and_aig[f] for f in chosen]
    n = [not_cell[f] for f in chosen]
    t = [and_cell[f] + not_cell[f] for f in chosen]
    and_const += (min(a) == max(a))
    and_spread[max(a) - min(a)] += 1
    not_spread[max(n) - min(n)] += 1
    tot_spread[max(t) - min(t)] += 1
    dist_and[min(a)] += 1
    for f in chosen[1:]:
        flips = min(w for g, w in B4 if act(rep, g) == f)
        via_rep = and_cell[rep] + not_cell[rep] + flips
        direct = and_cell[f] + not_cell[f]
        transform_cmp["대표+변환이 더 큼" if via_rep > direct else ("같음" if via_rep == direct else "대표+변환이 더 작음")] += 1
    rows.append((rep, len(orbits[label[rep]]), len(chosen), min(a), max(a), min(n), max(n), min(t), max(t)))

print("[AND 수] 뽑은 함수들의 AND 수가 모두 같은 궤도: %d / %d" % (and_const, len(sample)))
print("         궤도 안 (최대 - 최소) 의 분포:", dict(sorted(and_spread.items())))
print("         궤도별 최소 AND 수의 분포 (궤도 수):", dict(sorted(dist_and.items())))
print("[NOT 수] 궤도 안 (최대 - 최소) 의 분포:", dict(sorted(not_spread.items())))
print("[AND+NOT] 궤도 안 (최대 - 최소) 의 분포:", dict(sorted(tot_spread.items())))
print("[대표+변환 대 직접 합성] (AND+NOT, 대표가 아닌 함수 %d개):" % sum(transform_cmp.values()), dict(transform_cmp))
if CSV:
    with open(CSV, "w") as fh:
        fh.write("representative_hex,orbit_size,sampled,and_min,and_max,not_min,not_max,total_min,total_max\n")
        for r in rows:
            fh.write("%04x,%d,%d,%d,%d,%d,%d,%d,%d\n" % r)
    print("궤도별 표:", CSV)
