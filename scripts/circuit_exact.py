#!/usr/bin/env python3
"""
circuit_exact.py
================
4변수 부울 함수의 NPN 동치류 222개마다, 회로의 최소 크기를 ABC 의 정확 합성(SAT 풀이)으로 구한다
(PLAN.md 6절 "둘째 갈래"의 C6). 도구가 낸 값(circuit_metrics.py)이 최적값인지 견주는 데 쓴다.

기준
  and : 2입력 AND 수. 입력과 출력의 반전은 0으로 센다 (AIG 의 크기).     ABC 명령 `twoexact -a`
  l2  : 2입력 게이트의 수. 게이트의 종류는 가리지 않는다 (2입력 LUT 의 수).   ABC 명령 `twoexact`
  l3  : 3입력 LUT 의 수.                                                    ABC 명령 `lutexact -K 3`

이 세 기준에서는 입력 순열·입력 반전·출력 반전이 비용을 바꾸지 않으므로, 최솟값은 NPN 동치류마다 하나이다.
동치류의 대표로는 진리표가 가장 작은 함수를 쓴다.

방법: 대표 함수를 실제로 의존하는 변수만의 함수로 줄인 뒤, 노드 수 N 을 하한부터 하나씩 늘려 가며
"N 개로 만들 수 있는가"를 묻는다. 처음으로 만들 수 있다는 답이 나온 N 이 최솟값이다.
만들 수 있다는 답은 ABC 가 내놓은 회로를 이 스크립트가 다시 계산해 진리표가 같은지 확인한다.
만들 수 없다는 답은 SAT 풀이기(Glucose)의 답을 믿는다.
아주 작은 문제에서 ABC 가 답을 적지 않고 끝나는 경우가 있었는데, 그때는 이 스크립트가 모든 회로를 따져 답한다 (N 이 3 이하일 때).
시간 제한에 걸리면 그때까지 증명된 하한만 적는다 (상태 `lower`).

실행: python scripts/circuit_exact.py <abc 실행 파일> <기준> <출력 파일> [<호출 한 번의 시간 제한(초), 기본 600>] [<나눌 수> <나머지>]
      마지막 두 인자를 주면 대표들 가운데 (순번 % 나눌 수 == 나머지)인 것만 한다 (여러 개를 동시에 돌릴 때).
      출력 파일에 이미 끝난 대표가 있으면 건너뛰고, 하한만 적힌 대표는 그 하한부터 이어서 한다.
      환경 변수 CIRCUIT_EXACT_SKIP 에 16진수 진리표를 쉼표로 이어 주면 그 대표들은 건너뛴다 (어려운 것을 나중으로 미룰 때).
출력: "대표 진리표  의존 변수 수  기준  값  상태(exact|lower)  걸린 시간(초)  회로" 한 줄씩.
      상태가 lower 인 줄의 값은 "이보다 작은 N 으로는 만들 수 없다"까지만 증명된 하한이다.
      같은 대표가 여러 줄에 나오면 마지막 줄이 유효하다.
준비: ABC 를 소스에서 빌드한다 (https://github.com/berkeley-abc/abc, `make ABC_USE_NO_READLINE=1`).
"""
import itertools
import os
import re
import subprocess
import sys
import tempfile
import time

N = 16
FULL = 1 << N
WORKDIR = tempfile.mkdtemp(prefix="circuit_exact_")


def point_map(p, a):
    out = []
    for x in range(N):
        y = 0
        for i in range(4):
            y |= ((x >> p[i]) & 1) << i
        out.append(y ^ a)
    return out


MAPS = [point_map(p, a) for p in itertools.permutations(range(4)) for a in range(N)]


def act(f, t):
    h = 0
    for x in range(N):
        if (f >> x) & 1:
            h |= 1 << t[x]
    return h


def npn_representatives():
    seen = bytearray(FULL)
    reps = []
    for f in range(FULL):
        if seen[f]:
            continue
        reps.append(f)
        for t in MAPS:
            h = act(f, t)
            seen[h] = 1
            seen[h ^ (FULL - 1)] = 1
    return reps


def support(f):
    return [i for i in range(4) if any(((f >> x) & 1) != ((f >> (x ^ (1 << i))) & 1) for x in range(N))]


def compress(f, sup):
    """f 를 의존 변수만의 함수로 줄인 진리표."""
    g = 0
    for y in range(1 << len(sup)):
        x = 0
        for j, i in enumerate(sup):
            x |= ((y >> j) & 1) << i
        g |= ((f >> x) & 1) << y
    return g


GATE = re.compile(r"^([A-Z])\s*=\s*\d+'b([01]+)\(\s*([A-Za-z ]+?)\s*\)\s*$")


def simulate(lines, s):
    """ABC 가 출력한 회로(줄마다 '이름 = 진리표( 입력들 )')를 계산해 s 변수 진리표를 돌려준다.
    노드의 진리표에서 비트 번호는 첫째 입력을 가장 높은 자리로 한 이진수이다.
    ABC 는 출력 노드를 맨 먼저 적으므로, 계산은 아래 줄부터 한다."""
    gates = []
    for line in lines:
        m = GATE.match(line.strip())
        if m:
            gates.append((m.group(1), m.group(2), m.group(3).split()))
    if not gates:
        return None, ""
    out_name = gates[0][0]
    table = 0
    for y in range(1 << s):
        val = {chr(ord("a") + i): (y >> i) & 1 for i in range(s)}
        for name, tt, fanins in reversed(gates):
            idx = 0
            for fi in fanins:
                idx = (idx << 1) | val[fi]
            val[name] = int(tt[len(tt) - 1 - idx])
        table |= val[out_name] << y
    return table, ";".join("%s=%s(%s)" % (n, tt, "".join(fi)) for n, tt, fi in gates)


def ask(abc, metric, s, hexstr, n, limit):
    # Glucose(-g)가 빠르지만, 풀기도 전에 모순이 드러나는 작은 문제에서는 답을 적지 않고 끝나는 경우가 있었다.
    # 그때는 ABC 의 기본 풀이기로 다시 묻는다. ABC 는 풀이를 파일로도 적으므로 임시 폴더에서 돌린다.
    text = ""
    for solver in ("-g ", ""):
        if metric == "and":
            cmd = "twoexact -I %d -N %d -a %s%s" % (s, n, solver, hexstr)
        elif metric == "l2":
            cmd = "twoexact -I %d -N %d %s%s" % (s, n, solver, hexstr)
        else:
            cmd = "lutexact -N %d -M %d -K 3 %s%s" % (s, n, solver, hexstr)
        try:
            done = subprocess.run([abc, "-q", cmd], stdin=subprocess.DEVNULL, capture_output=True, text=True,
                                  timeout=limit, cwd=WORKDIR)
        except subprocess.TimeoutExpired:
            return "timeout", []
        text = done.stdout + done.stderr
        if "no solution" in text or "cannot be implemented" in text:
            return "unsat", []
        if "Realization" in text:
            return "sat", text.splitlines()
    # 두 풀이기 모두 답을 적지 않고 끝난 경우 (아주 작은 문제에서 생긴다): 직접 모두 따져 본다.
    if metric in ("and", "l2") and n <= 3:
        return ("sat-enumerated" if enumerate_small(metric, s, int(hexstr, 16), n) else "unsat"), []
    return "error: " + text.strip()[-200:], []


def enumerate_small(metric, s, g, n):
    """게이트 n 개(마지막 게이트가 출력)로 s 변수 함수 g 를 만들 수 있는지를 모든 경우를 따져 답한다. n 이 작을 때만 쓴다.
    게이트는 두 입력의 반전을 자유롭게 고르는 AND 이고, 기준이 l2 이면 XOR 도 쓴다. 출력의 반전은 자유이다."""
    size = 1 << s
    ones = (1 << size) - 1
    start = []
    for i in range(s):
        v = 0
        for y in range(size):
            v |= ((y >> i) & 1) << y
        start.append(v)

    def rec(signals, k):
        if k == n:
            return signals[-1] == g or (signals[-1] ^ ones) == g
        for i in range(len(signals)):
            for j in range(i + 1, len(signals)):
                a, b = signals[i], signals[j]
                cands = [a & b, a & (b ^ ones), (a ^ ones) & b, (a ^ ones) & (b ^ ones)]
                if metric == "l2":
                    cands.append(a ^ b)
                for v in cands:
                    if rec(signals + [v], k + 1):
                        return True
        return False

    return rec(start, 0)


def solve(abc, metric, f, limit, start=0):
    sup = support(f)
    s = len(sup)
    if s <= 1:
        return 0, "exact", ""
    if metric == "l3" and s <= 3:
        return 1, "exact", ""
    g = compress(f, sup)
    hexstr = "%0*x" % (max(1, (1 << s) // 4), g)
    n = max(start, (s - 1) if metric != "l3" else 2)
    while True:
        answer, lines = ask(abc, metric, s, hexstr, n, limit)
        if answer == "unsat":
            n += 1
            continue
        if answer == "sat-enumerated":
            return n, "exact", "(enumerated)"
        if answer == "sat":
            table, circuit = simulate(lines, s)
            if table != g:
                raise RuntimeError("ABC 가 낸 회로의 진리표가 다르다: %04x N=%d (%s)" % (f, n, circuit))
            return n, "exact", circuit
        if answer == "timeout":
            return n, "lower", ""
        raise RuntimeError("%04x N=%d: %s" % (f, n, answer))


def main():
    abc, metric, out_path = os.path.abspath(sys.argv[1]), sys.argv[2], sys.argv[3]
    limit = int(sys.argv[4]) if len(sys.argv) > 4 else 600
    mod, rem = (int(sys.argv[5]), int(sys.argv[6])) if len(sys.argv) > 6 else (1, 0)
    assert metric in ("and", "l2", "l3")
    done, lower = {}, {}
    if os.path.exists(out_path):
        for line in open(out_path):
            parts = line.split()
            if len(parts) >= 6 and parts[4] == "exact":
                done[int(parts[0], 16)] = True
            elif len(parts) >= 6 and parts[4] == "lower":      # 그 값보다 작은 N 으로는 안 된다는 것까지 증명되었다
                lower[int(parts[0], 16)] = max(lower.get(int(parts[0], 16), 0), int(parts[3]))
    reps = npn_representatives()
    assert len(reps) == 222
    skip = {int(x, 16) for x in os.environ.get("CIRCUIT_EXACT_SKIP", "").split(",") if x}
    with open(out_path, "a") as out:
        for i, f in enumerate(reps):
            if i % mod != rem or f in done or f in skip:
                continue
            t = time.time()
            value, status, circuit = solve(abc, metric, f, limit, lower.get(f, 0))
            out.write("%04x %d %s %d %s %.1f %s\n" % (f, len(support(f)), metric, value, status, time.time() - t, circuit))
            out.flush()


if __name__ == "__main__":
    main()
