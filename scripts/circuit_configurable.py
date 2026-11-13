#!/usr/bin/env python3
"""
circuit_configurable.py
=======================
설정값으로 4변수 함수 65,536개 가운데 하나를 고르는 회로 두 가지를 실제로 합성해 크기를 견준다
(탐색 노트 4.2절, PLAN.md 6절 "둘째 갈래"의 C4, 예상 P3).

  S1 (LUT)        : 설정 16비트 c. 출력은 c[x]. 16:1 멀티플렉서 하나이다.
  S2 (대표와 변환) : 설정 18비트 = 대표 번호 r(9비트, 궤도 402개) + 순열 p(5비트) + 반전 a(4비트).
                    입력 x 를 a 로 반전하고(XOR 4개), 스위치 5개짜리 망으로 순서를 바꾼 뒤(스위치마다 MUX 2개),
                    "대표 블록" U(r, w) = (r번 대표 함수의 w 에서의 값)에 넣는다.
  U 만            : S2 에서 반전과 순열을 회로 밖(배선, 이웃 셀)에 맡긴다고 볼 때 남는 부분. 설정 9비트.

합성은 circuit_metrics.py 와 같은 방식이다 (Yosys 안의 ABC, 흐름 B, 같은 셀 묶음, `dch -f` 뒤 `map -a`/`amap`, `if -K k -a`).
설정 비트는 회로의 입력으로 두고, 그것을 저장하는 기억 소자는 세지 않는다 (비트 수만 따로 적는다).
합성된 회로는 모든 입력 조합에서 다시 계산해 확인한다: S1 은 설정 65,536가지 전부, S2 는 쓰이는 설정 205,824가지
전부에서 모형과 같은 함수를 내는지, 그리고 S2 가 내는 함수를 모으면 65,536개 전부인지.

대표 블록 U 의 크기는 대표를 무엇으로 고르고 번호를 어떻게 매기느냐에 따라 달라진다. 세 가지를 재 본다.
  min    : 궤도에서 진리표가 가장 작은 함수를 대표로, 진리표 순으로 번호를 매긴다.
  random : 같은 대표에 번호를 무작위로 매긴다 (seed 1).
  weight : 같은 대표를 (1의 개수, 진리표) 순으로 번호를 매긴다.
쓰지 않는 번호(402~511)는 무엇이 나와도 되는 것으로 둔다.

실행: python scripts/circuit_configurable.py [작업 폴더]
준비: pip install yowasp-yosys
"""
import itertools
import os
import random
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import circuit_metrics as cm   # 셀 묶음과 최적화 흐름을 그대로 쓴다

N = 16
FULL = 1 << N

# 스위치 5개: 4입력 정렬망의 비교기 자리. 정렬망은 어떤 순서의 입력도 정렬하므로, 비교기를 자유 스위치로 바꾸면
# 24가지 순열이 모두 나온다 (아래 check_network 에서 확인).
SWITCHES = [(0, 1), (2, 3), (0, 2), (1, 3), (1, 2)]


def check_network():
    seen = set()
    for bits in range(32):
        w = [0, 1, 2, 3]
        for k, (i, j) in enumerate(SWITCHES):
            if (bits >> k) & 1:
                w[i], w[j] = w[j], w[i]
        seen.add(tuple(w))
    assert len(seen) == 24, "스위치 망이 모든 순열을 만들지 못한다"


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


def orbit_representatives():
    seen = bytearray(FULL)
    reps = []
    for f in range(FULL):
        if seen[f]:
            continue
        reps.append(f)
        for t in MAPS:
            seen[act(f, t)] = 1
    assert len(reps) == 402
    return reps


def verilog_s1():
    return "module top(input [15:0] c, input [3:0] x, output y);\n  assign y = c[x];\nendmodule\n"


def verilog_rom(order):
    lines = ["  reg [15:0] tt;", "  always @* begin", "    case (r)"]
    for i, f in enumerate(order):
        lines.append("      9'd%d: tt = 16'h%04x;" % (i, f))
    lines += ["      default: tt = 16'hxxxx;", "    endcase", "  end"]
    return "\n".join(lines) + "\n"


def verilog_u(order):
    return ("module top(input [8:0] r, input [3:0] w, output y);\n" + verilog_rom(order)
            + "  assign y = tt[w];\nendmodule\n")


def verilog_s2(order):
    v = ["module top(input [8:0] r, input [4:0] p, input [3:0] a, input [3:0] x, output y);",
         "  wire [3:0] s0 = x ^ a;"]
    for k, (i, j) in enumerate(SWITCHES):
        cur, nxt = "s%d" % k, "s%d" % (k + 1)
        v.append("  wire [3:0] %s;" % nxt)
        for b in range(4):
            if b == i:
                v.append("  assign %s[%d] = p[%d] ? %s[%d] : %s[%d];" % (nxt, b, k, cur, j, cur, i))
            elif b == j:
                v.append("  assign %s[%d] = p[%d] ? %s[%d] : %s[%d];" % (nxt, b, k, cur, i, cur, j))
            else:
                v.append("  assign %s[%d] = %s[%d];" % (nxt, b, cur, b))
    v.append("  wire [3:0] w = s%d;" % len(SWITCHES))
    return "\n".join(v) + "\n" + verilog_rom(order) + "  assign y = tt[w];\nendmodule\n"


def measure(d, name, verilog):
    """Verilog 한 모듈을 합성해 값들을 돌려준다 (흐름 B 만).
    Yosys 로 Verilog 를 AND 와 인버터의 회로(pre.aig)로 옮긴 다음, circuit_metrics.py 와 같은 방법으로
    Yosys 안의 ABC 에 스크립트를 넘겨 최적화한다 (opt.aig). 입력의 이름은 끝까지 따라간다."""
    work = os.path.join(d, name)
    os.makedirs(work, exist_ok=True)
    with open(os.path.join(work, "top.v"), "w") as fh:
        fh.write(verilog)
    with open(os.path.join(work, "pre.ys"), "w") as fh:
        fh.write("read_verilog top.v\nhierarchy -top top\nproc\nopt\nmemory\nopt\ntechmap\nopt\n"
                 "aigmap\nopt\nwrite_aiger -symbols pre.aig\n")
    subprocess.run(["yowasp-yosys", "-q", "-s", "pre.ys"], cwd=work, check=True)
    for lib_name, lib in cm.LIBS.items():
        cm.write_genlib(os.path.join(work, lib_name + ".genlib"), lib)
    cmds = ["backup", "read_aiger pre.aig", "strash"] + cm.FLOWS["B"] + ["write_aiger -s opt.aig"]
    for lib_name in cm.LIBS:
        cmds.append("read_library %s.genlib" % lib_name)
        for mapper, tag in (("map -a", "map"), ("amap", "amap")):
            cmds += ["read_aiger opt.aig", "dch -f", mapper, "write_blif %s.%s.blif" % (lib_name, tag)]
    for lut_name, k in list(cm.LUTS.items()) + [("l6", 6)]:
        cmds += ["read_aiger opt.aig", "dch -f", "if -K %d -a" % k, "write_blif %s.blif" % lut_name]
    cmds += ["restore", "strash"]
    with open(os.path.join(work, "s.abc"), "w") as fh:
        fh.write("\n".join(cmds) + "\n")
    with open(os.path.join(work, "d.v"), "w") as fh:
        fh.write("module dummy(input a, input b, output y); assign y = a & b; endmodule\n")
    with open(os.path.join(work, "s.ys"), "w") as fh:
        fh.write("read_verilog d.v\nproc\nopt\ntechmap\nabc -script s.abc\n")
    subprocess.run(["yowasp-yosys", "-q", "-s", "s.ys"], cwd=work, check=True)
    values = {}
    with open(os.path.join(work, "opt.aig"), "rb") as fh:
        values["and"] = int(fh.readline().split()[5])
    for lib_name, lib in cm.LIBS.items():
        areas = []
        for tag in ("map", "amap"):
            area = 0
            for line in open(os.path.join(work, "%s.%s.blif" % (lib_name, tag))):
                if line.startswith(".gate"):
                    cell = line.split()[1]
                    if cell not in ("zero", "one"):
                        area += lib[cell][0]
            areas.append(area)
        values[lib_name] = min(areas)
    for lut_name in list(cm.LUTS) + ["l6"]:
        values[lut_name] = sum(1 for line in open(os.path.join(work, lut_name + ".blif"))
                               if line.startswith(".names") and len(line.split()) >= 4)
    return values


def simulate_aig(path, position):
    """이진 AIGER 회로를 모든 입력 조합에 대해 한꺼번에 계산한다 (큰 정수의 비트 하나가 입력 조합 하나).
    position: 입력 이름 -> 입력 조합 번호에서 그 입력이 차지하는 비트 자리. 출력의 값을 큰 정수로 돌려준다."""
    with open(path, "rb") as fh:
        data = fh.read()
    end = data.index(b"\n")
    _, m, i, l, o, a = data[:end].split()
    i, l, o, a = int(i), int(l), int(o), int(a)
    assert l == 0 and o == 1
    pos = end + 1
    end = data.index(b"\n", pos)
    out_lit = int(data[pos:end])
    pos = end + 1

    def number():
        nonlocal pos
        x, shift = 0, 0
        while True:
            byte = data[pos]
            pos += 1
            x |= (byte & 0x7F) << shift
            if not byte & 0x80:
                return x
            shift += 7

    ands = []
    for k in range(a):
        lhs = 2 * (i + k + 1)
        r0 = lhs - number()
        r1 = r0 - number()
        ands.append((lhs, r0, r1))
    names = {}
    for ln in data[pos:].decode("ascii", "replace").split("\n"):
        if ln.startswith("i") and " " in ln and ln[1:ln.index(" ")].isdigit():
            k, name = ln[1:].split(" ", 1)
            names[int(k)] = name.strip()
    n = len(position)
    assert i == n and len(names) == n, "입력 수나 이름이 다르다"
    size = 1 << n
    ones = (1 << size) - 1
    val = {0: 0}
    for k in range(i):
        b = position[names[k]]
        block = ((1 << (1 << b)) - 1) << (1 << b)
        val[k + 1] = block * (ones // ((1 << (1 << (b + 1))) - 1))

    def get(lit):
        v = val[lit >> 1]
        return (v ^ ones) if (lit & 1) else v

    for lhs, r0, r1 in ands:
        val[lhs >> 1] = get(r0) & get(r1)
    return get(out_lit), size


def tables(out, size):
    """출력 값을 설정마다의 16비트 진리표로 나눈다 (x 가 입력 조합 번호의 아래 4비트일 때)."""
    import array
    arr = array.array("H")
    arr.frombytes(out.to_bytes(size // 8, "little"))
    return arr


def verify(d, orders):
    # S1: 설정 c 가 그대로 진리표여야 한다
    pos = {"x[%d]" % b: b for b in range(4)}
    pos.update({"c[%d]" % b: 4 + b for b in range(16)})
    tt = tables(*simulate_aig(os.path.join(d, "s1", "opt.aig"), pos))
    assert all(tt[c] == c for c in range(FULL))
    print("검증 S1: 합성된 회로가 설정 65,536가지에서 각각 그 진리표의 함수를 낸다.")
    switch_maps = []
    for bits in range(32):
        w = [0, 1, 2, 3]
        for k, (i, j) in enumerate(SWITCHES):
            if (bits >> k) & 1:
                w[i], w[j] = w[j], w[i]
        switch_maps.append(tuple(w))      # 망을 지난 뒤의 b번 선에는 처음의 w[b]번 선이 온다
    for tag, order in orders.items():
        pos = {"w[%d]" % b: b for b in range(4)}
        pos.update({"r[%d]" % b: 4 + b for b in range(9)})
        tt = tables(*simulate_aig(os.path.join(d, "u_" + tag, "opt.aig"), pos))
        assert all(tt[r] == order[r] for r in range(len(order)))
        pos = {"x[%d]" % b: b for b in range(4)}
        pos.update({"a[%d]" % b: 4 + b for b in range(4)})
        pos.update({"p[%d]" % b: 8 + b for b in range(5)})
        pos.update({"r[%d]" % b: 13 + b for b in range(9)})
        tt = tables(*simulate_aig(os.path.join(d, "s2_" + tag, "opt.aig"), pos))
        realized = set()
        for r, f in enumerate(order):
            for p in range(32):
                wm = switch_maps[p]
                for a in range(N):
                    g = 0
                    for x in range(N):
                        z = x ^ a
                        w = 0
                        for b in range(4):
                            w |= ((z >> wm[b]) & 1) << b
                        g |= ((f >> w) & 1) << x
                    assert tt[(r << 9) | (p << 4) | a] == g, "S2 의 회로가 모형과 다르다"
                    realized.add(g)
        assert len(realized) == FULL
        print("검증 번호 %s: U 는 대표 402개를 그대로 내고, S2 는 설정 %d가지에서 모형과 같은 함수를 내며 "
              "그 함수들을 모으면 65,536개 전부이다." % (tag, len(order) * 32 * 16))


def main():
    check_network()
    reps = orbit_representatives()
    orders = {
        "min": list(reps),
        "random": random.Random(1).sample(reps, len(reps)),
        "weight": sorted(reps, key=lambda f: (bin(f).count("1"), f)),
    }
    d = sys.argv[1] if len(sys.argv) > 1 else tempfile.mkdtemp(prefix="circuit_configurable_")
    rows = [("S1 (LUT)", 16, measure(d, "s1", verilog_s1()))]
    for tag, order in orders.items():
        rows.append(("U 만, 번호 %s" % tag, 9, measure(d, "u_" + tag, verilog_u(order))))
        rows.append(("S2, 번호 %s" % tag, 18, measure(d, "s2_" + tag, verilog_s2(order))))
    keys = ["and", "ai", "g2", "cmos", "l2", "l3", "l4", "l6"]
    print("%-18s %8s " % ("구조", "설정 비트") + " ".join("%7s" % k for k in keys))
    for name, bits, v in rows:
        print("%-18s %8d " % (name, bits) + " ".join("%7d" % v[k] for k in keys))
    print()
    verify(d, orders)
    print()
    print("S2 의 앞부분(반전과 순열)만의 크기: XOR 4개 + MUX 10개. AND 수로는 XOR 3, MUX 3 으로 쳐서 42.")
    print("S1 은 MUX 15개(AND 45개)이다.")


if __name__ == "__main__":
    main()
