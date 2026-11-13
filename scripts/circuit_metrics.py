#!/usr/bin/env python3
"""
circuit_metrics.py
==================
4변수 부울 함수 65,536개 전부를 Yosys 안의 ABC 로 합성해, 세 가지 비용 기준으로 잰 값을 기록한다
(PLAN.md 6절 "둘째 갈래"의 C3. 비용 기준은 저자가 2026-11-12에 "세 가지 모두"로 정했다).

최적화 흐름 (circuit_all_functions.py 와 같다)
  A : strash 뒤에 resyn2 에 해당하는 묶음을 세 번
  B : strash, dc2 뒤에 compress2rs 에 해당하는 순서를 두 번, 끝에 dc2

비용 기준 (흐름마다 일곱 가지 값)
  and  : 2입력 AND 수 (AIG 크기). 인버터는 0.
  ai   : 게이트 수. AND2 와 INV 만 쓰고 둘 다 1로 센다.
  g2   : 게이트 수. INV, AND2, OR2, NAND2, NOR2, XOR2, XNOR2 를 쓰고 모두 1로 센다.
  cmos : 표준 셀 면적의 대용값. 작은 CMOS 셀 묶음에서 트랜지스터 쌍의 수를 면적으로 삼는다 (아래 CMOS).
         실제 공정의 셀 라이브러리가 아니다.
  l2, l3, l4 : FPGA LUT 수. 각각 2입력, 3입력, 4입력 LUT 로 옮겼을 때의 LUT 수.

셀로 옮기는 단계는 ABC 의 `dch -f` 뒤에 면적 위주의 매퍼 둘(`map -a`, `amap`)을 각각 돌려 작은 쪽을 택한다.
(시험 삼아 300개쯤의 함수에서 비교했더니 라이브러리에 따라 나은 쪽이 달랐고, `dch -f` 를 빼면 대체로 나빠졌다.)
LUT 로 옮기는 단계는 `dch -f` 뒤에 `if -K k -a` 이다.
도구가 낸 값이며 최적값이 아닐 수 있다.

실행: python scripts/circuit_metrics.py <출력 파일> <묶음 번호> [<묶음 번호> ...]     (묶음은 0~15, 4,096개씩)
      환경 변수 CIRCUIT_ABC 에 따로 빌드한 ABC 실행 파일의 경로를 주면, Yosys 안의 ABC 대신 그것으로 같은 스크립트를 돌린다
      (도구의 판이 달라도 결과가 같은지 보는 데 쓴다).
출력: "묶음 흐름 기준 문자열" 한 줄씩. 문자열은 함수마다 16진수 두 자리.
      전부 모은 뒤에는 cat 으로 이어 붙여 gzip 으로 묶어 둔다 (notes/circuit_metrics.txt.gz).
준비: pip install yowasp-yosys
"""
import os
import subprocess
import sys
import tempfile

N = 16
RESYN2 = ["balance", "rewrite", "refactor", "balance", "rewrite", "rewrite -z", "balance",
          "refactor -z", "rewrite -z", "balance"]
COMPRESS2RS = ["balance -l", "resub -K 6 -l", "rewrite -l", "resub -K 6 -N 2 -l", "refactor -l",
               "resub -K 8 -l", "balance -l", "resub -K 8 -N 2 -l", "rewrite -l", "resub -K 10 -l",
               "rewrite -z -l", "resub -K 10 -N 2 -l", "balance -l", "resub -K 12 -l", "refactor -z -l",
               "resub -K 12 -N 2 -l", "rewrite -z -l", "balance -l"]
FLOWS = {"A": RESYN2 * 3, "B": ["dc2"] + COMPRESS2RS * 2 + ["dc2"]}

# 셀 이름 -> (면적, 논리식). 면적은 genlib 의 둘째 칸에 들어간다.
AI = {"inv": (1, "O=!a"), "and2": (1, "O=a*b")}
G2 = {"inv": (1, "O=!a"), "and2": (1, "O=a*b"), "or2": (1, "O=a+b"), "nand2": (1, "O=!(a*b)"),
      "nor2": (1, "O=!(a+b)"), "xor2": (1, "O=a*!b+!a*b"), "xnor2": (1, "O=a*b+!a*!b")}
CMOS = {"inv": (1, "O=!a"), "nand2": (2, "O=!(a*b)"), "nor2": (2, "O=!(a+b)"),
        "nand3": (3, "O=!(a*b*c)"), "nor3": (3, "O=!(a+b+c)"),
        "nand4": (4, "O=!(a*b*c*d)"), "nor4": (4, "O=!(a+b+c+d)"),
        "and2": (3, "O=a*b"), "or2": (3, "O=a+b"),
        "aoi21": (3, "O=!(a*b+c)"), "oai21": (3, "O=!((a+b)*c)"),
        "aoi22": (4, "O=!(a*b+c*d)"), "oai22": (4, "O=!((a+b)*(c+d))"),
        "xor2": (5, "O=a*!b+!a*b"), "xnor2": (5, "O=a*b+!a*!b")}
LIBS = {"ai": AI, "g2": G2, "cmos": CMOS}
LUTS = {"l2": 2, "l3": 3, "l4": 4}
METRICS = ["and", "ai", "g2", "cmos", "l2", "l3", "l4"]


def support(f):
    return [i for i in range(4) if any(((f >> x) & 1) != ((f >> (x ^ (1 << i))) & 1) for x in range(N))]


def degenerate_costs(f):
    """상수와 변수 하나짜리 함수의 값: 반전된 변수만 인버터 하나가 든다."""
    s = support(f)
    inv = 1 if (len(s) == 1 and not (f >> (1 << s[0])) & 1) else 0
    return {"and": 0, "ai": inv, "g2": inv, "cmos": inv, "l2": 0, "l3": 0, "l4": 0}


def write_genlib(path, lib):
    with open(path, "w") as fh:
        fh.write("GATE zero 0 O=CONST0;\nGATE one 0 O=CONST1;\n")
        for name, (area, expr) in lib.items():
            phase = "INV" if expr.startswith("O=!") and "+!" not in expr[3:] and "*!" not in expr[3:] else "UNKNOWN"
            fh.write("GATE %s %d %s; PIN * %s 1 999 1 0 1 0\n" % (name, area, expr, phase))


def run_chunk(chunk):
    return run_functions(list(range(chunk * 4096, (chunk + 1) * 4096)))


def run_functions(funcs):
    real = [f for f in funcs if len(support(f)) >= 2]
    values = {(fl, m): {} for fl in FLOWS for m in METRICS}
    for f in funcs:
        if len(support(f)) < 2:
            for fl in FLOWS:
                for m, v in degenerate_costs(f).items():
                    values[(fl, m)][f] = v
    with tempfile.TemporaryDirectory() as d:
        os.mkdir(os.path.join(d, "in"))
        os.mkdir(os.path.join(d, "out"))
        for f in real:
            with open(os.path.join(d, "in", "%04x.blif" % f), "w") as fh:
                fh.write(".model f\n.inputs x0 x1 x2 x3\n.outputs y\n.names x0 x1 x2 x3 y\n")
                for x in range(N):
                    if (f >> x) & 1:
                        fh.write("%d%d%d%d 1\n" % (x & 1, (x >> 1) & 1, (x >> 2) & 1, (x >> 3) & 1))
                fh.write(".end\n")
        for name, lib in LIBS.items():
            write_genlib(os.path.join(d, name + ".genlib"), lib)
        cmds = ["backup"]   # Yosys 가 넘겨준 회로를 잠시 치워 두었다가 끝에 되돌린다
        for fl, opt in FLOWS.items():
            for f in real:
                cmds += ["read_blif in/%04x.blif" % f, "strash"] + opt + ["write_aiger out/%04x.%s.aig" % (f, fl)]
        for name in LIBS:
            cmds.append("read_library %s.genlib" % name)
            for fl in FLOWS:
                for f in real:
                    for mapper, tag in (("map -a", "map"), ("amap", "amap")):
                        cmds += ["read_aiger out/%04x.%s.aig" % (f, fl), "dch -f", mapper,
                                 "write_blif out/%04x.%s.%s.%s.blif" % (f, fl, name, tag)]
        for name, k in LUTS.items():
            for fl in FLOWS:
                for f in real:
                    cmds += ["read_aiger out/%04x.%s.aig" % (f, fl), "dch -f", "if -K %d -a" % k,
                             "write_blif out/%04x.%s.%s.blif" % (f, fl, name)]
        native = os.environ.get("CIRCUIT_ABC")      # 따로 빌드한 ABC 실행 파일의 경로를 주면 Yosys 를 거치지 않는다
        if native:
            with open(os.path.join(d, "s.abc"), "w") as fh:
                fh.write("\n".join(cmds[1:]) + "\n")
            done = subprocess.run([os.path.abspath(native), "-q", "source s.abc"], cwd=d, stdin=subprocess.DEVNULL,
                                  stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        else:
            cmds += ["restore", "strash"]
            with open(os.path.join(d, "s.abc"), "w") as fh:
                fh.write("\n".join(cmds) + "\n")
            with open(os.path.join(d, "d.v"), "w") as fh:
                fh.write("module dummy(input a, input b, output y); assign y = a & b; endmodule\n")
            with open(os.path.join(d, "s.ys"), "w") as fh:
                fh.write("read_verilog d.v\nproc\nopt\ntechmap\nabc -script s.abc\n")
            done = subprocess.run(["yowasp-yosys", "-q", "-s", "s.ys"], cwd=d)
        expected = len(real) * len(FLOWS) * (1 + 2 * len(LIBS) + len(LUTS))
        if len(os.listdir(os.path.join(d, "out"))) != expected:
            raise RuntimeError("합성이 끝까지 가지 못했다 (종료 코드 %d)" % done.returncode)
        for fl in FLOWS:
            for f in real:
                with open(os.path.join(d, "out", "%04x.%s.aig" % (f, fl)), "rb") as fh:
                    values[(fl, "and")][f] = int(fh.readline().split()[5])
                for name, lib in LIBS.items():
                    areas = []
                    for tag in ("map", "amap"):
                        area = 0
                        for line in open(os.path.join(d, "out", "%04x.%s.%s.%s.blif" % (f, fl, name, tag))):
                            if line.startswith(".gate"):
                                cell = line.split()[1]
                                if cell not in ("zero", "one"):
                                    area += lib[cell][0]
                        areas.append(area)
                    values[(fl, name)][f] = min(areas)
                for name in LUTS:
                    luts = 0
                    for line in open(os.path.join(d, "out", "%04x.%s.%s.blif" % (f, fl, name))):
                        if line.startswith(".names") and len(line.split()) >= 4:   # 입력이 둘 이상인 노드
                            luts += 1
                    values[(fl, name)][f] = luts
    return funcs, values


if __name__ == "__main__":
    out_path, chunks = sys.argv[1], [int(a) for a in sys.argv[2:]]
    with open(out_path, "w") as out:
        for chunk in chunks:
            funcs, values = run_chunk(chunk)
            for fl in FLOWS:
                for m in METRICS:
                    assert all(values[(fl, m)][f] < 256 for f in funcs)
                    out.write("%d %s %s %s\n" % (chunk, fl, m, "".join("%02x" % values[(fl, m)][f] for f in funcs)))
            out.flush()
