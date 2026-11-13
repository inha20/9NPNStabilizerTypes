#!/usr/bin/env python3
"""
circuit_all_functions.py
========================
4변수 부울 함수 65,536개 전부를 Yosys 안의 ABC 로 합성해 2입력 AND 수(AIG 크기)를 기록한다.
도구가 낸 값이며 최적값이 아닐 수 있다. 흐름은 둘 가운데 고른다.
  resyn2x3    : strash 뒤에 balance/rewrite/refactor 묶음을 세 번 (circuit_pilot.py 의 "측정 1")
  compress2rs : strash, dc2 뒤에 ABC 의 compress2rs 에 해당하는 순서(resub 포함)를 두 번, 끝에 dc2

실행: python scripts/circuit_all_functions.py <흐름> <출력 파일> <묶음 번호> [<묶음 번호> ...]
      함수를 4,096개씩 16묶음(0~15)으로 나눈다. 한 묶음에 5분쯤 걸리므로 묶음을 나눠 여러 번 돌려도 되고,
      출력 파일들을 이어 붙이면 된다 (예: cat part*.txt > notes/circuit_and_counts.txt).
출력: "묶음번호 16진수문자열" 한 줄씩. 문자열의 i번째 글자가 함수 (묶음번호*4096 + i) 의 AND 수(16진수 한 자리).
준비: pip install yowasp-yosys
"""
import json
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
FLOWS = {"resyn2x3": ["strash"] + RESYN2 * 3,
         "compress2rs": ["strash", "dc2"] + COMPRESS2RS * 2 + ["dc2"]}
AIG_SCRIPT = FLOWS[sys.argv[1]]
out_path, chunks = sys.argv[2], [int(a) for a in sys.argv[3:]]


def support_size(f):
    return sum(1 for i in range(4)
               if any(((f >> x) & 1) != ((f >> (x ^ (1 << i))) & 1) for x in range(N)))


with open(out_path, "w") as out:
    for chunk in chunks:
        funcs = range(chunk * 4096, (chunk + 1) * 4096)
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, "f.v"), "w") as fh:
                for f in funcs:
                    fh.write("module f_%04x(input [3:0] x, output y); wire [15:0] t = 16'h%04x; assign y = t[x]; endmodule\n" % (f, f))
            with open(os.path.join(d, "s.abc"), "w") as fh:
                fh.write("\n".join(AIG_SCRIPT) + "\n")
            with open(os.path.join(d, "s.ys"), "w") as fh:
                fh.write("read_verilog f.v\nproc\nopt\ntechmap\nopt\nabc -script s.abc\nopt_clean\ntee -q -o stat.json stat -json\n")
            subprocess.run(["yowasp-yosys", "-q", "-s", "s.ys"], cwd=d, check=True)
            stat = json.load(open(os.path.join(d, "stat.json")))
        counts = {}
        for name, info in stat["modules"].items():
            counts[int(name.split("_")[-1], 16)] = info.get("num_cells_by_type", {}).get("$lut", 0)
        digits = []
        for f in funcs:
            a = 0 if support_size(f) <= 1 else counts[f]
            assert a < 16
            digits.append("%x" % a)
        out.write("%d %s\n" % (chunk, "".join(digits)))
        out.flush()
