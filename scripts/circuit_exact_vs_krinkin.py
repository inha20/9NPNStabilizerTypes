#!/usr/bin/env python3
"""
circuit_exact_vs_krinkin.py
===========================
이 폴더의 AND 정확 합성 결과(circuit_exact.py 의 `and` 기준)를 Krinkin 의 공개 자료와 동치류마다 견준다.
자료: https://github.com/krinkin/bounds  (data/npn4_opt_aig.csv, data/npn4_witnesses.json, MIT 라이선스)
논문: K. Krinkin, arXiv:2603.09379.  자료는 읽기만 하고 그 저장소의 코드는 실행하지 않는다.

하는 일
  1. 동치류 대표(진리표의 정수)가 두 쪽에서 같은 222개인지 본다.
  2. 이 폴더가 최솟값을 확정(exact)한 동치류마다 Krinkin 의 값과 같은지 본다. 다른 것이 있으면 오류로 멈춘다.
  3. 이 폴더가 확정하지 못한 동치류는 (증명된 하한 ~ 상한)을 적고 Krinkin 의 값을 곁에 둔다.
     상한은 ABC 로 흐름 B(circuit_metrics.py 와 같은 순서)를 돌려 얻은 회로를 AIGER 로 받아 다시 계산해
     진리표가 같음을 확인한 AND 수이다(ABC 실행 파일 경로를 주었을 때).
  4. Krinkin 의 witness 회로(10 AND)를 직접 계산해 진리표를 확인한다(자료 폴더를 주었을 때).
     그 저장소의 DRAT 하한 증명(k = 9)은 여기서 확인하지 않는다.

실행: python3 -I scripts/circuit_exact_vs_krinkin.py <ABC 실행 파일 또는 -> <Krinkin 저장소 폴더> <and 결과 파일>...
"""
import csv
import json
import os
import subprocess
import sys
import tempfile
from collections import Counter

N = 16
ONES = (1 << N) - 1

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from circuit_metrics import COMPRESS2RS  # noqa: E402


def load_mine(paths):
    """파일들에서 동치류마다 마지막 줄을 유효한 줄로 읽는다. 줄 모양: <16진수> <입력 수> and <값> <exact|lower> <초> <회로>"""
    res = {}
    for p in paths:
        for line in open(p):
            f = line.split()
            if len(f) < 5 or f[2] != "and":
                continue
            res[int(f[0], 16)] = (f[4], int(f[3]))
    return res


def load_krinkin(d):
    out = {}
    with open(os.path.join(d, "data", "npn4_opt_aig.csv")) as fh:
        for row in csv.DictReader(fh):
            out[int(row["npn_rep_dec"])] = (int(row["opt_aig"]), row["status"])
    return out


def input_words():
    """입력 x_j 의 진리표(행 t 에서 t 의 j번 비트)."""
    return [sum(((t >> j) & 1) << t for t in range(N)) for j in range(4)]


def eval_witness(gates, out_inv):
    """Krinkin 의 형식: 노드 1..4 는 입력 x0..x3, 게이트 i 는 노드 4+i, [a, pa, b, pb] 는 AND(v[a]^pa, v[b]^pb)."""
    v = {j + 1: w for j, w in enumerate(input_words())}
    node = 4
    for a, pa, b, pb in gates:
        node += 1
        v[node] = (v[a] ^ (ONES if pa else 0)) & (v[b] ^ (ONES if pb else 0))
    return v[node] ^ (ONES if out_inv else 0)


def parse_aiger(path):
    """ABC 가 쓴 이진 AIGER(aig)를 읽어 (AND 수, 출력의 진리표)를 돌려준다. 입력 4개, 출력 1개, 래치 없음."""
    data = open(path, "rb").read()
    nl = data.index(b"\n")
    hdr = data[:nl].split()
    assert hdr[0] == b"aig"
    M, I, L, O, A = (int(x) for x in hdr[1:6])
    assert I == 4 and L == 0 and O == 1 and M == I + A
    pos = nl + 1
    nl2 = data.index(b"\n", pos)
    out_lit = int(data[pos:nl2])
    pos = nl2 + 1

    def read_delta():
        nonlocal pos
        x = shift = 0
        while True:
            ch = data[pos]
            pos += 1
            x |= (ch & 0x7F) << shift
            if not ch & 0x80:
                return x
            shift += 7

    val = {0: 0}
    for j, w in enumerate(input_words()):
        val[2 * (j + 1)] = w

    def lit(l):
        w = val[l & ~1]
        return w ^ ONES if l & 1 else w

    for k in range(A):
        lhs = 2 * (I + 1 + k)
        d0 = read_delta()
        d1 = read_delta()
        r0 = lhs - d0
        r1 = r0 - d1
        val[lhs] = lit(r0) & lit(r1)
    return A, lit(out_lit)


def flow_b_upper(abc, f):
    """흐름 B 로 합성한 AIG 의 AND 수와, 그 회로를 다시 계산했을 때 진리표가 f 와 같은지."""
    with tempfile.TemporaryDirectory() as d:
        with open(os.path.join(d, "f.blif"), "w") as fh:
            fh.write(".model f\n.inputs x0 x1 x2 x3\n.outputs y\n.names x0 x1 x2 x3 y\n")
            for x in range(N):
                if (f >> x) & 1:
                    fh.write("%d%d%d%d 1\n" % (x & 1, (x >> 1) & 1, (x >> 2) & 1, (x >> 3) & 1))
            fh.write(".end\n")
        cmds = ["read_blif f.blif", "strash", "dc2"] + COMPRESS2RS * 2 + ["dc2", "write_aiger o.aig"]
        with open(os.path.join(d, "s.abc"), "w") as fh:
            fh.write("\n".join(cmds) + "\n")
        subprocess.run([os.path.abspath(abc), "-q", "source s.abc"], cwd=d, stdin=subprocess.DEVNULL,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        A, tt = parse_aiger(os.path.join(d, "o.aig"))
        return A, tt == f


def main():
    abc, kdir, files = sys.argv[1], sys.argv[2], sys.argv[3:]
    mine = load_mine(files)
    kr = load_krinkin(kdir)
    try:
        commit = subprocess.run(["git", "-C", kdir, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
    except OSError:
        commit = "?"
    print("Krinkin 자료: github.com/krinkin/bounds, 커밋 %s, 동치류 %d개, 상태 %s" % (commit, len(kr), sorted({s for _, s in kr.values()})))
    print("이 폴더의 결과: 동치류 %d개 (확정 %d개)" % (len(mine), sum(1 for s, _ in mine.values() if s == "exact")))
    ok = True

    same_reps = set(mine) == set(kr)
    print("1. 동치류 대표가 두 쪽에서 같은 집합인가: %s" % ("예 (222개)" if same_reps else "아니오"))
    ok &= same_reps

    exact = {h: v for h, (s, v) in mine.items() if s == "exact"}
    diff = [(h, v, kr[h][0]) for h, v in exact.items() if h in kr and v != kr[h][0]]
    agree = sum(1 for h, v in exact.items() if h in kr and v == kr[h][0])
    print("2. 이 폴더가 확정한 %d개 가운데 Krinkin 의 값과 같은 것: %d개, 다른 것: %d개 %s" % (
        len(exact), agree, len(diff), diff if diff else ""))
    ok &= not diff

    openc = sorted(h for h, (s, _) in mine.items() if s != "exact")
    print("3. 확정하지 못한 동치류 %d개 (대표, 증명된 하한, ABC 흐름 B 의 상한과 그 회로 확인, Krinkin 의 값):" % len(openc))
    for h in openc:
        lo = mine[h][1]
        kv = kr[h][0]
        if abc != "-":
            hi, good = flow_b_upper(abc, h)
            hi_txt = "%d (회로를 다시 계산해 진리표 %s)" % (hi, "일치" if good else "불일치")
            ok &= good
            inside = lo <= kv <= hi
        else:
            hi_txt = "계산하지 않음"
            inside = lo <= kv
        ok &= inside
        print("   %04x: 하한 %d, 상한 %s, Krinkin %d (%s)" % (h, lo, hi_txt, kv, "구간 안" if inside else "구간 밖"))

    wp = os.path.join(kdir, "data", "npn4_witnesses.json")
    if os.path.exists(wp):
        w = json.load(open(wp))
        print("4. Krinkin 의 witness 회로(자료의 입력 순서: x_j 는 행 번호의 j번 비트)를 직접 계산:")
        for e in w["witnesses"]:
            tt = eval_witness(e["gates"], e["out_inv"])
            good = tt == int(e["npn_rep_hex"], 16)
            print("   %s: AND %d개, 진리표 %04x -> %s" % (e["npn_rep_hex"], len(e["gates"]), tt, "일치" if good else "불일치"))
            ok &= good

    allv = Counter()
    for h in kr:
        allv[kr[h][0]] += 1
    print("Krinkin 의 값의 분포(동치류 수): %s" % dict(sorted(allv.items())))
    mv = Counter(mine[h][1] for h in mine if mine[h][0] == "exact")
    print("이 폴더가 확정한 값의 분포(동치류 수): %s" % dict(sorted(mv.items())))
    print("결과:", "불일치 없음" if ok else "불일치가 있다")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
