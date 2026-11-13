#!/usr/bin/env python3
"""
circuit_lut_flows.py
====================
LUT 수(l2, l3)에서 도구의 값이 최솟값보다 크게 나온 것이 약한 설정 탓인지 보려고, ABC 의 다른 LUT 매핑 순서
몇 가지를 NPN 동치류의 대표 222개에서 돌려 본다 (예비 실험에서 약한 기준선 때문에 생긴 수치를 결과로 읽었던 일의 재발 방지).

  F0 : 흐름 B 뒤에 `dch -f; if -K k -a`                      (circuit_metrics.py 가 쓰는 것)
  F1 : F0 뒤에 `mfs2` 두 번
  F2 : (`dch -f; if -K k -a; mfs2; strash`) 세 번 뒤에 다시 매핑
  F3 : `&dch -f; &if -K k -a; &satlut`                       (SAT 로 LUT 를 다시 묶는다)
  F4 : 선택 노드 없이 `if -K k -a` 두 번

대표의 값을 그 동치류의 모든 함수의 값으로 쳐서 65,536개의 합을 어림한다 (전수 측정이 아니다).
노드 수는 ABC 의 `print_stats` 가 적는 nd 를 쓴다.

실행: python scripts/circuit_lut_flows.py <abc 실행 파일> <k: 2 또는 3> <정확 합성 결과 파일>
"""
import re
import subprocess
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import circuit_exact as ce
import circuit_metrics as cm


def main():
    abc, k, exact_path = os.path.abspath(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    b = "strash; " + "; ".join(cm.FLOWS["B"])
    flows = {
        "F0 흐름 B; dch -f; if -a": b + "; dch -f; if -K %d -a; print_stats" % k,
        "F1 F0; mfs2; mfs2": b + "; dch -f; if -K %d -a; mfs2; mfs2; print_stats" % k,
        "F2 (dch; if; mfs2; strash) x3": b + "; " + "; ".join(["dch -f; if -K %d -a; mfs2; strash" % k] * 3)
                                         + "; dch -f; if -K %d -a; print_stats" % k,
        "F3 &dch; &if; &satlut": b + "; &get; &dch -f; &if -K %d -a; &satlut -N 128 -C 10000; &put; print_stats" % k,
        "F4 if -a; if -a": b + "; if -K %d -a; if -K %d -a; print_stats" % (k, k),
    }
    reps = ce.npn_representatives()
    size = {}
    for f in reps:
        members = set()
        for t in ce.MAPS:
            h = ce.act(f, t)
            members.add(h)
            members.add(h ^ 0xFFFF)
        size[f] = len(members)
    exact = {}
    for line in open(exact_path):
        p = line.split()
        if len(p) >= 5 and p[4] == "exact":
            exact[int(p[0], 16)] = int(p[3])
    total_exact = sum(exact[f] * size[f] for f in reps)
    print("k = %d.  최솟값의 합(65,536개): %d" % (k, total_exact))
    for name, flow in flows.items():
        total = reached = failed = 0
        for f in reps:
            if len(ce.support(f)) < 2:
                reached += 1          # 게이트가 필요 없는 함수 (값 0)
                continue
            out = subprocess.run([abc, "-q", "read_truth %04x; %s" % (f, flow)], stdin=subprocess.DEVNULL,
                                 capture_output=True, text=True).stdout
            m = re.findall(r"nd =\s*(\d+)", out)
            if not m:
                failed += 1
                continue
            v = int(m[-1])
            total += v * size[f]
            reached += v <= exact[f]
        print("  %-30s 합 %7d (최솟값의 합보다 %+6.2f%%), 최솟값에 닿은 동치류 %3d/222, 실패 %d" % (
            name, total, 100.0 * (total - total_exact) / total_exact, reached, failed))


if __name__ == "__main__":
    main()
