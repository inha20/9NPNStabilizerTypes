#!/usr/bin/env python3
"""
circuit_metrics_analyze.py
==========================
circuit_metrics.py 가 기록한 함수별 비용(도구의 결과)을 B_4-궤도(402개)와 NPN 동치류(222개)로 묶어 본다.
비용 기준 일곱 가지(and, ai, g2, cmos, l2, l3, l4)와 최적화 흐름 둘(A, B)마다 다음을 낸다.

  직접      : 함수마다 도구로 직접 합성한 값의 합.
  일정한 묶음: 묶음(궤도 또는 NPN 동치류) 안에서 값이 모두 같은 묶음의 수와 (최대 - 최소)의 분포.
  대표 하나  : 묶음마다 회로를 하나만 두고, 나머지 함수는 그 회로에 입력 순열(배선, 비용 0)과
               입력 반전(아래 인버터 비용)을 붙여 얻을 때의 합. 대표는 이 합이 가장 작아지도록 고른다.
  모두 허용  : 묶음의 모든 함수의 회로를 쓸 수 있을 때, 함수마다 "어떤 회로 + 변환" 가운데 가장 싼 값의 합.
               자기 회로를 그대로 쓰는 경우를 포함하므로 직접 합성한 값보다 클 수 없다.
  더 싸지는 함수: "모두 허용"의 값이 직접 합성한 값보다 작은 함수의 수.
               도구가 그 함수를 직접 합성한 결과가 최적이 아님을 보여 주는 예이다.

인버터 비용: and 와 LUT 기준(l2, l3, l4)에서는 0, 게이트·셀 기준(ai, g2, cmos)에서는 인버터 하나에 1.
  AIG 에서는 반전이 간선의 표시일 뿐이고, LUT 는 입력 반전을 표 안으로 흡수하기 때문이다.
NPN 동치류에서는 출력 반전도 변환으로 허용하며 같은 인버터 비용을 한 번 더 센다.

실행: python scripts/circuit_metrics_analyze.py notes/circuit_metrics.txt.gz [궤도별 CSV 출력 파일]
준비: pip install numpy
"""
import gzip
import itertools
import sys
from collections import Counter

import numpy as np

N = 16
FULL = 1 << N
FLOWS = ["A", "B"]
METRICS = ["and", "ai", "g2", "cmos", "l2", "l3", "l4"]
INV_COST = {"and": 0, "ai": 1, "g2": 1, "cmos": 1, "l2": 0, "l3": 0, "l4": 0}


def load(path):
    vals = {(fl, m): np.full(FULL, -1, dtype=np.int64) for fl in FLOWS for m in METRICS}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            chunk, fl, m, digits = line.split()
            base = int(chunk) * 4096
            vals[(fl, m)][base:base + 4096] = np.frombuffer(bytes.fromhex(digits), dtype=np.uint8)
    for key, v in vals.items():
        assert (v >= 0).all(), "빠진 값이 있다: %s %s" % key
    return vals


# ---- B_4 의 원소 384개: 변수 순열 p 와 반전 마스크 a. 점 x 를 (p 로 비트를 옮긴 것) xor a 로 보낸다.
F = np.arange(FULL, dtype=np.uint32)


def point_map(p, a):
    out = []
    for x in range(N):
        y = 0
        for i in range(4):
            y |= ((x >> p[i]) & 1) << i
        out.append(y ^ a)
    return out


def function_perm(t):
    """함수 f 를, 점 t(x) 에서의 값이 f(x) 인 함수로 보내는 대응을 65,536개 전체에 대해 한 번에 계산한다."""
    h = np.zeros(FULL, dtype=np.uint32)
    for x in range(N):
        h |= ((F >> x) & 1) << t[x]
    return h


ELEMENTS = []   # (반전한 입력의 수, 함수들의 순열)
for p in itertools.permutations(range(4)):
    for a in range(N):
        ELEMENTS.append((bin(a).count("1"), function_perm(point_map(p, a))))
ELEMENTS.sort(key=lambda e: e[0])
WEIGHTS = np.array([w for w, _ in ELEMENTS], dtype=np.int64)
PERMS = np.stack([p for _, p in ELEMENTS]).astype(np.uint16)          # (384, 65536)
COMPLEMENT = (F ^ (FULL - 1)).astype(np.uint16)

canon_b4 = PERMS.min(axis=0)                                           # 궤도의 가장 작은 진리표
canon_npn = np.minimum(canon_b4, canon_b4[COMPLEMENT])


def distance_sums(perms, weights):
    """S[g] = (g 의 묶음의 함수 h 마다, g 를 h 로 옮기는 변환 가운데 인버터가 가장 적은 것의 인버터 수)의 합."""
    cols = np.ascontiguousarray(perms.T)                               # (65536, 변환 수), 변환은 무게 순
    out = np.zeros(FULL, dtype=np.int64)
    for g in range(FULL):
        _, first = np.unique(cols[g], return_index=True)
        out[g] = weights[first].sum()
    return out


def analyse(c, w, canon, perms, weights, dist_sum):
    labels_rep, labels = np.unique(canon, return_inverse=True)
    k = len(labels_rep)
    size = np.bincount(labels, minlength=k)
    lo = np.full(k, 10 ** 9, dtype=np.int64)
    hi = np.zeros(k, dtype=np.int64)
    np.minimum.at(lo, labels, c)
    np.maximum.at(hi, labels, c)
    # 모두 허용
    reuse = c.copy()
    for wt, perm in zip(weights, perms):
        cand = c + w * wt
        np.minimum.at(reuse, perm, cand)
    # 대표 하나
    total_g = size[labels] * c + w * dist_sum
    best = np.full(k, 10 ** 12, dtype=np.int64)
    np.minimum.at(best, labels, total_g)
    return {
        "groups": k, "direct": int(c.sum()), "constant": int((lo == hi).sum()),
        "spread": dict(sorted(Counter((hi - lo).tolist()).items())),
        "single": int(best.sum()), "any": int(reuse.sum()), "cheaper": int((reuse < c).sum()),
        "lo": lo, "hi": hi, "labels_rep": labels_rep, "size": size,
    }


def pct(a, b):
    return 100.0 * (a - b) / b if b else 0.0


def main():
    vals = load(sys.argv[1])
    print("B_4-궤도 %d개, NPN 동치류 %d개" % (len(np.unique(canon_b4)), len(np.unique(canon_npn))))

    # NPN 쪽 변환: 출력 반전까지. 무게는 (입력 인버터 수 + 출력 인버터 수).
    perms_npn = np.concatenate([PERMS, COMPLEMENT[PERMS]])
    weights_npn = np.concatenate([WEIGHTS, WEIGHTS + 1])
    order = np.argsort(weights_npn, kind="stable")
    perms_npn, weights_npn = perms_npn[order], weights_npn[order]
    dist_b4 = distance_sums(PERMS, WEIGHTS)
    dist_npn = distance_sums(perms_npn, weights_npn)

    per_orbit = {}
    for title, canon, perms, weights, dist in (
            ("B_4-궤도 (입력 순열과 입력 반전)", canon_b4, PERMS, WEIGHTS, dist_b4),
            ("NPN 동치류 (출력 반전까지)", canon_npn, perms_npn, weights_npn, dist_npn)):
        print()
        print("=" * 100)
        print(title)
        print("=" * 100)
        print("%-4s %-5s %10s %9s %10s %8s %10s %8s %9s  %s" % (
            "흐름", "기준", "직접", "일정한 묶음", "대표 하나", "(%)", "모두 허용", "(%)", "더 싸지는 함수", "(최대 - 최소)의 분포"))
        for fl in FLOWS:
            for m in METRICS:
                r = analyse(vals[(fl, m)], INV_COST[m], canon, perms, weights, dist)
                print("%-4s %-5s %10d %6d/%-3d %10d %+8.2f %10d %+8.2f %9d  %s" % (
                    fl, m, r["direct"], r["constant"], r["groups"], r["single"], pct(r["single"], r["direct"]),
                    r["any"], pct(r["any"], r["direct"]), r["cheaper"], r["spread"]))
                if canon is canon_b4:
                    per_orbit[(fl, m)] = r

    print()
    print("=" * 100)
    print("흐름 A 와 B 의 비교 (함수 65,536개)")
    print("=" * 100)
    print("%-5s %10s %10s %10s %10s %10s %12s" % ("기준", "A 의 합", "B 의 합", "B < A", "B = A", "B > A", "둘 중 작은 값의 합"))
    for m in METRICS:
        a, b = vals[("A", m)], vals[("B", m)]
        print("%-5s %10d %10d %10d %10d %10d %12d" % (
            m, a.sum(), b.sum(), (b < a).sum(), (b == a).sum(), (b > a).sum(), np.minimum(a, b).sum()))

    print()
    print("값의 분포 (흐름 B, 함수 수)")
    for m in METRICS:
        print("  %-5s %s" % (m, dict(sorted(Counter(vals[("B", m)].tolist()).items()))))

    if len(sys.argv) > 2:
        with open(sys.argv[2], "w") as out:
            out.write("대표 진리표,궤도 크기," + ",".join(
                "%s_%s_최소,%s_%s_최대" % (fl, m, fl, m) for fl in FLOWS for m in METRICS) + "\n")
            any_r = per_orbit[("A", "and")]
            for i in range(any_r["groups"]):
                row = ["%04x" % any_r["labels_rep"][i], str(any_r["size"][i])]
                for fl in FLOWS:
                    for m in METRICS:
                        r = per_orbit[(fl, m)]
                        row += [str(r["lo"][i]), str(r["hi"][i])]
                out.write(",".join(row) + "\n")
        print()
        print("궤도별 표를 %s 에 적었다." % sys.argv[2])


if __name__ == "__main__":
    main()
