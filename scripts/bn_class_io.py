#!/usr/bin/env python3
"""
bn_class_io.py
==============
B_n 의 부분군 켤레류를 담은 표(나열 파일, 공개 표)를 읽는 도우미이다. 다른 스크립트가 불러 쓴다.

두 가지 줄 모양을 읽는다.
  나열 파일 (scripts/enumerate_bn_subgroups.g 의 출력):  G <탭> j <탭> 위수 <탭> k <탭> 생성원 수 <탭> 생성원 ...
  공개 표 (notes/stabilizer_types_n7_*.tsv):               j <탭> 위수 <탭> k <탭> 상태 <탭> 생성원 수 <탭> 생성원 ...
생성원 하나는 n + 1 개의 수이다: 꼭짓점 0 의 상, 꼭짓점 2^(i-1) 의 상 (i = 1..n).
꼭짓점은 0..2^n - 1 이고 x -> pi(x) xor c 꼴이다. '#' 로 시작하는 줄은 머리말이다.
파일 이름이 .gz 로 끝나면 압축 파일로 읽는다.

2026-11-14. Claude가 쓴 코드이고 저자 검토 전이다. 2026-11-13의 임시 스크립트(run_enum.py, enum_check.py 의 load_vp)를 정리한 것이다.
"""
import gzip
from collections import Counter

import numpy as np


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, encoding="utf-8")


def gens_to_vp(cube, vals, ngens):
    """생성원 ngens 개(각 n+1 개의 수)를 꼭짓점 순열(길이 2^n 배열)의 목록으로 바꾼다."""
    n, X = cube.n, cube.X
    out = []
    for t in range(ngens):
        im = vals[t * (n + 1):(t + 1) * (n + 1)]
        c = im[0]
        vp = np.full(cube.V, c, dtype=np.int64)
        for i in range(n):
            vp ^= ((X >> i) & 1) * (im[1 + i] ^ c)
        out.append(vp.astype(cube.dt))
    return out


def read_rows(path):
    """줄 하나씩 (j, 위수, k, 상태 또는 None, 생성원 수, 생성원 수들의 목록, 원래 줄)을 낸다. 생성원은 아직 순열로 바꾸지 않는다."""
    with open_text(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if p[0] == "G":
                yield int(p[1]), int(p[2]), int(p[3]), None, int(p[4]), [int(t) for t in p[5:]], line
            else:
                yield int(p[0]), int(p[1]), int(p[2]), p[3], int(p[4]), [int(t) for t in p[5:]], line


def read_classes(path, cube, parts=1, part=0, statuses=None):
    """read_rows 와 같지만 생성원을 꼭짓점 순열로 바꿔서 (j, 위수, k, 상태, 생성원 순열 목록, 원래 줄)을 낸다.
    parts, part: 류를 parts 개로 나눠 part 번째(0부터)만 낸다(처리 순서의 번호 mod parts).
    statuses: 주어지면 그 상태인 줄만 낸다(공개 표에서)."""
    for idx, (j, order, k, status, ng, vals, line) in enumerate(read_rows(path)):
        if idx % parts != part:
            continue
        if statuses is not None and status not in statuses:
            continue
        yield j, order, k, status, gens_to_vp(cube, vals, ng), line


def summarize(rows):
    """공개 표의 (j, 위수, k, 상태) 목록에서 개수를 센다. 반환 dict."""
    rows = list(rows)
    st = Counter(r[3] for r in rows)
    ab = [r for r in rows if r[3] == "absent"]
    closed = [r for r in rows if r[3] in ("present", "absent")]
    kmax = max((r[2] for r in ab), default=0)
    above = [r for r in closed if r[2] > kmax]
    return dict(
        classes=len(rows), closed=len(closed), present=st["present"], absent=st["absent"], notclosed=st["notclosed"],
        absent_k=sorted(Counter(r[2] for r in ab).items()), kmax=kmax,
        closed_above=len(above), closed_above_present=sum(1 for r in above if r[3] == "present"),
    )


def describe(s):
    """summarize 의 결과를 한국어 두 줄로 쓴다."""
    return (
        "류 %d개 = 닫힘 %d (나타남 %d + 안 나타남 %d) + 닫히지 않음 %d\n"
        "안 나타나는 류의 꼭짓점 궤도 수 k 분포 %s, 최대 %d; k > %d 인 닫힌 류 %d개 가운데 나타나는 것 %d개"
        % (s["classes"], s["closed"], s["present"], s["absent"], s["notclosed"],
           "{" + ", ".join("%d:%d" % kv for kv in s["absent_k"]) + "}", s["kmax"],
           s["kmax"], s["closed_above"], s["closed_above_present"])
    )
