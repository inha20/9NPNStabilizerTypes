#!/usr/bin/env python3
"""
merge_bn_class_tables.py
========================
나열 파일(모든 켤레류)과 classify_bn_classes.py 가 낸 .res 파일들(닫힌 군의 판정)을 합쳐 공개 표를 만든다.
.res 에 없는 류는 닫히지 않은 것(notclosed)이다. 합치기 전에 다음을 확인하고 어긋나면 멈춘다.
  - .res 의 번호가 나열 파일에 있고 겹치지 않으며, 위수와 꼭짓점 궤도 수 k 가 나열 파일과 같다.
  - 상태가 present 나 absent 뿐이다 (ERROR, unknown 이 있으면 멈춘다).
  - 선택: --expect 류,닫힘,나타남,안나타남 을 주면 개수가 같은지 본다 (n <= 6 은 표지 표 방법의 값과 견준다).

쓰는 것
  --classes 표.tsv[.gz]     모든 류: 번호, 위수, k, 상태, 생성원 수, 생성원들 (머리말 5줄)
  --absent  안나타남.tsv     상태가 absent 인 줄만

사용
  python3 -I scripts/merge_bn_class_tables.py N 나열파일 --res a.res b.res --classes 표.tsv.gz --absent 안나타남.tsv [--expect 55200,4598,4414,184]

2026-11-14. Claude가 쓴 코드이고 저자 검토 전이다. 2026-11-13에 n = 7 표를 만든 일회용 코드를 정리한 것이다.
"""
import argparse
import gzip
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bn_class_io import read_rows, summarize, describe   # noqa: E402


def header(n, total):
    nv = 1 << n
    return [
        "# n = %d: B_%d = Z_2^%d x| S_%d 의 부분군 켤레류 %s개 전부 (scripts/enumerate_bn_subgroups.g 로 나열, 켤레류마다 하나, OEIS A269890 의 개수와 같다)"
        % (n, n, n, n, format(total, ",")),
        "# 열: 번호, 위수, 꼭짓점 궤도 수 k, 상태, 생성원 수, 생성원들",
        "# 상태: present = 어떤 함수의 정확한 대칭군이 된다 / absent = 닫혀 있지만 어떤 함수의 대칭군도 못 된다 / notclosed = 닫혀 있지 않아 나타날 수 없다",
        "# 생성원 하나는 n + 1 = %d 개의 수: 꼭짓점 0 의 상, 꼭짓점 2^(i-1) 의 상 (i = 1..%d). 꼭짓점은 0..%d 이고 x -> pi(x) xor c 꼴이다 (scripts/export_bn_generators.g 와 같은 모양)"
        % (n + 1, n, nv - 1),
        "# 판정: scripts/np_cegar.py (present 는 대칭군이 정확히 이 군인 함수를 찾아 확인, absent 는 군 원소에서 나온 절의 UNSAT 으로 확인, notclosed 는 궤도 라벨의 대칭군이 더 큼)",
    ]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("n", type=int)
    ap.add_argument("enum_file")
    ap.add_argument("--res", nargs="+", required=True)
    ap.add_argument("--classes", required=True)
    ap.add_argument("--absent", required=True)
    ap.add_argument("--expect", default=None)
    a = ap.parse_args()

    enum = {}
    order_of_rows = []
    for j, order, k, _, ng, vals, _line in read_rows(a.enum_file):
        if j in enum:
            sys.exit("나열 파일에 번호 %d 가 두 번 있다" % j)
        enum[j] = (order, k, ng, vals)
        order_of_rows.append(j)
    status = {}
    for path in a.res:
        with open(path, encoding="utf-8") as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                j, order, k, st = int(p[0]), int(p[1]), int(p[2]), p[3]
                if st not in ("present", "absent"):
                    sys.exit("%s: 번호 %d 의 상태가 %s 이다" % (path, j, st))
                if j not in enum:
                    sys.exit("%s: 번호 %d 가 나열 파일에 없다" % (path, j))
                if j in status:
                    sys.exit("번호 %d 의 판정이 두 번 있다" % j)
                if (order, k) != enum[j][:2]:
                    sys.exit("번호 %d: 위수, k 가 나열 파일과 다르다 %s %s" % (j, (order, k), enum[j][:2]))
                status[j] = st

    hd = header(a.n, len(enum))
    rows = []
    out_all, out_abs = [], []
    for j in order_of_rows:
        order, k, ng, vals = enum[j]
        st = status.get(j, "notclosed")
        rows.append((j, order, k, st))
        line = "\t".join(str(x) for x in [j, order, k, st, ng] + vals)
        out_all.append(line)
        if st == "absent":
            out_abs.append(line)

    # 머리말: 안 나타남 표는 첫 줄만 다르다
    hd_abs = list(hd)
    hd_abs[0] = ("# n = %d: 닫혀 있지만 어떤 함수의 대칭군도 못 되는 B_%d 의 켤레류 %d개 (%s 에서 상태가 absent 인 줄만 옮김)"
                 % (a.n, a.n, len(out_abs), os.path.basename(a.classes)))

    text_all = "\n".join(hd + out_all) + "\n"
    text_abs = "\n".join(hd_abs + out_abs) + "\n"
    if a.classes.endswith(".gz"):
        with open(a.classes, "wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as g:   # 같은 내용이면 같은 바이트가 되도록
                g.write(text_all.encode("utf-8"))
    else:
        with open(a.classes, "w", encoding="utf-8") as f:
            f.write(text_all)
    with open(a.absent, "w", encoding="utf-8") as f:
        f.write(text_abs)

    s = summarize(rows)
    print("n=%d: %s" % (a.n, describe(s).replace("\n", "\n     ")))
    if a.expect:
        want = [int(t) for t in a.expect.split(",")]
        got = [s["classes"], s["closed"], s["present"], s["absent"]]
        if got != want:
            sys.exit("기대한 개수 (류, 닫힘, 나타남, 안 나타남) %s 와 다르다: %s" % (want, got))
        print("     기대한 개수(류, 닫힘, 나타남, 안 나타남) %s 와 같다" % want)


if __name__ == "__main__":
    main()
