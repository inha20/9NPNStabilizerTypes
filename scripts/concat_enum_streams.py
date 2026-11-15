#!/usr/bin/env python3
"""
concat_enum_streams.py
======================
scripts/enumerate_bn_subgroups.g 를 여러 프로세스로 나눠 돌린 출력 파일들(나열 파일)을 주어진 순서로 이어 붙이고
류 번호 j 를 1, 2, 3, ... 으로 다시 매긴다. 같은 류가 두 파일에 나오지 않는지(위수, k, 생성원이 같은 줄이 겹치는지)도 본다.

사용: python3 -I scripts/concat_enum_streams.py 합친파일.tsv 줄기1.tsv 줄기2.tsv ...
2026-11-13 의 n = 7 나열은 줄기 0, 1, 2, 3, A, B, C, D 의 순서로 이어 붙였다(합쳐서 55,200류, OEIS A269890 와 같다).

2026-11-14. Claude가 쓴 코드이고 저자 검토 전이다. 2026-11-13에 손으로 한 일을 스크립트로 옮겼다.
"""
import sys


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    out, ins = sys.argv[1], sys.argv[2:]
    seen = set()
    j = 0
    with open(out, "w", encoding="utf-8") as fo:
        for path in ins:
            with open(path, encoding="utf-8") as f:
                for line in f:
                    if not line.startswith("G\t"):
                        continue
                    p = line.rstrip("\n").split("\t")
                    key = tuple(p[2:])
                    if key in seen:
                        sys.exit("%s: 같은 류가 두 번 나온다: %s" % (path, key[:3]))
                    seen.add(key)
                    j += 1
                    fo.write("\t".join(["G", str(j)] + p[2:]) + "\n")
    print("%d류를 %s 에 썼다" % (j, out))


if __name__ == "__main__":
    main()
