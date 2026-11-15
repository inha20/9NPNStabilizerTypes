#!/usr/bin/env python3
"""
classify_bn_classes.py
======================
나열 파일(scripts/enumerate_bn_subgroups.g 의 출력)의 B_n 부분군 켤레류마다 scripts/np_cegar.py 의 appears 로 상태를 정한다.
  present   : 어떤 함수의 정확한 대칭군이 된다 (그런 함수를 찾아 대칭군을 정확히 세어 확인)
  absent    : 닫혀 있지만 어떤 함수의 대칭군도 못 된다 (군 원소에서 나온 절들의 모순으로 확인)
  notclosed : 닫혀 있지 않아 나타날 수 없다 (궤도 라벨 함수의 대칭군이 더 크다)

쓰는 것
  <접두사>.res     닫힌 군만: 번호 <탭> 위수 <탭> k <탭> 상태 <탭> 반복 수 <탭> 절 수
  <접두사>.absent  안 나타나는 군의 나열 줄 그대로

사용
  python3 -I scripts/classify_bn_classes.py N 나열파일 접두사 [--parts P --part I]
  --parts P --part I : 류를 P 몫으로 나눠 I 번째(0부터)만 처리한다. 코어가 둘이면 두 프로세스를 따로 돌린다
                       (예: --parts 2 --part 0 / --parts 2 --part 1). 결과 .res 파일들은 merge_bn_class_tables.py 가 합친다.
  n = 7(55,200류)은 코어 둘로 약 3분이 걸렸다(2026-11-13).

2026-11-14. Claude가 쓴 코드이고 저자 검토 전이다. 2026-11-13의 임시 스크립트 run_enum.py 를 정리한 것이다.
"""
import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import np_cegar as C            # noqa: E402
from bn_class_io import read_classes   # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("n", type=int)
    ap.add_argument("enum_file")
    ap.add_argument("prefix")
    ap.add_argument("--parts", type=int, default=1)
    ap.add_argument("--part", type=int, default=0)
    a = ap.parse_args()
    cube = C.Cube(a.n)
    t0 = time.time()
    cnt = closed = 0
    with open(a.prefix + ".res", "w") as fo, open(a.prefix + ".absent", "w") as fa:
        for j, order, k, _, gens, line in read_classes(a.enum_file, cube, a.parts, a.part):
            cnt += 1
            try:
                r = C.appears(cube, gens, order=order)
            except Exception as e:                    # 판정 실패는 숨기지 않고 표시한다 (병합 단계가 멈춘다)
                fo.write("%d\t%d\t%d\tERROR\t%s\n" % (j, order, k, str(e)[:60]))
                continue
            if r["status"] != "notclosed":
                closed += 1
                fo.write("%d\t%d\t%d\t%s\t%s\t%s\n" % (j, order, k, r["status"], r.get("iters"), r.get("clauses")))
                if r["status"] != "present":
                    fa.write(line)
            if cnt % 5000 == 0:
                fo.flush()
                fa.flush()
                print("part %d: %d classes, closed %d, %.0f s" % (a.part, cnt, closed, time.time() - t0), flush=True)
    print("done part %d/%d: %d classes, closed %d, %.0f s" % (a.part, a.parts, cnt, closed, time.time() - t0))


if __name__ == "__main__":
    main()
