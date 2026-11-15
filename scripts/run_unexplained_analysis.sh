#!/bin/bash
# run_unexplained_analysis.sh
# 닫혔지만 안 나타나는 류에 대해 (1) 정규화군 몫 판정법이 설명하지 못하는 류의 구조, (2) 병합 그래프 판정법의 적용 범위,
# (3) SAT 를 쓰지 않는 색칠 전수 확인, (4) 자유 색칠의 안정화군 구조(GAP)를 낸다.
# n = 7 은 공개 표를 쓰고, n = 2..6 은 scripts/run_bn_class_pipeline.sh 의 작업 폴더에 있는 표를 쓴다(그 폴더를 주었을 때만).
# 화면 출력은 notes/unexplained_output.txt 에 있다(n = 2..6, 7 을 모두 돌린 것).
# 판정법의 정의와 논증은 scripts/analyze_unexplained_classes.py 의 머리말에 있고, 노트는 notes/2026-11-12_stabilizer_types.md 3절이다.
#
# 사용 (저장소 맨 위 폴더에서):  GAP=gap bash scripts/run_unexplained_analysis.sh 작업폴더 [파이프라인작업폴더]
#   파이프라인작업폴더 = scripts/run_bn_class_pipeline.sh 에 준 작업 폴더(classes_n*.tsv 와 annot_n*.tsv 가 들어 있다). 없으면 n = 7 만 돈다.
#   필요: GAP 4.x (환경 변수 GAP 이 실행 파일, 기본 gap), python3 + numpy, scipy (python-sat 은 쓰지 않는다)
#   시간: 코어 하나에 약 2분(2026-11-15 측정; n = 7 의 파이썬 부분 약 100초, n = 2..6 과 GAP 부분은 합쳐서 몇 초), 메모리 최대 약 0.5GB.
#
# 2026-11-15. Claude가 쓴 코드이고 저자 검토 전이다.
set -e
W=${1:-/tmp/unexplained}
PIPE=${2:-}
GAP=${GAP:-gap}
mkdir -p "$W"

run_n() {   # run_n n 표 주석표
  local n=$1 T=$2 A=$3
  echo "######## n = $n"
  echo "=== 설명 못 하는 류: 모든 색칠의 안정화군, 극소 실현 분할, 쌍 병합 그래프"
  python3 -I scripts/analyze_unexplained_classes.py $n $T $A
  echo "=== 병합 그래프 판정법을 안 나타나는 류 전체에 적용"
  python3 -I scripts/analyze_unexplained_classes.py $n $T $A --pair-graph-all
  echo "=== SAT 를 쓰지 않는 확인: 안 나타나는 모든 류의 모든 2-색칠(보색 제외)의 안정화군이 H 보다 큰가"
  python3 -I scripts/analyze_unexplained_classes.py $n $T $A --verify-all --kmax 16
  echo "=== 자유 색칠의 안정화군 K 의 구조 (GAP)"
  python3 -I scripts/export_unexplained_for_gap.py $n $T $A "$W/unexplained_n$n.g"
  printf 'DATAFILE := "%s";\nRead("scripts/identify_unexplained_overgroups.g");\n' "$W/unexplained_n$n.g" > "$W/id_n$n.g"
  $GAP -A -q -o 2g "$W/id_n$n.g" < /dev/null
}

if [ -n "$PIPE" ]; then
  for n in 2 3 4 5 6; do
    run_n $n "$PIPE/classes_n$n.tsv" "$PIPE/annot_n$n.tsv"
  done
fi
run_n 7 notes/stabilizer_types_n7_classes.tsv.gz notes/stabilizer_types_n7_closed_absent_annotated.tsv
echo "=== 끝"
