#!/bin/bash
# run_bn_class_pipeline.sh
# B_n 부분군 켤레류의 나열 -> 판정 -> 표 -> 분석을 n = 2..6 에 처음부터 돌리고(표지 표 방법의 값과 견준다),
# n = 7 은 공개 표(notes/stabilizer_types_n7_classes.tsv.gz)를 분석한다. 마지막에 안 나타나는 군 가운데 꼭짓점 궤도가 가장 많은 군을
# 전수로 확인하고(n = 6, 7), 몫 군의 정체를 GAP 으로 확인한다. 화면 출력은 notes/check_bn_class_pipeline_output.txt 에 있다.
#
# 사용 (저장소 맨 위 폴더에서):  GAP=gap bash scripts/run_bn_class_pipeline.sh 작업폴더
#   필요: GAP 4.x (환경 변수 GAP 이 실행 파일, 기본 gap), python3 + numpy, scipy, python-sat(pysat)
#   시간: 코어 둘에 약 4분. 쓰는 곳: 작업 폴더(나열 파일, 판정, 표)와 notes/stabilizer_types_n7_closed_absent_annotated.tsv (n = 7 분석표).
# n = 7 의 나열부터 다시 하려면: scripts/enumerate_bn_stream.g 를 PARTS=8 등으로 조각 내어 돌리고,
#   scripts/concat_enum_streams.py 로 이어 붙이고, scripts/classify_bn_classes.py 와 scripts/merge_bn_class_tables.py 를 n = 2..6 과 같이 쓴다.
#   시간: 11-13~14 에는 정리 전의 임시 판으로 두 번에 나눠 돌렸다. notes/enumerate_bn_n7_log.txt 에 기록된 부분공간 류 148개의 시간 합이 코어 하나 기준 약 49분이고
#   가장 오래 걸린 류가 약 18분이다(중간에 멈춘 실행에서 잃은 시간은 기록에 없다). 정리한 판으로 n = 7 을 처음부터 다시 돌려 시간을 재지는 않았다.
#
# 2026-11-14. Claude가 쓴 코드이고 저자 검토 전이다.
set -e
W=${1:-/tmp/bn_pipeline}
GAP=${GAP:-gap}
mkdir -p "$W"
T7=notes/stabilizer_types_n7_classes.tsv.gz
# 표지 표 방법(scripts/stabilizer_types.g)이 낸 값: 류, 닫힘, 나타남, 안 나타남 (notes/2026-11-12_stabilizer_types.md 2.1절의 표)
declare -A EXP=( [2]=8,5,4,1 [3]=33,14,10,4 [4]=193,49,34,15 [5]=953,171,143,28 [6]=7440,874,790,84 )

for n in 2 3 4 5 6; do
  echo "=== n = $n: 나열, 판정, 표, 정규화군 몫 분석"
  cat > "$W/drv_$n.g" <<EOG
N_ENUM := $n; KMIN := 1; OUTFILE := "$W/enum_n$n.tsv";
Read("scripts/enumerate_bn_subgroups.g");
QUIT;
EOG
  $GAP -A -q -o 3g "$W/drv_$n.g" < /dev/null > "$W/gap_$n.log" 2>&1
  python3 -I scripts/classify_bn_classes.py $n "$W/enum_n$n.tsv" "$W/cls_n${n}a" --parts 2 --part 0 > /dev/null & p1=$!
  python3 -I scripts/classify_bn_classes.py $n "$W/enum_n$n.tsv" "$W/cls_n${n}b" --parts 2 --part 1 > /dev/null & p2=$!
  wait $p1; wait $p2
  python3 -I scripts/merge_bn_class_tables.py $n "$W/enum_n$n.tsv" --res "$W/cls_n${n}a.res" "$W/cls_n${n}b.res" \
      --classes "$W/classes_n$n.tsv" --absent "$W/absent_n$n.tsv" --expect "${EXP[$n]}"
  python3 -I scripts/analyze_bn_absent_classes.py $n "$W/classes_n$n.tsv" "$W/annot_n$n.tsv" | sed -n '3,7p'
done

echo "=== n = 7: 공개 표 $T7 의 정규화군 몫 분석"
python3 -I scripts/analyze_bn_absent_classes.py 7 $T7 notes/stabilizer_types_n7_closed_absent_annotated.tsv

echo "=== 안 나타나는 군 가운데 꼭짓점 궤도가 가장 많은 군을 전수로 확인 (n = 6)"
python3 -I scripts/check_bn_absent_certificate.py 6 "$W/classes_n6.tsv"
echo "=== 같은 확인 (n = 7)"
python3 -I scripts/check_bn_absent_certificate.py 7 $T7

for n in 6 7; do
  T="$W/classes_n6.tsv"; [ $n = 7 ] && T=$T7
  echo "=== 몫 군 N(H)/H 의 정체 (n = $n, 꼭짓점 궤도 12개 이상인 안 나타나는 류)"
  python3 -I scripts/analyze_bn_absent_classes.py $n $T --kmin 12 --gap-out "$W/q_n$n.g" > /dev/null
  printf 'QFILE := "%s";\nRead("scripts/identify_quotient_groups.g");\n' "$W/q_n$n.g" > "$W/id_n$n.g"
  $GAP -A -q -o 2g "$W/id_n$n.g" < /dev/null
done
echo "=== 끝"
