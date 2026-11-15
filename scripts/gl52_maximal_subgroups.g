# GL(5,2) 의 극대부분군 켤레류 (5.3절의 순수 평행이동군 논증이 쓴다).
# 기대: 5 류, 지수 31, 31 (점 안정화군, 초평면 안정화군), 155, 155 (직선, 평면 안정화군), 64512 (싱어 정규화군 31:5).
# 실행: gap -A -q -o 3g scripts/gl52_maximal_subgroups.g   (1초)
# 2026-11-15. Claude가 쓴 코드이고 저자 검토 전이다.
G := GL(5, 2);
iso := IsomorphismPermGroup(G);
P := Image(iso);
Print("|GL(5,2)| = ", Size(P), ", 차수(31개 0 아닌 벡터) ", NrMovedPoints(P), "\n");
mx := MaximalSubgroupClassReps(P);
Print("극대부분군 켤레류: ", Length(mx), "개\n");
info := List(mx, H -> [Index(P, H), Size(H), SortedList(List(Orbits(H, MovedPoints(P)), Length))]);
Sort(info);
for x in info do
  Print("  지수 = ", x[1], "  |H| = ", x[2], "  31개 점 위의 궤도 길이: ", x[3], "\n");
od;
QUIT;
