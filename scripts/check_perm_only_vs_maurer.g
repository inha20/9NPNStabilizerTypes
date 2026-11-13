# check_perm_only_vs_maurer.g
# 변수 순열만 허용했을 때(입력 반전 없음), 순열 대칭군 S(f) 의 켤레류별로 함수가 몇 개인지를 표지 표로 센다.
#
# 목적: Maurer (2015, SpringerPlus 4:418)가 4변수(그림 27)와 5변수(그림 28)에서 전수 계산으로 낸 표와,
#       R9에서 쓴 것과 같은 표지 표 계산이 같은 값을 내는지 대조한다. 같으면 R9의 계산법에 대한 독립된 확인이 된다.
#       n = 6 은 Maurer 의 표에 없는 값이므로 대조 대상이 아니고, 합이 2^64 인지만 본다.
#
# 방법: 대칭군 G = S_n 이 정점 2^n 개에 작용한다. 켤레류 j 의 군 H_j 가 고정하는 함수는 2^(H_j 의 정점 궤도 수)개이다.
#       이 고정점 벡터를 G 의 표지 표로 분해하면, 완전 대칭군이 H_j 와 켤레인 함수의 수가
#       (분해 성분) x |G| / |H_j| 로 나온다. (같은 식을 B_n 에 쓴 것이 stabilizer_types.g 이다.)
#
# 실행: gap -q scripts/check_perm_only_vs_maurer.g          (n = 4, 5, 6)
#       n 을 골라 돌리려면 파일 앞에 NS := [4]; 처럼 먼저 지정한다.
# 출력은 notes/check_perm_only_vs_maurer_output.txt 에 있다.
SizeScreen([4096, 24]);
CheckPermOnly := function(n)
  local vert, toPerm, gens, G, tom, m, fix, j, H, dec, res;
  vert := Tuples([0,1], n);
  toPerm := function(p) return PermList(List(vert, v -> Position(vert, Permuted(v, p)))); end;
  gens := List(GeneratorsOfGroup(SymmetricGroup(n)), toPerm);
  G := Group(gens);
  tom := TableOfMarks(G);
  m := Length(OrdersTom(tom));
  fix := [];
  for j in [1..m] do
    H := RepresentativeTom(tom, j);
    Add(fix, 2^Length(Orbits(H, [1..2^n])));
  od;
  dec := DecomposedFixedPointVector(tom, fix);
  res := [];
  for j in [1..m] do
    if dec[j] > 0 then
      H := RepresentativeTom(tom, j);
      Add(res, [dec[j] * Size(G) / Size(H), Size(H), dec[j], StructureDescription(H)]);
    fi;
  od;
  Sort(res, function(a, b) return a[1] > b[1]; end);
  return [m, res];
end;
if not IsBound(NS) then NS := [4, 5, 6]; fi;
for n in NS do
  r := CheckPermOnly(n);
  Print("n = ", n, ": S_n의 부분군 켤레류 ", r[1], "개 중 순열 대칭군으로 나타나는 것 ", Length(r[2]), "개, 함수 수의 합 ", Sum(r[2], x -> x[1]), " (2^(2^n) = ", 2^(2^n), ")\n");
  for x in r[2] do
    Print("   함수 ", x[1], "개  군 위수 ", x[2], "  궤도 수 ", x[3], "  ", x[4], "\n");
  od;
od;
QUIT;
