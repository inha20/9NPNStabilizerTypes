# identify_quotient_groups.g
# analyze_bn_absent_classes.py --gap-out 이 쓴 파일(QUOTIENTS)을 읽어, 몫 P = N(H)/H 가 궤도 집합 위에 일으키는 순열군이
# 어떤 군인지 알려진 군과 견준다. 사용(저장소 맨 위 폴더에서):
#   QFILE := "/tmp/q_n7.g"; Read("scripts/identify_quotient_groups.g");     # gap -q -A
# 견주는 군: W(D4) = C_2^3 : S_4 (짝수 개의 부호 바꿈과 좌표 순열; 8점 위의 군), C_2 x AGL(3,2) (위수 2688).
# 위수가 2000 이하이고 1024 가 아니면 SmallGroups 의 번호(IdGroup)도 쓴다.
# 2026-11-14. Claude가 쓴 코드이고 저자 검토 전이다.
if not IsBound(QFILE) then Error("QFILE 을 먼저 정하세요"); fi;
SizeScreen([2000, 24]);
Read(QFILE);
WD4 := Group((1,2)(3,4), (3,4)(5,6), (5,6)(7,8), (1,3)(2,4), (3,5)(4,6), (5,7)(6,8));
# AGL(3,2) = F_2^3 : GL(3,2), F_2^3 의 8개 벡터 위의 순열군으로 만든다
V3 := AsSSortedList(GF(2)^3);
AGL32 := Group(Concatenation(GeneratorsOfGroup(Action(GL(3,2), V3, OnRight)),
                             List(Basis(GF(2)^3), b -> PermList(List(V3, v -> Position(V3, v + b))))));
if Size(AGL32) <> 1344 then Error("AGL(3,2) 의 위수가 1344 가 아니다"); fi;
Z2xAGL32 := DirectProduct(CyclicGroup(IsPermGroup, 2), AGL32);
Print("W(D4): 위수 ", Size(WD4), ", IdGroup ", IdGroup(WD4), "\n");
Print("C2 x AGL(3,2): 위수 ", Size(Z2xAGL32), "\n");
for r in QUOTIENTS do
  P := Group(List(r.gens, PermList));
  Print("번호 ", r.j, ": |H|=", r.order, ", k=", r.k, ", |N(H)/H|=", r.NH, ", 순열군 P 의 위수 ", Size(P),
        ", 이행적인가 ", IsTransitive(P, [1..r.k]), ", 가해군인가 ", IsSolvableGroup(P));
  if Size(P) <> r.NH then Print("  [위수 불일치!]"); fi;
  if Size(P) <= 2000 and Size(P) <> 1024 then Print(", IdGroup ", IdGroup(P)); fi;
  if IsSolvableGroup(P) then Print(", 구조 ", StructureDescription(P)); fi;
  if Size(P) = Size(WD4) then Print(", W(D4) 와 동형인가 ", IsomorphismGroups(P, WD4) <> fail); fi;
  if Size(P) = Size(Z2xAGL32) then Print(", C2 x AGL(3,2) 와 동형인가 ", IsomorphismGroups(P, Z2xAGL32) <> fail); fi;
  Print("\n");
od;
QUIT;
