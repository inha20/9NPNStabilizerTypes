# identify_unexplained_overgroups.g
# export_unexplained_for_gap.py 가 쓴 파일을 읽어, 정규화군 몫 판정법이 설명하지 못하는 류마다
# 자유 색칠 c 의 안정화군 K = Stab_{B_n}(f_c) 의 구조를 적는다.
#   - K 의 위수가 파이썬이 센 값과 같은가, H <= K 인가, H 가 K 에서 정규가 아닌가(자유 색칠이면 K 와 N(H)의 교집합이 H 이므로 정규일 수 없다)를 확인한다.
#   - K 가 H 의 오른쪽 잉여류들에 일으키는 작용의 상(K/핵, 핵은 H 의 K 에서의 정규핵)의 StructureDescription 을 적는다.
#     지수 [K:H] 가 m 이면 이 작용은 m 개의 점 위의 추이적 작용이고 H/핵은 점 하나의 고정군이다.
# 사용: DATAFILE := "...g"; Read("scripts/identify_unexplained_overgroups.g");   (scripts/run_unexplained_analysis.sh 가 부른다)
#
# 2026-11-15. Claude가 쓴 코드이고 저자 검토 전이다.
Read(DATAFILE);
Bn := Group(List(BNGENS, PermList));
if Size(Bn) <> 2^NBITS * Factorial(NBITS) then Error("B_n 의 위수가 다르다"); fi;
key := fail; tally := []; pos := fail;     # 람다식이 쓰는 전역 변수를 미리 선언한다(경고를 피하려고)
for r in DATA do
  H := Group(List(r.gens, PermList));
  if Size(H) <> r.order then Error("H 의 위수가 다르다"); fi;
  Print("j=", r.j, " |H|=", r.order, " k=", r.k, " |N/H|=", r.NH, " 자유 색칠 ", Length(r.free), "개\n");
  tally := [];
  for f in r.free do
    K := Stabilizer(Bn, f.S, OnSets);
    if Size(K) <> f.sK then Error("K 의 위수가 다르다"); fi;
    if not IsSubgroup(K, H) then Error("H 가 K 에 들어 있지 않다"); fi;
    if IsNormal(K, H) then Error("H 가 K 에서 정규이다(자유 색칠이면 일어나지 않는다)"); fi;
    act := Action(K, RightCosets(K, H), OnRight);
    if Size(act) <= 5000 then
      sd := StructureDescription(act : short);
    else
      sd := Concatenation("위수 ", String(Size(act)));
    fi;
    key := [Index(K, H), Size(K) / Size(act), Size(act), sd];
    pos := PositionProperty(tally, t -> t[1] = key);
    if pos = fail then Add(tally, [key, 1]); else tally[pos][2] := tally[pos][2] + 1; fi;
  od;
  for t in tally do
    Print("   ", t[2], "개: [K:H]=", t[1][1], ", |핵|=", t[1][2], ", |K/핵|=", t[1][3], ", K/핵 = ", t[1][4], "\n");
  od;
od;
QUIT;
