# check_s6_outer_automorphism.g
# n = 6 에서, S_6 의 부분군 L 이 완전 매칭(변수 쌍 분할) 15개에 만드는 궤도 수가
# "L 을 S_6 의 외부 자기동형으로 옮긴 군이 변수 쌍 15개에 만드는 궤도 수"와 같은지 56개 켤레류 전부에서 본다.
#
# 외부 자기동형은 Sylvester 의 구성으로 만든다: 완전 매칭 15개 가운데 서로 쌍을 공유하지 않는 5개의 묶음
# (synthematic total)이 정확히 6개 있고, S_6 가 이 6개에 작용하는 방식이 외부 자기동형이다.
# 실행: gap -q scripts/check_s6_outer_automorphism.g
SizeScreen([4096, 24]);
G := SymmetricGroup(6);
duads := Combinations([1..6], 2);
synthemes := Orbit(G, Set([[1,2],[3,4],[5,6]]), OnSetsSets);
# total: 15개 쌍을 빠짐없이 한 번씩 덮는 매칭 5개
totals := Filtered(Combinations(synthemes, 5), t -> Length(Union(t)) = 15);
Print("쌍 ", Length(duads), "개, 완전 매칭 ", Length(synthemes), "개, 묶음(total) ", Length(totals), "개\n");
OnTotals := function(t, g) return Set(List(t, s -> OnSetsSets(s, g))); end;
totals := List(totals, Set);
theta := ActionHomomorphism(G, totals, OnTotals);
Print("묶음 6개 위의 작용이 S_6 전체인가: ", Size(Image(theta)) = 720,
      ",  맞바꿈 (1,2) 의 상의 순환 구조: ", CycleStructurePerm((1,2)^theta), "\n");
bad := 0;  moved := 0;
for C in ConjugacyClassesSubgroups(G) do
  L := Representative(C);
  a := Length(Orbits(L, synthemes, OnSetsSets));
  b := Length(Orbits(Image(theta, L), duads, OnSets));
  c := Length(Orbits(L, duads, OnSets));
  if a <> b then bad := bad + 1; fi;
  if a <> c then moved := moved + 1; fi;
od;
Print("켤레류 56개 가운데 (매칭 위의 궤도 수) <> (옮긴 군의 쌍 위의 궤도 수) 인 것: ", bad, "개\n");
Print("참고: 옮기지 않고 비교하면 (매칭 위의 궤도 수) <> (쌍 위의 궤도 수) 인 것: ", moved, "개\n");
QUIT;
