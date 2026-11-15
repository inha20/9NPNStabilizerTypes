# export_bn_generators.g
# B_n 의 부분군 켤레류 L_j (표지 표의 대표원)의 생성원을 부호 있는 순열 (sigma, c) 로 읽을 수 있는 수로 내보낸다.
# scripts/np_appearance.py 와 scripts/verify_npn_class_identities.py 가 표지 표 없이 부분군을 직접 다루는 데 쓴다.
#
# 실행 (n = 6 은 저장해 둔 작업 공간 `gap -L tom6.ws` 로 PRE 를 쓴다):
#   N_EXPORT := 6; OUTFILE := "/tmp/bn_gens_n6.tsv"; JLIST := [...]; Read("scripts/export_bn_generators.g");
# 출력:  G  j  위수  꼭짓점 궤도 수  생성원 수  g1 의 7개 수  g2 의 7개 수 ...
#   생성원 하나는 꼭짓점 0 의 상, 꼭짓점 2^(i-1) 의 상 (i = 1..n) 의 n+1 개 수 (모두 0 부터 센 꼭짓점 번호)다.
#   원소 (sigma, c) 는 x -> sigma(x) xor c 이므로 c = (0 의 상), sigma(e_i) = (e_i 의 상) xor c 로 되돌린다.

if not IsBound(N_EXPORT) then N_EXPORT := 4; fi;
if not IsBound(OUTFILE) then OUTFILE := "/tmp/bn_gens.tsv"; fi;
NO_MAIN := true;
Read("scripts/npn_pair_resolution.g");   # MakeBn 을 쓴다
Unbind(NO_MAIN);

ExportGenerators := function(n)
  local R, tom, G, m, j, L, gens, g, str, jlist, t0;
  t0 := Runtime();
  if IsBound(PRE) and PRE.n = n and IsBound(PRE.tom) then
    R := PRE;  tom := R.tom;
  else
    R := MakeBn(n);
    tom := TableOfMarks(R.G);
  fi;
  m := Length(MarksTom(tom));
  if IsBound(JLIST) then jlist := JLIST; else jlist := [1..m]; fi;
  PrintTo(OUTFILE, "# n = ", n, ", 부분군 켤레류 ", m, "개의 생성원\n");
  for j in jlist do
    L := RepresentativeTom(tom, j);
    if Size(L) = 1 then
      gens := [];
    else
      gens := SmallGeneratingSet(L);
    fi;
    str := "";
    for g in gens do
      str := Concatenation(str, "\t", String((1^g) - 1),
                           JoinStringsWithSeparator(List([1..n], i -> Concatenation("\t", String(((2^(i-1) + 1)^g) - 1))), ""));
    od;
    AppendTo(OUTFILE, "G\t", j, "\t", Size(L), "\t", Length(Orbits(L, [1..2^n])), "\t", Length(gens), str, "\n");
  od;
  Print("끝 (", Int((Runtime() - t0) / 1000), " 초)\n");
end;

ExportGenerators(N_EXPORT);
