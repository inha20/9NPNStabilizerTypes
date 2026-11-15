# export_bn_class_data.g
# B_n 의 켤레류와, B_n 의 부분군 켤레류 L_j 마다 "L_j 의 원소가 B_n 의 어느 켤레류에 몇 개 들어 있는가"를 내보낸다.
# scripts/verify_npn_class_identities.py 가 이것으로 G_n = B_n x Z_2 의 켤레류별 고정점 항등식을 세운다.
#
# 실행 (n <= 5):  N_EXPORT := 5; OUTFILE := "/tmp/bn_class_n5.tsv"; Read("scripts/export_bn_class_data.g");
# n = 6: 표지 표를 저장해 둔 작업 공간(tom6.ws, scripts/npn_pair_resolution.g 머리말 참고)을 `gap -L tom6.ws` 로 읽어 PRE 를 쓴다.
#        JLIST := [...]; 를 정하면 그 번호의 류만 내보낸다 (없으면 모두).
#
# 출력 (탭으로 나눈다):
#   B  b  류 크기  꼭짓점 위 순환 수  모든 순환이 짝수 길이인가(0/1)  양의 순환 길이들  음의 순환 길이들
#   L  j  위수  꼭짓점 궤도 수  aO  b1:개수1 b2:개수2 ...        (aO = H 가 L_j 와 켤레인 B_n-궤도 수, 번호 대조용)
# B_n 의 원소 (sigma, c): x -> sigma(x) xor c.  켤레류는 순환의 길이와 부호(순환 위 c 의 합의 홀짝)로 정해진다.

if not IsBound(N_EXPORT) then N_EXPORT := 4; fi;
if not IsBound(OUTFILE) then OUTFILE := "/tmp/bn_class_data.tsv"; fi;
NO_MAIN := true;
Read("scripts/npn_pair_resolution.g");   # MakeBn 을 쓴다
Unbind(NO_MAIN);

ExportClassData := function(n)
  local R, tom, G, m, bits, ClassKey, ccl, keys, sizes, info, i, c, g, key, pos, Orbs, jlist, j, L, cl, cnt, k,
        str, fixO, aO, t0, nb, rnd, key2, cyc;
  t0 := Runtime();
  if IsBound(PRE) and PRE.n = n and IsBound(PRE.tom) then
    R := PRE;  tom := R.tom;
    Print("n = ", n, ": 저장해 둔 표지 표를 쓴다\n");
  else
    R := MakeBn(n);
    tom := TableOfMarks(R.G);
  fi;
  G := UnderlyingGroup(tom);
  m := Length(MarksTom(tom));
  Print("n = ", n, ": |B_n| = ", Size(G), ", 부분군 켤레류 ", m, "개\n");

  bits := function(x)
    local l;
    l := CoefficientsQadic(x, 2);
    return Concatenation(l, ListWithIdenticalEntries(n - Length(l), 0));
  end;
  # g 의 켤레류 열쇠 [양의 순환 길이들, 음의 순환 길이들]
  ClassKey := function(g)
    local cb, pi, i, y, e, cycs, pos, neg, cy, s;
    cb := bits((1^g) - 1);
    pi := [];
    for i in [1..n] do
      y := bits(((2^(i-1) + 1)^g) - 1);
      e := List([1..n], k -> (y[k] + cb[k]) mod 2);
      if Sum(e) <> 1 then Error("B_n 의 원소가 아니다"); fi;
      Add(pi, Position(e, 1));
    od;
    cycs := Cycles(PermList(pi), [1..n]);
    pos := [];  neg := [];
    for cy in cycs do
      s := Sum(cy, k -> cb[k]) mod 2;
      if s = 0 then Add(pos, Length(cy)); else Add(neg, Length(cy)); fi;
    od;
    return [SortedList(pos), SortedList(neg)];
  end;

  ccl := ConjugacyClasses(G);
  keys := [];  sizes := [];  info := [];
  for c in ccl do
    g := Representative(c);
    key := ClassKey(g);
    if key in keys then Error("두 켤레류가 같은 열쇠를 가진다"); fi;
    Add(keys, key);  Add(sizes, Size(c));
    Orbs := Orbits(Group(g), [1..2^n]);
    Add(info, [Length(Orbs), ForAll(Orbs, o -> IsEvenInt(Length(o)))]);
    # 같은 켤레류의 임의 원소가 같은 열쇠를 갖는지 점검
    for i in [1..20] do
      if ClassKey(Random(c)) <> key then Error("열쇠가 켤레류를 정하지 않는다"); fi;
    od;
  od;
  nb := Length(keys);
  if Sum(sizes) <> Size(G) then Error("류 크기의 합이 |B_n| 이 아니다"); fi;
  Print("   B_n 의 켤레류 ", nb, "개\n");

  PrintTo(OUTFILE, "# n = ", n, ", |B_n| = ", Size(G), ", B_n 의 켤레류 ", nb, "개, 부분군 켤레류 ", m, "개\n");
  for i in [1..nb] do
    AppendTo(OUTFILE, "B\t", i, "\t", sizes[i], "\t", info[i][1], "\t", Number([info[i][2]], x -> x), "\t",
             keys[i][1], "\t", keys[i][2], "\n");
  od;

  fixO := List([1..m], j -> 2^Length(Orbits(RepresentativeTom(tom, j), [1..2^n])));
  aO := DecomposedFixedPointVector(tom, fixO);
  aO := Concatenation(aO, ListWithIdenticalEntries(m - Length(aO), 0));

  if IsBound(JLIST) then jlist := JLIST; else jlist := [1..m]; fi;
  for j in jlist do
    L := RepresentativeTom(tom, j);
    cnt := ListWithIdenticalEntries(nb, 0);
    for c in ConjugacyClasses(L) do
      k := Position(keys, ClassKey(Representative(c)));
      cnt[k] := cnt[k] + Size(c);
    od;
    if Sum(cnt) <> Size(L) then Error("L 의 원소 수가 맞지 않는다"); fi;
    str := JoinStringsWithSeparator(List(Filtered([1..nb], k -> cnt[k] > 0), k -> Concatenation(String(k), ":", String(cnt[k]))), " ");
    AppendTo(OUTFILE, "L\t", j, "\t", Size(L), "\t", Length(Orbits(L, [1..2^n])), "\t", aO[j], "\t", str, "\n");
    if j mod 500 = 0 then Print("   ... ", j, " / ", m, " (", Int((Runtime() - t0) / 1000), " 초)\n"); fi;
  od;
  Print("   끝 (", Int((Runtime() - t0) / 1000), " 초)\n");
end;

ExportClassData(N_EXPORT);
