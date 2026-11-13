# npn_stabilizer_types.g
# 부울 함수의 NPN 안정화군(입력 반전, 변수 순열, 출력 반전이 이루는 B_n x Z_2 안에서)의 종류를 모두 적는다.
#
#   B_n x Z_2 는 (g, e) 가 함수 f 를 (g, e) f = (f o g^-1) xor e 로 바꾸며 작용한다 (e = 1 이면 출력을 뒤집는다).
#   K(f) = { (g, e) : (g, e) f = f } 를 f 의 NPN 안정화군이라 한다. 두 갈래가 나온다.
#     - H(f) = { g : f o g = f }                          (e = 0, 대칭. stabilizer_types.g 의 H(f))
#     - A(f) = { g : f o g = not f }                      (e = 1, 반대칭)
#   A(f) 는 비어 있거나 H(f) 의 잉여류이다. Kbar(f) = H(f) 합 A(f) 는 B_n 의 부분군이고 [Kbar : H] <= 2 이다.
#   Kbar 는 쌍 {f, not f} 를 보존하는 g 전체, 곧 "출력 부호를 무시한 대칭군"이다.
#   K(f) 는 Kbar 에서 Z_2 로 가는 준동형 phi 의 그래프이다 (phi 의 핵 = H, phi = 0 이면 반대칭이 없다).
#   반대칭 원소는 꼭짓점을 하나도 고정하지 못하므로(f(gx) = not f(x) 이면 gx != x), Kbar != H 이면 Kbar 는 고정 꼭짓점이 없다.
#
# 고정점 수. B_n 의 부분군 L 과 준동형 phi: L -> Z_2 에 대해 그래프 부분군 K = {(g, phi(g))} 가 고정하는 함수는
#   f(gx) = f(x) xor phi(g)   (모든 g in L)
# 을 만족하는 f 이다. 꼭짓점의 L-궤도 하나에서 f 는 한 점의 값으로 정해지고, 그 점의 안정화군이 phi 의 핵 안에 있어야만 값이 이어진다.
# 그러므로
#   고정 함수 수 = 2^(L 의 꼭짓점 궤도 수)   (phi = 0 이거나, phi != 0 이고 핵의 꼭짓점 궤도 수가 L 의 두 배일 때)
#              = 0                           (그렇지 않을 때)
# (phi != 0 이면 L 의 궤도마다 핵의 궤도가 정확히 둘로 갈라져야 한다). K 가 (1, 1) 을 품으면 f = not f 가 되어 고정 함수가 없다.
#
# 두 가지 계산을 한다.
#  방법 A (n <= NA_MAX): B_n x Z_2 의 표지 표. K 의 켤레류마다 NPN 궤도 수와 함수 수를 얻는다. 완전한 표이다.
#  방법 B (모든 n, n = 6 포함): B_n 의 표지 표만 쓴다. B_n 이 쌍 {f, not f} 의 집합 Y 에 작용한다고 보고,
#      L 이 쌍 하나를 보존하려면 f 가 어떤 phi 에 대해 위 조건을 만족하면 된다. phi 를 모두 더하고 쌍마다 f 가 둘이므로
#      Y 에서 L 의 고정점 수 = 2^(oL - 1) * (1 + ng),  oL = L 의 꼭짓점 궤도 수,
#      ng = 지수 2 인 부분군 N 가운데 N 의 꼭짓점 궤도 수가 2 oL 인 것의 수.
#      표지 표 분해로 Kbar 가 L 과 켤레인 NPN 류의 수를 얻는다. (Kbar 만 알고 phi 는 모른다.)
#
# 실행: gap -A -q -o 6g scripts/npn_stabilizer_types.g
#       N_MIN, N_MAX, NA_MAX, OUTDIR 을 먼저 정의하면 바꿀 수 있다. 기본은 N_MIN = 1, N_MAX = 5, NA_MAX = 5, OUTDIR = "notes".
#       n = 6 은 N_MIN := 6; N_MAX := 6; NA_MAX := 0; 으로 방법 B 만 돌린다.
# 출력: notes/npn_stabilizer_types_n<n>.tsv       (방법 A: K 의 켤레류 가운데 나타나는 것)
#       notes/npn_projected_types_n<n>.tsv        (방법 B: Kbar 의 켤레류 가운데 나타나는 것)
#       notes/npn_projected_types_n<n>_absent.tsv (방법 B: 나타나지 않는 것)
#       켤레류 번호는 실행할 때마다 달라질 수 있으므로 번호로 군을 가리키지 말 것.

SizeScreen([4096, 24]);
if not IsBound(N_MIN) then N_MIN := 1; fi;
if not IsBound(N_MAX) then N_MAX := 5; fi;
if not IsBound(NA_MAX) then NA_MAX := 5; fi;
if not IsBound(OUTDIR) then OUTDIR := "notes"; fi;
LoadPackage("smallgrp");

Tab := "\t";
Join := function(strs, sep)
  local s, r, first;
  r := "";  first := true;
  for s in strs do
    if not first then r := Concatenation(r, sep); fi;
    r := Concatenation(r, s);
    first := false;
  od;
  return r;
end;
CompactLengths := function(lens)   # [1,1,2,2,2] -> "1^2 2^3"
  local vals, parts, v;
  vals := Set(lens);
  parts := List(vals, v -> Concatenation(String(v), "^", String(Number(lens, x -> x = v))));
  return Join(parts, " ");
end;
B2I := function(b) if b then return 1; else return 0; fi; end;
Struct := function(H)
  if Size(H) = 1 then return "1"; fi;
  return StructureDescription(H);
end;

# 지수 2 인 부분군을 모두 돌려준다 (모두 정규부분군이다).
# L/<L^2> 는 원소 위수가 2 이하인 아벨군이므로 그 극대부분군이 지수 2 인 부분군의 전부이다.
IndexTwoSubgroups := function(L)
  local gens, sq, N, hom, V;
  if Size(L) = 1 then return []; fi;
  gens := GeneratorsOfGroup(L);
  sq := Concatenation(List(gens, g -> g^2), Concatenation(List(gens, g -> List(gens, h -> Comm(g, h)))));
  N := NormalClosure(L, Subgroup(L, sq));
  if Index(L, N) = 1 then return []; fi;
  hom := NaturalHomomorphismByNormalSubgroup(L, N);
  V := Image(hom);
  return List(MaximalSubgroups(V), M -> PreImage(hom, M));
end;

NPNRun := function(n)
  local pts, bit, negs, perms, s, G, gsize, dirs, phi, Sn, SnClasses, SnPairOrbs, keyList, ties, p, q,
        SnTieSig, baseLabs, labs, labs2, SnLabel, tom, m, fixO, aO, fixY, ngood, decY, j, L, oL, ng, subs,
        classKey, keyOf, ukeys, buckets, ClassIdx, tr, vorb, outA, outB, outBabs, i, totalfun, totalcls, nocc,
        swap, GA, tomA, mA, fixA, K, Ho, Kb, Hb, swapin, oK, oH, decA, ord, hord, kbord, hbord, d, dK, varorb,
        twisted, fixv, vorblensK, antisym, FindFunction, S, repTT, Hbar, Kbar, idxK, idxH, sumB, sumR9, bad, tmp,
        nclsocc, nclsB, nanti, nantiB, hbordK;

  pts := [0 .. 2^n - 1];
  bit := function(x, i) return QuoInt(x, 2^(i-1)) mod 2; end;
  negs := List([1..n], i -> PermList(List(pts, x -> (x + (1 - 2*bit(x,i)) * 2^(i-1)) + 1)));
  perms := [];
  for s in GeneratorsOfGroup(SymmetricGroup(n)) do
    Add(perms, PermList(List(pts, x -> Sum([1..n], i -> bit(x, i^s) * 2^(i-1)) + 1)));
  od;
  G := Group(Concatenation(negs, perms));
  gsize := Size(G);
  dirs := List([1..n], i -> Set(List(Filtered(pts, x -> bit(x,i) = 0), x -> Set([x+1, x+2^(i-1)+1]))));
  phi := ActionHomomorphism(G, dirs, OnSetsSets);
  Sn := Image(phi);

  # S_n 의 부분군 켤레류의 이름 (stabilizer_types.g 와 같다)
  SnClasses := List(ConjugacyClassesSubgroups(Sn), Representative);
  SnPairOrbs := function(L) return Length(Orbits(L, Combinations([1..n], 2), OnSets)); end;
  keyList := List(SnClasses, L -> [Size(L), SortedList(List(Orbits(L, [1..n]), Length)), Struct(L), SnPairOrbs(L)]);
  SnTieSig := function(L)
    return String(Collected(List(AsList(L), g -> SortedList(CycleLengths(g, [1..n])))));
  end;
  baseLabs := List(keyList, k -> Concatenation(String(k[1]), " ", String(k[2]), " ", k[3]));
  labs := [];
  for p in [1..Length(keyList)] do
    labs[p] := baseLabs[p];
    if Number(baseLabs, x -> x = baseLabs[p]) > 1 then
      labs[p] := Concatenation(labs[p], " 쌍궤도 ", String(keyList[p][4]));
    fi;
  od;
  labs2 := [];
  for p in [1..Length(keyList)] do
    labs2[p] := labs[p];
    if Number(labs, x -> x = labs[p]) > 1 then
      labs2[p] := Concatenation(labs[p], " ", SnTieSig(SnClasses[p]));
    fi;
  od;
  if Length(Set(labs2)) < Length(labs2) then
    Print("오류: S_", n, " 의 켤레류 이름이 아직 겹친다\n");
  fi;
  SnLabel := function(L)
    local k;
    k := First([1..Length(SnClasses)], k -> Size(SnClasses[k]) = Size(L) and IsConjugate(Sn, SnClasses[k], L));
    return labs2[k];
  end;

  tr := function(h) return (1^h) - 1; end;   # h 가 꼭짓점 0 을 보내는 곳 = h 의 평행이동 성분

  # ---------------------------------------------------------------- B_n 의 표지 표
  tom := TableOfMarks(G);
  m := Length(MarksTom(tom));
  fixO := List([1..m], j -> 2^Length(Orbits(RepresentativeTom(tom, j), [1..2^n])));
  aO := DecomposedFixedPointVector(tom, fixO);       # H(f) 가 켤레류 j 인 B_n-궤도 수 (stabilizer_types.g 의 a[j])

  # B_n 의 부분군을 켤레류 번호로 바꾸는 함수
  keyOf := function(U)
    return [Size(U), SortedList(List(Orbits(U, [1..2^n]), Length)), Size(Image(phi, U))];
  end;
  classKey := List([1..m], j -> keyOf(RepresentativeTom(tom, j)));
  ukeys := Set(classKey);
  buckets := List(ukeys, k -> Filtered([1..m], j -> classKey[j] = k));
  ClassIdx := function(U)
    local key, pos, k;
    key := keyOf(U);
    pos := PositionSorted(ukeys, key);
    if pos > Length(ukeys) or ukeys[pos] <> key then return fail; fi;
    for k in buckets[pos] do
      if IsConjugate(G, U, RepresentativeTom(tom, k)) then return k; fi;
    od;
    return fail;
  end;

  # ---------------------------------------------------------------- 방법 B: 쌍 {f, not f} 위의 B_n
  fixY := [];  ngood := [];
  for j in [1..m] do
    L := RepresentativeTom(tom, j);
    oL := Length(Orbits(L, [1..2^n]));
    subs := IndexTwoSubgroups(L);
    ng := Number(subs, H -> Length(Orbits(H, [1..2^n])) = 2 * oL);
    Add(ngood, ng);
    Add(fixY, 2^(oL - 1) * (1 + ng));
  od;
  decY := DecomposedFixedPointVector(tom, fixY);
  decY := Concatenation(decY, ListWithIdenticalEntries(m - Length(decY), 0));   # 분해 함수는 끝의 0 을 잘라 낸다
  if ForAny(decY, x -> not IsInt(x) or x < 0) then Print("   !!! 방법 B 의 분해에 음수나 분수가 나왔다\n"); fi;

  outB := Concatenation(OUTDIR, "/npn_projected_types_n", String(n), ".tsv");
  outBabs := Concatenation(OUTDIR, "/npn_projected_types_n", String(n), "_absent.tsv");
  PrintTo(outB, "# n = ", n, ", |B_n| = ", gsize, ", B_n 의 켤레류 ", m, "개 가운데 Kbar(f) = {g : f o g = f 또는 not f} 로 나타나는 것\n");
  AppendTo(outB, "# 열: 위수 |Kbar|, 지수(쌍 하나의 B_n-궤도 크기), 함수 수 (쌍마다 함수 둘), NPN 류 수, Kbar 를 H(f) 로 갖는 B_n-궤도 수(참고: stabilizer_types.g 의 값),\n");
  AppendTo(outB, "#     ng (지수 2 인 부분군 가운데 반대칭이 가능한 것의 수; 0 이면 이 Kbar 의 NPN 류는 모두 반대칭이 없다),\n");
  AppendTo(outB, "#     Kbar 의 변수 위 궤도 길이, 평행이동 부분의 차원 d, 꼭짓점 위 궤도 수, 고정 꼭짓점 수, 구조\n");
  PrintTo(outBabs, "# n = ", n, ", Kbar 로 나타나지 않는 B_n 의 켤레류\n");
  AppendTo(outBabs, "# 열: 위수, 변수 위 궤도 길이, 평행이동 부분의 차원 d, 꼭짓점 위 궤도 수, 고정 꼭짓점 수, ng, H(f) 로는 나타나는 B_n-궤도 수\n");

  totalfun := 0;  totalcls := 0;  nocc := 0;
  for j in [1..m] do
    L := RepresentativeTom(tom, j);
    ord := Size(L);
    vorb := Orbits(L, [1..2^n]);
    hbord := Size(Image(phi, L));
    d := LogInt(ord / hbord, 2);
    varorb := SortedList(List(Orbits(Image(phi, L), [1..n]), Length));
    fixv := Number(vorb, o -> Length(o) = 1);
    if decY[j] <> 0 then
      nocc := nocc + 1;
      totalcls := totalcls + decY[j];
      totalfun := totalfun + 2 * decY[j] * gsize / ord;
      AppendTo(outB, ord, Tab, gsize / ord, Tab, 2 * decY[j] * gsize / ord, Tab, decY[j], Tab, aO[j], Tab, ngood[j], Tab,
               varorb, Tab, d, Tab, Length(vorb), Tab, fixv, Tab, Struct(L), "\n");
    else
      AppendTo(outBabs, ord, Tab, varorb, Tab, d, Tab, Length(vorb), Tab, fixv, Tab, ngood[j], Tab, aO[j], "\n");
    fi;
  od;
  Print("n = ", n, " [방법 B] |B_n| = ", gsize, ", 켤레류 ", m, "개, Kbar 로 나타나는 것 ", nocc, "개\n");
  Print("   NPN 류 수 ", totalcls, ", 함수 수 합 ", totalfun, " (2^(2^n) = ", 2^(2^n), ")\n");
  if totalfun <> 2^(2^n) then Print("   !!! 함수 수의 합이 맞지 않는다\n"); fi;
  nclsB := totalcls;
  # NPN 류 수 = (B_n-궤도 수 + 자기 여집합 궤도 수) / 2, 자기 여집합 궤도 수 = 2 (NPN 류 수) - (B_n-궤도 수)
  Print("   B_n-궤도 수 ", Sum(aO), ", 자기 여집합 B_n-궤도 수 = 2 NPN - NP = ", 2 * totalcls - Sum(aO), "\n");

  # ---------------------------------------------------------------- 방법 A: B_n x Z_2 의 표지 표
  if n <= NA_MAX then
    swap := MappingPermListList([2^n + 1, 2^n + 2], [2^n + 2, 2^n + 1]);
    GA := Group(Concatenation(GeneratorsOfGroup(G), [swap]));
    tomA := TableOfMarks(GA);
    mA := Length(MarksTom(tomA));
    fixA := [];
    for i in [1..mA] do
      K := RepresentativeTom(tomA, i);
      if swap in K then
        Add(fixA, 0);
      else
        Ho := Stabilizer(K, 2^n + 1);          # e = 0 인 원소들 (K 가 점 2^n+1 을 움직이면 e = 1)
        Kb := Subgroup(G, List(GeneratorsOfGroup(K), g -> RestrictedPerm(g, [1..2^n])));
        Hb := Subgroup(G, List(GeneratorsOfGroup(Ho), g -> RestrictedPerm(g, [1..2^n])));
        oK := Length(Orbits(Kb, [1..2^n]));
        if Size(Ho) = Size(K) then
          Add(fixA, 2^oK);
        elif Length(Orbits(Hb, [1..2^n])) = 2 * oK then
          Add(fixA, 2^oK);
        else
          Add(fixA, 0);
        fi;
      fi;
    od;
    decA := DecomposedFixedPointVector(tomA, fixA);
    decA := Concatenation(decA, ListWithIdenticalEntries(mA - Length(decA), 0));
    if ForAny(decA, x -> not IsInt(x) or x < 0) then Print("   !!! 방법 A 의 분해에 음수나 분수가 나왔다\n"); fi;

    FindFunction := function(K, Kb, Hb, antisym)
      # K 가 정확히 NPN 안정화군인 함수를 찾는다. 점 집합 S = { x : f(x) = 1 } 로 돌려준다.
      local orbs, kk, t, S, mask, act, tryS, gensD, Kd, orbsD, pairs, choice, o, key, i, x;
      act := function(S, el)
        local T;
        T := OnSets(S, el);
        if (2^n + 1)^el <> 2^n + 1 then return Difference([1..2^n], T); fi;
        return T;
      end;
      tryS := function(S)
        return Size(Stabilizer(GA, S, act)) = Size(K);
      end;
      if not antisym then
        orbs := Orbits(Kb, [1..2^n]);
        kk := Length(orbs);
        for t in [1..3000] do
          S := Union(Filtered(orbs, o -> Random([0,1]) = 1));
          if tryS(S) then return S; fi;
        od;
        if kk <= 16 then
          for mask in [0 .. 2^kk - 1] do
            S := Union(List(Filtered([1..kk], i -> QuoInt(mask, 2^(i-1)) mod 2 = 1), i -> orbs[i]));
            if tryS(S) then return S; fi;
          od;
        fi;
        return fail;
      fi;
      # 반대칭이 있는 경우: K 가 (정점, 비트) 쌍 위에 (x, b) -> (x^g, b xor e(g)) 로 작용한다고 보고 한 쪽 궤도를 고른다.
      gensD := List(GeneratorsOfGroup(K), function(g)
                  local e;
                  e := B2I((2^n + 1)^g <> 2^n + 1);
                  return PermList(Concatenation(
                    List([1..2^n], x -> x^g + 2^n * ((0 + e) mod 2)),
                    List([1..2^n], x -> x^g + 2^n * ((1 + e) mod 2))));
                end);
      Kd := Group(gensD, ());
      orbsD := Orbits(Kd, [1 .. 2^(n+1)]);
      # 궤도마다 (정점 쪽 투영의 최솟값) 을 키로 묶으면 둘씩 짝이 된다
      pairs := [];
      for o in orbsD do
        key := Minimum(List(o, x -> ((x - 1) mod 2^n) + 1));
        if not IsBound(pairs[key]) then pairs[key] := []; fi;
        Add(pairs[key], o);
      od;
      pairs := Compacted(pairs);
      if ForAny(pairs, x -> Length(x) <> 2) then return fail; fi;
      kk := Length(pairs);
      for t in [1..3000] do
        S := [];
        for o in pairs do
          choice := o[Random([1, 2])];
          UniteSet(S, Set(List(Filtered(choice, x -> x > 2^n), x -> x - 2^n)));
        od;
        S := Set(S);
        if tryS(S) then return S; fi;
      od;
      if kk <= 16 then
        for mask in [0 .. 2^kk - 1] do
          S := [];
          for i in [1..kk] do
            choice := pairs[i][1 + (QuoInt(mask, 2^(i-1)) mod 2)];
            UniteSet(S, Set(List(Filtered(choice, x -> x > 2^n), x -> x - 2^n)));
          od;
          if tryS(Set(S)) then return Set(S); fi;
        od;
      fi;
      return fail;
    end;

    outA := Concatenation(OUTDIR, "/npn_stabilizer_types_n", String(n), ".tsv");
    PrintTo(outA, "# n = ", n, ", |B_n x Z_2| = ", 2 * gsize, ", B_n x Z_2 의 켤레류 ", mA, "개 가운데 NPN 안정화군 K(f) 로 나타나는 것\n");
    AppendTo(outA, "# 열: 1 위수 |K|, 2 지수(함수 하나의 NPN 궤도 크기), 3 함수 수, 4 NPN 류 수, 5 반대칭 유무(1 이면 K 가 e = 1 인 원소를 가짐, |Kbar| = 2 |H|),\n");
    AppendTo(outA, "#     6 |Kbar|, 7 |H|, 8 H 의 순열 부분의 위수, 9 Kbar 의 순열 부분의 위수, 10 H 의 순열 부분의 변수 위 궤도 길이, 11 Kbar 의 순열 부분의 변수 위 궤도 길이,\n");
    AppendTo(outA, "#     12 H 의 순열 부분의 S_n 켤레류, 13 H 의 평행이동 부분의 차원 d (|W| = 2^d), 14 Kbar 의 평행이동 부분의 차원, 15 H 의 꼭짓점 위 궤도 수, 16 Kbar 의 꼭짓점 위 궤도 수,\n");
    AppendTo(outA, "#     17 H 의 고정 꼭짓점 수, 18 H 의 유형(U 고정 꼭짓점 있음, T0 꼬임, T1 평행이동 있음), 19 구조(K), 20 구조(Kbar),\n");
    AppendTo(outA, "#     21 B_n 켤레류 번호(Kbar 의 것과 H 의 것; 실행마다 달라질 수 있음), 22 함수 하나(진리표 16진수)\n");
    totalfun := 0;  totalcls := 0;  nocc := 0;  nanti := 0;
    sumB := List([1..m], j -> 0);
    sumR9 := List([1..m], j -> 0);
    for i in [1..mA] do
      if decA[i] = 0 then continue; fi;
      K := RepresentativeTom(tomA, i);
      ord := Size(K);
      Ho := Stabilizer(K, 2^n + 1);
      Kb := Subgroup(G, List(GeneratorsOfGroup(K), g -> RestrictedPerm(g, [1..2^n])));
      Hb := Subgroup(G, List(GeneratorsOfGroup(Ho), g -> RestrictedPerm(g, [1..2^n])));
      antisym := (Size(Ho) <> ord);
      nocc := nocc + 1;
      totalcls := totalcls + decA[i];
      totalfun := totalfun + decA[i] * 2 * gsize / ord;
      if antisym then nanti := nanti + decA[i]; fi;
      idxK := ClassIdx(Kb);
      idxH := ClassIdx(Hb);
      sumB[idxK] := sumB[idxK] + decA[i];
      sumR9[idxH] := sumR9[idxH] + decA[i] * 2 * gsize / ord;
      vorb := Orbits(Hb, [1..2^n]);
      fixv := Number(vorb, o -> Length(o) = 1);
      hbord := Size(Image(phi, Hb));
      d := LogInt(Size(Hb) / hbord, 2);
      if fixv > 0 then twisted := "U"; elif d = 0 then twisted := "T0"; else twisted := "T1"; fi;
      S := FindFunction(K, Kb, Hb, antisym);
      if S = fail then repTT := "-"; else repTT := LowercaseString(HexStringInt(Sum(S, p -> 2^(p-1)))); fi;
      hbordK := Size(Image(phi, Kb));
      AppendTo(outA, ord, Tab, 2 * gsize / ord, Tab, decA[i] * 2 * gsize / ord, Tab, decA[i], Tab, B2I(antisym), Tab, Size(Kb), Tab, Size(Hb), Tab,
               hbord, Tab, hbordK, Tab, SortedList(List(Orbits(Image(phi, Hb), [1..n]), Length)), Tab,
               SortedList(List(Orbits(Image(phi, Kb), [1..n]), Length)), Tab, SnLabel(Image(phi, Hb)), Tab,
               d, Tab, LogInt(Size(Kb) / hbordK, 2), Tab, Length(vorb), Tab, Length(Orbits(Kb, [1..2^n])), Tab, fixv, Tab, twisted, Tab,
               Struct(K), Tab, Struct(Kb), Tab, idxK, " ", idxH, Tab, repTT, "\n");
    od;
    Print("n = ", n, " [방법 A] |B_n x Z_2| = ", 2 * gsize, ", 켤레류 ", mA, "개, 나타나는 것 ", nocc, "개\n");
    Print("   NPN 류 수 ", totalcls, " (반대칭이 있는 류 ", nanti, "), 함수 수 합 ", totalfun, " (2^(2^n) = ", 2^(2^n), ")\n");
    if totalfun <> 2^(2^n) then Print("   !!! 함수 수의 합이 맞지 않는다\n"); fi;
    # 방법 B 와의 대조: Kbar 가 켤레류 j 인 NPN 류 수
    bad := 0;
    for j in [1..m] do
      if sumB[j] <> decY[j] then bad := bad + 1; fi;
    od;
    Print("   방법 B 와 대조 (Kbar 별 NPN 류 수): 다른 켤레류 ", bad, "개\n");
    # stabilizer_types.g 의 표(H 별 함수 수)와 대조: H(f) 가 켤레류 j 인 함수 수 = aO[j] |B_n| / |L_j|
    bad := 0;
    for j in [1..m] do
      if sumR9[j] <> aO[j] * gsize / Size(RepresentativeTom(tom, j)) then bad := bad + 1; fi;
    od;
    Print("   H 별 함수 수를 B_n 표(aO)와 대조: 다른 켤레류 ", bad, "개\n");
    Print("   자기 여집합 NPN 류(반대칭 있는 류) ", nanti, ", 2 NPN - NP = ", 2 * totalcls - Sum(aO), "\n");
  fi;
end;

for n in [N_MIN .. N_MAX] do
  NPNRun(n);
od;
QUIT;
