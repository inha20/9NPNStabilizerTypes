# stabilizer_types.g
# 부울 함수의 안정화군(입력 반전과 변수 순열이 이루는 B_n 안에서)의 종류를 모두 적는다.
#
#   H(f) = { g in B_n : f(g x) = f(x) 모든 x }.
#   B_n 의 부분군 켤레류 하나하나에 대해:
#     - 어떤 함수의 H(f) 로 나타나는가, 나타난다면 그 함수는 몇 개이고 B_n-궤도는 몇 개인가
#       (표지 표의 분해: 켤레류 j 의 함수 수 = a[j] * |B_n| / |H_j|),
#     - 군의 구조: 위수, 변수 위의 궤도 길이(순열 부분 Hbar 의), 평행이동 부분 W = H 와 Z_2^n 의 교집합의 차원,
#       꼭짓점 2^n 개 위의 궤도 수와 고정 꼭짓점 수(고정 꼭짓점이 없으면 "꼬였다"),
#     - 나타나지 않는다면 왜 그런가: 궤도 닫힘(H 의 모든 꼭짓점 궤도를 각각 보존하는 군 Hhat 가 H 와 같은가).
#       닫히지 않았으면 H 는 Hhat 보다 작아서 어떤 함수도 H 만 안정화할 수 없다.
#       닫혔는데도 나타나지 않으면 부울 함수가 값을 둘밖에 갖지 못해서 생기는 장애이다.
#
# 군 H 가 B_n 안에서 고정 꼭짓점을 가지면 H 는 S_n 의 부분군과 B_n 안에서 켤레이다(고정점을 0 으로 옮기면 된다).
# 꼬인 군(고정 꼭짓점이 없음)은 순열만으로는 나타낼 수 없는 대칭이다.
#
# 실행: gap -A -q scripts/stabilizer_types.g          (n = 1..5)
#       N_MAX := 6; 을 먼저 정의하면 6 까지 (표지 표 계산이 수 분 걸린다)
# 출력: notes/stabilizer_types_n<n>.tsv (나타나는 켤레류), notes/stabilizer_types_n<n>_absent.tsv (나타나지 않는 켤레류)
#       탭으로 구분한 표이다. 켤레류 번호는 실행할 때마다 달라질 수 있으므로 번호로 군을 가리키지 말 것.

SizeScreen([4096, 24]);   # 긴 줄을 중간에 끊지 않게 한다
if not IsBound(N_MIN) then N_MIN := 1; fi;
if not IsBound(N_MAX) then N_MAX := 5; fi;
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

Run := function(n)
  local pts, bit, negs, perms, s, G, gsize, dirs, phi, Sn, tom, m, fix, a, j, H, Hb, ord, hbord, dW,
        vorb, fixv, twisted, Khat, closed, varorb, i, vorblens, Hbi, lift, tr, sig, orbs, S,
        nocc, nclosed, nclosedabsent, totalfun, totalorb, SnLabel, SnClasses, SnPairOrbs,
        twlist, repTT, FindFunction, outocc, outabs, ties, keyList, p, q,
        SnTieSig, baseLabs, labs, labs2;

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

  # S_n 의 부분군 켤레류에 결정적인 이름을 붙인다: (위수, 변수 위 궤도 길이, 구조, 2-부분집합 위 궤도 수)
  SnClasses := List(ConjugacyClassesSubgroups(Sn), Representative);
  SnPairOrbs := function(L) return Length(Orbits(L, Combinations([1..n], 2), OnSets)); end;
  keyList := List(SnClasses, L -> [Size(L), SortedList(List(Orbits(L, [1..n]), Length)), Struct(L), SnPairOrbs(L)]);
  ties := 0;
  for p in [1..Length(keyList)] do
    for q in [p+1..Length(keyList)] do
      if keyList[p] = keyList[q] then ties := ties + 1; fi;
    od;
  od;
    # 이름: "위수 변수 위 궤도 길이 구조". 겹치면 2-부분집합 위 궤도 수를, 그래도 겹치면 원소들의 순환형 분포를 덧붙인다
  # (S_6 에는 바깥 자기동형으로 서로 옮겨지는 켤레류의 짝이 있어 앞의 세 값만으로는 구별되지 않는다).
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

  tom := TableOfMarks(G);
  m := Length(MarksTom(tom));
  fix := List([1..m], j -> 2^Length(Orbits(RepresentativeTom(tom, j), [1..2^n])));
  a := DecomposedFixedPointVector(tom, fix);

  FindFunction := function(H)
    local orbs, t, S, mask, kk;
    orbs := Orbits(H, [1..2^n]);
    kk := Length(orbs);
    for t in [1..3000] do
      S := Union(Filtered(orbs, o -> Random([0,1]) = 1));
      if Size(Stabilizer(G, S, OnSets)) = Size(H) then return S; fi;
    od;
    if kk <= 16 then
      for mask in [0 .. 2^kk - 1] do
        S := Union(List(Filtered([1..kk], i -> QuoInt(mask, 2^(i-1)) mod 2 = 1), i -> orbs[i]));
        if Size(Stabilizer(G, S, OnSets)) = Size(H) then return S; fi;
      od;
    fi;
    return fail;
  end;

  tr := function(h) return (1^h) - 1; end;   # h 가 꼭짓점 0 을 보내는 곳 = h 의 평행이동 성분

  outocc := Concatenation(OUTDIR, "/stabilizer_types_n", String(n), ".tsv");
  outabs := Concatenation(OUTDIR, "/stabilizer_types_n", String(n), "_absent.tsv");
  PrintTo(outocc, "# n = ", n, ", |B_n| = ", gsize, ", B_n 의 켤레류 ", m, "개 가운데 안정화군으로 나타나는 것\n");
  AppendTo(outocc, "# 열: 위수, 지수(함수 하나의 B_n-궤도 크기), 함수 수, B_n-궤도 수, Hbar 의 위수, Hbar 의 변수 위 궤도 길이, Hbar 의 S_n 켤레류,\n");
  AppendTo(outocc, "#     평행이동 부분의 차원 d (|W| = 2^d), 꼭짓점 위 궤도 수, 고정 꼭짓점 수, 유형(U 고정 꼭짓점 있음, T0 꼬임, T1 평행이동 있음),\n");
  AppendTo(outocc, "#     닫힘, 구조, 꼬임 자료(T0 일 때 변수 궤도마다 (궤도 길이, 점 안정화군 위수, 반전 여부)), 꼭짓점 궤도 길이, 함수 하나(진리표 16진수)\n");
  PrintTo(outabs, "# n = ", n, ", 안정화군으로 나타나지 않는 B_n 의 켤레류\n");
  AppendTo(outabs, "# 열: 위수, Hbar 의 위수, Hbar 의 변수 위 궤도 길이, 평행이동 부분의 차원 d, 꼭짓점 위 궤도 수, 고정 꼭짓점 수, 닫힘(1 이면 H = Hhat), Hhat 의 위수\n");

  nocc := 0;  nclosed := 0;  nclosedabsent := 0;  totalfun := 0;  totalorb := 0;
  for j in [1..m] do
    H := RepresentativeTom(tom, j);
    ord := Size(H);
    Hb := Image(phi, H);
    hbord := Size(Hb);
    dW := LogInt(ord / hbord, 2);
    varorb := SortedList(List(Orbits(Hb, [1..n]), Length));
    vorb := Orbits(H, [1..2^n]);
    vorblens := SortedList(List(vorb, Length));
    fixv := Number(vorb, o -> Length(o) = 1);
    Khat := Stabilizer(G, List(vorb, o -> Set(o)), OnTuplesSets);
    closed := (Size(Khat) = ord);
    if closed then nclosed := nclosed + 1; fi;
    if a[j] <> 0 then
      nocc := nocc + 1;
      totalfun := totalfun + a[j] * gsize / ord;
      totalorb := totalorb + a[j];
      if fixv > 0 then twisted := "U"; elif dW = 0 then twisted := "T0"; else twisted := "T1"; fi;
      # 꼬임 자료
      twlist := [];
      if twisted = "T0" then
        lift := RestrictedMapping(phi, H);
        for orbs in Orbits(Hb, [1..n]) do
          i := orbs[1];
          Hbi := Stabilizer(Hb, i);
          sig := false;
          for s in GeneratorsOfGroup(Hbi) do
            if bit(tr(PreImagesRepresentative(lift, s)), i) = 1 then sig := true; fi;
          od;
          Add(twlist, Concatenation("(", String(Length(orbs)), ",", String(Size(Hbi)), ",", String(B2I(sig)), ")"));
        od;
      fi;
      S := FindFunction(H);
      if S = fail then repTT := "-"; else repTT := LowercaseString(HexStringInt(Sum(S, p -> 2^(p-1)))); fi;
      AppendTo(outocc, ord, Tab, gsize / ord, Tab, a[j] * gsize / ord, Tab, a[j], Tab, hbord, Tab, varorb, Tab, SnLabel(Hb), Tab,
               dW, Tab, Length(vorb), Tab, fixv, Tab, twisted, Tab, B2I(closed), Tab, Struct(H), Tab,
               Join(twlist, " "), Tab, CompactLengths(vorblens), Tab, repTT, "\n");
    else
      if closed then nclosedabsent := nclosedabsent + 1; fi;
      AppendTo(outabs, ord, Tab, hbord, Tab, varorb, Tab, dW, Tab, Length(vorb), Tab, fixv, Tab, B2I(closed), Tab, Size(Khat), "\n");
    fi;
  od;
  Print("n = ", n, ": |B_n| = ", gsize, ", 켤레류 ", m, "개, 나타나는 것 ", nocc, "개, 안 나타나는 것 ", m - nocc, "개\n");
  Print("   닫힌 켤레류(H = Hhat) ", nclosed, "개 가운데 나타나지 않는 것 ", nclosedabsent, "개\n");
  Print("   함수 수 합 ", totalfun, " (2^(2^n) = ", 2^(2^n), "), B_n-궤도 수 합 ", totalorb, "\n");
  if totalfun <> 2^(2^n) then Print("   !!! 함수 수의 합이 맞지 않는다\n"); fi;
end;

for n in [N_MIN .. N_MAX] do
  Run(n);
od;
QUIT;
