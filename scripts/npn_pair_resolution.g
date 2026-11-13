# npn_pair_resolution.g
# 쌍 (Kbar, H) 의 분해와 그 세분에 필요한 자료를 GAP 으로 내보낸다. 풀이는 scripts/solve_npn_pair_resolution.py 가 한다.
#
# 배경. 함수 f 의 NPN 안정화군 K(f) 는 B_n 안의 군 Kbar(f) 에서 Z_2 로 가는 준동형 phi_f 의 그래프이다 (핵 = H(f)).
#   B_n 의 켤레류 j 마다 적는 값은 다음과 같다.
#     aO[j]    H(f) 가 j 와 켤레인 B_n-궤도(NP 류)의 수
#     decY[j]  Kbar(f) 가 j 와 켤레인 NPN 류의 수 (쌍 {f, not f} 위의 B_n 표지 표 분해, 고정점 수 2^(o-1) (1 + ng))
#   B_n 에서 Z_2 로 가는 준동형은 0 과 psi1 = (변수 순열의 부호), psi2 = (반전한 변수 수의 홀짝), psi3 = psi1 + psi2 이다.
#   psi 를 하나 고르면 B_n 이 함수 집합에 g * f = (f o g^-1) xor psi(g) 로 작용한다 ("꼬인 작용"). 이 작용의 안정화군은
#   {g in Kbar(f) : phi_f(g) = psi(g)} 이고, 고정 함수 수는 psi 가 L 의 꼭짓점 안정화군에서 모두 0 이면 2^o(L), 아니면 0 이다.
#   그래서 같은 표지 표로 꼬인 작용의 궤도 분해 aOpsi[j] (안정화군이 j 인 궤도의 수)도 얻는다.
#   NPN 류 하나가 (Kbar, phi) 일 때 꼬인 궤도로 기여하는 것은
#       phi = psi|Kbar 이면 안정화군 Kbar 인 궤도 둘, 아니면 안정화군 ker(phi + psi|Kbar) 인 궤도 하나
#   이다 (psi = 0 이면 원래의 aO 이다). phi 마다 이 정보 (서명)를 적으면 풀이 프로그램이 연립방정식을 세운다.
#
# 실행: N_MIN := 1; N_MAX := 5; NA_MAX := 5; OUTDIR := "..."; Read("scripts/npn_pair_resolution.g");
#       n = 6 은 N_MIN := 6; N_MAX := 6; NA_MAX := 0; 으로 돌린다 (`gap -A -q -o 6g 파일.g < /dev/null`).
#       n = 6 의 표지 표는 7분쯤 걸린다. 최상위에서 표지 표를 만들어 작업 공간으로 저장해 두고 이어 돌릴 수도 있다 (SaveWorkspace 는
#       최상위에서만 되므로 표준 입력으로 읽힌다):
#         NO_MAIN := true; Read("scripts/npn_pair_resolution.g");
#         PRE := MakeBn(6); PRE.tom := TableOfMarks(PRE.G); SaveWorkspace("tom6.ws"); QUIT;
#       그 뒤 `gap -L tom6.ws -A -q 파일.g < /dev/null` 의 파일에서 Unbind(NO_MAIN); N_MIN := 6; N_MAX := 6; ... 로 읽으면 PRE 를 쓴다.
# 켤레류 번호 찾기: 부분군 U 가 켤레류 j 의 대표원 L_j 안에 있으면 U 의 류는 SubsTom(tom)[j] (L_j 의 부분군 류들) 안에 있으므로
#   그 가운데 불변량(위수, 꼭짓점 궤도 길이, S_n 속 상의 크기)이 같은 류만 본다. 후보가 하나뿐이면 켤레 검사 없이 정해진다.
# 출력: <OUTDIR>/npn_pair_data_n<n>.tsv   (켤레류 번호는 실행마다 달라질 수 있으므로 같은 파일 안에서만 쓸 것)
#   C 줄: j  위수  꼭짓점 궤도 수  aO  decY  ng  validN(쉼표로 나눈 켤레류 번호, 없으면 -)  aOpsi1  aOpsi2  aOpsi3
#   P 줄: j  kerClass  t1  t2  t3        (decY[j] > 0 인 j 의 유효한 phi 마다 한 줄. phi = 0 이면 kerClass = 0.
#                                         t_k = -1 이면 phi = psi_k|Kbar, 아니면 ker(phi + psi_k|Kbar) 의 켤레류 번호)
#   Q 줄: j  유효한 phi != 0 의 N_B(L)-궤도 수  (decY[j] > 0 인 j)
#   D 줄: j  위수  변수 위 궤도 길이  평행이동 부분의 차원 d  꼭짓점 궤도 수  고정 꼭짓점 수  구조  (aO > 0 또는 decY > 0 인 j)
#   A 줄: idxK  idxH  반대칭(0/1)  NPN 류 수  kerClass  t1  t2  t3     (방법 A 의 참값. n <= NA_MAX 일 때만)

SizeScreen([4096, 24]);
if not IsBound(N_MIN) then N_MIN := 1; fi;
if not IsBound(N_MAX) then N_MAX := 5; fi;
if not IsBound(NA_MAX) then NA_MAX := 5; fi;
if not IsBound(OUTDIR) then OUTDIR := "notes"; fi;

Tab := "\t";
B2I := function(b) if b then return 1; else return 0; fi; end;

# 지수 2 인 부분군을 모두 돌려준다 (npn_stabilizer_types.g 와 같다).
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

# L 의 원소에 대한 술어 p 가 L 에서 Z_2 로 가는 준동형의 핵 {g : p(g)} 을 정할 때 그 핵을 돌려준다 (슈라이어 생성원).
KernelOfPredicate := function(L, p)
  local gens, g0, s, sg, K;
  gens := GeneratorsOfGroup(L);
  g0 := First(gens, s -> not p(s));
  if g0 = fail then return L; fi;
  sg := [];
  for s in gens do
    if p(s) then
      Add(sg, s);
      Add(sg, g0 * s * g0^-1);
    else
      Add(sg, s * g0^-1);
      Add(sg, g0 * s);
    fi;
  od;
  K := Subgroup(L, sg);
  if Index(L, K) <> 2 then Error("핵의 지수가 2 가 아니다"); fi;
  return K;
end;

# B_n = Z_2^n x| S_n 을 2^n 개의 꼭짓점 위의 순열군으로 만들고, 방향 위의 작용 phi 와 준동형 psi_1..3 을 돌려준다.
MakeBn := function(n)
  local pts, bit, negs, perms, s, G, dirs, phi, psis, i, k, K, Hb;
  pts := [0 .. 2^n - 1];
  bit := function(x, i) return QuoInt(x, 2^(i-1)) mod 2; end;
  negs := List([1..n], i -> PermList(List(pts, x -> (x + (1 - 2*bit(x,i)) * 2^(i-1)) + 1)));
  perms := [];
  for s in GeneratorsOfGroup(SymmetricGroup(n)) do
    Add(perms, PermList(List(pts, x -> Sum([1..n], i -> bit(x, i^s) * 2^(i-1)) + 1)));
  od;
  G := Group(Concatenation(negs, perms));
  dirs := List([1..n], i -> Set(List(Filtered(pts, x -> bit(x,i) = 0), x -> Set([x+1, x+2^(i-1)+1]))));
  phi := ActionHomomorphism(G, dirs, OnSetsSets);

  # B_n 에서 Z_2 로 가는 0 이 아닌 준동형 세 개 (값은 0 또는 1)
  psis := [
    function(g) if SignPerm(Image(phi, g)) = 1 then return 0; else return 1; fi; end,
    function(g) return Sum(CoefficientsQadic((1^g) - 1, 2)) mod 2; end,
    function(g)
      return (Sum(CoefficientsQadic((1^g) - 1, 2)) + B2I(SignPerm(Image(phi, g)) <> 1)) mod 2;
    end ];
  # 준동형인지 임의의 원소로 점검한다
  for i in [1 .. 200] do
    K := PseudoRandom(G);  Hb := PseudoRandom(G);
    for k in [1..3] do
      if (psis[k](K * Hb) - psis[k](K) - psis[k](Hb)) mod 2 <> 0 then Error("psi 가 준동형이 아니다"); fi;
    od;
  od;

  return rec(n := n, G := G, phi := phi, psis := psis);
end;

PairRun := function(n)
  local R, G, gsize, phi, tom, m, fixO, aO, keyOf, classKey, ukeys, buckets, ClassIdx, subsOf, ClassIdxIn,
        fixY, ngood, decY, validIdx, validSubs, j, L, oL, subs, vN, out, swap, GA, tomA, mA, fixA, i, K, Ho, Kb, Hb, oK,
        decA, antisym, t0, nA, nAanti, cnt, str, psis, PsiValid, fixP, decP, k, Sig, N, sigs, Kps, jK, jH,
        hbord, varorb, vorb, fixv, StructStr, Nrm, orbsN;

  t0 := Runtime();
  if IsBound(PRE) and PRE.n = n and IsBound(PRE.tom) then
    R := PRE;
    tom := R.tom;
    Print("n = ", n, ": 저장해 둔 표지 표를 쓴다\n");
  else
    R := MakeBn(n);
    tom := TableOfMarks(R.G);
  fi;
  G := UnderlyingGroup(tom);
  gsize := Size(G);
  phi := R.phi;
  psis := R.psis;
  m := Length(MarksTom(tom));
  Print("n = ", n, ": B_n 의 표지 표 ", m, "개 켤레류, ", Int((Runtime() - t0) / 1000), " 초\n");
  fixO := List([1..m], j -> 2^Length(Orbits(RepresentativeTom(tom, j), [1..2^n])));
  aO := DecomposedFixedPointVector(tom, fixO);
  aO := Concatenation(aO, ListWithIdenticalEntries(m - Length(aO), 0));

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

  # U 가 켤레류 j 의 대표원 안의 부분군일 때 U 의 켤레류 번호 (후보는 SubsTom 의 j 번째 목록)
  subsOf := SubsTom(tom);
  ClassIdxIn := function(U, j)
    local key, cands, k;
    key := keyOf(U);
    cands := Filtered(subsOf[j], k -> classKey[k] = key);
    if Length(cands) = 1 then return cands[1]; fi;
    for k in cands do
      if IsConjugate(G, U, RepresentativeTom(tom, k)) then return k; fi;
    od;
    return fail;
  end;

  # psi 가 U 의 모든 꼭짓점 안정화군에서 0 인가
  PsiValid := function(U, psi)
    local o, x;
    for o in Orbits(U, [1..2^n]) do
      x := o[1];
      if ForAny(GeneratorsOfGroup(Stabilizer(U, x)), g -> psi(g) <> 0) then return false; fi;
    od;
    return true;
  end;

  fixY := [];  ngood := [];  validIdx := [];  validSubs := [];
  for j in [1..m] do
    L := RepresentativeTom(tom, j);
    oL := Length(Orbits(L, [1..2^n]));
    subs := IndexTwoSubgroups(L);
    vN := Filtered(subs, H -> Length(Orbits(H, [1..2^n])) = 2 * oL);
    Add(ngood, Length(vN));
    Add(fixY, 2^(oL - 1) * (1 + Length(vN)));
    Add(validSubs, vN);
    Add(validIdx, Set(List(vN, N -> ClassIdxIn(N, j))));
    if fail in validIdx[j] then Print("   !!! 켤레류를 찾지 못한 부분군이 있다 (j = ", j, ")\n"); fi;
    if j mod 1000 = 0 then Print("   ... 켤레류 ", j, " / ", m, " (", Int((Runtime() - t0) / 1000), " 초)\n"); fi;
  od;
  Print("   지수 2 부분군과 켤레류 번호: ", Int((Runtime() - t0) / 1000), " 초\n");
  decY := DecomposedFixedPointVector(tom, fixY);
  decY := Concatenation(decY, ListWithIdenticalEntries(m - Length(decY), 0));
  if ForAny(decY, x -> not IsInt(x) or x < 0) then Print("   !!! 분해에 음수나 분수가 나왔다\n"); fi;

  # 꼬인 작용의 분해
  decP := [];
  for k in [1..3] do
    fixP := List([1..m], j -> function() local U; U := RepresentativeTom(tom, j);
                   if PsiValid(U, psis[k]) then return 2^Length(Orbits(U, [1..2^n])); else return 0; fi; end());
    decP[k] := DecomposedFixedPointVector(tom, fixP);
    decP[k] := Concatenation(decP[k], ListWithIdenticalEntries(m - Length(decP[k]), 0));
    if ForAny(decP[k], x -> not IsInt(x) or x < 0) then Print("   !!! 꼬인 작용(", k, ")의 분해에 음수나 분수가 나왔다\n"); fi;
  od;
  Print("   꼬인 작용의 분해: ", Int((Runtime() - t0) / 1000), " 초\n");

  StructStr := function(U)
    if Size(U) = 1 then return "1"; fi;
    return StructureDescription(U);
  end;

  # (L, N) 에서 서명을 구한다. N = L 이면 phi = 0. jL 은 L 의 켤레류 번호, Kps 는 psi_k 의 L 위 핵 세 개.
  Sig := function(L, N, jL, Kps)
    local kc, t, k, Kp, N3;
    if Size(N) = Size(L) then kc := 0; else kc := ClassIdxIn(N, jL); fi;
    t := [];
    for k in [1..3] do
      Kp := Kps[k];                                              # psi_k 의 L 위 핵 (psi_k|L = 0 이면 L 자신)
      if Size(Kp) = Size(L) then
        # psi_k|L = 0: phi = 0 이면 같다, 아니면 ker(phi + 0) = N
        if kc = 0 then Add(t, -1); else Add(t, kc); fi;
      elif Size(N) = Size(L) then
        Add(t, ClassIdxIn(Kp, jL));                              # phi = 0 이고 psi|L != 0
      elif N = Kp then
        Add(t, -1);                                              # phi = psi|L
      else
        N3 := KernelOfPredicate(L, g -> (g in N) = (g in Kp));  # phi + psi|L 의 핵
        Add(t, ClassIdxIn(N3, jL));
      fi;
    od;
    return Concatenation([kc], t);
  end;

  out := Concatenation(OUTDIR, "/npn_pair_data_n", String(n), ".tsv");
  PrintTo(out, "# n = ", n, ", |B_n| = ", gsize, ", B_n 의 켤레류 ", m, "개 (번호는 이 실행 안에서만 뜻이 있다)\n");
  AppendTo(out, "# C 줄: j, 위수, 꼭짓점 궤도 수, aO, decY, ng, validN(켤레류 번호들), aOpsi1, aOpsi2, aOpsi3\n");
  AppendTo(out, "# P 줄: j, kerClass(phi=0 이면 0), t1, t2, t3 (decY[j] > 0 인 j 의 유효한 phi 마다)\n");
  AppendTo(out, "# A 줄: idxK, idxH, 반대칭(0/1), NPN 류 수, kerClass, t1, t2, t3 (B_n x Z_2 의 표지 표; n <= NA_MAX 일 때만)\n");
  for j in [1..m] do
    L := RepresentativeTom(tom, j);
    if Length(validIdx[j]) = 0 then str := "-"; else str := JoinStringsWithSeparator(List(validIdx[j], String), ","); fi;
    AppendTo(out, "C", Tab, j, Tab, Size(L), Tab, Length(Orbits(L, [1..2^n])), Tab, aO[j], Tab, decY[j], Tab, ngood[j], Tab, str,
             Tab, decP[1][j], Tab, decP[2][j], Tab, decP[3][j], "\n");
  od;
  for j in [1..m] do
    if decY[j] = 0 then continue; fi;
    L := RepresentativeTom(tom, j);
    Kps := List([1..3], k -> KernelOfPredicate(L, g -> psis[k](g) = 0));
    sigs := Sig(L, L, j, Kps);
    AppendTo(out, "P", Tab, j, Tab, sigs[1], Tab, sigs[2], Tab, sigs[3], Tab, sigs[4], "\n");
    for N in validSubs[j] do
      sigs := Sig(L, N, j, Kps);
      AppendTo(out, "P", Tab, j, Tab, sigs[1], Tab, sigs[2], Tab, sigs[3], Tab, sigs[4], "\n");
    od;
  od;
  Print("   서명: ", Int((Runtime() - t0) / 1000), " 초\n");
  # 서명이 K 의 켤레류를 구분하는지 보려고, N_B(L) 이 유효한 phi 들(핵 N 들)에 일으키는 궤도의 수를 적는다
  AppendTo(out, "# Q 줄: j, 유효한 phi != 0 의 N_B(L)-궤도 수 (phi = 0 은 따로 하나)\n");
  for j in [1..m] do
    if decY[j] = 0 then continue; fi;
    L := RepresentativeTom(tom, j);
    if Length(validSubs[j]) = 0 then
      orbsN := [];
    else
      Nrm := Normalizer(G, L);
      orbsN := Orbits(Nrm, validSubs[j], OnPoints);
    fi;
    AppendTo(out, "Q", Tab, j, Tab, Length(orbsN), "\n");
  od;
  Print("   정규화군 궤도: ", Int((Runtime() - t0) / 1000), " 초\n");
  AppendTo(out, "# D 줄: j, 위수, 변수 위 궤도 길이, 평행이동 부분의 차원 d, 꼭짓점 궤도 수, 고정 꼭짓점 수, 구조 (aO > 0 또는 decY > 0 인 j 만)\n");
  for j in [1..m] do
    if aO[j] = 0 and decY[j] = 0 then continue; fi;
    L := RepresentativeTom(tom, j);
    hbord := Size(Image(phi, L));
    vorb := Orbits(L, [1..2^n]);
    fixv := Number(vorb, o -> Length(o) = 1);
    varorb := SortedList(List(Orbits(Image(phi, L), [1..n]), Length));
    AppendTo(out, "D", Tab, j, Tab, Size(L), Tab, varorb, Tab, LogInt(Size(L) / hbord, 2), Tab, Length(vorb), Tab, fixv, Tab, StructStr(L), "\n");
  od;
  Print("   구조 설명: ", Int((Runtime() - t0) / 1000), " 초\n");

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
        Ho := Stabilizer(K, 2^n + 1);
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
    nA := 0;  nAanti := 0;
    for i in [1..mA] do
      if decA[i] = 0 then continue; fi;
      K := RepresentativeTom(tomA, i);
      Ho := Stabilizer(K, 2^n + 1);
      Kb := Subgroup(G, List(GeneratorsOfGroup(K), g -> RestrictedPerm(g, [1..2^n])));
      Hb := Subgroup(G, List(GeneratorsOfGroup(Ho), g -> RestrictedPerm(g, [1..2^n])));
      antisym := (Size(Ho) <> Size(K));
      jK := ClassIdx(Kb);
      jH := ClassIdx(Hb);
      Kps := List([1..3], k -> KernelOfPredicate(Kb, g -> psis[k](g) = 0));
      sigs := Sig(Kb, Hb, jK, Kps);
      AppendTo(out, "A", Tab, jK, Tab, jH, Tab, B2I(antisym), Tab, decA[i], Tab,
               sigs[1], Tab, sigs[2], Tab, sigs[3], Tab, sigs[4], "\n");
      nA := nA + decA[i];
      if antisym then nAanti := nAanti + decA[i]; fi;
    od;
    Print("   방법 A: NPN 류 ", nA, ", 반대칭이 있는 류 ", nAanti, "\n");
  fi;
  Print("   끝 (", Int((Runtime() - t0) / 1000), " 초)\n");
end;

if not IsBound(NO_MAIN) then
  for n in [N_MIN .. N_MAX] do
    PairRun(n);
  od;
fi;
# (QUIT 은 if 안에서 쓸 수 없다. `< /dev/null` 로 돌리면 입력이 끝나는 대로 GAP 이 끝난다.)
