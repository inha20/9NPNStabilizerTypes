# enum_fix_prep.g
# B_6 의 켤레류 j 마다, 유효한 준동형 phi (phi = 0 포함)마다 "L 이 phi 로 고정하는 함수"의 매개변수를 내보낸다.
#   f(x) = b[궤도(x)] xor 부호(x),  부호(x) = phi(t) (t in L, x0^t = x, x0 은 궤도 대표),  b 는 궤도마다 자유로운 비트.
# 이 함수들 가운데 안정화군이 정확히 L 인 것의 수 E(L, phi) 는 C 프로그램(sample_npn_stabilizers.c enumfix)이 센다.
# 그러면 L 을 Kbar 로 가지고 서명이 s 인 NPN 류의 수는  |L| / (2 |N_B(L)|) * sum_{phi: 서명 s} E(L, phi)  이다
#   (K = {(g, phi(g))} 의 G_n 안 정규화군은 N_B(L, phi) x Z_2 이고, 함수가 정확히 K 인 것을 E(K) 개로 세면 류의 수는 E(K) |K| / |N_G(K)|;
#    phi 를 N_B(L)-궤도 전체에서 더하면 정규화군의 지수와 궤도 크기가 맞물려 위 식이 된다).
#
# 실행(저장소 맨 위 폴더에서): 먼저 표지 표를 작업 공간으로 저장해 둔다(npn_pair_resolution.g 머리말).
#   gap -L tom6.ws -A -q 파일.g < /dev/null   (파일.g: JOBS := [j 번호들]; OUTFILE := "..."; Read("scripts/enum_fix_prep.g");)
# 켤레류 번호는 저장한 표지 표의 것이다 (npn_pair_data_n6.tsv 의 j 와 같다).

SizeScreen([4096, 24]);
Unbind(NO_MAIN);
NO_MAIN := true;
Read("scripts/npn_pair_resolution.g");   # IndexTwoSubgroups, KernelOfPredicate 를 쓴다 (NO_MAIN 이라 PairRun 은 돌지 않는다)

EnumFixPrep := function(tom, psis, phi, n, jobs, outfile)
  local G, m, keyOf, classKey, subsOf, ClassIdxIn, Sig, j, L, oL, vN, Kps, phis, nrm, idx, N, sg, orb, sgn, k, x0, x, t,
        ints, Join;
  G := UnderlyingGroup(tom);
  m := Length(MarksTom(tom));
  keyOf := function(U)
    return [Size(U), SortedList(List(Orbits(U, [1..2^n]), Length)), Size(Image(phi, U))];
  end;
  classKey := List([1..m], j -> keyOf(RepresentativeTom(tom, j)));
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
  Sig := function(L, N, jL, Kps)
    local kc, t, k, Kp, N3;
    if Size(N) = Size(L) then kc := 0; else kc := ClassIdxIn(N, jL); fi;
    t := [];
    for k in [1..3] do
      Kp := Kps[k];
      if Size(Kp) = Size(L) then
        if kc = 0 then Add(t, -1); else Add(t, kc); fi;
      elif Size(N) = Size(L) then
        Add(t, ClassIdxIn(Kp, jL));
      elif N = Kp then
        Add(t, -1);
      else
        N3 := KernelOfPredicate(L, g -> (g in N) = (g in Kp));
        Add(t, ClassIdxIn(N3, jL));
      fi;
    od;
    return Concatenation([kc], t);
  end;
  Join := function(l) return JoinStringsWithSeparator(List(l, String), " "); end;

  PrintTo(outfile, "");
  for j in jobs do
    L := RepresentativeTom(tom, j);
    oL := Orbits(L, [1..2^n]);
    vN := Filtered(IndexTwoSubgroups(L), N -> Length(Orbits(N, [1..2^n])) = 2 * Length(oL));
    Kps := List([1..3], k -> KernelOfPredicate(L, g -> psis[k](g) = 0));
    phis := Concatenation([L], vN);
    nrm := Size(Normalizer(G, L));
    AppendTo(outfile, "J ", j, " ", Size(L), " ", nrm, " ", Length(oL), " ", Length(phis), "\n");
    for idx in [1 .. Length(phis)] do
      N := phis[idx];
      sg := Sig(L, N, j, Kps);
      orb := List([1 .. 2^n], x -> 0);
      sgn := List([1 .. 2^n], x -> 0);
      for k in [1 .. Length(oL)] do
        x0 := oL[k][1];
        for x in oL[k] do
          t := RepresentativeAction(L, x0, x);
          orb[x] := k - 1;
          if Size(N) = Size(L) or t in N then sgn[x] := 0; else sgn[x] := 1; fi;
        od;
      od;
      AppendTo(outfile, "F ", idx, " ", Join(sg), "\n");
      AppendTo(outfile, "V ", Join(orb), "\n");
      AppendTo(outfile, "S ", Join(sgn), "\n");
    od;
  od;
end;

EnumFixPrep(PRE.tom, PRE.psis, PRE.phi, PRE.n, JOBS, OUTFILE);
Print("작업 파일을 썼다: ", OUTFILE, "\n");
