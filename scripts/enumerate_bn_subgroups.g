# enumerate_bn_subgroups.g
# B_n = Z_2^n x| S_n 의 모든 부분군 켤레류를 "표지 표 없이" 구조로 나열한다.
#   부분군 H 마다 C = H ∩ (평행이동군 N), Sigma = H 의 S_n 속 상 이 정해지고, H/C 는 N/C 의 보군(complement)이다.
#   (1) S_n 의 작용으로 N 의 부분공간 C 의 류를 모두 구하고
#   (2) 각 C 의 정규화군 Aut(C) <= S_n 의 부분군 류 Sigma 마다
#   (3) (N x| Sigma)/C 에서 N/C 의 보군 류를 모두 구해 H 로 되돌린다.
# 같은 켤레류가 두 번 나올 수는 있으나 빠지는 류는 없다 (모든 부분군은 이 세 단계 데이터로 정해지므로).
# 실행: N_ENUM := 6; OUTFILE := "/tmp/enum_n6.tsv"; KMIN := 1; Read("scripts/enumerate_bn_subgroups.g");
# 출력: G j 위수 꼭짓점 궤도 수 생성원 수 (생성원마다 n+1 개 수) -- export_bn_generators.g 와 같은 모양
if not IsBound(N_ENUM) then N_ENUM := 4; fi;
if not IsBound(OUTFILE) then OUTFILE := "/tmp/enum_bn.tsv"; fi;
if not IsBound(KMIN) then KMIN := 1; fi;
if not IsBound(FUSE) then FUSE := true; fi;
if not IsBound(REP_STEP) then REP_STEP := 1; fi;   # 부분공간 류를 나눠서 여러 프로세스로 돌릴 때
if not IsBound(REP_OFF) then REP_OFF := 0; fi;
if not IsBound(LOGFILE) then LOGFILE := false; fi;       # 진행 기록(선택): 부분공간 류마다 한 줄, 20초 넘게 걸린 Sigma 마다 한 줄
if not IsBound(STATEFILE) then STATEFILE := false; fi;   # 다시 시작용 상태(선택): 끝낸 부분공간 류의 다음 번호와 지금까지의 개수
if not IsBound(RESUME) then RESUME := false; fi;          # true 이고 STATEFILE 이 있으면 거기서 이어서 시작한다 (OUTFILE 에 덧붙인다)
NO_MAIN := true;
Read("scripts/npn_pair_resolution.g");   # MakeBn
Unbind(NO_MAIN);

EnumerateBnSubgroups := function(n)
  local R, G, pts, bit, negs, N, S, vec, elt, Rref, KeyOf, ActKey, allkeys, V, sub, key, seen, i, orb, x, y, s,
        gensS, reps, repcls, rep, Aut, sigmas, Sig, GS, Cb, Cg, hom, Q, Nq, cc, K, H, gens, str, g, cnt, j, t0, cntC,
        k, posn, lst, Cb2, Xor2, Nrm, GN, Hs, inv, keep, a, b, reject, invOf, nclasses, dupes, ri, MovedDist,
        buf, startri, st, tri, cnt0, tu, nc, si, edges, StrongInv, strong;
  t0 := Runtime();
  MovedDist := function(H)
    if Size(H) <= 3000 then return Collected(List(Elements(H), h -> NrMovedPoints(h))); fi;
    return 0;
  end;
  R := MakeBn(n);  G := R.G;
  pts := [0 .. 2^n - 1];
  bit := function(x, i) return QuoInt(x, 2^(i-1)) mod 2; end;
  negs := List([1..n], i -> PermList(List(pts, x -> (x + (1 - 2*bit(x,i)) * 2^(i-1)) + 1)));
  N := Group(negs);
  # 켤레 검사 전에 걸러내는 불변량용: 큐브의 변(꼭짓점 쌍)들
  edges := Set(Concatenation(List([1..n], i -> List([1..2^n], p -> Set([p, p^negs[i]])))));
  StrongInv := function(H)   # B_n 안의 켤레로 바뀌지 않는 값: 변 위 궤도 길이 분포, 원소별 순환 구조 분포
    local c;
    c := fail;
    if Size(H) <= 3000 then c := Collected(List(Elements(H), CycleStructurePerm)); fi;
    return [Collected(List(Orbits(H, edges, OnSets), Length)), c];
  end;
  S := Stabilizer(G, 1);
  vec := g -> (1^g) - 1;
  elt := function(c)
    local r, i;
    r := One(G);
    for i in [1..n] do if bit(c, i) = 1 then r := r * negs[i]; fi; od;
    return r;
  end;
  Xor2 := function(a, b)
    local r, i;
    r := 0;
    for i in [0 .. n-1] do
      if (QuoInt(a, 2^i) mod 2) <> (QuoInt(b, 2^i) mod 2) then r := r + 2^i; fi;
    od;
    return r;
  end;
  # F_2 위 정수 벡터들의 기약 사다리꼴 (선두 비트 내림차순)
  Rref := function(vs)
    local basis, v, b, lead, idx;
    basis := [];
    for v in vs do
      for b in basis do
        if bit(v, LogInt(b, 2) + 1) = 1 then v := Xor2(v, b); fi;
      od;
      if v <> 0 then
        lead := LogInt(v, 2) + 1;
        for idx in [1..Length(basis)] do
          if bit(basis[idx], lead) = 1 then basis[idx] := Xor2(basis[idx], v); fi;
        od;
        Add(basis, v);
        Sort(basis, function(a, b) return a > b; end);
      fi;
    od;
    return basis;
  end;
  ActKey := function(key, s) return Rref(List(key, c -> vec(elt(c)^s))); end;
  # 모든 부분공간
  V := GF(2)^n;
  allkeys := [];
  for sub in Subspaces(V) do
    Add(allkeys, Rref(List(BasisVectors(Basis(sub)), v -> Sum([1..n], i -> IntFFE(v[i]) * 2^(i-1)))));
  od;
  Sort(allkeys);
  Print("부분공간 ", Length(allkeys), "개\n");
  gensS := GeneratorsOfGroup(S);
  seen := BlistList([1..Length(allkeys)], []);
  reps := [];
  for i in [1..Length(allkeys)] do
    if not seen[i] then
      Add(reps, i);
      orb := [i];  seen[i] := true;
      j := 1;
      while j <= Length(orb) do
        for s in gensS do
          y := PositionSorted(allkeys, ActKey(allkeys[orb[j]], s));
          if not seen[y] then seen[y] := true; Add(orb, y); fi;
        od;
        j := j + 1;
      od;
    fi;
  od;
  Print("S_n 아래 부분공간 류 ", Length(reps), "개 (", Int((Runtime() - t0) / 1000), " 초)\n");
  startri := 1;  cnt := 0;
  if RESUME and STATEFILE <> false and IsExistingFile(STATEFILE) then
    st := ReadAsFunction(STATEFILE)();
    startri := st[1];  cnt := st[2];
    Print("이어서 시작: 부분공간 류 ", startri, " 번부터, 지금까지 ", cnt, "개\n");
  else
    PrintTo(OUTFILE, "# n = ", n, ", 구조로 나열한 부분군 (KMIN = ", KMIN, ")\n");
  fi;
  for ri in [1..Length(reps)] do
    if ri < startri then continue; fi;
    if (ri - 1) mod REP_STEP <> REP_OFF then continue; fi;
    tri := Runtime();  cnt0 := cnt;  buf := "";   # 부분공간 류 하나의 결과는 끝난 뒤 한꺼번에 쓴다 (중간에 멈춰도 파일에는 끝난 류만 있다)
    rep := reps[ri];
    key := allkeys[rep];
    Aut := Stabilizer(S, key, ActKey);
    sigmas := List(ConjugacyClassesSubgroups(Aut), Representative);
    si := 0;
    for Sig in sigmas do
      si := si + 1;  tu := Runtime();
      GS := Group(Concatenation(negs, GeneratorsOfGroup(Sig)));
      Cb := List(key, elt);
      if Length(Cb) = 0 then Cg := TrivialSubgroup(GS); else Cg := Subgroup(GS, Cb); fi;
      hom := NaturalHomomorphismByNormalSubgroup(GS, Cg);
      Q := Image(hom);
      Nq := Image(hom, N);
      cc := ComplementClassesRepresentatives(Q, Nq);
      nc := Length(cc);
      Hs := List(cc, K -> PreImage(hom, K));
      if Length(Hs) > 1 and FUSE then
        # 같은 (C, Sigma) 에서 나온 보군 류들은 N_{Aut(C)}(Sigma) 의 켤레로 합쳐질 수 있다
        Nrm := Normalizer(Aut, Sig);
        GN := Group(Concatenation(negs, GeneratorsOfGroup(Nrm)));
        inv := List(Hs, H -> [Length(Orbits(H, [1..2^n])), Collected(List(Orbits(H, [1..2^n]), Length)), MovedDist(H)]);
        keep := [];
        strong := List(Hs, h -> fail);   # 강한 불변량은 약한 불변량이 같을 때에만, 군마다 한 번씩 계산한다
        for a in [1..Length(Hs)] do
          reject := false;
          for b in keep do
            if inv[a] = inv[b] then
              if strong[a] = fail then strong[a] := StrongInv(Hs[a]); fi;
              if strong[b] = fail then strong[b] := StrongInv(Hs[b]); fi;
              if strong[a] = strong[b] and IsConjugate(GN, Hs[a], Hs[b]) then reject := true; break; fi;
            fi;
          od;
          if not reject then Add(keep, a); fi;
        od;
        Hs := Hs{keep};
      fi;
      for H in Hs do
        k := Length(Orbits(H, [1..2^n]));
        if k < KMIN then continue; fi;
        if Size(H) = 1 then gens := []; else gens := SmallGeneratingSet(H); fi;
        str := "";
        for g in gens do
          str := Concatenation(str, "\t", String((1^g) - 1),
                 JoinStringsWithSeparator(List([1..n], i -> Concatenation("\t", String(((2^(i-1) + 1)^g) - 1))), ""));
        od;
        cnt := cnt + 1;
        buf := Concatenation(buf, "G\t", String(cnt), "\t", String(Size(H)), "\t", String(k), "\t", String(Length(gens)), str, "\n");
      od;
      if LOGFILE <> false and Runtime() - tu > 20000 then
        AppendTo(LOGFILE, "  Sigma ", si, "/", Length(sigmas), " |Sigma|=", Size(Sig), " 보군 류 ", nc, " -> ", Length(Hs), "개, ",
                 Int((Runtime() - tu) / 1000), " 초\n");
      fi;
    od;
    AppendTo(OUTFILE, buf);
    if STATEFILE <> false then PrintTo(STATEFILE, "return [", ri + 1, ", ", cnt, "];\n"); fi;
    if LOGFILE <> false then
      AppendTo(LOGFILE, "rep ", ri, "/", Length(reps), " dim C=", Length(key), " |Aut|=", Size(Aut), " Sigma ", Length(sigmas),
               "개, 이번 ", cnt - cnt0, "개, 합 ", cnt, "개, ", Int((Runtime() - tri) / 1000), " 초\n");
    fi;
  od;
  Print("끝: 부분군 ", cnt, "개 (", Int((Runtime() - t0) / 1000), " 초)\n");
end;

EnumerateBnSubgroups(N_ENUM);
