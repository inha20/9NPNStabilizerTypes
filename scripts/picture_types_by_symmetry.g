# picture_types_by_symmetry.g
# 그림 종류 수 T 를 대칭의 종류별로 나눠 본다 (모델 A).
#   S_n 의 부분군 켤레류 L 마다:
#   (1) 위수, 구조, 변수 위의 궤도 길이, M_n 위의 궤도 수 (= 그 L 이 주는 T)
#   (2) "입력 반전을 허용한 안정화군의 순열 부분" Hbar 가 L 인 함수의 수와 B_n-궤도 수
#   (3) "변수 순열만 본 대칭군" G(f) 가 정확히 L 인 함수의 수
#       (Clote-Kranakis, Kisielewicz 의 뜻으로 L 이 부울 함수의 대칭군이 되는가)
# 실행: gap -q scripts/picture_types_by_symmetry.g     (n = 4, 5; N_MAX := 6; 을 먼저 정의하면 6까지)

SizeScreen([4096, 24]);   # 긴 줄을 중간에 끊지 않게 한다
if not IsBound(N_MAX) then N_MAX := 5; fi;

BySymmetry := function(n)
  local pts, bit, negs, perms, G, s, dirs, phi, Sn, p0, M, tom, m, fix, a, j, H, Hbar,
        cls, reps, k, T, cnt, orb, size, idx, inv, key, poss, real, i,
        S, tomS, mS, fixS, aS, cntS, L, onlyHbar, onlyG;
  pts := [0 .. 2^n - 1];
  bit := function(x, i) return QuoInt(x, 2^(i-1)) mod 2; end;
  negs := List([1..n], i -> PermList(List(pts, x -> (x + (1 - 2*bit(x,i)) * 2^(i-1)) + 1)));
  perms := [];
  for s in GeneratorsOfGroup(SymmetricGroup(n)) do
    Add(perms, PermList(List(pts, x -> Sum([1..n], i -> bit(x, i^s) * 2^(i-1)) + 1)));
  od;
  G := Group(Concatenation(negs, perms));
  size := Size(G);
  dirs := List([1..n], i -> Set(List(Filtered(pts, x -> bit(x,i) = 0), x -> Set([x+1, x+2^(i-1)+1]))));
  phi := ActionHomomorphism(G, dirs, OnSetsSets);
  Sn := Image(phi);
  p0 := Set(List([1..QuoInt(n,2)], k -> [2*k-1, 2*k]));
  if n mod 2 = 1 then AddSet(p0, [n]); fi;
  M := Orbit(Sn, p0, OnSetsSets);
  cls := ConjugacyClassesSubgroups(Sn);
  reps := List(cls, Representative);
  cnt := List(reps, x -> 0);  orb := List(reps, x -> 0);
  tom := TableOfMarks(G);
  m := Length(MarksTom(tom));
  fix := List([1..m], j -> 2^Length(Orbits(RepresentativeTom(tom, j), [1..2^n])));
  a := DecomposedFixedPointVector(tom, fix);
  for j in [1..m] do
    if a[j] <> 0 then
      H := RepresentativeTom(tom, j);
      Hbar := Image(phi, H);
      idx := First([1..Length(reps)], k -> Size(reps[k]) = Size(Hbar) and IsConjugate(Sn, reps[k], Hbar));
      cnt[idx] := cnt[idx] + a[j] * size / Size(H);
      orb[idx] := orb[idx] + a[j];
    fi;
  od;
  # 변수 순열만 본 대칭군: S = S_n 이 2^n 개 점에 작용
  S := Group(perms);
  tomS := TableOfMarks(S);
  mS := Length(MarksTom(tomS));
  fixS := List([1..mS], j -> 2^Length(Orbits(RepresentativeTom(tomS, j), [1..2^n])));
  aS := DecomposedFixedPointVector(tomS, fixS);
  cntS := List(reps, x -> 0);
  for j in [1..mS] do
    if aS[j] <> 0 then
      L := Image(phi, RepresentativeTom(tomS, j));
      idx := First([1..Length(reps)], k -> Size(reps[k]) = Size(L) and IsConjugate(Sn, reps[k], L));
      cntS[idx] := cntS[idx] + aS[j] * Size(S) / Size(L);
    fi;
  od;
  Print("n = ", n, ": S_n 의 부분군 켤레류 ", Length(reps), "개, |M_n| = ", Length(M), "\n");
  Print("  번호 | 위수 | 구조 | 변수 위 궤도 길이 | T | Hbar = L 인 함수 수 | 그 B_n-궤도 수 | G(f) = L 인 함수 수\n");
  poss := [];  real := [];
  for k in [1..Length(reps)] do
    T := Length(Orbits(reps[k], M, OnSetsSets));
    AddSet(poss, T);
    if cnt[k] > 0 then AddSet(real, T); fi;
    Print("  ", k, " | ", Size(reps[k]), " | ", StructureDescription(reps[k]), " | ",
          SortedList(List(Orbits(reps[k], [1..n]), Length)), " | ", T, " | ", cnt[k], " | ", orb[k], " | ", cntS[k], "\n");
  od;
  Print("  부분군이 줄 수 있는 T 값: ", poss, "\n");
  Print("  실제로 나타나는 T 값  : ", real, "\n");
  Print("  Hbar 로 나타나지 않는 켤레류 번호: ", Filtered([1..Length(reps)], k -> cnt[k] = 0), "\n");
  Print("  G(f) 로 나타나지 않는 켤레류 번호: ", Filtered([1..Length(reps)], k -> cntS[k] = 0), "\n");
  onlyHbar := Filtered([1..Length(reps)], k -> cnt[k] > 0 and cntS[k] = 0);
  onlyG := Filtered([1..Length(reps)], k -> cnt[k] = 0 and cntS[k] > 0);
  Print("  Hbar 로는 나타나지만 G(f) 로는 나타나지 않는 것: ", onlyHbar, "\n");
  Print("  G(f) 로는 나타나지만 Hbar 로는 나타나지 않는 것: ", onlyG, "\n");
  Print("  함수 수 합 확인: Hbar 쪽 ", Sum(cnt), ", G(f) 쪽 ", Sum(cntS), " (2^(2^n) = ", 2^(2^n), ")\n\n");
end;

if not IsBound(N_MIN) then N_MIN := 3; fi;
for n in [N_MIN..N_MAX] do BySymmetry(n); od;
QUIT;
