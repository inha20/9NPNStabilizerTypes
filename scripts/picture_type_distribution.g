# picture_type_distribution.g
# n변수 부울 함수 2^(2^n)개를 '그림 종류 수' T(f)별로 정확히 센다 (모델 A: 4x4 지도를 쌓는 규약).
#
# 방법: B_n = Z_2^n : S_n 의 부분군 켤레류 전부에 대해 표지 표(table of marks)를 만들고,
#   부분군 K가 고정하는 함수의 수 2^(K의 점 궤도 수)를 표지 표로 분해해
#   "안정화군이 정확히 H(의 켤레)"인 B_n-궤도의 수 a_H 를 얻는다.
#   그 궤도의 함수들은 T(f) = (H의 S_n-상이 변수 쌍 분할 집합 M_n 위에 만드는 궤도 수) 를 갖는다
#   (7번 저장소 GeneralizationNote.md 의 명제 A2; n=4 에서는 7번 논문 정리 5).
#
# 실행: gap -q scripts/picture_type_distribution.g      (n = 3, 4, 5)
#       N_MAX := 6; 을 먼저 정의하면 n = 6 까지 시도한다.

SizeScreen([4096, 24]);   # 긴 줄을 중간에 끊지 않게 한다
if not IsBound(N_MAX) then N_MAX := 5; fi;

PictureTypeDistribution := function(n)
  local pts, bit, negs, perms, G, i, s, dirs, phi, Sn, p0, M, tom, m, fix, a, j, H,
        Hbar, T, dist, orbs, total, nf, triv, size, t0;
  t0 := Runtime();
  pts := [0 .. 2^n - 1];
  bit := function(x, i) return QuoInt(x, 2^(i-1)) mod 2; end;   # i = 1..n
  # 입력 반전: x -> x xor e_i
  negs := List([1..n], i -> PermList(List(pts, x -> (x + (1 - 2*bit(x,i)) * 2^(i-1)) + 1)));
  # 변수 순열 s: 새 점의 i번째 좌표 = 옛 점의 (i^s)번째 좌표
  perms := [];
  for s in GeneratorsOfGroup(SymmetricGroup(n)) do
    Add(perms, PermList(List(pts, x -> Sum([1..n], i -> bit(x, i^s) * 2^(i-1)) + 1)));
  od;
  G := Group(Concatenation(negs, perms));
  size := Size(G);
  # phi: B_n -> S_n, 간선의 n개 방향에 대한 작용
  dirs := List([1..n], i -> Set(List(Filtered(pts, x -> bit(x,i) = 0), x -> Set([x+1, x+2^(i-1)+1]))));
  phi := ActionHomomorphism(G, dirs, OnSetsSets);
  Sn := Image(phi);
  # 변수 쌍 분할 (n 홀수이면 외톨이 하나 포함)
  p0 := Set(List([1..QuoInt(n,2)], k -> [2*k-1, 2*k]));
  if n mod 2 = 1 then AddSet(p0, [n]); fi;
  M := Orbit(Sn, p0, OnSetsSets);
  tom := TableOfMarks(G);
  m := Length(MarksTom(tom));
  fix := List([1..m], j -> 2^Length(Orbits(RepresentativeTom(tom, j), [1..2^n])));
  a := DecomposedFixedPointVector(tom, fix);
  dist := [];  orbs := [];  total := 0;  nf := 0;  triv := 0;
  for j in [1..m] do
    if a[j] <> 0 then
      H := RepresentativeTom(tom, j);
      Hbar := Image(phi, H);
      T := Length(Orbits(Hbar, M, OnSetsSets));
      if not IsBound(dist[T]) then dist[T] := 0; orbs[T] := 0; fi;
      dist[T] := dist[T] + a[j] * size / Size(H);
      orbs[T] := orbs[T] + a[j];
      total := total + a[j];
      nf := nf + a[j] * size / Size(H);
      if Size(H) = 1 then triv := a[j] * size; fi;
    fi;
  od;
  Print("n = ", n, ": |B_n| = ", size, ", |M_n| = ", Length(M), ", 부분군 켤레류 ", m, "개, ",
        Runtime() - t0, " ms\n");
  Print("  B_n-궤도 수 합 = ", total, ", 함수 수 합 = ", nf, " (2^(2^n) = ", 2^(2^n), ")\n");
  Print("  안정화군이 자명한 함수 = ", triv, "\n");
  for T in [1..Length(dist)] do
    if IsBound(dist[T]) then
      Print("  T = ", T, ": 함수 ", dist[T], "개, B_n-궤도 ", orbs[T], "개\n");
    fi;
  od;
end;

for n in [3..N_MAX] do PictureTypeDistribution(n); od;
QUIT;
