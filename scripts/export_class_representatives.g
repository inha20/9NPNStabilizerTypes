# export_class_representatives.g
# picture_type_distribution.g 의 계산이 기대는 두 명제를 다른 길로 확인하고, 대조용 함수를 내보낸다.
#
# (1) 명제 A1: 격자 대칭군 W 를 지도의 기하(축의 한 칸 이동과 뒤집기, 축 맞바꿈, 층 뒤집기)에서
#     직접 만들어, "쌍 분할 p0 를 보존하는 부호 순열 전체"와 같은 군인지 본다.
# (2) 명제 A2: 안정화군 켤레류 전부에서, 기하로 만든 W 로 센 이중잉여류의 수 |H\B/W| 와
#     "H 의 S_n-상이 M_n 에 만드는 궤도 수"가 같은지 본다.
# (3) 켤레류마다 그 군을 정확히 안정화군으로 갖는 함수를 하나 찾아 진리표를 내보낸다.
#     check_pictures_by_definition.py 가 그 함수의 그림을 직접 비교해 T 를 다시 센다.
#
# 실행: gap -q scripts/export_class_representatives.g                 (n = 5, 켤레류 전부)
#       N := 6; SAMPLE := 250; 을 먼저 정의하면 n = 6 에서 표본 250개만 내보낸다 ((1),(2)는 전부 확인).

SizeScreen([4096, 24]);   # 긴 줄을 중간에 끊지 않게 한다
if not IsBound(N) then N := 5; fi;
if not IsBound(SAMPLE) then SAMPLE := 0; fi;
if not IsBound(OUT) then OUT := Concatenation("notes/class_representatives_n", String(N), ".txt"); fi;

n := N;
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
Sn := Image(phi);
k := QuoInt(n, 2);
odd := n mod 2;
p0 := Set(List([1..k], a -> [2*a-1, 2*a]));
if odd = 1 then AddSet(p0, [n]); fi;
M := Orbit(Sn, p0, OnSetsSets);

# ---- (1) 기하에서 만든 격자 대칭군 ----
gray := [[0,0],[0,1],[1,1],[1,0]];
cells := Tuples([0..3], k);
if odd = 1 then
  cells := Concatenation(List(cells, c -> [Concatenation(c, [0]), Concatenation(c, [1])]));
fi;
PointOfCell := function(c)
  local x, a;
  x := 0;
  for a in [1..k] do
    x := x + gray[c[a]+1][1] * 2^(2*a-2) + gray[c[a]+1][2] * 2^(2*a-1);
  od;
  if odd = 1 then x := x + c[k+1] * 2^(n-1); fi;
  return x + 1;
end;
CellPerm := function(g)
  local img, c;
  img := [];
  for c in cells do img[PointOfCell(c)] := PointOfCell(g(c)); od;
  return PermList(img);
end;
wgens := [];
Add(wgens, CellPerm(function(c) local d; d := ShallowCopy(c); d[1] := (d[1] + 1) mod 4; return d; end));
Add(wgens, CellPerm(function(c) local d; d := ShallowCopy(c); d[1] := (-d[1]) mod 4; return d; end));
for ax in [1..k-1] do
  Add(wgens, CellPerm(function(c) local d, t; d := ShallowCopy(c); t := d[ax]; d[ax] := d[ax+1]; d[ax+1] := t; return d; end));
od;
if odd = 1 then
  Add(wgens, CellPerm(function(c) local d; d := ShallowCopy(c); d[k+1] := 1 - d[k+1]; return d; end));
fi;
Wgeo := Group(wgens);
Wst := PreImage(phi, Stabilizer(Sn, p0, OnSetsSets));
Print("n = ", n, "\n");
Print("[A1] |W(기하)| = ", Size(Wgeo), ", 기대값 8^k k! (홀수면 2배) = ", 8^k * Factorial(k) * 2^odd, "\n");
Print("[A1] W(기하) <= B_n: ", IsSubgroup(G, Wgeo), ",  W(기하) = N : Stab(p0): ", Wgeo = Wst,
      ",  [B_n : W] = ", Index(G, Wgeo), " = |M_n| = ", Length(M), "\n");

# ---- (2) 모든 안정화군 켤레류에서 이중잉여류 수와 궤도 수 비교 ----
tom := TableOfMarks(G);
m := Length(MarksTom(tom));
fix := List([1..m], j -> 2^Length(Orbits(RepresentativeTom(tom, j), [1..2^n])));
a := DecomposedFixedPointVector(tom, fix);
cos := RightCosets(G, Wgeo);
hom2 := ActionHomomorphism(G, cos, OnRight);
realized := Filtered([1..m], j -> a[j] <> 0);
Tpred := [];  Tdc := [];  bad := 0;
for j in realized do
  H := RepresentativeTom(tom, j);
  Tpred[j] := Length(Orbits(Image(phi, H), M, OnSetsSets));
  Tdc[j] := Length(Orbits(Image(hom2, H), [1..Length(cos)]));
  if Tpred[j] <> Tdc[j] then bad := bad + 1; fi;
od;
Print("[A2] 안정화군으로 나타나는 켤레류 ", Length(realized), "개 (전체 ", m, "개) 가운데 ",
      "이중잉여류 수와 궤도 수가 다른 것: ", bad, "개\n");

# ---- (3) 켤레류마다 함수 하나 ----
FindFunction := function(H, tries)
  local orbs, t, S, mask, kk;
  orbs := Orbits(H, [1..2^n]);
  kk := Length(orbs);
  for t in [1..tries] do
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

if SAMPLE = 0 or SAMPLE >= Length(realized) then
  chosen := realized;
else
  # 대칭의 종류가 고르게 들어가도록: (T, 상의 위수, 변수 위 궤도 길이)마다 하나씩 먼저 넣고 나머지는 무작위
  chosen := [];  seen := [];
  for j in realized do
    H := Image(phi, RepresentativeTom(tom, j));
    key := [Tpred[j], Size(H), SortedList(List(Orbits(H, [1..n]), Length))];
    if not key in seen then Add(seen, key); Add(chosen, j); fi;
  od;
  rest := Difference(realized, chosen);
  while Length(chosen) < SAMPLE and Length(rest) > 0 do
    j := Random(rest);  RemoveSet(rest, j);  Add(chosen, j);
  od;
  Sort(chosen);
fi;

PrintTo(OUT, "# n = ", n, "\n");
AppendTo(OUT, "# 열: 켤레류 번호, 안정화군의 위수, 이 유형의 B_n-궤도 수, 예측한 T, 진리표\n");
AppendTo(OUT, "# 진리표는 십진수이다. 비트 x 가 f(x) 이고, x 의 i번째 비트(1이 최하위)가 변수 i 의 값이다.\n");
fails := 0;
for j in chosen do
  H := RepresentativeTom(tom, j);
  S := FindFunction(H, 400);
  if S = fail then
    fails := fails + 1;
  else
    AppendTo(OUT, j, " ", Size(H), " ", a[j], " ", Tpred[j], " ", Sum(S, p -> 2^(p-1)), "\n");
  fi;
od;
Print("[내보내기] ", Length(chosen), "개 켤레류 가운데 함수를 찾지 못한 것: ", fails, "개 -> ", OUT, "\n");
QUIT;
