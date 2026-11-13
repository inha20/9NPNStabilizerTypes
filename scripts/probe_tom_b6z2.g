# probe_tom_b6z2.g
# 한도를 둔 시험: B_6 x Z_2 (꼭짓점 64개와 출력 반전이 맞바꾸는 두 점, 모두 66점 위의 순열군, 위수 92,160)의
# 표지 표를 GAP 이 만들 수 있는가. 방법 A(npn_stabilizer_types.g)를 n = 6 에 쓰려면 이 표가 필요하다.
# 실행: timeout 1800 gap -A -q -o 6g scripts/probe_tom_b6z2.g < /dev/null > probe.log 2>&1
# 끝나면 "done: classes = ..." 한 줄이 로그에 남는다. 30분 안에 끝나지 않았으면 "start" 한 줄만 남는다.
n := 6;
pts := [0 .. 2^n - 1];
bit := function(x, i) return QuoInt(x, 2^(i-1)) mod 2; end;
negs := List([1..n], i -> PermList(List(pts, x -> (x + (1 - 2*bit(x,i)) * 2^(i-1)) + 1)));
perms := [];
for s in GeneratorsOfGroup(SymmetricGroup(n)) do
  Add(perms, PermList(List(pts, x -> Sum([1..n], i -> bit(x, i^s) * 2^(i-1)) + 1)));
od;
G := Group(Concatenation(negs, perms));
swap := MappingPermListList([2^n + 1, 2^n + 2], [2^n + 2, 2^n + 1]);
GA := Group(Concatenation(GeneratorsOfGroup(G), [swap]));
Print("start: |G| = ", Size(G), ", |GA| = ", Size(GA), "\n");
t0 := Runtime();
tomA := TableOfMarks(GA);
Print("done: classes = ", Length(MarksTom(tomA)), ", ms = ", Runtime() - t0, "\n");
QUIT;
