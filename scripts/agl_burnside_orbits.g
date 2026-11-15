# AGL(m,2) 와 GL(m,2) 가 m 변수 부울 함수 2^(2^m) 개 위에서 가지는 궤도 수 (번사이드 보조정리).
# 기대값: AGL(m,2) 는 m=2: 5, m=3: 10, m=4: 32, m=5: 382 (scripts/agl52_orbits.c 의 전수 열거와 견준다).
#         GL(m,2) 는 점 0 을 고정하므로 점 0 의 값을 자유로 두고 센다 (scripts/agl52_orbits.c -DLINEAR_ONLY 와 견준다).
# 점은 0..2^m-1 (m 비트 벡터), 군은 평행이동 x->x^1, 좌표 순환, 좌표 교환 (0 1), 전단 x1->x1+x0 이 만든다.
# 실행: gap -A -q -o 3g scripts/agl_burnside_orbits.g   (1초)
# 2026-11-15. Claude가 쓴 코드이고 저자 검토 전이다.
xor := function(a, b)
  local r, i;
  r := 0;
  for i in [0..11] do
    if (QuoInt(a, 2^i) + QuoInt(b, 2^i)) mod 2 = 1 then r := r + 2^i; fi;
  od;
  return r;
end;
runm := function(m)
  local N, trans, cycle, swap, tv, G, glorder, cc, tot, L;
  N := 2^m;
  trans := PermList(List([0..N-1], x -> xor(x, 1) + 1));
  cycle := PermList(List([0..N-1], function(x)
              local y, i;
              y := 0;
              for i in [0..m-1] do
                if QuoInt(x, 2^i) mod 2 = 1 then y := y + 2^((i+1) mod m); fi;
              od;
              return y + 1; end));
  swap := PermList(List([0..N-1], function(x)
              local b0, b1;
              b0 := x mod 2; b1 := QuoInt(x, 2) mod 2;
              return x - b0 - 2*b1 + 2*b0 + b1 + 1; end));
  tv := PermList(List([0..N-1], x -> xor(x, 2*(x mod 2)) + 1));
  G := Group(trans, cycle, swap, tv);
  glorder := Product([0..m-1], i -> 2^m - 2^i);
  Print("m=", m, " |G|=", Size(G), " (2^m*|GL(m,2)| = ", N*glorder, ")\n");
  cc := ConjugacyClasses(G);
  tot := Sum(cc, c -> Size(c) * 2^(Length(Cycles(Representative(c), [1..N]))));
  Print("   AGL(", m, ",2): 켤레류 ", Length(cc), "개, 함수 위의 궤도 수 = ", tot / Size(G), "\n");
  L := Group(cycle, swap, tv);
  cc := ConjugacyClasses(L);
  tot := Sum(cc, c -> Size(c) * 2^(Length(Cycles(Representative(c), [1..N]))));
  Print("   GL(", m, ",2) (위수 ", Size(L), "): 켤레류 ", Length(cc), "개, 함수 위의 궤도 수 = ", tot / Size(L), "\n");
end;
for m in [2, 3, 4, 5] do runm(m); od;
QUIT;
