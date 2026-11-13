# goursat_class_count_run.g
# scripts/goursat_class_count.g 로 n = 1..6 에서 G_n = B_n x Z_2 의 부분군 켤레류 수를 낸다.
#
# 실행 (저장소 맨 위에서):
#   gap -q -o 6g -L tom6.ws scripts/goursat_class_count_run.g > notes/goursat_class_count_output.txt
#   tom6.ws 는 scripts/npn_pair_resolution.g 머리말에 적힌 대로 만든 B_6 의 표지 표 작업 공간이다 (레코드 PRE 가 있어야 n = 6 을 돈다).
#   작업 공간이 없으면 n = 1..5 만 돈다 (gap -q -o 6g scripts/goursat_class_count_run.g, 몇 초).
#   (작업 공간을 읽을 때 GAP 이 찍는 '#I' 한 줄은 출력 파일에서 지웠다.)
# 점검:
#   n = 1..5 : 같은 수를 G_n 의 표지 표로 먼저 구한 값 5, 27, 131, 994, 6382 와 견주고, B_n 의 켤레류 수를 OEIS A269890 의 값
#              2, 8, 33, 193, 953 과 견준다.
#   n = 6    : 켤레류 7440 개(A269890)를 모두 돈 뒤, 일부 류(U/W 의 차원 d >= 5 인 류 전부와 위수 3000 이하의 나머지에서 무작위
#              300개)의 궤도 수를 극대부분군을 직접 나열하는 다른 방법으로 다시 센다 (지수 2 부분군의 개수가 2^d - 1 인지,
#              그 개수들의 N_G(U)-궤도 수가 같은지).

Read("scripts/goursat_class_count.g");

expectedTotal := [5, 27, 131, 994, 6382];
expectedClasses := [2, 8, 33, 193, 953];

for n in [1..5] do
    G := WreathProduct(CyclicGroup(IsPermGroup, 2), SymmetricGroup(n));
    reps := List(ConjugacyClassesSubgroups(G), Representative);
    r := GoursatCount(G, reps);
    Print("n=", n, " |B_n|=", Size(G), " classes(B_n)=", r.C, " (A269890: ", expectedClasses[n], ")",
          " nonzero_orbits=", r.nonzeroOrbits, " total=", r.total, " (tom: ", expectedTotal[n], ")",
          " maxd=", r.maxd, " dist=", r.dist,
          " ", r.C = expectedClasses[n] and r.total = expectedTotal[n], "\n");
od;

if IsBound(PRE) and IsRecord(PRE) and PRE.n = 6 then
    G := UnderlyingGroup(PRE.tom);
    C := Length(MarksTom(PRE.tom));
    sumNonzero := 0; maxd := 0; dist := []; perClass := [];
    for i in [1..C] do
        U := RepresentativeTom(PRE.tom, i);
        r := OrbitsOnNonzeroHom(G, U);
        perClass[i] := [Size(U), r[1], r[2]];
        sumNonzero := sumNonzero + r[2];
        if r[1] > maxd then maxd := r[1]; fi;
        if not IsBound(dist[r[1] + 1]) then dist[r[1] + 1] := 0; fi;
        dist[r[1] + 1] := dist[r[1] + 1] + 1;
    od;
    Print("n=6 |B_n|=", Size(G), " classes(B_n)=", C, " (A269890: 7440)",
          " nonzero_orbits=", sumNonzero, " total=", 2 * C + sumNonzero,
          " maxd=", maxd, " dist=", dist, "\n");
    # 다른 방법으로 일부를 다시 센다.
    big := Filtered([1..C], i -> perClass[i][2] >= 5);
    small := Filtered([1..C], i -> perClass[i][2] < 5 and perClass[i][1] <= 3000);
    Reset(GlobalMersenneTwister, 20261113);
    pick := Set(Concatenation(big, List([1..300], k -> Random(small))));
    bad := 0;
    for i in pick do
        U := RepresentativeTom(PRE.tom, i);
        N := Normalizer(G, U);
        idx2 := Filtered(MaximalSubgroups(U), K -> Index(U, K) = 2);
        if Length(idx2) <> 2^perClass[i][2] - 1 then
            Print("MISMATCH count at class ", i, "\n"); bad := bad + 1;
        elif Length(Orbits(N, idx2, function(K, g) return K^g; end)) <> perClass[i][3] then
            Print("MISMATCH orbits at class ", i, "\n"); bad := bad + 1;
        fi;
    od;
    Print("n=6 check: ", Length(pick), " classes (d>=5: ", Length(big), ") recounted by maximal subgroups, mismatches=", bad, "\n");
else
    Print("n=6 skipped: no PRE workspace (see header)\n");
fi;
QUIT;
