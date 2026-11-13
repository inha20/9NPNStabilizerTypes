# goursat_class_count.g
# G_n = B_n x Z_2 의 부분군 켤레류 수를 Goursat 보조정리로 센다. B_n 의 부분군 켤레류 대표만 쓰고 G_n 의 표지 표는 쓰지 않는다.
#
# 배경. A = B_n 일 때 A x Z_2 의 부분군 L 은 둘 중 하나이다.
#   (i)  L 이 {1} x Z_2 를 포함한다: L = U x Z_2 (U 는 L 의 A 위 사영). A-켤레류 [U] 마다 하나이다.
#   (ii) L 이 {1} x Z_2 를 포함하지 않는다: L 은 어떤 준동형 phi: U -> Z_2 의 그래프 {(u, phi(u))} 이다 (U 는 사영).
# (a, c) 로 켤레하면 (U, phi) 는 (aUa^-1, phi o conj(a^-1)) 로 간다. 그래서 (ii)의 켤레류는 [U] 마다 N_A(U)-궤도의 수
# (Hom(U, Z_2) 위, 영 사상 포함)이다. 영 사상은 늘 궤도 하나이므로
#     (A x Z_2 의 부분군 켤레류 수) = 2 C + sum_[U] (Hom(U, Z_2) 의 영 아닌 원소 위 N_A(U)-궤도 수),   C = A 의 켤레류 수.
# Hom(U, Z_2) = Hom(U/W, Z_2), W = [U, U] <u^2 : u in U> 이고 U/W = F_2^d 는 U 의 가장 큰 기본 아벨 2-몫이다.
# N_A(U) 는 U/W 에 GL(d, 2) 의 부분군으로 작용한다. GL(d, 2) 의 부분군은 V 와 V* 위에서 궤도 수가 같다(원소마다 고정점 수가
# 2^dim ker(g - 1) 로 같으므로 Burnside 에 의해). 그래서 궤도 수는 U/W 의 영 아닌 벡터 위에서 센다.
#
# 쓰는 법: Read("scripts/goursat_class_count.g");  뒤에 GoursatCount(A, reps) (reps 는 A 의 부분군 켤레류 대표 목록).
#   돌려주는 값은 레코드 rec(C, nonzeroOrbits, total, maxd, dist). dist[d+1] 은 U/W 의 차원이 d 인 류의 수이다.
#   돌리는 쪽은 scripts/goursat_class_count_run.g 이다.

# [d, 궤도 수] 를 돌려준다. d = dim U/W, 궤도 수 = N_G(U) 가 Hom(U, Z_2) 의 영 아닌 원소에 일으키는 궤도의 수.
OrbitsOnNonzeroHom := function(G, U)
    local W, ind, d, N, hom, Q, bas, perms, g, rows, i, k, u, ex, v, w, p, j;
    # W = U' 와 모든 u^2 이 만드는 U 의 정규부분군 (U/W 가 가장 큰 기본 아벨 2-몫이 되도록)
    W := NormalClosure(U, ClosureGroup(DerivedSubgroup(U),
             List(GeneratorsOfGroup(U), x -> x^2)));
    ind := Index(U, W);
    d := LogInt(ind, 2);
    if ind <> 2^d then Error("U/W is not an elementary abelian 2-group"); fi;
    if d = 0 then return [0, 0]; fi;
    N := Normalizer(G, U);
    hom := NaturalHomomorphismByNormalSubgroup(U, W);
    Q := Image(hom);
    bas := IndependentGeneratorsOfAbelianGroup(Q);
    if Length(bas) <> d then Error("rank mismatch"); fi;
    # N 의 생성원 g 마다 Q 위의 F_2-선형 사상을 (행 = 기저의 상) 구해, 영 아닌 벡터 1..2^d-1 의 순열로 바꾼다.
    perms := [];
    for g in GeneratorsOfGroup(N) do
        rows := [];
        for i in [1..d] do
            u := PreImagesRepresentative(hom, bas[i]);
            ex := IndependentGeneratorExponents(Q, Image(hom, u^g));
            Add(rows, List(ex, x -> x mod 2));
        od;
        p := [];
        for k in [1..2^d - 1] do
            v := List([1..d], i -> Int(k / 2^(i-1)) mod 2);
            w := List([1..d], j -> (Sum([1..d], i -> v[i] * rows[i][j])) mod 2);
            Add(p, Sum([1..d], i -> w[i] * 2^(i-1)));
        od;
        Add(perms, PermList(p));
    od;
    if ForAny(perms, x -> x = fail) then Error("not a permutation"); fi;
    return [d, Length(Orbits(Group(Concatenation(perms, [()])), [1..2^d - 1]))];
end;

# A 의 켤레류 대표 목록 reps 에서 A x Z_2 의 부분군 켤레류 수를 센다.
GoursatCount := function(A, reps)
    local C, sumNonzero, maxd, dist, U, r;
    C := Length(reps);
    sumNonzero := 0; maxd := 0; dist := [];
    for U in reps do
        r := OrbitsOnNonzeroHom(A, U);
        sumNonzero := sumNonzero + r[2];
        if r[1] > maxd then maxd := r[1]; fi;
        if not IsBound(dist[r[1] + 1]) then dist[r[1] + 1] := 0; fi;
        dist[r[1] + 1] := dist[r[1] + 1] + 1;
    od;
    return rec(C := C, nonzeroOrbits := sumNonzero, total := 2 * C + sumNonzero,
               maxd := maxd, dist := dist);
end;
