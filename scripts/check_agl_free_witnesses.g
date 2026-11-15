# check_agl_free_witnesses.g
# ==========================
# 노트 3절 "일반 n (2)" 의 순수 평행이동군 정리에서 쓰는 "자유로운 함수가 있다" 는 말을 GAP 으로 따로 확인한다.
# (scripts/pure_translation_absent.py 가 같은 말을 파이썬으로 이미 확인했고, 이 스크립트는 다른 코드로 한 번 더 한다.)
#
# 함수 g : F_2^m -> F_2 를 값이 1 인 점들의 집합 S 로 적고, 군 G 가 F_2^m 의 점 2^m 개 위에 순열로 작용할 때
# g 의 안정화군 = Stabilizer(G, S, OnSets). "자유롭다" = 이 안정화군이 자명하다.
#
#  A. m = 5: GL(5,2) 의 극대부분군 켤레류 M 마다 F_2^5 x| M (평행이동과 M 으로 만든 군, 위수 32 |M|) 에 자유로운 함수가 있다.
#     (그래서 S 를 제자리에 두는 GL(5,2) 의 진부분군 A_S 가 어떤 극대부분군의 켤레에 들어 있으면 F_2^5 x| A_S 에도 자유로운 함수가 있다.)
#  B. m = 6: AGL(6,2) 에 자유로운 함수가 있다. (m = 7, 8 은 GAP 의 일반 집합 안정화군 계산이 너무 느려서 하지 않고,
#     scripts/check_agl_free_witnesses_nauty.py 가 nauty 로 m = 7 까지, 파이썬의 정확한 되짚기가 m = 8 까지 확인한다.)
#
# 함수는 고정된 씨앗의 무작위 부분집합에서 안정화군이 자명한 것이 나올 때까지 고른다.
# (GAP 의 극대부분군 대표원이 실행마다 다를 수 있어서 극대부분군 쪽은 어떤 함수를 찾았는지는 적지 않는다.)
#
# 실행: gap -A -q -o 3g scripts/check_agl_free_witnesses.g < /dev/null       (수 초)
# 2026-11-15. Claude가 쓴 코드이고 저자 검토 전이다.

SizeScreen([250, 24]);

truth_table := function(vecs, S)
    local x, tt, j, i, v;
    tt := 0;
    for j in S do
        v := vecs[j];
        x := Sum([1..Length(v)], i -> IntFFE(v[i]) * 2^(i-1));
        tt := tt + 2^x;
    od;
    return tt;
end;

# 평행이동 생성원 (기본벡터 e_i 만큼)
translations_of := function(vecs, m)
    local i, e, res;
    res := [];
    for i in [1..m] do
        e := ShallowCopy(vecs[1]);
        e[i] := One(GF(2));
        Add(res, PermList(List(vecs, v -> Position(vecs, v + e))));
    od;
    return res;
end;

findfree := function(rs, G, N, maxtries)
    local tries, S, st;
    for tries in [1..maxtries] do
        S := Filtered([1..N], i -> Random(rs, [0, 1]) = 1);
        st := Size(Stabilizer(G, S, OnSets));
        if st = 1 then return [tries, S]; fi;
    od;
    return fail;
end;

partA := function()
    local vecs, GLg, hom, P, cls, tr, rs, i, M, GM, res, ok;
    Print("=== A. m = 5: GL(5,2) 의 극대부분군 켤레류마다 F_2^5 x| M 에 자유로운 함수가 있다 ===\n");
    vecs := Elements(GF(2)^5);
    GLg := GL(5, 2);
    hom := ActionHomomorphism(GLg, vecs, OnRight);
    P := Image(hom);
    tr := translations_of(vecs, 5);
    cls := MaximalSubgroupClassReps(P);
    # GAP 이 고르는 대표원과 나열 순서는 실행마다 다를 수 있으므로 (|M|, 32개 점 위의 궤도 길이) 로 정렬해서 적는다
    cls := ShallowCopy(cls);
    Sort(cls, function(A, B)
        return [Size(A), SortedList(List(Orbits(A, [1..32]), Length))] < [Size(B), SortedList(List(Orbits(B, [1..32]), Length))];
    end);
    Print("  |GL(5,2)| = ", Size(P), ", 극대부분군 켤레류 ", Length(cls), "개\n");
    rs := RandomSource(IsMersenneTwister, 20261115);
    ok := true;
    for M in cls do
        GM := Group(Concatenation(tr, GeneratorsOfGroup(M)));
        res := findfree(rs, GM, 32, 2000);
        if res = fail then
            Print("  지수 ", Size(P) / Size(M), " |M| = ", Size(M), ": 찾지 못함\n");
            ok := false;
        else
            # (대표원이 실행마다 다를 수 있어 어떤 함수를 찾았는지는 적지 않는다. 찾은 함수의 안정화군 위수가 1 임은 위에서 GAP 이 계산했다.)
            Print("  지수 ", Size(P) / Size(M), "  |M| = ", Size(M), "  M 의 32개 점 위의 궤도 길이 ", SortedList(List(Orbits(M, [1..32]), Length)),
                  "  |F_2^5 x| M| = ", Size(GM), ": 자유로운 함수가 있다 (안정화군 위수 1)\n");
        fi;
    od;
    return ok;
end;

partB := function(m)
    local vecs, GLm, gens, perms, tr, G, rs, res, glorder;
    vecs := Elements(GF(2)^m);
    GLm := GL(m, 2);
    perms := List(GeneratorsOfGroup(GLm), A -> Permutation(A, vecs, OnRight));
    tr := translations_of(vecs, m);
    G := Group(Concatenation(tr, perms));
    glorder := Product([0..m-1], i -> 2^m - 2^i);
    rs := RandomSource(IsMersenneTwister, 20261115 + m);
    res := findfree(rs, G, 2^m, 200);
    if Size(G) <> 2^m * glorder then Print("  m=", m, ": 군의 위수가 AGL(m,2) 와 다르다\n"); return false; fi;
    if res = fail then
        Print("  m=", m, "  |AGL(m,2)| = ", Size(G), ": 찾지 못함\n");
        return false;
    fi;
    Print("  m=", m, "  |AGL(m,2)| = ", Size(G), ": 자유로운 함수가 있다 (안정화군 위수 1, 값이 1 인 점 ", Length(res[2]), "개)\n");
    return true;
end;

okA := partA();
Print("\n=== B. m = 6: AGL(m,2) 에 자유로운 함수가 있다 ===\n");
okB := ForAll([6], partB);
if okA and okB then Print("\n결과: 모두 찾았다\n"); else Print("\n결과: 못 찾은 경우가 있다\n"); fi;
QUIT_GAP(0);
