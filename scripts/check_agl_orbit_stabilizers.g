# check_agl_orbit_stabilizers.g
# ===========================
# scripts/agl52_orbits.c 의 결과를 GAP 으로 따로 점검한다.
#
# agl52_orbits.c 는 AGL(5,2) (위수 319,979,520) 와 GL(5,2) (위수 9,999,360, -DLINEAR_ONLY) 가 5변수 부울 함수 2^32 개 위에서
# 하는 작용의 궤도를 너비우선탐색으로 전부 나열하고 궤도마다 안정화군 위수(= |G| / 궤도 크기)를 적는다.
# 이 스크립트는 그 출력 파일(notes/agl52_orbits_output.txt, notes/gl52_orbits_output.txt)의 궤도마다
#   (1) 대표 함수의 지지집합(값이 1인 점들)의 안정화군을 GAP 의 Stabilizer(G, 집합, OnSets) 로 구해 위수를 C 의 값과 견주고,
#   (2) 궤도 크기의 합이 2^32 인지 보고,
#   (3) 안정화군이 자명한(자유) 궤도의 수를 센다.
# 군은 C 프로그램의 생성원과 같은 순열로 만들고 위수를 먼저 확인한다 (생성원이 군 전체를 만드는지의 점검).
#
# 점 x (0..31) 는 GAP 에서 x+1. 함수 f 는 32 비트 진리표(비트 x = 점 x 의 값). 지지집합 = {x+1 : 비트 x 가 1}.
# GL(5,2) 는 점 0 을 고정하므로 점 0 의 값은 안정화군에 영향을 주지 않는다 (그래서 자유 궤도 수가 짝수이다).
#
# 실행 (저장소 루트에서): gap -A -q -o 3g scripts/check_agl_orbit_stabilizers.g
# 입력 경로를 바꾸려면 CHECK_DIR := "경로";  를 먼저 두고 Read 한다 (기본값 "notes").
#
# 2026-11-15. Claude가 쓴 코드이고 저자 검토 전이다.

if not IsBound(CHECK_DIR) then CHECK_DIR := "notes"; fi;
SizeScreen([250, 24]);

# 점 0..31 위의 순열 생성원 (agl52_orbits.c 의 build_tables 와 같은 정의)
bitsof := function(x) return List([0..4], i -> Int(x / 2^i) mod 2); end;
frombits := function(b) return Sum([0..4], i -> b[i+1] * 2^i); end;
xor := function(a, b) return Sum([0..4], i -> ((Int(a / 2^i) + Int(b / 2^i)) mod 2) * 2^i); end;

gen_translation := PermList(List([0..31], x -> xor(x, 1) + 1));
gen_cycle := PermList(List([0..31], function(x)
    local b, y, i;
    b := bitsof(x); y := [0,0,0,0,0];
    for i in [0..4] do y[((i+1) mod 5) + 1] := b[i+1]; od;
    return frombits(y) + 1; end));
gen_swap := PermList(List([0..31], function(x)
    local b; b := bitsof(x); return frombits([b[2], b[1], b[3], b[4], b[5]]) + 1; end));
gen_shear := PermList(List([0..31], function(x)
    local b; b := bitsof(x); return frombits([b[1], (b[2] + b[1]) mod 2, b[3], b[4], b[5]]) + 1; end));

AGL52 := Group(gen_translation, gen_cycle, gen_swap, gen_shear);
GL52 := Group(gen_cycle, gen_swap, gen_shear);
Print("|<4 개 생성원>| = ", Size(AGL52), "  (AGL(5,2) 의 위수 32 * 9999360 = ", 32 * 9999360, ")\n");
Print("|<3 개 생성원>| = ", Size(GL52), "  (GL(5,2) 의 위수 9999360)\n");
if Size(AGL52) <> 32 * 9999360 or Size(GL52) <> 9999360 then Print("생성원이 군 전체를 만들지 않는다\n"); QUIT_GAP(1); fi;

# 출력 파일 한 개를 점검한다
check_file := function(path, G, label)
    local f, line, parts, rep, stab, size, tot, bad, n, free, support, sz, line2;
    f := InputTextFile(path);
    if f = fail then Print("파일을 못 연다: ", path, "\n"); QUIT_GAP(1); fi;
    tot := 0; bad := 0; n := 0; free := 0;
    line := ReadLine(f);
    while line <> fail do
        if Length(line) > 6 and line{[1..6]} = "orbit " then
            line2 := Filtered(line, c -> c <> '\n');
            parts := SplitString(line2, " ");
            # orbit <i> rep=XXXXXXXX size=S stab=ST weight=W
            rep := IntHexString(parts[3]{[5..Length(parts[3])]});
            size := Int(parts[4]{[6..Length(parts[4])]});
            stab := Int(parts[5]{[6..Length(parts[5])]});
            support := Filtered([0..31], x -> (QuoInt(rep, 2^x) mod 2) = 1);
            sz := Size(Stabilizer(G, List(support, x -> x + 1), OnSets));
            if sz <> stab or size * stab <> Size(G) then
                bad := bad + 1;
                Print("  불일치: ", line2, "  GAP 의 안정화군 위수 ", sz, "\n");
            fi;
            n := n + 1; tot := tot + size;
            if stab = 1 then free := free + 1; fi;
        fi;
        line := ReadLine(f);
    od;
    CloseStream(f);
    Print(label, ": 궤도 ", n, "개, 크기의 합 ", tot, " (2^32 = ", 2^32, ": ", tot = 2^32, "), 자유 궤도 ", free, "개, 안정화군 위수가 어긋난 궤도 ", bad, "개\n");
    return [n, tot, free, bad];
end;

r1 := check_file(Concatenation(CHECK_DIR, "/agl52_orbits_output.txt"), AGL52, "AGL(5,2)");
r2 := check_file(Concatenation(CHECK_DIR, "/gl52_orbits_output.txt"), GL52, "GL(5,2) ");
if r1[4] = 0 and r2[4] = 0 and r1[2] = 2^32 and r2[2] = 2^32 then
    Print("결과: 모두 일치\n");
else
    Print("결과: 불일치가 있다\n");
fi;
QUIT_GAP(0);
