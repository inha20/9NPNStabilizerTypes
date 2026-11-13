#!/usr/bin/env python3
"""
circuit_lut4_exact_n5.py
========================
5변수 함수를 4입력 LUT 몇 개로 만들 수 있는지의 정확한 최솟값을 구하고, 표본 합성(circuit_sample_orbits.py)과 ABC 의 LUT 매핑
순서들이 그 값에서 얼마나 떨어져 있는지 잰다 (PLAN.md 6절 C8 (가)를 5변수에서 정확한 기준과 함께 하는 것).

정확한 값을 손으로 구할 수 있는 까닭
  5개 변수에 모두 의존하는 함수는 4입력 LUT 하나로는 안 되고, 어느 한 변수 x_i 로 갈라(Shannon 전개) 두 부분함수를 각각 LUT4 하나로,
  합치는 데 LUT3 하나를 쓰면 3개로 된다. 그러므로 최솟값은 2 또는 3이다.
  2개로 되려면 마지막 LUT 의 입력 하나가 앞 LUT 의 출력이어야 하므로 f(x) = g(h(x_A), x_B') 이고 |A| <= 4, |B'| <= 3, A 와 B' 의 합집합이
  5개 변수 전부이다. 입력을 더해서 늘려도 되므로 |A| = 4, |B'| = 3 이라 해도 된다. 곧 A = {x_i 가 아닌 네 변수}, B' = {x_i} + S (S 는 A 의
  두 변수)이다. S 의 값(네 가지)마다 A 에서 S 를 뺀 두 변수의 값(네 가지)을 행, x_i 의 값(둘)을 열로 하는 4x2 표를 만들면,
  각 표에서 서로 다른 행이 둘 이하일 때만 h 와 g 를 정할 수 있다. 필요한 까닭은 h 가 0 또는 1 뿐이라 같은 s 에서 행의 종류가
  둘을 넘을 수 없기 때문이고, 충분한 까닭은 행의 종류 둘에 h = 0, 1 을 매기고 g 를 그 표 그대로 정하면 되기 때문이다.
  모든 (i, S) 30가지 가운데 하나라도 이 조건을 만족하면 최솟값은 2, 아니면 3이다.
  이 판정은 ABC 의 `lutexact -N 5 -M 2 -K 4` (SAT 로 LUT 두 개의 회로를 찾는다)와 견주어 확인한다 (selftest).

실행:  python scripts/circuit_lut4_exact_n5.py selftest <abc 실행 파일> [<함수 수>]
       python scripts/circuit_lut4_exact_n5.py gap <표본 출력 파일> [<표본 출력 파일> ...]
       python scripts/circuit_lut4_exact_n5.py flows <abc 실행 파일> <표본 출력 파일>
       python scripts/circuit_lut4_exact_n5.py relabel <abc 실행 파일> <표본 출력 파일> [<흐름 이름의 앞머리, 기본 F2>]
  gap     : 표본 출력(circuit_sample_orbits.py run 의 결과, l4 열)에서 함수마다 변환 K개의 LUT4 수가 정확한 값과 얼마나 다른지 낸다.
  flows   : 같은 함수 목록(앞의 FLOW_M 개, 기본 200)에서 ABC 의 LUT 매핑 순서(F0~F6)와 Shannon 전개를 변환 없이 한 번씩 돌려 정확한
            값과 견준다. 순서를 바꿔 가며 돌린 최솟값을 변환 K개를 돌린 최솟값(같은 횟수)과도 견준다.
  relabel : 한 순서(기본 F2)에서도 변환(입력 이름 바꾸기)의 이득이 남는지 본다. 앞의 RELABEL_M 개(기본 100)에 대해 변환 RELABEL_K 개
            (기본은 표본 파일과 같은 수)를 만들어 그 순서로 돌린다. 앞의 16개는 표본 파일과 같은 변환이다. 순서는 F0~F6 의 앞머리로 고른다.
  환경 변수: FLOW_M (flows 가 보는 함수 수, 기본 200; F6 은 한 함수에 2초쯤 걸려 앞의 50개만), RELABEL_M, RELABEL_K.
"""
import itertools
import os
import random
import subprocess
import sys
import tempfile
import time

N = 5
SIZE = 1 << N


def support(f):
    return [i for i in range(N) if any(((f >> x) & 1) != ((f >> (x ^ (1 << i))) & 1) for x in range(SIZE))]


def exact_lut4(f):
    """f: 5변수 진리표(32비트 정수, 비트 x 가 f(x), x 의 i번째 비트가 변수 i).
    입력이 둘 이상인 LUT4 의 수의 최솟값. 의존하는 변수가 1개 이하이면 0, 4개 이하이면 1, 5개이면 2 또는 3."""
    sup = support(f)
    if len(sup) <= 1:
        return 0
    if len(sup) <= 4:
        return 1
    for i in range(N):
        others = [v for v in range(N) if v != i]
        for s_vars in itertools.combinations(others, 2):
            t_vars = [v for v in others if v not in s_vars]
            ok = True
            for s in itertools.product((0, 1), repeat=2):
                base = sum(b << v for b, v in zip(s, s_vars))
                rows = set()
                for t in itertools.product((0, 1), repeat=2):
                    x0 = base + sum(b << v for b, v in zip(t, t_vars))
                    rows.add((((f >> x0) & 1), ((f >> (x0 | (1 << i))) & 1)))
                if len(rows) > 2:
                    ok = False
                    break
            if ok:
                return 2
    return 3


def abc_two_luts(abc, f):
    """ABC 의 lutexact 로 LUT4 둘로 되는지 (True/False)."""
    r = subprocess.run([abc, "-q", "lutexact -N 5 -M 2 -K 4 %08x" % f], capture_output=True, text=True, timeout=300,
                       stdin=subprocess.DEVNULL)
    out = r.stdout + r.stderr
    if "no solution" in out:
        return False
    if "Solution" in out or "solution" in out or "Verification" in out or "LUT" in out:
        return True
    raise RuntimeError("lutexact 출력을 읽지 못했다: " + out[:300])


def selftest(abc, count):
    rng = random.Random(4242)
    funcs = []
    while len(funcs) < count:
        f = rng.getrandbits(32)
        if len(support(f)) == 5:
            funcs.append(f)
    # 2 가 나오는 함수도 충분히 있도록, 손으로 만든 2-LUT 함수도 넣는다
    for _ in range(count // 4):
        h = rng.getrandbits(16)
        g = rng.getrandbits(16)
        i = rng.randrange(N)
        others = [v for v in range(N) if v != i]
        s_vars = rng.sample(others, 2)
        t_vars = [v for v in others if v not in s_vars]
        f = 0
        for x in range(SIZE):
            hv = (h >> sum(((x >> v) & 1) << k for k, v in enumerate(s_vars + t_vars))) & 1
            gv = (g >> (hv | (((x >> i) & 1) << 1) | (sum(((x >> v) & 1) << k for k, v in enumerate(s_vars)) << 2))) & 1
            f |= gv << x
        if len(support(f)) == 5:
            funcs.append(f)
    agree = twos = 0
    for f in funcs:
        mine = exact_lut4(f) == 2
        theirs = abc_two_luts(abc, f)
        if mine != theirs:
            print("불일치: %08x  내 판정 %s, lutexact %s" % (f, mine, theirs))
            return 1
        agree += 1
        twos += mine
    print("selftest: 함수 %d개(무작위 %d + 2-LUT 로 만든 %d)에서 판정이 모두 일치, 그 가운데 LUT4 둘로 되는 것 %d개" % (
        agree, count, len(funcs) - count, twos))
    return 0


def read_samples(paths):
    rows = []
    for p in paths:
        for line in open(p):
            if line.startswith("#") or not line.strip():
                continue
            c = line.split()
            l4 = c[-1]
            vals, canon = l4.split(";")
            rows.append((int(c[0]), int(c[1], 16), [int(v) for v in vals.split(",")], int(canon)))
    return rows


def exp_min_of_k(vals, k):
    """K개 가운데 무작위로 k개를 고를 때 최솟값의 기댓값 (circuit_sample_orbits.py 와 같은 식)."""
    from math import comb
    v = sorted(vals)
    K = len(v)
    return sum(v[i] * comb(K - 1 - i, k - 1) for i in range(K - k + 1)) / comb(K, k)


def gap(paths):
    rows = read_samples(paths)
    ex = [exact_lut4(f) for _, f, _, _ in rows]
    m = len(rows)
    n2 = sum(1 for e in ex if e == 2)
    print("함수 %d개 (모두 5변수에 의존). 정확한 최솟값: LUT4 둘로 되는 것 %d개 (%.1f%%), 셋이 필요한 것 %d개. 평균 %.4f" % (
        m, n2, 100.0 * n2 / m, m - n2, sum(ex) / m))
    mean_ex = sum(ex) / m
    bad = 0
    for (idx, f, vals, canon), e in zip(rows, ex):
        if min(vals) < e or canon < e:
            bad += 1
    print("정확한 값보다 작은 합성 결과가 있는 함수: %d개 (0 이어야 한다)" % bad)
    K = len(rows[0][2])
    print("변환 %d개 (자기 자신 포함)의 LUT4 수를 정확한 값과 견줌" % K)
    ks = [k for k in (1, 2, 4, 8, 16, 32, 64, 128, 256) if k <= K]
    for k in ks:
        mean = sum(exp_min_of_k(vals, k) for _, _, vals, _ in rows) / m
        print("  k = %3d: 최솟값의 평균 %.3f  정확한 값보다 %+.3f (%+.1f%%)" % (k, mean, mean - mean_ex, 100 * (mean - mean_ex) / mean_ex))
    ident = sum(vals[0] for _, _, vals, _ in rows) / m
    canon = sum(canon for _, _, _, canon in rows) / m
    print("  자기 자신만: %.3f (%+.1f%%), 표준형만: %.3f (%+.1f%%)" % (ident, 100 * (ident - mean_ex) / mean_ex, canon, 100 * (canon - mean_ex) / mean_ex))
    reach = sum(1 for (_, _, vals, _), e in zip(rows, ex) if min(vals) == e)
    print("  변환 %d개 가운데 최솟값이 정확한 값과 같은 함수: %d개 (%.1f%%)" % (K, reach, 100.0 * reach / m))
    # Shannon 전개는 항상 LUT4 셋 이하로 만든다 (최솟값 2 또는 3 의 상한). 도구의 값이 그보다 큰 함수는 얼마나 되나
    over1 = sum(1 for _, _, vals, _ in rows if vals[0] > 3)
    overk = sum(1 for _, _, vals, _ in rows if min(vals) > 3)
    print("  Shannon 전개(셋 이하)보다 큰 함수: 자기 자신만의 값에서 %d개 (%.1f%%), 변환 %d개의 최솟값에서 %d개 (%.1f%%)" % (
        over1, 100.0 * over1 / m, K, overk, 100.0 * overk / m))
    for e in (2, 3):
        sel = [(vals, c) for (_, _, vals, c), x in zip(rows, ex) if x == e]
        if sel:
            print("  정확한 값이 %d인 함수 %d개: 자기 자신만의 평균 %.3f, 변환 %d개 최솟값의 평균 %.3f" % (
                e, len(sel), sum(v[0] for v, _ in sel) / len(sel), K, sum(min(v) for v, _ in sel) / len(sel)))
    return 0


def blif_of(f):
    names = " ".join("x%d" % i for i in range(N))
    lines = [".model f", ".inputs " + names, ".outputs y", ".names " + names + " y"]
    for x in range(SIZE):
        if (f >> x) & 1:
            lines.append("".join(str((x >> i) & 1) for i in range(N)) + " 1")
    lines.append(".end")
    return "\n".join(lines) + "\n"


FLOWS = {
    "F0 흐름 B; dch -f; if -a": "strash; {B}; dch -f; if -K 4 -a",
    "F1 F0; mfs2; mfs2": "strash; {B}; dch -f; if -K 4 -a; mfs2; mfs2",
    "F2 (dch; if; mfs2; strash) x3": "strash; {B}; " + "; ".join(["dch -f; if -K 4 -a; mfs2; strash"] * 3) + "; dch -f; if -K 4 -a",
    "F3 &dch; &if; &satlut": "strash; {B}; &get; &dch -f; &if -K 4 -a; &satlut -N 128 -C 10000; &put",
    "F4 if -a; if -a": "strash; {B}; if -K 4 -a; if -K 4 -a",
    "F5 흐름 없이 if -a": "strash; if -K 4 -a",
    "F6 흐름 없이 &deepsyn -T 2; &if; &satlut": "strash; &get; &deepsyn -T 2; &put; &get; &dch -f; &if -K 4 -a; &satlut -N 128 -C 10000; &put",
}
F6_CAP = 50          # F6 은 한 함수에 2초쯤 걸리므로 앞의 50개에서만 돌린다


def flow_command(prefix):
    """FLOWS 에서 이름이 prefix 로 시작하는 흐름의 ABC 명령 (흐름 B 를 채워 넣은 것)."""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import circuit_sample_orbits as cso
    b = "; ".join(["dc2"] + cso.COMPRESS2RS * 2 + ["dc2"])
    for name, cmd in FLOWS.items():
        if name.startswith(prefix):
            return name, cmd.replace("{B}", b)
    raise KeyError(prefix)


def blif_truth_table(text, n=N):
    """ABC 가 쓴 BLIF(입력 n개, 출력 하나, .names 노드들)를 모든 입력에서 계산해 진리표 정수(비트 x = 입력 x 에서의 값)를 낸다.
    입력 i 는 .inputs 줄의 i번째 이름이다 (ABC 는 &get/&put 을 거치면 이름을 pi0.. 로 바꾸지만 순서는 지킨다)."""
    inputs, outputs, nodes = None, None, []
    for line in text.replace("\\\n", " ").splitlines():
        t = line.split()
        if not t or t[0].startswith("#"):
            continue
        if t[0] == ".inputs":
            inputs = t[1:]
        elif t[0] == ".outputs":
            outputs = t[1:]
        elif t[0] == ".names":
            nodes.append((t[1:-1], t[-1], []))
        elif t[0] == ".end":
            break
        elif t[0] == ".model":
            continue
        elif t[0].startswith("."):
            raise ValueError("읽지 못하는 BLIF 줄: " + line)
        else:
            nodes[-1][2].append((t[0], t[1]) if len(t) == 2 else ("", t[0]))
    if inputs is None or outputs is None or len(inputs) != n or len(outputs) != 1:
        raise ValueError("입력 %d개, 출력 1개가 아니다" % n)
    index = {name: i for i, name in enumerate(inputs)}
    by_out = {out: (ins, cubes) for ins, out, cubes in nodes}

    def value(name, x, memo):
        if name in memo:
            return memo[name]
        if name in index:
            v = (x >> index[name]) & 1
        else:
            ins, cubes = by_out[name]
            vals = [value(a, x, memo) for a in ins]
            if not cubes:
                v = 0
            else:
                phase = int(cubes[0][1])
                hit = any(all(c == "-" or int(c) == b for c, b in zip(cube, vals)) for cube, _ in cubes)
                v = phase if hit else 1 - phase
        memo[name] = v
        return v

    return sum(value(outputs[0], x, {}) << x for x in range(1 << n))


def run_flow(abc, cmd, funcs, chunk=50):
    """funcs: 진리표(정수)의 목록. 각각에 ABC 명령 cmd 를 돌려 LUT 수(입력이 둘 이상인 .names 의 수)의 목록, 함수당 ABC 실행 시간(초),
    결과가 원래 함수와 다른 함수의 수(BLIF 를 모든 입력에서 계산해 견준다)를 낸다."""
    counts = []
    bad = 0
    t_abc = 0.0
    for start in range(0, len(funcs), chunk):
        part = funcs[start:start + chunk]
        with tempfile.TemporaryDirectory() as d:
            script = []
            for j, f in enumerate(part):
                with open(os.path.join(d, "%d.blif" % j), "w") as fh:
                    fh.write(blif_of(f))
                script += ["read_blif %d.blif" % j, cmd, "write_blif %d.out.blif" % j]
            with open(os.path.join(d, "s.abc"), "w") as fh:
                fh.write("\n".join(script) + "\n")
            t0 = time.perf_counter()
            r = subprocess.run([os.path.abspath(abc), "-q", "source s.abc"], cwd=d, stdin=subprocess.DEVNULL,
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=3000)
            t_abc += time.perf_counter() - t0
            for j, f in enumerate(part):
                out_path = os.path.join(d, "%d.out.blif" % j)
                if not os.path.exists(out_path):
                    raise RuntimeError("ABC 가 결과를 쓰지 못했다 (반환값 %s, 묶음 %d, 함수 %08x)" % (r.returncode, start // chunk, f))
                with open(out_path) as fh:
                    text = fh.read()
                luts = sum(1 for line in text.splitlines() if line.startswith(".names") and len(line.split()) >= 4)
                counts.append(luts)
                if blif_truth_table(text) != f:
                    bad += 1
    return counts, t_abc / len(funcs), bad


def flows(abc, path):
    rows = read_samples([path])
    rows = rows[:int(os.environ.get("FLOW_M", "200"))]
    funcs = [f for _, f, _, _ in rows]
    ex = [exact_lut4(f) for f in funcs]
    m = len(rows)
    mean_ex = sum(ex) / m
    print("함수 %d개. 정확한 최솟값의 평균 %.4f (LUT4 둘로 되는 함수 %d개)" % (m, mean_ex, sum(1 for e in ex if e == 2)), flush=True)
    print("  Shannon 전개(LUT4 셋) 평균 3.000 (정확한 값보다 %+.1f%%)" % (100 * (3 - mean_ex) / mean_ex), flush=True)
    counts, times = {}, {}
    for name in FLOWS:
        cap = min(m, F6_CAP) if name.startswith("F6") else m
        _, cmd = flow_command(name[:2])
        c, t, wrong = run_flow(abc, cmd, funcs[:cap])
        counts[name], times[name] = c, t
        me = sum(ex[:cap]) / cap
        mean = sum(c) / cap
        reach = sum(1 for x, e in zip(c, ex) if x == e)
        below = sum(1 for x, e in zip(c, ex) if x < e)
        print("  %-40s 함수 %3d개: 평균 %.3f (정확한 값 %.3f 보다 %+.1f%%), 정확한 값에 닿은 함수 %d개 (%.1f%%), 함수당 ABC %.3f초%s" % (
            name, cap, mean, me, 100 * (mean - me) / me, reach, 100.0 * reach / cap, t,
            "  [정확한 값보다 작은 것 %d개: 오류]" % below if below else ""), flush=True)
        print("      검증: 결과 %d개를 BLIF 계산으로 원래 함수와 견줌, 다른 것 %d개%s" % (cap, wrong, "" if wrong == 0 else "  [오류]"), flush=True)
    # 같은 함수에서 F0 과 견준 짝지은 차이 (함수마다 LUT 수의 차이, 평균과 표준오차)
    f0 = [n for n in FLOWS if n.startswith("F0")][0]
    print("  F0 과 견준 짝지은 차이 (함수마다 [그 흐름의 LUT 수 - F0 의 LUT 수], 평균 ± 표준오차):", flush=True)
    for name in FLOWS:
        if name == f0:
            continue
        cap = len(counts[name])
        d = [counts[name][j] - counts[f0][j] for j in range(cap)]
        md = sum(d) / cap
        se = (sum((x - md) ** 2 for x in d) / (cap - 1) / cap) ** 0.5
        print("    %s 함수 %3d개: %+.3f ± %.3f" % (name[:2], cap, md, se), flush=True)
    # 표본 파일의 l4 열과 F0 이 같은 흐름인지 (자기 자신 하나)
    diff = sum(1 for (_, _, vals, _), c in zip(rows, counts[f0]) if vals[0] != c)
    print("  검사: F0 의 값이 표본 파일의 (변환하지 않은) l4 값과 다른 함수 %d개 (0 이어야 한다)" % diff, flush=True)
    # 흐름을 바꿔 가며 돌린 최솟값 대 변환을 바꿔 가며 돌린 최솟값 (같은 횟수)
    cheap = [n for n in FLOWS if not n.startswith("F6")]
    best = [min(counts[n][j] for n in cheap) for j in range(m)]
    kk = len(cheap)
    mb = sum(best) / m
    print("같은 횟수(%d번)로 견주기, 함수 %d개" % (kk, m), flush=True)
    print("  흐름 F0~F%d 를 한 번씩 돌린 최솟값: 평균 %.3f (정확한 값보다 %+.1f%%), 정확한 값에 닿은 함수 %d개 (%.1f%%)" % (
        kk - 1, mb, 100 * (mb - mean_ex) / mean_ex, sum(1 for b, e in zip(best, ex) if b == e),
        100.0 * sum(1 for b, e in zip(best, ex) if b == e) / m), flush=True)
    mr = sum(exp_min_of_k(vals, kk) for _, _, vals, _ in rows) / m
    print("  F0 을 변환 %d개(무작위)에 돌린 최솟값의 기댓값: 평균 %.3f (정확한 값보다 %+.1f%%)" % (kk, mr, 100 * (mr - mean_ex) / mean_ex), flush=True)
    both = [min(min(vals), b) for (_, _, vals, _), b in zip(rows, best)]
    mx = sum(both) / m
    print("  둘을 모두 쓴 최솟값(변환 16개의 F0, 그리고 F0~F%d 를 변환 없이): 평균 %.3f (정확한 값보다 %+.1f%%)" % (
        kk - 1, mx, 100 * (mx - mean_ex) / mean_ex), flush=True)
    print("  흐름별로 최솟값을 가진 흐름 수(같으면 모두 센다): " + ", ".join(
        "%s %d" % (n[:2], sum(1 for j in range(m) if counts[n][j] == best[j])) for n in cheap), flush=True)
    print("# 함수별 값: 번호 진리표 정확한값 " + " ".join(n[:2] for n in FLOWS))
    for j in range(m):
        print("%d %08x %d %s" % (rows[j][0], funcs[j], ex[j], " ".join(str(counts[n][j]) if j < len(counts[n]) else "-" for n in FLOWS)))
    return 0


def relabel(abc, path, prefix="F2"):
    """더 강한 순서에서도 변환(입력 이름 바꾸기)의 이득이 남는가: 앞의 RELABEL_M 개 함수에 대해 표본 생성과 같은 변환 16개를 만들어
    그 순서(기본 F2)를 돌리고, F0 의 표본 결과(변환 16개)와 견준다."""
    import random as _r
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import circuit_sample_orbits as cso
    rows = read_samples([path])
    rows = rows[:int(os.environ.get("RELABEL_M", "100"))]
    header = open(path).readline()
    seed = int(header.split("seed=")[1].split()[0])
    space = cso.Space(N)
    group_n = len(space.maps)
    name, cmd = flow_command(prefix)
    ex = [exact_lut4(f) for _, f, _, _ in rows]
    m = len(rows)
    mean_ex = sum(ex) / m
    K0 = len(rows[0][2])                                  # 표본 파일의 변환 수 (F0 의 값이 있는 만큼)
    K = int(os.environ.get("RELABEL_K", str(K0)))          # 이 흐름에 쓸 변환 수. 앞의 K0 개는 표본 파일과 같은 변환이다
    flat = []
    for idx, f, _, _ in rows:
        vrng = _r.Random(seed * 1000003 + idx)
        variants = [f] + [space.act(f, vrng.randrange(group_n), vrng.getrandbits(1)) for _ in range(K - 1)]
        flat += variants
    c, t, wrong = run_flow(abc, cmd, flat)
    per = [c[i * K:(i + 1) * K] for i in range(m)]
    # 검사: 다른 흐름의 값은 F0 과 같을 필요가 없으므로 정확한 값보다 작은 것이 없는지만 본다.
    bad = sum(1 for p, e in zip(per, ex) if min(p) < e)
    # F0 이면 앞의 K0 개 변환이 표본 파일과 같은 변환이므로 표본기가 낸 l4 값과 하나씩 견줄 수 있다
    # (표본기는 결과를 원래 함수와 견주지 않으므로, BLIF 계산으로 검증한 이 실행이 표본기 값을 검증하는 셈이다).
    sample_check = None
    if name.startswith("F0"):
        kk = min(K, K0)
        sample_check = (kk * m, sum(1 for p, (_, _, vals, _) in zip(per, rows) for a, b in zip(p[:kk], vals[:kk]) if a != b))
    print("함수 %d개 x 변환 %d개, 흐름 %s (함수당 ABC %.3f초). 정확한 값의 평균 %.4f. 정확한 값보다 작은 결과가 있는 함수: %d개 (0 이어야 한다)" % (
        m, K, name, t, mean_ex, bad), flush=True)
    print("  검증: 변환한 함수 %d개의 결과를 모두 BLIF 계산으로 원래 함수와 견줌, 다른 것 %d개%s" % (
        len(flat), wrong, "" if wrong == 0 else "  [오류]"), flush=True)
    if sample_check is not None:
        print("  검사: 표본 파일의 l4 값 %d개(함수마다 앞의 %d개 변환)와 다른 값 %d개 (0 이어야 한다)" % (
            sample_check[0], sample_check[0] // m, sample_check[1]), flush=True)
    print("            %-30s %-30s" % ("F0 (표본, 변환 %d개)" % K0, "%s (변환 %d개)" % (name[:2], K)))
    for k in [k for k in (1, 2, 4, 6, 8, 16, 32, 64, 128, 256) if k <= K]:
        b = sum(exp_min_of_k(p, k) for p in per) / m
        if k <= K0:
            a = sum(exp_min_of_k(vals, k) for _, _, vals, _ in rows) / m
            left = "평균 %.3f (%+.1f%%)" % (a, 100 * (a - mean_ex) / mean_ex)
        else:
            left = "-"
        print("  k = %3d    %-30s 평균 %.3f (%+.1f%%)" % (k, left, b, 100 * (b - mean_ex) / mean_ex), flush=True)
    r0 = sum(1 for (_, _, vals, _), e in zip(rows, ex) if min(vals) == e)
    r1 = sum(1 for p, e in zip(per, ex) if min(p) == e)
    print("  변환을 모두 쓴 최솟값이 정확한 값과 같은 함수: F0 (변환 %d개) %d개 (%.1f%%), %s (변환 %d개) %d개 (%.1f%%)" % (
        K0, r0, 100.0 * r0 / m, name[:2], K, r1, 100.0 * r1 / m), flush=True)
    const = sum(1 for p in per if len(set(p)) == 1)
    const0 = sum(1 for _, _, vals, _ in rows if len(set(vals)) == 1)
    print("  변환의 값이 모두 같은 함수: F0 (변환 %d개) %d개, %s (변환 %d개) %d개" % (K0, const0, name[:2], K, const), flush=True)
    if K <= 16:
        print("# 함수별 값: 번호 진리표 정확한값 변환 %d개의 %s 값(쉼표)" % (K, name[:2]))
        for (idx, f, _, _), e, p in zip(rows, ex, per):
            print("%d %08x %d %s" % (idx, f, e, ",".join(str(v) for v in p)))
    else:
        print("# 함수별 값: 번호 진리표 정확한값 변환 %d개의 %s 값의 최솟값 평균" % (K, name[:2]))
        for (idx, f, _, _), e, p in zip(rows, ex, per):
            print("%d %08x %d %d %.3f" % (idx, f, e, min(p), sum(p) / len(p)))
    return 0


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    cmd = sys.argv[1]
    if cmd == "selftest":
        sys.exit(selftest(os.path.abspath(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 40))
    if cmd == "gap":
        sys.exit(gap(sys.argv[2:]))
    if cmd == "flows":
        sys.exit(flows(sys.argv[2], sys.argv[3]))
    if cmd == "relabel":
        sys.exit(relabel(sys.argv[2], sys.argv[3], sys.argv[4] if len(sys.argv) > 4 else "F2"))
    sys.exit(__doc__)
