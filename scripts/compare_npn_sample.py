#!/usr/bin/env python3
"""compare_npn_sample.py -- 무작위 표본(sample_npn_stabilizers.c)과 K 표의 함수 수를 견준다.

사용:  python3 -I scripts/compare_npn_sample.py <K 표 tsv> <n> <표본 출력 파일>... [--zmax 4.5] [--top 25]

K 표는 두 가지 형식을 읽는다.
  (가) enumerate_npn_stabilizers.c 의 전수 열거 결과(n <= 5, notes/npn_stabilizer_types_n<n>.tsv),
  (나) make_npn_k_table.py 가 쓴 표 (K 의 종류별; n = 6 도 있음).
표본 프로그램이 세는 '종류'는 (반대칭 유무, H 의 [위수, 순열 부분의 위수, 평행이동 차원 d, 변수 위 궤도 유형, 꼭짓점 궤도 수, 고정 꼭짓점 수],
Kbar 의 같은 여섯 값)이다. K 표의 행을 같은 열쇠로 묶어 함수 수를 더하고, 표본에서 센 개수와 견준다.
함수 하나가 열쇠 k 에 속할 확률은 (K 표의 함수 수 합)/2^(2^n) 이고, 표본 개수는 이항분포를 따른다.
기대 개수가 20 이상인 열쇠만 z 값으로 시험하고, 표본에 나왔는데 표에 없는 열쇠가 있으면 오류로 센다.
"""
import sys, math, re
from collections import defaultdict


def pat(s):
    xs = sorted(int(t) for t in re.findall(r"\d+", s))
    return int("".join(map(str, xs)))


def read_table(path, n):
    rows = []
    for ln in open(path, encoding="utf-8"):
        if ln.startswith("#") or not ln.strip():
            continue
        c = ln.rstrip("\n").split("\t")
        if c[9].startswith("["):
            # (가) 전수 열거 형식: 열 1..22
            order_k, funcs, npn, anti = int(c[0]), int(c[2]), int(c[3]), int(c[4])
            kbar, h = int(c[5]), int(c[6])
            hperm, kperm = int(c[7]), int(c[8])
            hpat, kpat = pat(c[9]), pat(c[10])
            hd, kd = int(c[12]), int(c[13])
            hvo, kvo = int(c[14]), int(c[15])
            hfix = int(c[16])
            kfix = 0 if anti else hfix
        else:
            # (나) make_npn_k_table.py 형식: 열 1..19
            order_k, funcs, npn, anti = int(c[0]), int(c[2]), int(c[3]), int(c[4])
            kbar, h = int(c[5]), int(c[6])
            kpat, kd, kvo, kfix = pat(c[7]), int(c[8]), int(c[9]), int(c[10])
            hpat, hd, hvo, hfix = pat(c[12]), int(c[13]), int(c[14]), int(c[15])
            hperm, kperm = h >> hd, kbar >> kd
        assert hperm == h >> hd and kperm == kbar >> kd, (hperm, h, hd, kperm, kbar, kd)
        key = (anti, (h, hperm, hd, hpat, hvo, hfix), (kbar, kperm, kd, kpat, kvo, kfix))
        rows.append((key, funcs, npn))
    return rows


def read_samples(paths):
    cnt = defaultdict(int)
    total = 0
    for p in paths:
        for ln in open(p, encoding="utf-8"):
            if not ln.startswith("npn "):
                continue
            t = [int(x) for x in ln.split()[1:]]
            anti = t[0]
            h = tuple(t[1:7])
            k = tuple(t[7:13])
            c = t[13]
            cnt[(anti, h, k)] += c
            total += c
    return cnt, total


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    opts = {a.split("=")[0]: a.split("=")[1] for a in sys.argv[1:] if a.startswith("--") and "=" in a}
    zmax = float(opts.get("--zmax", 4.5))
    top = int(opts.get("--top", 25))
    path, n = args[0], int(args[1])
    samples = args[2:]
    total_funcs = 2 ** (2 ** n)
    rows = read_table(path, n)
    exp = defaultdict(int)
    for key, funcs, npn in rows:
        exp[key] += funcs
    tot_exp = sum(exp.values())
    cnt, N = read_samples(samples)
    print(f"K 표: {len(rows)}행, 열쇠 {len(exp)}개, 함수 수 합 {tot_exp} (2^(2^n) = {total_funcs}, {'같다' if tot_exp == total_funcs else '다르다'})")
    print(f"표본: {N}개, 열쇠 {len(cnt)}개")
    extra = [k for k in cnt if k not in exp]
    for k in extra:
        print("표에 없는 열쇠가 표본에 나왔다:", k, cnt[k])
    tested = 0
    chi2 = 0.0
    zs = []
    for key, f in exp.items():
        p = f / total_funcs
        mu = N * p
        x = cnt.get(key, 0)
        if mu >= 20:
            z = (x - mu) / math.sqrt(mu * (1 - p))
            tested += 1
            chi2 += z * z
            zs.append((abs(z), z, key, mu, x))
    zs.sort(reverse=True)
    bad = [t for t in zs if t[0] > zmax]
    print(f"시험한 열쇠(기대 개수 >= 20): {tested}개, |z| 최대 {zs[0][0]:.2f}" if zs else "시험한 열쇠 없음")
    if tested:
        print(f"   sum z^2 = {chi2:.1f} (자유도 {tested}; 평균 {tested}, 표준편차 {math.sqrt(2 * tested):.1f}), |z| > {zmax}: {len(bad)}개")
    # 기대 개수가 작은데 표본에 나온 열쇠(=표에서 확률이 작은 열쇠가 많이 나온 경우)도 본다
    surprise = []
    for key, x in cnt.items():
        if key in exp:
            mu = N * exp[key] / total_funcs
            if mu < 20 and x > 0:
                # 포아송 꼬리 확률 P(X >= x)
                if mu == 0:
                    pv = 0.0
                else:
                    pv = 1.0 - sum(math.exp(-mu) * mu ** i / math.factorial(i) for i in range(x))
                surprise.append((pv, key, mu, x))
    surprise.sort()
    small_bad = [s for s in surprise if s[0] < 1e-6]
    print(f"기대 개수가 20 미만인데 표본에 나온 열쇠: {len(surprise)}개, 포아송 꼬리 확률 < 1e-6: {len(small_bad)}개")
    for pv, key, mu, x in small_bad[:10]:
        print("   ", key, f"기대 {mu:.3f}, 관측 {x}, P = {pv:.2e}")
    print(f"\n기대 개수가 큰 열쇠 {top}개 (열쇠: 반대칭, H[위수, 순열위수, d, 궤도유형, 꼭짓점궤도, 고정꼭짓점], Kbar[...])")
    for key, f in sorted(exp.items(), key=lambda t: -t[1])[:top]:
        p = f / total_funcs
        mu = N * p
        x = cnt.get(key, 0)
        z = (x - mu) / math.sqrt(mu * (1 - p)) if mu > 0 else float("nan")
        print(f"   기대 {mu:14.1f}  관측 {x:12d}  z {z:6.2f}   {key}")
    ok = (not extra) and (not bad) and (not small_bad) and (tot_exp == total_funcs)
    print("\n결과:", "일치 (모든 시험을 통과)" if ok else "불일치가 있다")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
