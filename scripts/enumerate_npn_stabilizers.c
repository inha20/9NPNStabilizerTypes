/*
 * enumerate_npn_stabilizers.c
 * ===========================
 * n변수 부울 함수 2^(2^n)개를 하나도 빼지 않고 훑어서 (n = 3, 4, 5), 함수마다 NPN 안정화군을 직접 구하고
 * 군의 종류별 함수 수를 센다. enumerate_full_stabilizers.c 에 반대칭(출력 반전을 겸한 대칭)을 덧붙인 것이다.
 *
 *   H(f) = { (pi, c) : f(pi(x xor c)) = f(x)      모든 x }    대칭      (입력 반전 c, 변수 순열 pi)
 *   A(f) = { (pi, c) : f(pi(x xor c)) = not f(x)  모든 x }    반대칭
 *   Kbar(f) = H(f) 합 A(f).  A(f) 는 비어 있거나 H(f) 의 잉여류이다 (프로그램이 크기를 점검한다).
 *
 * 반대칭이 있으려면 f 의 무게(f = 1 인 점의 수)가 2^(n-1) 이어야 하므로 그런 함수만 반대칭을 찾는다.
 * 군의 종류는 H 와 Kbar 각각에 대해 여섯 값 (위수, 순열 부분의 위수, 평행이동 부분의 차원 d,
 * 순열 부분의 변수 위 궤도 유형, 꼭짓점 위 궤도 수, 고정 꼭짓점 수)을 적고 반대칭 유무를 더한다.
 * GAP 의 표지 표 계산(npn_stabilizer_types.g)과 독립인 검증이다. compare 는 scripts/check_npn_bruteforce.py 가 한다.
 *
 * 거르기. 변수 순열 p 가 대칭이 되려면 d1, d2 (enumerate_full_stabilizers.c 참고)가 맞아야 한다.
 * 반대칭의 경우 f 가 균형이면 d1 은 f 와 not f 에서 같고 (모서리 방향마다 11 쌍 수 = 00 쌍 수), d2 는
 * f 와 not f 의 값을 견주어야 한다. f o P_p 의 i 방향은 f 의 (p 의 역순열)(i) 방향이므로 역순열을 쓴다.
 * 거르기가 반대칭을 놓치지 않는지는 일정 간격으로 뽑은 함수에서 모든 순열을 거르지 않고 시험해 비교한다 (다르면 멈춘다).
 *
 * f 와 not f 는 군이 같으므로 최상위 비트가 0 인 절반만 훑는다 (합칠 때 2배 한다).
 *
 * 컴파일: gcc -O2 -march=native -DNV=5 -o enumnpn5 scripts/enumerate_npn_stabilizers.c
 * 실행:   ./enumnpn5 0 2 > npn_part0.txt &  ./enumnpn5 1 2 > npn_part1.txt     (범위를 둘로 나눠 두 코어에서)
 *         python3 -I scripts/check_npn_bruteforce.py 5 notes/npn_stabilizer_types_n5.tsv npn_part0.txt npn_part1.txt
 * n = 3, 4 는 -DNV=3, -DNV=4 로 컴파일해 ./enumnpn4 0 1 로 돌린다 (바로 끝난다).
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef NV
#define NV 5
#endif
#define NP (1 << NV) /* 입력점의 수 = 진리표의 비트 수 */
#define NB (NP / 8)  /* 진리표의 바이트 수 */
#define MAXPERM 120
#define NPAIRS2 (NV * (NV - 1) / 2)

typedef uint32_t tt_t;

static int nperm;
static int perm[MAXPERM][NV];
static int invp[MAXPERM][NV];      /* invp[p] = perm[p] 의 역순열 */
static tt_t tab[MAXPERM][NB][256]; /* tab[p][b][v]: f 의 b번째 바이트가 v 일 때 f o P_p 에 생기는 비트들 */
static tt_t lowmask[NV];           /* x_i = 0 인 점들 */
static int pairidx[NV][NV];
static int ncand[1 << NPAIRS2];
static unsigned char cand[1 << NPAIRS2][MAXPERM];

#define FULL (NP >= 32 ? (tt_t)0xFFFFFFFFu : (tt_t)(((uint64_t)1 << NP) - 1))
#define MAXELEM (2 * MAXPERM * NP)

typedef struct {
    int order, hbar_order, d, pattern, vorbits, fixv;
} gkey_t;
typedef struct {
    int antisym;
    gkey_t h, k;
    uint64_t count;
} rec_t;
static rec_t recs[16384];
static int nrecs;
static uint64_t trivial_count;

static void gen_perms(int *cur, int *used, int depth) {
    if (depth == NV) {
        memcpy(perm[nperm++], cur, sizeof(int) * NV);
        return;
    }
    for (int v = 0; v < NV; v++)
        if (!used[v]) {
            used[v] = 1;
            cur[depth] = v;
            gen_perms(cur, used, depth + 1);
            used[v] = 0;
        }
}

static void init(void) {
    int cur[NV], used[NV] = {0};
    gen_perms(cur, used, 0); /* 사전식 순서: perm[0] 이 항등 순열 */
    for (int i = 0; i < NV; i++) {
        lowmask[i] = 0;
        for (int x = 0; x < NP; x++)
            if (!((x >> i) & 1)) lowmask[i] |= (tt_t)1 << x;
    }
    int k = 0;
    for (int i = 0; i < NV; i++)
        for (int j = i + 1; j < NV; j++) pairidx[i][j] = pairidx[j][i] = k++;
    for (int p = 0; p < nperm; p++) {
        int pinv[NP];
        for (int i = 0; i < NV; i++) invp[p][perm[p][i]] = i;
        for (int x = 0; x < NP; x++) { /* P_p(x) 의 i번째 좌표 = x 의 perm[p][i]번째 좌표 */
            int y = 0;
            for (int i = 0; i < NV; i++) y |= ((x >> perm[p][i]) & 1) << i;
            pinv[y] = x;
        }
        /* (f o P_p)(x) = f(P_p(x)) 이므로, f 의 비트 y 는 f o P_p 의 비트 pinv[y] 로 간다 */
        for (int b = 0; b < NB; b++)
            for (int v = 0; v < 256; v++) {
                tt_t r = 0;
                for (int j = 0; j < 8; j++)
                    if ((v >> j) & 1) r |= (tt_t)1 << pinv[8 * b + j];
                tab[p][b][v] = r;
            }
    }
    for (int mask = 0; mask < (1 << NPAIRS2); mask++) {
        ncand[mask] = 0;
        for (int p = 1; p < nperm; p++) {
            int ok = 1;
            for (int i = 0; i < NV; i++)
                if (perm[p][i] != i && !((mask >> pairidx[i][perm[p][i]]) & 1)) ok = 0;
            if (ok) cand[mask][ncand[mask]++] = (unsigned char)p;
        }
    }
}

static inline tt_t apply_perm(tt_t f, int p) {
    tt_t g = 0;
    for (int b = 0; b < NB; b++) g |= tab[p][b][(f >> (8 * b)) & 255];
    return g;
}

/* 결과의 비트 x = g 의 비트 (x xor c) */
static inline tt_t translate(tt_t g, int c) {
    for (int i = 0; i < NV; i++)
        if ((c >> i) & 1) g = ((g >> (1 << i)) & lowmask[i]) | ((g & lowmask[i]) << (1 << i));
    return g;
}

/* 어떤 c 에 대해 g(x xor c) = f(x) 가 모든 x 에서 성립하는가 */
static int is_translate(tt_t g, tt_t f) {
    if (f == 0) return g == 0;
    int a = __builtin_ctz(f);
    tt_t s = g;
    while (s) {
        int b = __builtin_ctz(s);
        s &= s - 1;
        if (translate(g, a ^ b) == f) return 1;
    }
    return 0;
}

/* g(x xor c) = f(x) 가 모든 x 에서 성립하는 c 를 모두 cs 에 모으고 개수를 돌려준다. */
static int all_translates(tt_t g, tt_t f, int *cs) {
    int n = 0;
    if (f == 0 || f == FULL) {
        if (g != f) return 0;
        for (int c = 0; c < NP; c++) cs[n++] = c;
        return n;
    }
    int a = __builtin_ctz(f);
    tt_t s = g;
    while (s) {
        int b = __builtin_ctz(s);
        s &= s - 1;
        if (translate(g, a ^ b) == f) cs[n++] = a ^ b;
    }
    return n;
}

static void compute_d1(tt_t f, int *d1) {
    for (int i = 0; i < NV; i++) d1[i] = __builtin_popcount(f & (f >> (1 << i)) & lowmask[i]);
}

static void compute_d2(tt_t f, int *d2) {
    for (int i = 0; i < NV; i++)
        for (int j = i + 1; j < NV; j++) {
            tt_t m = lowmask[i] & lowmask[j];
            d2[pairidx[i][j]] = __builtin_popcount(f & (f >> ((1 << i) + (1 << j))) & m) +
                                __builtin_popcount((f >> (1 << i)) & (f >> (1 << j)) & m);
        }
}

static int find(int *par, int a) {
    while (par[a] != a) a = par[a];
    return a;
}

/* 꼭짓점 x 가 (pi, c) 에서 가는 곳: pi(x xor c). pi(y) 의 i번째 좌표 = y 의 perm[p][i]번째 좌표 */
static inline int vertex_image(int p, int c, int x) {
    int y = x ^ c, r = 0;
    for (int i = 0; i < NV; i++) r |= ((y >> perm[p][i]) & 1) << i;
    return r;
}

/* 원소 (ep[k], ec[k]) k < ne 들이 이루는 군의 종류 */
static void compute_gkey(const int *ep, const int *ec, int ne, gkey_t *out) {
    int par[NP];
    int seenp[MAXPERM] = {0};
    int nperms = 0, ntrans = 0;
    for (int i = 0; i < NV; i++) par[i] = i;
    for (int k = 0; k < ne; k++) {
        if (!seenp[ep[k]]) { seenp[ep[k]] = 1; nperms++; }
        if (ep[k] == 0) ntrans++;
        for (int i = 0; i < NV; i++) {
            int a = find(par, i), b = find(par, perm[ep[k]][i]);
            if (a != b) par[a] = b;
        }
    }
    int len[NV] = {0}, sorted[NV], ns = 0;
    for (int i = 0; i < NV; i++) len[find(par, i)]++;
    for (int i = 0; i < NV; i++)
        if (len[i]) sorted[ns++] = len[i];
    for (int i = 0; i < ns; i++)
        for (int j = i + 1; j < ns; j++)
            if (sorted[j] < sorted[i]) {
                int t = sorted[i];
                sorted[i] = sorted[j];
                sorted[j] = t;
            }
    int pattern = 0;
    for (int i = 0; i < ns; i++) pattern = pattern * 10 + sorted[i];
    for (int x = 0; x < NP; x++) par[x] = x;
    int fixed[NP];
    for (int x = 0; x < NP; x++) fixed[x] = 1;
    for (int k = 0; k < ne; k++)
        for (int x = 0; x < NP; x++) {
            int y = vertex_image(ep[k], ec[k], x);
            if (y != x) fixed[x] = 0;
            int a = find(par, x), b = find(par, y);
            if (a != b) par[a] = b;
        }
    int vorbits = 0, fixv = 0;
    for (int x = 0; x < NP; x++) {
        if (find(par, x) == x) vorbits++;
        fixv += fixed[x];
    }
    if (ntrans & (ntrans - 1)) {
        fprintf(stderr, "TRANSLATION COUNT NOT A POWER OF TWO\n");
        exit(1);
    }
    out->order = ne;
    out->hbar_order = nperms;
    out->d = __builtin_ctz(ntrans);
    out->pattern = pattern;
    out->vorbits = vorbits;
    out->fixv = fixv;
}

static int same_key(const gkey_t *a, const gkey_t *b) {
    return a->order == b->order && a->hbar_order == b->hbar_order && a->d == b->d && a->pattern == b->pattern &&
           a->vorbits == b->vorbits && a->fixv == b->fixv;
}

static void add_rec(int antisym, const gkey_t *h, const gkey_t *k) {
    for (int i = 0; i < nrecs; i++)
        if (recs[i].antisym == antisym && same_key(&recs[i].h, h) && same_key(&recs[i].k, k)) {
            recs[i].count++;
            return;
        }
    rec_t *r = &recs[nrecs++];
    r->antisym = antisym;
    r->h = *h;
    r->k = *k;
    r->count = 1;
}

/* f 의 반대칭 (p, c) 를 모두 찾는다. use_filter 가 0 이면 거르지 않고 모든 순열을 시험한다. 개수를 돌려준다. */
static int find_anti(tt_t f, unsigned mask, int use_filter, int *ap, int *ac) {
    tt_t nf = ~f & FULL;
    int na = 0, cs[NP];
    int d2f[NPAIRS2], d2nf[NPAIRS2];
    if (use_filter) {
        compute_d2(f, d2f);
        compute_d2(nf, d2nf);
    }
    int count = use_filter ? ncand[mask] + 1 : nperm;
    for (int k = 0; k < count; k++) {
        int p = use_filter ? (k == 0 ? 0 : cand[mask][k - 1]) : k;
        if (use_filter && p != 0) {
            int ok = 1;
            for (int i = 0; i < NV && ok; i++)
                for (int j = i + 1; j < NV; j++)
                    if (d2f[pairidx[invp[p][i]][invp[p][j]]] != d2nf[pairidx[i][j]]) {
                        ok = 0;
                        break;
                    }
            if (!ok) continue;
        }
        int nc = all_translates(apply_perm(f, p), nf, cs);
        for (int w = 0; w < nc; w++) {
            ap[na] = p;
            ac[na++] = cs[w];
        }
    }
    return na;
}

int main(int argc, char **argv) {
    if (argc != 3) {
        fprintf(stderr, "usage: %s part parts\n", argv[0]);
        return 2;
    }
    uint64_t part = strtoull(argv[1], NULL, 10), parts = strtoull(argv[2], NULL, 10);
    init();
    uint64_t half = (uint64_t)1 << (NP - 1);
    uint64_t lo = half / parts * part, hi = (part + 1 == parts) ? half : half / parts * (part + 1);
    const uint64_t stride = (NV >= 5) ? 4099 : 1;
    uint64_t checked = 0, checked_anti = 0, balanced = 0, nanti_funcs = 0;
    static int ep[MAXELEM], ec[MAXELEM], ap[MAXELEM], ac[MAXELEM], ap2[MAXELEM], ac2[MAXELEM], kp[MAXELEM], kc[MAXELEM];

    for (uint64_t ff = lo; ff < hi; ff++) {
        tt_t f = (tt_t)ff;
        int d1[NV];
        unsigned mask = 0;
        compute_d1(f, d1);
        for (int i = 0; i < NV; i++)
            for (int j = i + 1; j < NV; j++)
                if (d1[i] == d1[j]) mask |= 1u << pairidx[i][j];

        unsigned char hbar[MAXPERM];
        int nh = 0;
        if (ncand[mask]) {
            int d2[NPAIRS2];
            compute_d2(f, d2);
            for (int k = 0; k < ncand[mask]; k++) {
                int p = cand[mask][k], ok = 1;
                for (int i = 0; i < NV && ok; i++)
                    for (int j = i + 1; j < NV; j++)
                        if (d2[pairidx[perm[p][i]][perm[p][j]]] != d2[pairidx[i][j]]) {
                            ok = 0;
                            break;
                        }
                if (!ok) continue;
                tt_t g = apply_perm(f, p);
                if (is_translate(g, f)) hbar[nh++] = (unsigned char)p;
            }
        }

        if (ff % stride == 0) { /* 대칭: 거르지 않고 다시 구해 비교 */
            int nh2 = 0, same = 1;
            unsigned char hbar2[MAXPERM] = {0};
            for (int p = 1; p < nperm; p++) {
                tt_t g = apply_perm(f, p);
                if (is_translate(g, f)) hbar2[nh2++] = (unsigned char)p;
            }
            if (nh2 != nh) same = 0;
            for (int k = 0; same && k < nh; k++)
                if (hbar[k] != hbar2[k]) same = 0;
            if (!same) {
                fprintf(stderr, "SELF-CHECK (SYMMETRY) FAILED at f = %llu\n", (unsigned long long)ff);
                return 1;
            }
            checked++;
        }

        /* 평행이동 부분 W: f(x xor c) = f(x). 무게가 홀수이면 c = 0 뿐이다 */
        int wcs[NP], nW = 1;
        wcs[0] = 0;
        if (!(__builtin_popcount(f) & 1)) nW = all_translates(f, f, wcs);

        /* 반대칭: 무게가 2^(n-1) 인 함수만 */
        int na = 0;
        if (__builtin_popcount(f) == NP / 2) {
            balanced++;
            na = find_anti(f, mask, 1, ap, ac);
            if (ff % stride == 0) {
                int na2 = find_anti(f, mask, 0, ap2, ac2), same = (na == na2);
                for (int k = 0; same && k < na; k++)
                    if (ap[k] != ap2[k] || ac[k] != ac2[k]) same = 0;
                if (!same) {
                    fprintf(stderr, "SELF-CHECK (ANTISYMMETRY) FAILED at f = %llu\n", (unsigned long long)ff);
                    return 1;
                }
                checked_anti++;
            }
        }

        if (nh == 0 && nW == 1 && na == 0) {
            trivial_count++;
            continue;
        }
        int ne = 0;
        for (int w = 0; w < nW; w++) { ep[ne] = 0; ec[ne++] = wcs[w]; }
        for (int k = 0; k < nh; k++) {
            int p = hbar[k], cs[NP];
            int nc = all_translates(apply_perm(f, p), f, cs);
            if (nc != nW) { fprintf(stderr, "COSET SIZE MISMATCH at f = %llu\n", (unsigned long long)ff); return 1; }
            for (int w = 0; w < nc; w++) { ep[ne] = p; ec[ne++] = cs[w]; }
        }
        gkey_t kh, kk;
        compute_gkey(ep, ec, ne, &kh);
        if (na) {
            if (na != ne) { fprintf(stderr, "ANTISYMMETRY COSET SIZE MISMATCH at f = %llu: |H| %d |A| %d\n", (unsigned long long)ff, ne, na); return 1; }
            memcpy(kp, ep, sizeof(int) * ne);
            memcpy(kc, ec, sizeof(int) * ne);
            for (int k = 0; k < na; k++) { kp[ne + k] = ap[k]; kc[ne + k] = ac[k]; }
            compute_gkey(kp, kc, ne + na, &kk);
            nanti_funcs++;
        } else {
            kk = kh;
        }
        add_rec(na > 0, &kh, &kk);
    }

    printf("# NV %d part %llu of %llu range %llu %llu selfchecked %llu selfchecked_anti %llu balanced %llu\n", NV,
           (unsigned long long)part, (unsigned long long)parts, (unsigned long long)lo, (unsigned long long)hi,
           (unsigned long long)checked, (unsigned long long)checked_anti, (unsigned long long)balanced);
    printf("# npn antisym | H: order hbar_order d orbit_lengths vertex_orbits fixed_vertices | Kbar: same six | count   (count is for the half with top bit 0)\n");
    int allones = 0;
    for (int i = 0; i < NV; i++) allones = allones * 10 + 1;
    printf("npn 0 1 1 0 %d %d %d 1 1 0 %d %d %d %llu\n", allones, NP, NP, allones, NP, NP, (unsigned long long)trivial_count);
    for (int i = 0; i < nrecs; i++)
        printf("npn %d %d %d %d %d %d %d %d %d %d %d %d %d %llu\n", recs[i].antisym, recs[i].h.order, recs[i].h.hbar_order,
               recs[i].h.d, recs[i].h.pattern, recs[i].h.vorbits, recs[i].h.fixv, recs[i].k.order, recs[i].k.hbar_order,
               recs[i].k.d, recs[i].k.pattern, recs[i].k.vorbits, recs[i].k.fixv, (unsigned long long)recs[i].count);
    return 0;
}
