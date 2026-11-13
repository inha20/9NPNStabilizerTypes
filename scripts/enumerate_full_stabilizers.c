/*
 * enumerate_full_stabilizers.c
 * ============================
 * n변수 부울 함수 2^(2^n)개를 하나도 빼지 않고 훑어서 (n = 3, 4, 5), 함수마다 안정화군 전체
 *
 *   H(f) = { (pi, c) : f(pi(x xor c)) = f(x) 가 모든 x 에서 성립 }      (입력 반전 c 와 변수 순열 pi)
 *
 * 를 직접 구하고, 군의 종류별 함수 수를 센다. 군의 종류는 다음 여섯 값의 묶음으로 적는다.
 *   위수 |H|, 순열 부분 Hbar 의 위수, 평행이동 부분 W 의 차원 d (|W| = 2^d, f 의 주기들),
 *   Hbar 의 변수 위 궤도 길이, H 가 꼭짓점 2^n 개에 만드는 궤도 수, H 의 고정 꼭짓점 수.
 * (고정 꼭짓점이 있으면 H 는 S_n 의 부분군과 켤레이다. 없으면 "꼬인" 대칭이다.)
 *
 * GAP 의 표지 표 계산(stabilizer_types.g)과 독립인 검증이다. 여기서는 "이 (pi, c) 가 이 함수의 대칭인가"를
 * 함수마다 직접 시험할 뿐이다. enumerate_all_functions.c 에서 순열 부분만 구하던 것을 전체 군으로 넓혔다.
 *
 * 빠르게 하려고, 대칭이 될 수 없는 순열을 두 가지 필요조건으로 먼저 걸러 낸다(enumerate_all_functions.c 와 같다).
 *   d1[i]    = { x, x xor e_i } 가 둘 다 f = 1 인 쌍의 수
 *   d2[i][j] = { x, x xor e_i xor e_j } 가 둘 다 f = 1 인 쌍의 수
 * 걸러지지 않은 순열만 정확히 시험한다. 거르는 단계가 대칭을 놓치지 않는지는 일정 간격으로 뽑은 함수에서
 * 모든 순열을 거르지 않고 시험한 결과와 비교해 확인한다 (다르면 프로그램이 멈춘다).
 *
 * f 와 not f 는 안정화군이 같으므로 최상위 비트가 0 인 절반만 훑는다 (합칠 때 2배 한다).
 *
 * 컴파일: gcc -O2 -march=native -DNV=5 -o enumfull5 scripts/enumerate_full_stabilizers.c
 * 실행:   ./enumfull5 0 2 > full_part0.txt &  ./enumfull5 1 2 > full_part1.txt     (범위를 둘로 나눠 두 코어에서)
 *         python scripts/compare_full_enumeration.py notes/stabilizer_types_n5.tsv full_part0.txt full_part1.txt
 * n = 3, 4 는 -DNV=3, -DNV=4 로 컴파일해 ./enumfull4 0 1 로 돌린다 (바로 끝난다).
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
#define MAXPAIRING 15
#define NPAIRS2 (NV * (NV - 1) / 2)

typedef uint32_t tt_t;

static int nperm;
static int perm[MAXPERM][NV];
static tt_t tab[MAXPERM][NB][256]; /* tab[p][b][v]: f 의 b번째 바이트가 v 일 때 f o P_p 에 생기는 비트들 */
static tt_t lowmask[NV];           /* x_i = 0 인 점들 */
static int npairing;
static int mate[MAXPAIRING][NV];   /* mate[m][i] = 분할 m 에서 i 의 짝 (외톨이는 자기 자신) */
static int pairing_action[MAXPERM][MAXPAIRING];
static int pairidx[NV][NV];
static int ncand[1 << NPAIRS2];
static unsigned char cand[1 << NPAIRS2][MAXPERM];

#define FULL (NP >= 32 ? (tt_t)0xFFFFFFFFu : (tt_t)(((uint64_t)1 << NP) - 1))
#define MAXELEM (2 * MAXPERM * NP)

typedef struct {
    int order, hbar_order, d, pattern, vorbits, fixv, pairorbits;
    uint64_t count;
} rec_t;
static rec_t rec_full[4096];
static int nrec_full;
static uint64_t trivial_full;

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

static void gen_pairings(int *m, int single_used) {
    int i = 0;
    while (i < NV && m[i] >= 0) i++;
    if (i == NV) {
        memcpy(mate[npairing++], m, sizeof(int) * NV);
        return;
    }
    if ((NV & 1) && !single_used) {
        m[i] = i;
        gen_pairings(m, 1);
        m[i] = -1;
    }
    for (int j = i + 1; j < NV; j++)
        if (m[j] < 0) {
            m[i] = j;
            m[j] = i;
            gen_pairings(m, single_used);
            m[i] = m[j] = -1;
        }
}

static void init(void) {
    int cur[NV], used[NV] = {0}, m[NV];
    gen_perms(cur, used, 0); /* 사전식 순서: perm[0] 이 항등 순열 */
    for (int i = 0; i < NV; i++) {
        m[i] = -1;
        lowmask[i] = 0;
        for (int x = 0; x < NP; x++)
            if (!((x >> i) & 1)) lowmask[i] |= (tt_t)1 << x;
    }
    gen_pairings(m, 0);
    int k = 0;
    for (int i = 0; i < NV; i++)
        for (int j = i + 1; j < NV; j++) pairidx[i][j] = pairidx[j][i] = k++;
    for (int p = 0; p < nperm; p++) {
        int pinv[NP];
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
        for (int a = 0; a < npairing; a++) {
            int img[NV];
            for (int i = 0; i < NV; i++) img[perm[p][i]] = perm[p][mate[a][i]];
            int found = -1;
            for (int c = 0; c < npairing; c++)
                if (memcmp(img, mate[c], sizeof(int) * NV) == 0) found = c;
            pairing_action[p][a] = found;
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

static int find(int *par, int a) {
    while (par[a] != a) a = par[a];
    return a;
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

/* 꼭짓점 x 가 대칭 (pi, c) 에서 가는 곳: pi(x xor c). pi(y) 의 i번째 좌표 = y 의 perm[p][i]번째 좌표 */
static inline int vertex_image(int p, int c, int x) {
    int y = x ^ c, r = 0;
    for (int i = 0; i < NV; i++) r |= ((y >> perm[p][i]) & 1) << i;
    return r;
}

static void classify_full(const int *ep, const int *ec, int ne, int hbar_order, int nW) {
    int par[NP > MAXPAIRING ? NP : MAXPAIRING];
    for (int i = 0; i < NV; i++) par[i] = i;
    for (int k = 0; k < ne; k++)
        for (int i = 0; i < NV; i++) {
            int a = find(par, i), b = find(par, perm[ep[k]][i]);
            if (a != b) par[a] = b;
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
    /* 쌍 분할 위의 궤도 수 */
    for (int i = 0; i < npairing; i++) par[i] = i;
    for (int k = 0; k < ne; k++)
        for (int i = 0; i < npairing; i++) {
            int a = find(par, i), b = find(par, pairing_action[ep[k]][i]);
            if (a != b) par[a] = b;
        }
    int pairorbits = 0;
    for (int i = 0; i < npairing; i++)
        if (find(par, i) == i) pairorbits++;
    /* 꼭짓점 위의 궤도 수와 고정 꼭짓점 수 */
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
    int d = __builtin_ctz(nW);
    for (int k = 0; k < nrec_full; k++) {
        rec_t *r = &rec_full[k];
        if (r->order == ne && r->hbar_order == hbar_order && r->d == d && r->pattern == pattern &&
            r->vorbits == vorbits && r->fixv == fixv && r->pairorbits == pairorbits) {
            r->count++;
            return;
        }
    }
    rec_t *r = &rec_full[nrec_full++];
    r->order = ne; r->hbar_order = hbar_order; r->d = d; r->pattern = pattern;
    r->vorbits = vorbits; r->fixv = fixv; r->pairorbits = pairorbits; r->count = 1;
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
    uint64_t checked = 0;

    for (uint64_t ff = lo; ff < hi; ff++) {
        tt_t f = (tt_t)ff;
        int d1[NV];
        unsigned mask = 0;
        for (int i = 0; i < NV; i++) d1[i] = __builtin_popcount(f & (f >> (1 << i)) & lowmask[i]);
        for (int i = 0; i < NV; i++)
            for (int j = i + 1; j < NV; j++)
                if (d1[i] == d1[j]) mask |= 1u << pairidx[i][j];

        unsigned char hbar[MAXPERM];
        int nh = 0;
        if (ncand[mask]) {
            int d2[NPAIRS2];
            for (int i = 0; i < NV; i++)
                for (int j = i + 1; j < NV; j++) {
                    tt_t m = lowmask[i] & lowmask[j];
                    d2[pairidx[i][j]] = __builtin_popcount(f & (f >> ((1 << i) + (1 << j))) & m) +
                                        __builtin_popcount((f >> (1 << i)) & (f >> (1 << j)) & m);
                }
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

        if (ff % stride == 0) { /* 거르지 않고 다시 구해 비교 */
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
                fprintf(stderr, "SELF-CHECK FAILED at f = %llu\n", (unsigned long long)ff);
                return 1;
            }
            checked++;
        }

        /* 평행이동 부분 W: f(x xor c) = f(x). 무게가 홀수이면 c = 0 뿐이다 (c != 0 은 f = 1 인 점을 둘씩 짝짓는다) */
        int wcs[NP], nW = 1;
        wcs[0] = 0;
        if (!(__builtin_popcount(f) & 1)) nW = all_translates(f, f, wcs);
        if (nh == 0 && nW == 1) {
            trivial_full++;
        } else {
            static int ep[MAXELEM], ec[MAXELEM];
            int ne = 0;
            for (int w = 0; w < nW; w++) { ep[ne] = 0; ec[ne++] = wcs[w]; }
            for (int k = 0; k < nh; k++) {
                int p = hbar[k], cs[NP];
                int nc = all_translates(apply_perm(f, p), f, cs);
                if (nc != nW) { fprintf(stderr, "COSET SIZE MISMATCH at f = %llu\n", (unsigned long long)ff); return 1; }
                for (int w = 0; w < nc; w++) { ep[ne] = p; ec[ne++] = cs[w]; }
            }
            classify_full(ep, ec, ne, nh + 1, nW);
        }
    }

    printf("# NV %d part %llu of %llu range %llu %llu selfchecked %llu\n", NV, (unsigned long long)part,
           (unsigned long long)parts, (unsigned long long)lo, (unsigned long long)hi, (unsigned long long)checked);
    printf("# full order hbar_order d orbit_lengths vertex_orbits fixed_vertices pairing_orbits count   (count is for the half with top bit 0)\n");
    int allones = 0;
    for (int i = 0; i < NV; i++) allones = allones * 10 + 1;
    printf("full 1 1 0 %d %d %d %d %llu\n", allones, NP, NP, npairing, (unsigned long long)trivial_full);
    for (int k = 0; k < nrec_full; k++)
        printf("full %d %d %d %d %d %d %d %llu\n", rec_full[k].order, rec_full[k].hbar_order, rec_full[k].d,
               rec_full[k].pattern, rec_full[k].vorbits, rec_full[k].fixv, rec_full[k].pairorbits,
               (unsigned long long)rec_full[k].count);
    return 0;
}
