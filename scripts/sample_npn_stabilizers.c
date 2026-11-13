/*
 * sample_npn_stabilizers.c
 * ========================
 * 변수 NV 개(기본 6)의 부울 함수 가운데 균일하게 뽑은 함수들의 NPN 안정화군을 직접 구해 군의 종류별 개수를 센다.
 * enumerate_npn_stabilizers.c (n <= 5 의 전수 열거)와 같은 계산을 64비트 진리표로 옮기고 전수 열거 대신 무작위 표본을 쓴다.
 * 표지 표(GAP)와 풀이 계산과 독립인 점검이다. 표본이 작아서 드문 종류는 점검하지 못한다.
 *
 *   H(f)  = { (pi, c) : f(pi(x xor c)) = f(x) }       A(f) = { (pi, c) : f(pi(x xor c)) = not f(x) }       Kbar = H 합 A
 *   군의 종류는 H 와 Kbar 각각에 대해 (위수, 순열 부분의 위수, 평행이동 부분의 차원 d, 순열 부분의 변수 위 궤도 유형,
 *   꼭짓점 위 궤도 수, 고정 꼭짓점 수)이고 반대칭 유무를 더한다.
 *
 * 컴파일: gcc -O2 -march=native -DNV=6 -o sampnpn6 scripts/sample_npn_stabilizers.c      (NV=5 이면 2^32 개 함수, 참값과 견줄 수 있다)
 * 실행:   ./sampnpn6 <표본 수> <씨앗>  > out.txt          (두 코어에서 씨앗을 달리해 돌린다)
 *         ./sampnpn6 selftest <개수> <씨앗>               (대칭인 함수를 일부러 만들어 걸러내기를 모든 순열 시험과 견준다)
 *         ./sampnpn6 enumfix <작업 파일> [검산 간격]       (군 L 과 준동형 phi 마다 L 이 phi 로 고정하는 함수 2^o 개를 모두 열거해
 *                                                          안정화군이 정확히 L 인 것의 수를 센다. 작업 파일은 scripts/enum_fix_prep.g 가 만든다)
 *         ./sampnpn6 randfix <작업 파일> <표본 수> <씨앗>     (enumfix 와 같은 작업 파일에서 (L, phi) 마다 L 이 phi 로 고정하는 함수를 균일하게 뽑아
 *                                                          안정화군이 정확히 L 인 비율 E / 2^o 를 어림한다. o 가 커서 열거가 안 되는 군용)
 *         python3 -I scripts/compare_npn_sample.py <K 표 tsv> <n> out1.txt out2.txt
 *
 * 내부 점검: (1) 걸러내기(1차·2차 불변량)로 구한 H, A 가 모든 순열을 시험한 결과와 같은지를
 *   대칭이 있다고 나온 표본 전부와 200003 번째마다 한 표본에서 견준다. (2) selftest 는 한 원소로 고정되는(또는 그 원소가
 *   f 를 반전시키는) 함수를 만들어 같은 견줌을 한다. (3) |A| = |H| 와 평행이동 잉여류 크기를 확인한다.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifndef NV
#define NV 6
#endif
#define NP (1 << NV)
#define NB (NP / 8)
#define MAXPERM 720
#define NPAIRS2 (NV * (NV - 1) / 2)

typedef uint64_t tt_t;
#define FULL (~(tt_t)0 >> (64 - NP))
#define MAXELEM (2 * MAXPERM * NP)

static int nperm;
static int perm[MAXPERM][NV];
static int invp[MAXPERM][NV];
static tt_t tab[MAXPERM][NB][256];
static tt_t lowmask[NV];
static int pairidx[NV][NV];
static int ncand[1 << NPAIRS2];
static uint16_t cand[1 << NPAIRS2][MAXPERM];

typedef struct {
    int order, hbar_order, d, pattern, vorbits, fixv;
} gkey_t;
typedef struct {
    int antisym;
    gkey_t h, k;
    uint64_t count;
} rec_t;
static rec_t recs[65536];
static int nrecs;
static uint64_t trivial_count;

typedef struct {
    int nh;
    uint16_t hbar[MAXPERM];
    int nW;
    int wcs[NP];
    int na;
    int ap[MAXELEM], ac[MAXELEM];
} res_t;
static res_t RF, RB; /* 걸러내기로 구한 것, 모든 순열을 시험해 구한 것 */

static uint64_t rng_s[4];
static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static inline uint64_t rng_next(void) { /* xoshiro256** */
    uint64_t r = rotl(rng_s[1] * 5, 7) * 9, t = rng_s[1] << 17;
    rng_s[2] ^= rng_s[0];
    rng_s[3] ^= rng_s[1];
    rng_s[1] ^= rng_s[2];
    rng_s[0] ^= rng_s[3];
    rng_s[2] ^= t;
    rng_s[3] = rotl(rng_s[3], 45);
    return r;
}
static void rng_seed(uint64_t seed) { /* splitmix64 */
    for (int i = 0; i < 4; i++) {
        uint64_t z = (seed += 0x9e3779b97f4a7c15ULL);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        rng_s[i] = z ^ (z >> 31);
    }
}

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
    gen_perms(cur, used, 0);
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
        for (int x = 0; x < NP; x++) {
            int y = 0;
            for (int i = 0; i < NV; i++) y |= ((x >> perm[p][i]) & 1) << i;
            pinv[y] = x;
        }
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
            if (ok) cand[mask][ncand[mask]++] = (uint16_t)p;
        }
    }
}

static inline tt_t apply_perm(tt_t f, int p) {
    tt_t g = 0;
    for (int b = 0; b < NB; b++) g |= tab[p][b][(f >> (8 * b)) & 255];
    return g;
}

/* 진리표 g 의 평행이동: 결과(x) = g(x xor c) */
static inline tt_t translate(tt_t g, int c) {
    for (int i = 0; i < NV; i++)
        if ((c >> i) & 1) g = ((g >> (1 << i)) & lowmask[i]) | ((g & lowmask[i]) << (1 << i));
    return g;
}

static int is_translate(tt_t g, tt_t f) {
    if (f == 0) return g == 0;
    int a = __builtin_ctzll(f);
    tt_t s = g;
    while (s) {
        int b = __builtin_ctzll(s);
        s &= s - 1;
        if (translate(g, a ^ b) == f) return 1;
    }
    return 0;
}

/* translate(g, c) == f 인 c 를 모두 cs 에 쓰고 개수를 돌려준다 */
static int all_translates(tt_t g, tt_t f, int *cs) {
    int n = 0;
    if (f == 0 || f == FULL) {
        if (g != f) return 0;
        for (int c = 0; c < NP; c++) cs[n++] = c;
        return n;
    }
    int a = __builtin_ctzll(f);
    tt_t s = g;
    while (s) {
        int b = __builtin_ctzll(s);
        s &= s - 1;
        if (translate(g, a ^ b) == f) cs[n++] = a ^ b;
    }
    return n;
}

static void compute_d1(tt_t f, int *d1) {
    for (int i = 0; i < NV; i++) d1[i] = __builtin_popcountll(f & (f >> (1 << i)) & lowmask[i]);
}

static void compute_d2(tt_t f, int *d2) {
    for (int i = 0; i < NV; i++)
        for (int j = i + 1; j < NV; j++) {
            tt_t m = lowmask[i] & lowmask[j];
            d2[pairidx[i][j]] = __builtin_popcountll(f & (f >> ((1 << i) + (1 << j))) & m) +
                                __builtin_popcountll((f >> (1 << i)) & (f >> (1 << j)) & m);
        }
}

static int find(int *par, int a) {
    while (par[a] != a) a = par[a];
    return a;
}

/* 원소 (p, c) 가 꼭짓점 x 를 보내는 곳: 순열(x xor c) */
static inline int vertex_image(int p, int c, int x) {
    int y = x ^ c, r = 0;
    for (int i = 0; i < NV; i++) r |= ((y >> perm[p][i]) & 1) << i;
    return r;
}

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
            if (sorted[j] < sorted[i]) { int t = sorted[i]; sorted[i] = sorted[j]; sorted[j] = t; }
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
    if (ntrans & (ntrans - 1)) { fprintf(stderr, "TRANSLATION COUNT NOT A POWER OF TWO\n"); exit(1); }
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
        if (recs[i].antisym == antisym && same_key(&recs[i].h, h) && same_key(&recs[i].k, k)) { recs[i].count++; return; }
    if (nrecs >= 65536) { fprintf(stderr, "TOO MANY KEYS\n"); exit(1); }
    rec_t *r = &recs[nrecs++];
    r->antisym = antisym;
    r->h = *h;
    r->k = *k;
    r->count = 1;
}

/* f∘(p, c) = not f 인 (p, c) 를 모두 구한다. use_filter 이면 1·2차 불변량으로 순열을 거른다 */
static int find_anti(tt_t f, unsigned mask, int use_filter, int *ap, int *ac) {
    tt_t nf = ~f & FULL;
    int na = 0, cs[NP];
    int d2f[NPAIRS2], d2nf[NPAIRS2];
    if (use_filter) { compute_d2(f, d2f); compute_d2(nf, d2nf); }
    int count = use_filter ? ncand[mask] + 1 : nperm;
    for (int k = 0; k < count; k++) {
        int p = use_filter ? (k == 0 ? 0 : cand[mask][k - 1]) : k;
        if (use_filter && p != 0) {
            int ok = 1;
            for (int i = 0; i < NV && ok; i++)
                for (int j = i + 1; j < NV; j++)
                    if (d2f[pairidx[invp[p][i]][invp[p][j]]] != d2nf[pairidx[i][j]]) { ok = 0; break; }
            if (!ok) continue;
        }
        int nc = all_translates(apply_perm(f, p), nf, cs);
        for (int w = 0; w < nc; w++) { ap[na] = p; ac[na++] = cs[w]; }
    }
    return na;
}

/* f 의 H 와 A 를 구한다 (순열 부분 목록, 평행이동 부분, 반대칭 원소) */
static void analyze(tt_t f, int use_filter, res_t *r) {
    int d1[NV];
    unsigned mask = 0;
    compute_d1(f, d1);
    for (int i = 0; i < NV; i++)
        for (int j = i + 1; j < NV; j++)
            if (d1[i] == d1[j]) mask |= 1u << pairidx[i][j];
    r->nh = 0;
    if (use_filter) {
        if (ncand[mask]) {
            int d2[NPAIRS2];
            compute_d2(f, d2);
            for (int k = 0; k < ncand[mask]; k++) {
                int p = cand[mask][k], ok = 1;
                for (int i = 0; i < NV && ok; i++)
                    for (int j = i + 1; j < NV; j++)
                        if (d2[pairidx[perm[p][i]][perm[p][j]]] != d2[pairidx[i][j]]) { ok = 0; break; }
                if (!ok) continue;
                if (is_translate(apply_perm(f, p), f)) r->hbar[r->nh++] = (uint16_t)p;
            }
        }
    } else {
        for (int p = 1; p < nperm; p++)
            if (is_translate(apply_perm(f, p), f)) r->hbar[r->nh++] = (uint16_t)p;
    }
    r->nW = 1;
    r->wcs[0] = 0;
    if (!(__builtin_popcountll(f) & 1)) r->nW = all_translates(f, f, r->wcs);
    r->na = 0;
    if (__builtin_popcountll(f) == NP / 2) r->na = find_anti(f, mask, use_filter, r->ap, r->ac);
}

static int same_res(const res_t *a, const res_t *b) {
    if (a->nh != b->nh || a->nW != b->nW || a->na != b->na) return 0;
    if (memcmp(a->hbar, b->hbar, sizeof(uint16_t) * a->nh)) return 0;
    if (memcmp(a->wcs, b->wcs, sizeof(int) * a->nW)) return 0;
    if (memcmp(a->ap, b->ap, sizeof(int) * a->na)) return 0;
    if (memcmp(a->ac, b->ac, sizeof(int) * a->na)) return 0;
    return 1;
}

static int is_positive(const res_t *r) { return r->nh > 0 || r->nW > 1 || r->na > 0; }

/* H 의 원소 목록 (perm 번호, 평행이동) 을 만든다. 원소 수를 돌려준다 */
static int build_h(tt_t f, const res_t *r, int *ep, int *ec) {
    int ne = 0;
    for (int w = 0; w < r->nW; w++) { ep[ne] = 0; ec[ne++] = r->wcs[w]; }
    for (int k = 0; k < r->nh; k++) {
        int p = r->hbar[k], cs[NP];
        int nc = all_translates(apply_perm(f, p), f, cs);
        if (nc != r->nW) { fprintf(stderr, "COSET SIZE MISMATCH\n"); exit(1); }
        for (int w = 0; w < nc; w++) { ep[ne] = p; ec[ne++] = cs[w]; }
    }
    return ne;
}

static int selftest(int argc, char **argv) {
    if (argc != 4) { fprintf(stderr, "usage: %s selftest count seed\n", argv[0]); return 2; }
    uint64_t cnt = strtoull(argv[2], NULL, 10);
    rng_seed(strtoull(argv[3], NULL, 10));
    init();
    uint64_t made_h = 0, made_a = 0, pos_h = 0, pos_a = 0, bad = 0, noanti = 0;
    for (uint64_t it = 0; it < cnt; it++) {
        int p = (int)(rng_next() % nperm), c = (int)(rng_next() % NP);
        int img[NP];
        for (int x = 0; x < NP; x++) img[x] = vertex_image(p, c, x);
        tt_t f = 0;
        int mode = (int)(it & 1); /* 0: g 로 고정, 1: g 가 반전시킴 */
        if (mode == 1) {
            int alleven = 1;
            for (int x = 0; x < NP && alleven; x++) {
                int len = 0, y = x;
                do { y = img[y]; len++; } while (y != x);
                if (len & 1) alleven = 0;
            }
            if (!alleven) { noanti++; continue; }
            int seen[NP] = {0};
            for (int x = 0; x < NP; x++) {
                if (seen[x]) continue;
                int v = (int)(rng_next() & 1), y = x;
                do { seen[y] = 1; if (v) f |= (tt_t)1 << y; v ^= 1; y = img[y]; } while (y != x);
            }
            made_a++;
        } else {
            int seen[NP] = {0};
            for (int x = 0; x < NP; x++) {
                if (seen[x]) continue;
                int v = (int)(rng_next() & 1), y = x;
                do { seen[y] = 1; if (v) f |= (tt_t)1 << y; y = img[y]; } while (y != x);
            }
            made_h++;
        }
        analyze(f, 1, &RF);
        analyze(f, 0, &RB);
        if (!same_res(&RF, &RB)) {
            bad++;
            if (bad <= 5) fprintf(stderr, "SELFTEST MISMATCH mode %d f=%016llx p=%d c=%d: filter nh %d nW %d na %d, brute nh %d nW %d na %d\n", mode,
                                  (unsigned long long)f, p, c, RF.nh, RF.nW, RF.na, RB.nh, RB.nW, RB.na);
        }
        if (mode == 1 && RB.na == 0) { bad++; fprintf(stderr, "SELFTEST: antisymmetric element not found\n"); }
        if (mode == 0 && !(RB.nh > 0 || RB.nW > 1 || (p == 0 && c == 0))) { bad++; fprintf(stderr, "SELFTEST: symmetric element not found\n"); }
        if (RB.nh > 0 || RB.nW > 1) pos_h++;
        if (RB.na > 0) pos_a++;
    }
    printf("selftest NV=%d: H-type %llu (nontrivial H %llu), A-type %llu (antisymmetric found %llu, skipped %llu), mismatches %llu\n", NV,
           (unsigned long long)made_h, (unsigned long long)pos_h, (unsigned long long)made_a, (unsigned long long)pos_a,
           (unsigned long long)noanti, (unsigned long long)bad);
    return bad ? 1 : 0;
}

/* 작업 파일의 한 줄에서 정수를 모두 읽는다 (GAP 이 찍은 [ 1, 2, ... ] 꼴을 견딘다) */
static int parse_ints(const char *line, int *out, int maxn) {
    int n = 0;
    const char *p = line;
    while (*p) {
        if ((*p >= '0' && *p <= '9') || (*p == '-' && p[1] >= '0' && p[1] <= '9')) {
            char *e;
            long v = strtol(p, &e, 10);
            if (n < maxn) out[n++] = (int)v;
            p = e;
        } else p++;
    }
    return n;
}

/* 작업 파일 형식 (줄마다):
 *   J <j> <|L|> <|N_B(L)|> <꼭짓점 궤도 수 o> <phi 개수>
 *   F <번호> <서명 4개>
 *   V <꼭짓점 0..NP-1 의 궤도 번호 (0 부터)>
 *   S <꼭짓점 0..NP-1 의 부호 비트>
 * 부호가 s(x)이고 궤도 번호가 o(x)일 때 f(x) = b[o(x)] xor s(x) (b 는 o 비트)가 L 이 phi 로 고정하는 모든 함수이다. */
static int enumfix(const char *path, uint64_t chk_stride, uint64_t nsamp) {
    FILE *fp = fopen(path, "r");
    if (!fp) { perror("fopen"); return 2; }
    char *line = malloc(1 << 20);
    int vals[2 * NP + 16];
    int curj = 0, lsize = 0, nrm = 0, onum = 0, nphi = 0, idx = 0, sig[4] = {0};
    int orb[NP], have_v = 0;
    while (fgets(line, 1 << 20, fp)) {
        if (line[0] == 'J') {
            int n = parse_ints(line + 1, vals, 16);
            if (n < 5) { fprintf(stderr, "bad J line\n"); return 2; }
            curj = vals[0]; lsize = vals[1]; nrm = vals[2]; onum = vals[3]; nphi = vals[4];
        } else if (line[0] == 'F') {
            int n = parse_ints(line + 1, vals, 16);
            if (n < 5) { fprintf(stderr, "bad F line\n"); return 2; }
            idx = vals[0];
            for (int k = 0; k < 4; k++) sig[k] = vals[1 + k];
            have_v = 0;
        } else if (line[0] == 'V') {
            if (parse_ints(line + 1, orb, NP) != NP) { fprintf(stderr, "bad V line\n"); return 2; }
            have_v = 1;
        } else if (line[0] == 'S' && have_v) {
            int sgn[NP];
            if (parse_ints(line + 1, sgn, NP) != NP) { fprintf(stderr, "bad S line\n"); return 2; }
            tt_t masks[NP] = {0}, S = 0;
            for (int x = 0; x < NP; x++) {
                if (orb[x] < 0 || orb[x] >= onum) { fprintf(stderr, "orbit number out of range\n"); return 2; }
                masks[orb[x]] |= (tt_t)1 << x;
                if (sgn[x]) S |= (tt_t)1 << x;
            }
            if (nsamp > 0) {   /* 무작위 모드: 열거하지 않고 nsamp 개를 균일하게 뽑는다 */
                uint64_t exact = 0;
                for (uint64_t it = 0; it < nsamp; it++) {
                    tt_t f = S;
                    uint64_t b = rng_next();
                    for (int i = 0; i < onum; i++)
                        if ((b >> i) & 1) f ^= masks[i];
                    analyze(f, 1, &RF);
                    int size = RF.nW * (1 + RF.nh) + RF.na;
                    if (size < lsize) { fprintf(stderr, "STABILIZER SMALLER THAN L: j %d phi %d f=%016llx size %d |L| %d\n", curj, idx, (unsigned long long)f, size, lsize); return 1; }
                    if (size == lsize) exact++;
                }
                printf("Q %d %d %d %d %d %d %d %llu %llu %d %d\n", curj, idx, sig[0], sig[1], sig[2], sig[3], onum, (unsigned long long)exact,
                       (unsigned long long)nsamp, lsize, nrm);
                fflush(stdout);
                have_v = 0;
                continue;
            }
            if (onum > 40) { printf("R %d %d %d %d %d %d %d %d SKIP\n", curj, idx, sig[0], sig[1], sig[2], sig[3], onum, 0); fflush(stdout); continue; }
            uint64_t total = (uint64_t)1 << onum, exact = 0, nchk = 0;
            tt_t f = S;
            for (uint64_t i = 0; i < total; i++) {
                analyze(f, 1, &RF);
                int size = RF.nW * (1 + RF.nh) + RF.na;
                if (size < lsize) { fprintf(stderr, "STABILIZER SMALLER THAN L: j %d phi %d f=%016llx size %d |L| %d\n", curj, idx, (unsigned long long)f, size, lsize); return 1; }
                if (size == lsize) exact++;
                if (chk_stride && i % chk_stride == 0) {
                    analyze(f, 0, &RB);
                    if (!same_res(&RF, &RB)) { fprintf(stderr, "SELF-CHECK FAILED j %d phi %d f=%016llx\n", curj, idx, (unsigned long long)f); return 1; }
                    nchk++;
                }
                if (i + 1 < total) f ^= masks[__builtin_ctzll(i + 1)];
            }
            printf("R %d %d %d %d %d %d %d %llu %llu %d %d %llu\n", curj, idx, sig[0], sig[1], sig[2], sig[3], onum, (unsigned long long)exact,
                   (unsigned long long)total, lsize, nrm, (unsigned long long)nchk);
            fflush(stdout);
            have_v = 0;
        }
    }
    (void)nphi;
    free(line);
    return 0;
}

int main(int argc, char **argv) {
    if (argc >= 2 && !strcmp(argv[1], "selftest")) return selftest(argc, argv);
    if (argc >= 3 && !strcmp(argv[1], "enumfix")) { init(); return enumfix(argv[2], argc >= 4 ? strtoull(argv[3], NULL, 10) : 0, 0); }
    if (argc >= 5 && !strcmp(argv[1], "randfix")) { init(); rng_seed(strtoull(argv[4], NULL, 10)); return enumfix(argv[2], 0, strtoull(argv[3], NULL, 10)); }
    if (argc != 3) { fprintf(stderr, "usage: %s samples seed | selftest count seed\n", argv[0]); return 2; }
    uint64_t nsamp = strtoull(argv[1], NULL, 10), seed = strtoull(argv[2], NULL, 10);
    init();
    rng_seed(seed);
    const uint64_t stride = 200003;
    uint64_t checked_stride = 0, checked_pos = 0, balanced = 0;
    static int ep[MAXELEM], ec[MAXELEM], kp[MAXELEM], kc[MAXELEM];
    clock_t t0 = clock();

    for (uint64_t it = 0; it < nsamp; it++) {
        tt_t f = rng_next() & FULL;
        analyze(f, 1, &RF);
        int pos = is_positive(&RF);
        if (pos || it % stride == 0) {
            analyze(f, 0, &RB);
            if (!same_res(&RF, &RB)) {
                fprintf(stderr, "SELF-CHECK FAILED at sample %llu f=%016llx: filter nh %d nW %d na %d, brute nh %d nW %d na %d\n", (unsigned long long)it,
                        (unsigned long long)f, RF.nh, RF.nW, RF.na, RB.nh, RB.nW, RB.na);
                return 1;
            }
            if (pos) checked_pos++; else checked_stride++;
        }
        if (__builtin_popcountll(f) == NP / 2) balanced++;
        if (!pos) { trivial_count++; continue; }
        int ne = build_h(f, &RF, ep, ec);
        gkey_t kh, kk;
        compute_gkey(ep, ec, ne, &kh);
        if (RF.na) {
            if (RF.na != ne) { fprintf(stderr, "ANTISYMMETRY COSET SIZE MISMATCH at sample %llu: |H| %d |A| %d\n", (unsigned long long)it, ne, RF.na); return 1; }
            memcpy(kp, ep, sizeof(int) * ne);
            memcpy(kc, ec, sizeof(int) * ne);
            for (int k = 0; k < RF.na; k++) { kp[ne + k] = RF.ap[k]; kc[ne + k] = RF.ac[k]; }
            compute_gkey(kp, kc, ne + RF.na, &kk);
        } else {
            kk = kh;
        }
        add_rec(RF.na > 0, &kh, &kk);
    }

    printf("# NV %d samples %llu seed %llu checked_positive %llu checked_random %llu balanced %llu seconds %.0f\n", NV, (unsigned long long)nsamp,
           (unsigned long long)seed, (unsigned long long)checked_pos, (unsigned long long)checked_stride, (unsigned long long)balanced,
           (double)(clock() - t0) / CLOCKS_PER_SEC);
    printf("# npn antisym | H: order hbar_order d orbit_lengths vertex_orbits fixed_vertices | Kbar: same six | count\n");
    int allones = 0;
    for (int i = 0; i < NV; i++) allones = allones * 10 + 1;
    printf("npn 0 1 1 0 %d %d %d 1 1 0 %d %d %d %llu\n", allones, NP, NP, allones, NP, NP, (unsigned long long)trivial_count);
    for (int i = 0; i < nrecs; i++)
        printf("npn %d %d %d %d %d %d %d %d %d %d %d %d %d %llu\n", recs[i].antisym, recs[i].h.order, recs[i].h.hbar_order,
               recs[i].h.d, recs[i].h.pattern, recs[i].h.vorbits, recs[i].h.fixv, recs[i].k.order, recs[i].k.hbar_order,
               recs[i].k.d, recs[i].k.pattern, recs[i].k.vorbits, recs[i].k.fixv, (unsigned long long)recs[i].count);
    return 0;
}
