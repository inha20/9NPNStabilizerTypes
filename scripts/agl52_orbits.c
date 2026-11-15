/*
 * agl52_orbits.c
 * ==============
 * AGL(5,2) (위수 32 * |GL(5,2)| = 319,979,520) 가 5변수 부울 함수 2^32 개 위에 작용할 때의 궤도를 전부 나열한다.
 * 함수 = 32 비트 진리표.  점 x (0..31) 의 값이 비트 x.  군 원소 g 는 점 순열이고 함수 f 를 g.f (비트 x 를 g(x) 로 옮김) 로 보낸다.
 * 생성원 4 개: 평행이동 x -> x^1, 좌표 순환 (x0..x4) -> (x4,x0,..,x3), 좌표 교환 (0 1), 전단 x1 -> x1 + x0.
 *   (순열행렬들과 기본행렬 하나가 GL(5,2) 를 만들고, 평행이동 e0 와 S5 가 모든 평행이동을 만든다)
 * 각 궤도마다 너비우선탐색으로 크기를 센다. 안정화군 위수 = |G| / 크기 가 정수여야 한다 (검산).
 * 출력: 궤도마다 (대표 = 궤도의 최솟값, 크기, 안정화군 위수, 무게) 와 요약 (궤도 수, 크기 합 = 2^32, 안정화군 위수별 궤도 수).
 *
 * -DLINEAR_ONLY 를 주면 평행이동 생성원을 빼고 GL(5,2) (위수 9,999,360, 점 0 을 고정) 가 같은 2^32 개 함수 위에서 하는 작용의 궤도를 나열한다.
 *   (점 0 의 값은 군이 건드리지 않는 자유로운 값이므로, 이것은 GL(5,2) 가 0 이 아닌 점 31개 위의 함수 2^31 개에 하는 작용의 궤도를 값 하나씩 두 번 센 것이다.)
 *
 * 컴파일: gcc -O2 -march=native -o agl52_orbits agl52_orbits.c
 *         gcc -O2 -march=native -DLINEAR_ONLY -o gl52_orbits agl52_orbits.c
 * 실행:   ./agl52_orbits > agl52_orbits_output.txt   (약 5 분, 메모리 약 2 GB)
 *
 * 2026-11-15. Claude가 쓴 코드이고 저자 검토 전이다.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#ifdef LINEAR_ONLY
#define NGEN 3
#define GOFF 1
#define GORDER 9999360ULL
#else
#define NGEN 4
#define GOFF 0
#define GORDER 319979520ULL
#endif

static uint32_t TAB[NGEN][4][256];

static void build_tables(void) {
    int p[4][32];
    for (int x = 0; x < 32; x++) {
        p[0][x] = x ^ 1;
        /* 좌표 순환: 비트 i 가 비트 (i+1)%5 로 */
        int y = 0;
        for (int i = 0; i < 5; i++) y |= ((x >> i) & 1) << ((i + 1) % 5);
        p[1][x] = y;
        /* 좌표 교환 (0 1) */
        int b0 = x & 1, b1 = (x >> 1) & 1;
        p[2][x] = (x & ~3) | (b0 << 1) | b1;
        /* 전단: x1 += x0 */
        p[3][x] = x ^ ((x & 1) << 1);
    }
    for (int g = 0; g < NGEN; g++) {
        const int *pg = p[g + GOFF];
        /* pg 가 순열인지 확인 */
        int seen[32] = {0};
        for (int x = 0; x < 32; x++) seen[pg[x]]++;
        for (int x = 0; x < 32; x++) if (seen[x] != 1) { fprintf(stderr, "순열이 아님\n"); exit(1); }
        for (int b = 0; b < 4; b++)
            for (int v = 0; v < 256; v++) {
                uint32_t w = 0;
                for (int j = 0; j < 8; j++) if ((v >> j) & 1) w |= 1u << pg[8 * b + j];
                TAB[g][b][v] = w;
            }
    }
}

static inline uint32_t apply(int g, uint32_t f) {
    return TAB[g][0][f & 255] | TAB[g][1][(f >> 8) & 255] | TAB[g][2][(f >> 16) & 255] | TAB[g][3][f >> 24];
}

int main(void) {
    build_tables();
    uint64_t *vis = calloc((1ULL << 32) / 64, sizeof(uint64_t));
    const uint64_t QCAP = 330000000ULL;
    uint32_t *queue = malloc(QCAP * sizeof(uint32_t));
    if (!vis || !queue) { fprintf(stderr, "메모리 부족\n"); return 1; }
    uint64_t norbits = 0, total = 0, maxsize = 0;
    uint64_t hist_n[64], hist_stab[64];  /* 안정화군 위수별 */
    int nh = 0;
    for (uint64_t f0 = 0; f0 < (1ULL << 32); f0++) {
        if ((vis[f0 >> 6] >> (f0 & 63)) & 1) continue;
        uint64_t head = 0, tail = 0;
        queue[tail++] = (uint32_t)f0;
        vis[f0 >> 6] |= 1ULL << (f0 & 63);
        while (head < tail) {
            uint32_t f = queue[head++];
            for (int g = 0; g < NGEN; g++) {
                uint32_t h = apply(g, f);
                uint64_t w = h >> 6, bit = 1ULL << (h & 63);
                if (!(vis[w] & bit)) {
                    vis[w] |= bit;
                    if (tail >= QCAP) { fprintf(stderr, "큐 용량 초과\n"); return 2; }
                    queue[tail++] = h;
                }
            }
        }
        uint64_t size = tail;
        if (GORDER % size != 0) { fprintf(stderr, "크기 %llu 가 |G| 의 약수가 아님\n", (unsigned long long)size); return 3; }
        uint64_t stab = GORDER / size;
        int wt = __builtin_popcount((uint32_t)f0);
        printf("orbit %llu rep=%08llx size=%llu stab=%llu weight=%d\n", (unsigned long long)norbits, (unsigned long long)f0,
               (unsigned long long)size, (unsigned long long)stab, wt);
        norbits++;
        total += size;
        if (size > maxsize) maxsize = size;
        int found = -1;
        for (int i = 0; i < nh; i++) if (hist_stab[i] == stab) { found = i; break; }
        if (found < 0) { found = nh++; hist_stab[found] = stab; hist_n[found] = 0; }
        hist_n[found]++;
        if ((norbits & 15) == 0) fflush(stdout);
    }
    printf("SUMMARY orbits=%llu total=%llu (2^32=%llu) maxsize=%llu |G|=%llu\n", (unsigned long long)norbits,
           (unsigned long long)total, (unsigned long long)(1ULL << 32), (unsigned long long)maxsize, (unsigned long long)GORDER);
    /* 안정화군 위수 오름차순으로 출력 */
    for (int i = 0; i < nh; i++)
        for (int j = i + 1; j < nh; j++)
            if (hist_stab[j] < hist_stab[i]) {
                uint64_t t = hist_stab[i]; hist_stab[i] = hist_stab[j]; hist_stab[j] = t;
                t = hist_n[i]; hist_n[i] = hist_n[j]; hist_n[j] = t;
            }
    for (int i = 0; i < nh; i++)
        printf("HIST stab=%llu orbits=%llu\n", (unsigned long long)hist_stab[i], (unsigned long long)hist_n[i]);
    printf("FREE_ORBITS %llu\n", (unsigned long long)((nh > 0 && hist_stab[0] == 1) ? hist_n[0] : 0));
    return 0;
}
