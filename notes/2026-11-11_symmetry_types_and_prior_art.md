# 기록: 어떤 대칭이 어떤 그림 종류 수를 주는가, 그리고 선행 연구 1차 조사

2026-11-11. `2026-11-11_picture_types_n5_n6.md`의 분포를 대칭의 종류별로 나눠 보고, 관련 문헌을 처음으로 찾아본 기록이다.

## 1. 대칭 종류별 분해

함수 $f$의 $B_n$-안정화군이 변수 순열로 보이는 모습을 $\bar H \le S_n$이라 하자(입력 반전을 잊은 상). 그림 종류 수는 $\bar H$가 변수 쌍 분할 집합 $\mathcal M_n$에 만드는 궤도 수이다(명제 A2). $S_n$의 부분군 켤레류마다 그 궤도 수와, 그 켤레류가 실제로 어떤 함수의 $\bar H$로 나타나는지를 셌다. 전체 표는 `picture_types_by_symmetry_output.txt`, 스크립트는 `scripts/picture_types_by_symmetry.g`이다.

**(가) 빠진 값은 부울 함수 때문이 아니라 $S_n$ 때문이다.** $S_5$의 부분군 19종이 $\mathcal M_5$에 만드는 궤도 수는 $\{1,2,3,4,5,6,9,15\}$뿐이고, $S_6$의 56종은 $\{1,\dots,9,11,15\}$뿐이다. 10, 12, 13, 14개의 궤도를 만드는 부분군이 애초에 없다. 가장 작은 비자명 대칭인 변수 두 개의 맞바꿈이 이미 15를 9로 줄인다.

**(나) 부분군이 줄 수 있는 값은 모두 실제로 나타난다.** $n = 4, 5, 6$ 모두에서 "부분군이 줄 수 있는 $T$ 값의 집합"과 "함수에서 실제로 나타나는 $T$ 값의 집합"이 같다. 7번 논문 따름정리 6($n=4$에서 1, 2, 3)의 $n = 5, 6$ 판이다.

**(다) 그러나 모든 부분군이 나타나지는 않는다.** 어떤 함수의 $\bar H$로도 나타나지 않는 켤레류:

| $n$ | 켤레류 수 | 나타나지 않는 것 |
|:---:|:---:|:---|
| 4 | 11 중 4 | $C_3$, $C_4$, 추이적인 $C_2 \times C_2$, $A_4$ |
| 5 | 19 중 6 | $C_3$, $C_5$, $C_6$, 궤도 길이 [2,3]인 $S_3$, $A_4$, $A_5$ |
| 6 | 56 중 15 | 궤도 길이 [1,1,1,3]인 $C_3$, $C_5$, $C_3 \times C_3$, $A_4$ [1,1,4], $A_5$ [1,5], $A_6$ 등 (전체는 출력 파일) |

$n=4$에서 $C_3$과 $A_4$가 나타나지 않는다는 것이 7번 논문 따름정리 6(세 분할 위의 상의 위수는 1, 2, 6뿐)과 같은 내용이다. $S_4$의 부분군 가운데 세 분할 위의 상이 위수 3인 것은 이 둘뿐이기 때문이다. $C_4$와 추이적인 $C_2 \times C_2$가 나타나지 않는다는 것은 따름정리 6에 들어 있지 않은 정보이다.

*(2026-11-11 정정: 처음에는 이 문단에 "이 넷의 $S_3$-상이 위수 3이거나, 나타나지 않기 때문이다"라고 뜻이 불분명하게 적었다.)*

**(라) 각 값의 주된 출처** ($n = 5$): $T = 9$는 위수 2인 대칭 두 종류, 곧 변수 두 개의 맞바꿈(301,664,160개)과 두 쌍을 한꺼번에 맞바꾸는 것(46,811,520개)에서 나오고, $T = 6$은 $C_2 \times C_2$, $T = 5$는 $C_4$ 하나뿐(42,720개), $T = 4$는 $S_3$과 $D_8$, $T = 1$은 $S_5$ 전체 대칭(964개)뿐이다. $T = 5$가 드문 까닭은 $C_3$이 나타나지 않아 $C_4$만 남기 때문이다. $T = 3$은 $D_{10}$과 $D_{12}$, $T = 2$는 $C_5 \rtimes C_4$와 $S_4$에서 나온다.

*(2026-11-11 정정: 처음에는 "$T = 9$는 변수 맞바꿈 하나"라고 적었으나, 두 쌍을 한꺼번에 맞바꾸는 대칭도 9를 준다.)*

**(마) 안정화군으로 나타나는 $B_n$의 부분군은 일부이다.** $B_5$의 부분군 켤레류 953개 가운데 어떤 함수의 안정화군으로 나타나는 것은 143개, $B_6$에서는 7,440개 가운데 790개이다(`scripts/export_class_representatives.g`의 출력).

**(바) $n = 6$과 $S_6$의 외부 자기동형.** $n=6$의 표에서 변수 두 개의 맞바꿈은 $T=9$, 세 쌍을 한꺼번에 맞바꾸는 것은 $T=11$, 3-순환은 $T=5$, 3-순환 둘의 곱은 $T=7$을 준다. 이 값들은 각각 "세 쌍 맞바꿈", "맞바꿈 하나", "3-순환 둘의 곱", "3-순환"이 변수 쌍 15개에 만드는 궤도 수와 같다(손으로 확인: $(15+3)/2$, $(15+7)/2$, $15/3$, $(15+2\cdot3)/3$). 완전 매칭 15개 위의 작용이 쌍 15개 위의 작용을 $S_6$의 외부 자기동형으로 비튼 것이라는 고전적 사실(Sylvester의 synthemes와 duads)과 맞는 모습이다. 처음에는 다섯 줄에서만 확인했고, 2026-11-12에 56개 켤레류 전체에서 확인했다(`scripts/check_s6_outer_automorphism.g`: 다른 켤레류 0개. 옮기지 않고 비교하면 36개가 다르다).

## 1-2. 입력 반전을 허용하면 달라지는 것 (2026-11-12 추가)

선행 문헌(2절)은 변수 순열만으로 본 대칭군 $G(f) = \{\pi : f(\pi x) = f(x)\}$를 다룬다. 이 폴더의 $\bar H(f)$는 입력 반전까지 허용한 대칭의 순열 부분이라서 $G(f) \subseteq \bar H(f)$이고 더 클 수 있다. 두 뜻에서 각각 어떤 부분군이 나타나는지를 같은 방법(표지 표)으로 계산했다. 표는 `picture_types_by_symmetry_output.txt`의 마지막 열이고, 5변수 이하는 $2^{32}$개 전수 열거로도 같은 값이 나왔다(`enumeration_comparison_output.txt`).

| $n$ | 켤레류 | $G(f)$로 나타나지 않는 것 | $\bar H$로 나타나지 않는 것 | $\bar H$로만 나타나는 것 (함수 수) |
|:---:|:---:|:---:|:---:|:---|
| 3 | 4 | 1 | 2 | 없음 |
| 4 | 11 | 4 | 4 | 없음 (두 목록이 같다) |
| 5 | 19 | 8 | 6 | $C_4$ [1,4] (42,720), $C_5 \rtimes C_4$ [5] (768) |
| 6 | 56 | 21 | 15 | $C_4$ [1,1,4] (776,845,440), $C_4 \times C_2$ [2,4] (20,668,320), $D_8$ [2,4] (20,466,720), $C_5 \rtimes C_4$ [1,5] (297,216), $S_4$ [6] (15,360), $S_5$ [6] (3,456) |

대괄호는 변수 위의 궤도 길이이다.

- **4변수에서는 차이가 없다.** 두 뜻에서 나타나지 않는 부분군이 똑같이 $C_3$, $C_4$, 추이적인 $C_2 \times C_2$, $A_4$이다.
- **5변수부터 차이가 난다.** 변수 순열만으로는 어떤 함수의 대칭군도 되지 못하는 군이 입력 반전을 허용하면 나타난다. 5변수의 $C_4$가 그 예이고, 5변수에서 그림 종류 수가 5인 함수 42,720개는 모두 이 경우이다. 곧 $T = 5$는 입력 반전이 섞인 대칭 때문에만 생기는 값이다.
- **반대 방향은 없다.** $n = 4, 5, 6$에서 $G(f)$로 나타나는 부분군은 모두 $\bar H$로도 나타난다. $n=3$에서는 자명군이 예외이다: 대칭이 전혀 없는 3변수 함수는 변수 순열만 보면 96개 있지만 입력 반전을 허용하면 하나도 없다.
- 켤레류 번호는 실행할 때마다 달라질 수 있어서, 군은 (구조, 궤도 길이)로 가리켰다.

이 차이가 문헌에 있는지는 아직 모른다(2절).

## 2. 선행 연구 1차 조사 (웹 검색 7회)

**부울 함수의 대칭군 — 가장 가까운 선행 연구.** $S_n$의 어떤 부분군이 부울 함수의 대칭군(변수 순열만, 입력 반전 없음)이 되는가는 Clote–Kranakis (1991)가 시작하고 Kisielewicz (1998, *J. Algebra*), Grech–Kisielewicz (2014, *European J. Combin.*; arXiv:1303.1290) 등이 이어 온 주제이다. 검색 결과의 초록·서론에서 확인한 것:

- 교대군 $A_n$은 어떤 부울 함수의 대칭군도 아니다(Clote–Kranakis).
- 정칙(regular) 순열군은 $C_3$, $C_4$, $C_5$, $C_2 \times C_2$ 네 예외를 빼고 모두 부울 함수의 대칭군이다(Grech).
- 아벨 순열군은 알려진 예외를 빼고 모두 대칭군이다(Grech–Kisielewicz).

**본문에서 확인한 것 (2026-11-12, Grech–Kisielewicz 2014의 앞 3분의 2).**

- 대칭군의 정의는 변수 순열뿐이다. 입력 반전을 허용하는 경우는 읽은 범위에 없다.
- 함수를 대칭군별로 센 표는 읽은 범위에 없다.
- $C_3$, $C_4$, $C_5$에 다른 군을 직합한 $C_i \oplus H$도 대칭군이 되지 못한다. 고정점을 붙인 $C_4$(5변수의 $C_4$ [1,4])가 $G(f)$로 나타나지 않는다는 위 1-2절의 계산과 맞는다.
- $(1\,2\,3)(4\,5\,6)$으로 생성되는 순환군은 대칭군이 된다. 6변수의 $C_3$ [3,3]이 $G(f)$로 나타난다는 계산(273,853,440개)과 맞는다.
- 주정리의 예외 목록이 적힌 뒷부분은 읽지 못했다(도구의 사용 한도). Maurer (2015)는 이 절을 쓸 때는 읽지 못했고 나중에 앞부분을 읽었다(4절).

**Horváth 외 (arXiv:1210.1015)의 앞 4분의 1에서 확인한 것 (2026-11-12).**

- 이 논문도 변수 순열만 다룬다. 입력 반전은 나오지 않고, 함수를 대칭군별로 센 것도 없다.
- 원시(primitive) 순열군 가운데 부울 함수의 대칭군이 되지 못하는 것으로 4변수의 $A_4$, 5변수의 $C_5$·$AGL(1,5)$·$A_5$, 6변수의 $PGL(2,5)$·$A_6$을 든다. 11변수 이상에서는 교대군뿐이라고 한다.
- 1-2절의 계산과 맞는다. 위수 20인 $C_5 \rtimes C_4$가 $AGL(1,5)$이고, 6변수의 추이적인 $S_5$가 $PGL(2,5)$이다. 그리고 이 둘은 입력 반전을 허용하면 나타난다(각각 768개, 3,456개 함수). $C_5$, $A_5$, $A_6$은 입력 반전을 허용해도 나타나지 않는다.

위 (다)의 목록($C_3$, $C_4$, $C_5$, 추이적 $C_2 \times C_2$, $A_n$)은 이 예외들과 겹친다. 따라서 **(다)는 이 문헌의 결과에서 따라 나오거나 그 변형일 가능성이 크다.** 다른 점은 이쪽이 입력 반전을 허용한 안정화군의 상을 본다는 것인데, 그 설정을 다룬 문헌이 있는지는 아직 찾지 못했다. 논문들의 본문은 읽지 않았다.

**5·6변수 지도의 규약.** Kuphaldt의 교재(*Lessons in Electric Circuits* 4권 8장)는 5·6변수 지도를 두 방식으로 설명한다. 전체 Gray 부호로 번호를 매긴 반사형(reflection)은 옛 교재의 방식이고, 4×4 지도를 겹쳐 쌓는 겹침형(overlay)이 지금 더 흔하다고 적는다. 이 폴더의 계산은 겹침형(보충 노트의 모델 A)에 해당한다. 교재 한 권의 서술이므로 일반적 사실로 쓰려면 교재를 더 확인해야 한다.

**논리 합성 쪽의 대칭 연구.** 회로 설계 자동화 분야는 부울 함수의 대칭을 오래 다뤄 왔다. Kravets–Sakallah (2000, ICCAD, 526–532쪽)의 *Generalized Symmetries in Boolean Functions*가 있고(처음에는 서지만 확인했고 기술보고서 본문의 앞부분은 4절에서 읽었다), Katebi–Sakallah–Markov (2013, ICCAD)는 기존 알고리즘 대부분이 변수 두 개의 맞바꿈 같은 고전적 대칭만 찾는다고 적는다. Maurer (2015, *SpringerPlus*)는 입력 수 1~18에 대해 대칭의 수와 대칭 유형의 수를 표로 싣고, 부분 대칭과 전체 대칭은 입력이 4개 이상이면 유형의 절반에 못 미친다고 적는다. 대칭의 종류를 순열군으로 분류하는 틀은 이 분야에 이미 있다는 뜻이다. 본문은 처음에는 읽지 않았고 앞부분을 4절에서 읽었다.

**세는 쪽의 문헌.** NP·NPN 동치류의 수는 Zhang 외 (2019, *Theory of Computing Systems*)가 10변수까지 계산했다. 부울 함수는 초입방체 꼭짓점의 부분집합이므로, 초입방체 대칭 아래 0/1-다면체의 동치류를 세는 문헌(arXiv:1101.0410 등)도 같은 대상을 다룬다. Filmus는 5변수까지의 NPN 동치류 대표원 목록을 공개해 두었다.

**찾지 못한 것.** 변수 배치별 그림의 수를 5·6변수에서 센 문헌, 부울 함수를 $B_n$-안정화군의 종류별로 센 표(표지 표 방식), 입력 반전을 허용한 안정화군의 상으로 어떤 순열군이 나타나는지를 다룬 문헌. 검색 7회로 "없다"고 할 수는 없다.

**OEIS (2026-11-12 저녁에 확인함).** 이 폴더의 수열(예: $W_n$-궤도 수 6, 34, 805, 17107114; 안정화군이 자명한 함수 수 22656, 3926492160; 그런 $B_n$-궤도 수 59, 1022524)은 처음에는 OEIS에 있는지 확인하지 못했다(검색 페이지가 Claude의 웹 도구로 열리지 않았고, 숫자를 웹 검색해서는 나오지 않았다). 같은 날 저녁에 Claude가 내장 브라우저로 oeis.org를 검색했고 위 수열 모두 일치하는 것이 없었다. 검색한 수열의 전체 목록은 `2026-11-12_stabilizer_types.md` 4절에 있다.

## 3. 이것이 뜻하는 것

- 그림 종류 수의 "가능한 값" 문제는 부울 함수가 아니라 순열군의 문제로 환원된다((가), (나)). 새로움이 있다면 값의 목록보다 **개수 분포**와 **$B_n$-안정화군의 상으로 나타나는 부분군의 분류** 쪽이다.
- 1-2절의 차이(입력 반전을 허용하면 더 나타나는 군)가 이 연구에서 문헌과 가장 분명하게 갈리는 지점이다. 이미 알려진 것인지 확인하는 일이 선행 연구 정독의 핵심 질문이 된다.
- 그 분류는 Kisielewicz 계열 문헌과 직접 맞닿아 있다. 원고를 쓰기 전에 Kisielewicz (1998)와 Grech–Kisielewicz (2014)의 본문을 읽고, 입력 반전을 허용한 경우가 이미 다뤄졌는지 확인해야 한다.

## 4. 추가 검색: 입력 반전을 허용한 대칭군 (2026-11-12)

검색 두 번(이진 부호의 자기동형군, 부호 있는 순열로 본 부울 함수의 대칭군). 검색 결과의 일부 문장만 읽었다.

- **Slepian (1953, *On the number of symmetry types of Boolean functions of n variables*, Canad. J. Math.).** 변수의 순열과 일부 변수의 보수 취하기가 이루는 군(위수 $2^n n!$, 초팔면체군과 동형)의 궤도를 "대칭 유형"이라 부르고 그 수를 센다. 보고하는 값은 $n = 3, \dots, 6$에서 22, 402, 1,228,158, 400,507,806,843,728이고(PDF를 도구로 요약해 읽음), 이 폴더의 $B_n$-궤도 수와 같다(OEIS A000616의 출처). 안정화군은 궤도의 크기를 구할 때만 쓰이고, 안정화군의 종류를 분류하거나 함수를 안정화군별로 세지는 않는다고 한다. 곧 이 폴더가 하는 "안정화군별로 센다"는 Slepian이 한 일이 아니다.
- **Katebi–Sakallah–Markov (ICCAD 2013).** 이 논문은 대칭을 "입출력의 순열(반전 가능)"로 소개하고, 4:1 멀티플렉서가 입출력의 순열·반전 아래 16개의 대칭을 가진다는 예(그림 1)를 든다. 그러나 제안한 탐지 알고리즘은 순열만 다루고 반전은 앞으로의 과제로 남긴다(PDF를 도구로 요약해 읽은 결과와 검색 결과의 문장). 대칭군의 어떤 부분군이 나타나는지의 분류는 없다. 처음에 "반전을 포함한 대칭을 탐지한다"고 적었던 것은 틀렸으므로 앞서 고쳤고, 이번에 "정의는 반전을 포함하고 알고리즘은 순열만"으로 다시 다듬었다.
- **Kravets–Sakallah (ICCAD 2000, 미시간대 기술보고서 CSE-TR-420-00).** 보고서 12만 자 가운데 앞의 10만 자를 도구로 요약해 읽었다. 이 글의 대칭은 "변수 순열 가운데 보수 취하기를 허용한 것"이어서 입력 반전을 포함한다(5절에서 입력에 위상을 주는 확장을 다룬다). 그러나 "대칭군"이라는 말은 순열로 자유롭게 바꿀 수 있는 변수 묶음을 가리키고, $B_n$ 안의 완전한 안정화군으로 정의하지는 않는다고 한다. 분할 계층으로 순열 집합의 일부를 암묵적으로 나타내며 완전하다고 주장하지도 않는다. 탐지 알고리즘(1차 대칭은 이차 시간의 분할, 고차는 공인자 행렬 비교)과 벤치마크 표가 내용이고, 어떤 부분군이 나타나는지의 분류나 목록은 없다. 예는 $f = a + b'$가 위상을 준 쌍 $\{a', b\}$에서 대칭이라는 정도의 작은 것이다. 곧 입력 반전을 포함한 순열 대칭은 로직 합성에서 이미 다뤄진 대상이며, 이 폴더의 R5가 새롭다면 "탐지"가 아니라 "어떤 군이 나타나고 함수가 몇 개인가"의 분류 쪽이다.
- **Maurer (2015, *SpringerPlus*, "A universal symmetry detection algorithm").** 본문 12만 2천 자 가운데 앞의 10만 자를 도구로 요약해 읽었다. 순열군 하나를 통째로 검사한다. 진리표의 입력 벡터를 그 군의 궤도로 나누고, 함수가 궤도마다 한 값(또는 그 보수)을 갖는지 본다. 입력 반전은 "다상(multi-phase) 대칭"으로 다룬다. 보통의 대칭에서 일부 입력을 뒤집은 것이며, 구현은 각 궤도에 상수 벡터를 XOR하는 것이다. 이것은 순열군 $G$를 평행이동 $t$로 공액한 군 $tGt^{-1}$로 보이고(Claude의 읽기이며 본문에서 확인하지는 않았다), 5변수 $C_4$ 예처럼 반전이 순열과 엮여 있어 공액으로 풀리지 않는 군과는 다르다(그 예에서 계산으로 확인, `2026-11-12_example_c4_n5.md` 5절). 요약 도구는 이 논문이 틀을 GF(2) 위의 행렬군으로 넓힌다고도 전했으므로, 그 틀에서 평행이동이 순열과 엮인 군이 이미 다뤄졌을 가능성은 배제하지 못한다. 논문은 4변수의 65,536개 함수 전부를 대칭 유형($S_n$의 부분군의 공액류)별로 세고 5변수도 같은 식으로 센다고 한다(수치는 그림 1에 있어 텍스트로 얻지 못했다). 6변수 이상은 전수가 어려워 부분군의 목록을 쓴다. 정확한 대칭군이 될 수 없는 "금지된" 군(Klein 4군, 교대군 등)을 든다. 곧 **순열만 쓴 대칭 유형별 함수 수는 4·5변수에서 이미 있다.** 이 폴더가 새롭다면 반전까지 넣은 안정화군 쪽이다.
- **고전 스위칭 이론의 "군 불변성" (1956~1965).** McCluskey (1956, *Bell System Technical Journal* 35(6), 1445–1453쪽, MIT 박사 학위 논문에서 나옴)의 *Detection of Group Invariance or Total Symmetry of a Boolean Function*은 변수의 순열, 일부 변수의 프라임(반전), 또는 둘을 함께 해도 함수가 바뀌지 않는 연산이 무엇인지를 주어진 함수에서 가려내는 방법을 보인다(서지와 도입부 몇 문장만 읽었고 본문은 읽지 못했다). 곧 $B_n$ 안에서 함수의 안정화군을 구하는 문제는 1956년에 이미 세워져 있었다. Golomb (1959, *IRE Trans. Circuit Theory* 6(5), 176–186쪽, *On the Classification of Boolean Functions*)의 초록은 순열과 보수 취하기로 서로 옮겨지는 함수를 같은 "대칭류"로 묶고, 대칭류의 수를 세는 방법과 두 함수가 같은 류인지 가리는 완전한 불변량을 다루며 그 불변량이 류의 크기와 대칭 성질에 관한 정보도 준다고 적는다(초록만 읽음). Choudhury–Basu (1961, *Indian J. Physics* 36, 31–42쪽)는 McCluskey와 거의 같은 제목의 논문이나 내용은 모른다(서지만 확인). Sheng (1965, *IEEE Trans. Electronic Computers* EC-14, 924–926쪽)의 초록은 반전을 섞은 변수에 대해서도 완전 대칭을 검출하는 방법을 낸다고 적는다. 이 글들이 안정화군의 종류를 분류하거나 함수를 안정화군별로 세는지는 초록으로는 알 수 없다. 이 글들은 검출과 불변량이 중심이다. BSTJ의 McCluskey 논문은 공개된 PDF가 있어 저자가 직접 읽을 수 있을 것이다(이 작업 공간에서는 그 PDF를 내려받지 못했다). Urquhart (2009, *Bull. Symbolic Logic* 15(3), 273–299쪽, *Enumerating Types of Boolean Functions*)는 순열과 보수 취하기의 군 아래에서 함수의 유형을 세는 문제의 1870년대부터 1970년대까지의 역사를 다루는 서베이로, 이 문헌들의 소재를 찾는 길잡이가 될 수 있다(초록만 확인했고 본문은 접근이 막혀 읽지 못했다).
- **Forbes (1971, 일리노이대 Coordinated Science Laboratory 보고서 UILU-ENG 71-2229, R-526).** 제목은 *Maximal Groups of Permutations and Complementations of the Independent Variables Which Leave a Set of Boolean Functions Invariant*이다. 변수의 순열과 보수 취하기(곧 $B_n$의 부분군)가 부울 함수의 집합을 불변으로 두는 경우를 다루는 것으로 읽히므로, 이 폴더의 R5와 주제가 맞닿을 가능성이 있다. **본문은 읽지 못했다.** 일리노이대 IDEALS 저장소의 기록(https://hdl.handle.net/2142/75037)에는 서지 정보만 있고 초록과 PDF가 없으며, 열람이 그 대학 구성원으로 제한되어 있다. 따라서 단일 함수의 안정화군을 분류하는지, 함수 집합의 불변군만 다루는지, 함수의 수를 세는지는 알 수 없다. 도서관 간 상호대차가 유일한 접근 경로로 보인다. 이 보고서가 R5를 이미 다뤘는지는 구한 뒤에 판단해야 하며, 그전까지 "R5는 새롭다"고 쓰면 안 된다.
- **선형 대칭.** Clausen의 결과로 보이는 서지(*Information Processing Letters*, 1992)가 일반선형군 GL(n,2)의 작용 아래 거의 모든 n변수 부울 함수의 안정화군이 자명하다고 말한다. $B_n$은 아핀군 AGL(n,2)의 부분군이므로 같은 방향의 결과이며, 이 폴더의 R7(대칭이 없는 함수가 대부분)과 닿아 있다. 서지만 확인했다.
- **Grech–Kisielewicz (arXiv:1303.1290, 1910.11806).** $G(f)$의 정의는 순열만 쓴다(초록과 첫 문장에서 확인).
- **이진 부호의 자기동형군.** 검색에서 나온 논문들(순환 부호, 선형 부호)은 순열만 쓴다. 평행이동(입력 반전에 해당)까지 넣은 비선형 부호의 자기동형군과의 관계는 찾지 못했다.

입력 반전을 허용한 안정화군의 분류(R5)가 이미 알려져 있다는 증거는 이번에도 찾지 못했다. 이것이 "없다"는 증명은 아니다. 입력 반전을 포함한 군을 쓴 문헌은 Slepian(유형의 수), McCluskey(1956)와 Kravets–Sakallah(탐지), Golomb(1959, 불변량과 류의 크기)이고, Katebi 외는 정의에만 반전을 넣었다. 로직 합성 쪽이 순열만 다룬다고 한 앞의 문장은 Kravets–Sakallah에 대해서는 틀렸으므로 고쳤다. Maurer 2015는 반전을 다상 대칭(순열군의 평행이동 공액)으로 다룬다. 가장 가까워 보이는 것은 제목만 확인한 Forbes (1971)이며, 이 보고서를 구하기 전에는 R5의 새로움을 주장할 수 없다.

## 출처

- Grech, Kisielewicz, *Symmetry groups of boolean functions*: https://arxiv.org/abs/1303.1290
- Grech, Kisielewicz, *Symmetry groups of boolean functions: simple groups*: https://arxiv.org/abs/1910.11806
- Grech, *Regular symmetric groups of boolean functions* (초록): https://www.bohrium.com/paper-details/regular-symmetric-groups-of-boolean-functions/811771028370882561-992
- Kuphaldt, *Larger 5 & 6-variable Karnaugh Maps*: https://www.allaboutcircuits.com/textbook/digital/chpt-8/larger-5-6-variable-karnaugh-maps/
- Kravets, Sakallah, *Generalized Symmetries in Boolean Functions* (서지): https://researchr.org/publication/KravetsS00/bibtex
- Kravets, Sakallah, *Generalized Symmetries in Boolean Functions* (기술보고서 CSE-TR-420-00, 앞의 10만 자): https://www.eecs.umich.edu/techreports/cse/00/CSE-TR-420-00.pdf
- Katebi, Sakallah, Markov, *Generalized Boolean Symmetries Through Nested Partition Refinement*: https://web.eecs.umich.edu/~imarkov/pubs/conf/iccad13-sym.pdf
- Maurer, *A universal symmetry detection algorithm* (앞의 10만 자를 요약해 읽음): https://link.springer.com/10.1186/s40064-015-1156-7
- Zhang 외, *A Group Algebraic Approach to NPN Classification of Boolean Functions*: https://pdxscholar.library.pdx.edu/ece_fac/505
- 0/1-다면체의 동치류 열거: https://arxiv.org/abs/1101.0410
- Horváth, Makay, Pöschel, Waldhauser, *Invariance groups of finite functions and orbit equivalence of permutation groups*: https://arxiv.org/abs/1210.1015
- Filmus, *NPN equivalence classes of Boolean functions*: https://yuvalfilmus.cs.technion.ac.il/publications/manuscripts/?id=787
- Slepian, *On the number of symmetry types of Boolean functions of n variables* (Canad. J. Math., 1953): https://cambridge.org/core/services/aop-cambridge-core/content/view/D47B82FA1F7C3712D4367AD8DAB5793D/S0008414X00043844a.pdf/on-the-number-of-symmetry-types-of-boolean-functions-of-n-variables.pdf
- Urquhart, *Enumerating Types of Boolean Functions* (Bull. Symbolic Logic 15(3), 2009; 초록만 확인): https://www.projecteuclid.org/euclid.bsl/1246453975
- McCluskey, *Detection of Group Invariance or Total Symmetry of a Boolean Function* (BSTJ 35(6), 1956; 서지와 도입부 일부): https://www.nokia.com/bell-labs/publications-and-media/publications/detection-of-group-invariance-or-total-symmetry-of-a-boolean-function
- Golomb, *On the Classification of Boolean Functions* (IRE Trans. Circuit Theory, 1959; 초록): https://manuscript.isc.ac/Inventory/49/1201073.htm
- Sheng, *Detection of Totally Symmetric Boolean Functions* (IEEE TEC, 1965; 초록): https://manuscript.isc.ac/Inventory/49/886985.htm
- Forbes, *Maximal Groups of Permutations and Complementations of the Independent Variables Which Leave a Set of Boolean Functions Invariant* (서지만 확인, 본문은 열람 제한): https://hdl.handle.net/2142/75037
- Clausen의 서지(GL(n,2)의 작용): https://www.bohrium.com/paper-details/symmetry-groups-of-boolean-functions-and-constructions-of-permutation-groups/812453836965806082-914
