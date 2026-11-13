#!/usr/bin/env python3
"""노트의 수치(네 자리 이상 또는 소수)가 출력 파일 어딘가에 나타나는지 본다. 일치 여부의 필요조건만 본다.

사용: python scripts/check_note_numbers.py 노트.md 출력파일... (.gz 도 읽는다)
검사하는 줄은 "|", "- ", "**", "1." 꼴로 시작하는 줄이고, 코드 조각(`...`)과 수식($...$)은 건너뛴다.
출력에 그대로 없는 수(절 번호, 자리를 줄인 값, 계산해 낸 값)도 걸러 낸다. 그런 것은 사람이 따로 본다.
"""
import glob, gzip, re, sys, os
note = sys.argv[1]
outs = sys.argv[2:]
def norm_text(t):
    t = t.replace("−", "-")
    t = re.sub(r"(?<=\d),(?=\d{3}(?!\d))", "", t)   # 천 단위 쉼표 제거
    return t
hay = ""
for p in outs:
    if p.endswith(".gz"):
        hay += gzip.open(p, "rt", errors="ignore").read() + "\n"
    else:
        hay += open(p, errors="ignore").read() + "\n"
hay = norm_text(hay)
tok = re.compile(r"(?<![\w.])[+-]?\d+(?:\.\d+)?")
hayset = set(m.group(0).lstrip("+") for m in tok.finditer(hay))
# 부호 없는 값과 소수 둘째 자리 반올림 쪽도 맞출 수 있게 절댓값도 넣는다
hayset |= {h.lstrip("-") for h in hayset}
text = open(note, encoding="utf-8").read()
miss = []
checked = 0
for ln, line in enumerate(text.splitlines(), 1):
    if not (line.startswith("|") or line.startswith("- ") or line.startswith("**") or re.match(r"\d+\.", line)):
        continue
    l = norm_text(line)
    l = re.sub(r"`[^`]*`", " ", l)           # 코드 조각 제외
    l = re.sub(r"\$[^$]*\$", " ", l)         # 수식 제외
    for m in tok.finditer(l):
        s = m.group(0).lstrip("+")
        digits = s.lstrip("-").replace(".", "")
        if len(digits) < 4 and "." not in s:
            continue
        if re.fullmatch(r"20\d\d", s):      # 연도
            continue
        checked += 1
        if s not in hayset and s.lstrip("-") not in hayset:
            miss.append((ln, s, line.strip()[:110]))
print("검사한 수치 %d개, 출력 파일에서 찾지 못한 것 %d개" % (checked, len(miss)))
for ln, s, l in miss:
    print("  줄 %d: %s    | %s" % (ln, s, l))
