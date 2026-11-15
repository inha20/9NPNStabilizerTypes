# enumerate_bn_stream.g
# scripts/enumerate_bn_subgroups.g 를 부분공간 류를 PARTS 몫으로 나눠 한 몫(PART)씩 돌리는 구동 파일이다. n = 7 처럼 오래 걸릴 때 쓴다.
# 사용 (저장소 맨 위 폴더에서, 인수는 환경 변수로 준다):
#   N=7 PARTS=8 PART=0 OUTDIR=/tmp/enum7 gap -A -q -o 3g scripts/enumerate_bn_stream.g < /dev/null
#   PART 를 0..PARTS-1 로 바꿔 몫마다 프로세스를 돌린다(코어 수만큼 동시에). 몫마다 OUTDIR 에 out_<PART>.tsv(나열), log_<PART>.txt(진행 기록),
#   state_<PART>.g(다시 시작용 상태)가 생긴다. 도중에 멈추면 같은 명령을 다시 쓰면 끝낸 부분공간 류 다음부터 이어서 한다
#   (처음 시작할 때는 state 파일이 없어도 된다). 다 끝나면 scripts/concat_enum_streams.py 로 out_*.tsv 를 이어 붙여 번호를 다시 매긴다.
# 시간: n = 7 은 PARTS=8 에 코어 4개로 약 20분이 걸렸다(2026-11-13, 두 번에 나눠 돌렸다). n <= 6 은 PARTS=1 로 1분 안쪽이다.
# 2026-11-14. Claude가 쓴 코드이고 저자 검토 전이다. 2026-11-13 의 임시 구동 파일 rf_stream.g 를 정리한 것이다.
SE := GAPInfo.SystemEnvironment;
N_ENUM := Int(SE.N);
REP_STEP := Int(SE.PARTS);
REP_OFF := Int(SE.PART);
KMIN := 1;
RESUME := true;
OUTFILE := Concatenation(SE.OUTDIR, "/out_", SE.PART, ".tsv");
LOGFILE := Concatenation(SE.OUTDIR, "/log_", SE.PART, ".txt");
STATEFILE := Concatenation(SE.OUTDIR, "/state_", SE.PART, ".g");
Read("scripts/enumerate_bn_subgroups.g");
QUIT;
