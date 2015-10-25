#####################
## R 소개
#####################
### 도움말의 활용
help(plot) 		# 그래프 함수에 대한 도움말
?rnorm		# 정규 분포에 대한 난수 발생 도움말
help.start()	# 도움말에 대한 메인 화면
help.search("graph")	# 기능중에 해당 단어를 포함하는 대상 찾기
apropos("stat")	# R 객체중에 "stat"란 명칭을 부분적이라도 포함하는 대상 찾기

### R의 종료
q()

### R 기능에 대한 demo
demo()
demo(is.things)
demo(Japanese)
demo(lm.glm)

### R 기능 설명에서 예제 실행
example(Japanese)
example(plot)
example(rnorm)

### 데이터 셋의 활용
mtcars
mtcars$mpg
# 데이터 셋(mtcars)를 명시적으로 고정하기
attach(mtcars)
mpg
cyl
# 데이터 셋(mtcars)를 작업 공간에서 분리하기
detach(mtcars)

### 디렉토리 관련 명령
getwd() 			#작업 디렉토리 알아내기
# 디렉토리 설정시 "\" 대신에 "/" 사용 필요함
setwd("C:/work/R/DataSets") 		# 작업 디렉토리 설정하기
list.files()  		#작업 디렉토리의 파일 목록 보기

### 기본데이터 타입
# 수치형(numeric)
a=1
mode(a)
length(a)
# 문자열(character)
a="happy"
mode(a)
length(a)
# 논리형(logical)
a=TRUE
mode(a)
length(a)
# 리스트(List)
a=list(1,2,"A",3,4)
mode(a)
length(a)

## 기본 데이터 타입에 대한 연산
# 산술 연산(Arithmetic) - 수치형 데이터에 대하여 수행
1+1 	# 더하기
1-1 	# 빼기
2*2	# 곱하기
2/2	# 나누기
2^2	# 제곱(2의 제곱)
2^3	# 제곱(3의 제곱)
10%/%3 # 몫 구하기(정수)
10%%3	# 나머지 구하기(정수)

# 비교 연산 - 결과가 참(TRUE) 혹은 거짓(FALSE)
1 < 2		# 보다 작다(lesser than)
1 > 2 	# 보다 크다(greater than)
1 <= 2 	# 작거나 같다(lesser than or equal to)
1 >= 2	# 크거나 같다(greater than or equal to)
1 == 2	# 같다(equal)
1 != 2	# 같지 않다(different)

# 논리 연산 - 결과가 참(TRUE) 혹은 거짓(FALSE)
!FALSE		# 아님(logical NOT)
T & T			# 논리곱(logical AND)
T & F
T && T
T && F
T | F			# 논리합(logical OR)
T || F		
xor(T, T)		# 배타합(exclusive OR)
xor(F, F)
xor(T, F)

### 기본 데이터 타입의 복수개의 대상을 벡터로 만들 수 있음
# 벡터(Vector)
a=c(1,2,"A",3,4)	# 벡터의 요소는 모두 동일한 데이터 타입이여야 함(옆의 예제는 문자로 변환됨)
mode(a)
length(a)

### 매트릭스(Matrix)
m = matrix(1:12, 3,4)	# 열(ROW)3개, 행(Column) 4개
m
mode(m)
length(m)

### 배열(Array)
a = array(1, 10)		# 1을 길이가 10개까지 생성
a
a = array(1:3, 10)	# 1에서 3을 반복해서 길이가 10개까지
a

### 팩터(Factor) : 분류하기 위한 용도로 사용
trt = factor(rep(c("C", "T"), c(3,4))) # C를 3번, T를 4번 반복
str(trt)	# 변수의 구조(Structure) 보여 주기
summary(trt)

### 데이터 프레임(Data Frame) : 테이블 형태로 만들기
n = c(2,3,5)			# 숫자
s = c("aa", "bb", "cc")		# 문자
b = c(TRUE, FALSE, TRUE)	# 논리
df = data.frame(n,s,b)		# 데이터프레임
df
mode(df)				# 데이터 타입 알아내기
df[1:2,]				# 부분집합
fix(df)				# 데이터프레임 내용 수정하기

### 파일 읽어 데이터 프레임에 넣기
getwd()
setwd("c:/work/R/DataSets")
myData = read.table("R_Tutorial_Data.txt", header=TRUE, sep="\t")

### 데이터 프레임을 파일로 쓰고 다시 읽기
write.table(iris, file="c:/work/R/DataSets/iris.csv", sep=",", row.names=F, col.names=T, fileEncoding="UTF-8")
df = read.table(file="c:/work/R/DataSets/iris.csv", sep=",", header=T, encoding="UTF-8")

### 통계함수의 사용
attach(iris)
ssl = Sepal.Length[Species=="setosa"]
summary(ssl)		# 요약
mean(ssl)			# 평균 구하기
sd(ssl)			# 표준 편차
par(mfrow = c(1, 1))	# 그래프 그린 화면의 배치(한면)
hist(ssl, breaks=12)	# 히스토그램 그리기
plot(density(ssl), col="red")		# 분포도(확율 밀도) 그리기
boxplot(ssl, col="yellow")		# 박스 플랏 그리기

# 확률 밀도 함수
rnorm(10, mean=0, sd=1)	# 표준정규 분표, 평균 0, 표준편차 1에 대한 난수 10개 발생
dnorm(1.64, mean=0, sd=1)	# 표준정규 분표(평균 0, 표준편차 1)에서 1.64에 대한 확율밀도
pnorm(c(-1.64, 1.64, 2.56), mean=0, sd=1)	# 표준정규 분표(평균 0, 표준편차 1)에서 3개의 변수에 대한 누적밀도

# 통계적 추론 사례(mtcars 데이터셋 사용)
attach(mtcars)
cor(mpg, cyl)	# 상관계수 : 마력과 실린더 수(당근 비례함)
t.test(mpg, cyl)	# t 검증 수행 : 대립가설로 "두 집단의 평균의 차는 영이 아니다"에 대하여 95% 신뢰수준으로 검증(p-Value가 0.05 보다 작으면 기각)
var.test(mpg, cyl)
myModel = lm(mpg ~ cyl) # 선형회귀 분석
par(mfrow=c(2,2))
plot(myModel)

## 회귀 분석
attach(mtcars)
par(mfrow=c(1,1))
myModel = lm(mpg ~ cyl) # 선형회귀 분석
plot(cyl, mpg, col="red")	# 분포도 그리기
abline(myModel)		# 회귀 모델에 대한 직성 그리기

### ANOVA
attach(mtcars)
myANOVA=aov(mpg ~ cyl)
par(mfrow=c(2,2))
plot(myANOVA)
