#----------------
# <I-5>
#----------------
mean(abs(rnorm(100))) 
rnorn(10)
setwd("d:\\Work")
pdf("aa.pdf")
hist(rnorm(100))
dev.off()
rnorm(10) # 표준정규 분포하에서 생성
data() # R에서 기본 제공하는 데이터
#----------------
# <I-8>
#----------------
help.start() # 도움말(공식 사이트)
help(seq) # seq 함수에 대한 도움말(수열 생성)
?seq # help(seq) 동일한 기능
help(abs)
RSiteSearch("lm") # Linear Model에 대하여

# History 기능
history()
setwd("C:\\Work")
savehistory(file = "history_150515.log")
loadhistory(file = "history_150515.log")

#----------------
# <I-8>
#----------------
help(rnorm) 
rnorm(10) # random number 생성( mean = 0, sd = 1)
mean(abs(rnorm(100))) # 
hist(rnorm(10)) # Histgram 표시
#----------
getwd() # Work Directory 가져옴
setwd("c:/Work") # Work Directory 설정
dir.create("c:/Work/Rtraining")
getwd()
setwd("c:/Work/Rtraining")
data() # R에서 기본적으로 제공하는 샘플데이터
BOD # R data에서 제공하는 생물학적 산소 요구량
BJsales # Sales Data with Leading Indicator
help(BJsales) 
options() # 현재 작업 환경의 Option들

#----------------
# <I-11>
#----------------
mtcars
lm(mpg~wt, data=mtcars) # mtcars 예제 데이터의 wt(무게)가 출력(mpg)에 영향까지 선형 분석
fit <- lm(mpg~wt, data=mtcars) # 선형 분석 결과를 fit 변수에 대입
str(fit) # 분석과정를 구조체로 변환
help(str) # fit 라는 변수의 저장된 값 표시


#----------------
# <I-13>
#----------------
install.packages("vcd") # Visualizing Categorical Data
help(pakage="vcd") # 데이터를 저장하고 있는 팩키지
library(vcd) # 라이브러리 명려을 통한 팩키지 사용
help(Arthritis) # 관절염 임상 데이터
Arthritis
example(Arthritis) # 사용예
