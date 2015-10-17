#세계 인구 데이터셋(2015)
library("wpp2015")
 
#평균수명
#e0F:여성, e0M:남성
#_supplemental, proj, proj80l, proj80u, proj95l, proj95u
#1950-2015 5년주기의 남성의 평균 수명
data(e0M)
data(e0F)
subset(e0M, e0M[, 1] == "Republic of Korea")

e0M[e0M$country == 'Japan',]
e0Mckj=e0M[e0M$country == 'Japan'|e0M$country == 'Republic of Korea'|e0M$country == 'China',c(1,15)]
ggplot(e0Mckj, aes(x=country, y=e0Mckj$2010-2015))+geom_point()

#한국 정보

e0M[,1]
subset(e0M, e0M[, 1] == "Republic of Korea")
 
#2015년이후 예상 평균 수명:e0Mproj
data(e0Mproj)
#한국 정보
subset(e0Mproj, e0Mproj[, 1] == "Republic of Korea")
 
#1950-2100 5년주기 이주자 수:migration
data(migration)
#한국 정보
subset(migration, migration[, 2] == "Republic of Korea")

#1950-2100 5년 주기 사망율
#mxF:여성, mxM:남성
data(mxM)
#한국 정보(나이별 정보)
subset(mxM, mxM[, 1] == "Republic of Korea"|mxM[, 1] == "Japan"|mxM[, 1] == "China")
 
#1950-2100 5년 주기 15~49세의 여성 출산할 활율
data(percentASFR)
#한국
subset(percentASFR, percentASFR[, 1] == "Republic of Korea"|percentASFR[, 1] == "Japan"|percentASFR[, 1] == "China")
 
#인구수
#pop:남녀합계, popF:여성, popM:남성
#proj, proj80l, proj80u, proj95l, proj95u, projHigh, projLow
data(pop)
#한국 정보
subset(pop, pop[, 2] == "Republic of Korea"|pop[, 2] == "Japan"|pop[, 2] == "China")

#1950-2100 5년간 남성과 여성 비율
data(sexRatio)
#한국
subset(sexRatio, sexRatio[, 1] == "Republic of Korea"|sexRatio[, 1] == "Japan"|sexRatio[, 1] == "China")
country country_code 1950-1955 1955-1960 1960-1965 1965-1970 1970-1975 1975-1980 1980-1985 1985-1990
Japan          392     1.056     1.056     1.056     1.056     1.056     1.056     1.056 
