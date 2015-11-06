# JAGS(Just Another Gibbs Sampler) 설치
# http://sourceforge.net/projects/mcmc-jags/files/latest/download

#install.packages("R2jags",dependencies=T)
#install.packages("rjags")
# 아래 메세지가 발생하고 실행 불가(Windows 7 32bit,R 3-2-2)
#  Package which is only available in source form,and may need
#  compilation of C/C++/Fortran: ‘rjags’
#  These will not be installed

#library(R2jags)
#library(rjags)

# Windows 7 32bit,R 3-2-2에서 JAGS 설치 문제가 있어 WinBUGS 사용함
#######################################################
# 12.1 Psychophysical functions
#######################################################
# 워크 스페이스의 변수 초기화
rm(list=ls()) 

# sets working directories:
setwd("D:/work/bayesian_cognitive_code/CaseStudies/PsychophysicalFunctions")
library(R2WinBUGS)
# WinBUGS 설치 디렉토리
bugsdir = "C:/Program Files/WinBUGS14"
# WinBUG가 미설치 된 경우 아래 사이트에서 받기
# http://www.mrc-bsu.cam.ac.uk/software/bugs/the-bugs-project-winbugs/

# 자극(소리)의 지속 시간
# 8개의 대상(사람(?))에 대한 28번의 시도(interval(간격) 주기)
data_x =
"200,220,240,260,270,275,280,285,290,295,300,305,310,330,335,340,345,350,355,360,365,370,375,380,385,390,400,NA
200,220,240,260,270,280,285,290,295,300,305,310,315,320,325,330,340,350,355,360,365,370,380,400,NA,NA,NA,NA
200,220,225,230,235,240,245,250,255,260,265,270,275,280,285,290,295,300,305,310,315,320,330,340,360,380,400,NA
200,220,240,260,270,275,280,285,290,295,300,305,310,315,320,325,330,335,340,360,380,400,NA,NA,NA,NA,NA,NA
200,220,240,250,255,260,265,270,275,280,285,290,295,300,305,310,315,320,325,330,335,340,360,380,400,NA,NA,NA
200,220,240,260,265,270,275,280,285,290,295,300,305,310,315,320,325,330,335,340,345,350,355,360,380,400,NA,NA
180,200,210,220,230,235,240,245,250,255,260,265,270,305,310,315,320,325,330,335,340,345,350,355,360,380,400,420
200,210,215,220,225,230,235,240,245,250,260,265,270,275,280,285,290,295,300,305,310,315,320,330,340,360,380,400"
x       = matrix(NA,nrow=8,ncol=28)
#x[]    = as.matrix(read.table("data_x.txt",sep="\t"))
x[]     = as.matrix(read.table(file=textConnection(data_x),sep=","))

# 자극(소리)의 차이 식별
# 8개의 대상(사람(?))에 대한 28번의 시도(interval(간격) 주기)에 대한 식별
data_n=
"6,6,6,6,6,4,16,19,17,9,9,11,5,5,6,24,13,16,4,9,6,8,6,11,2,2,8,NA
6,6,6,6,6,12,12,11,13,10,6,10,15,19,11,5,9,7,8,19,16,12,7,8,NA,NA,NA,NA
8,15,4,9,12,20,7,11,5,7,3,3,10,16,11,9,7,10,14,15,8,6,6,6,6,6,6,NA
6,6,6,6,13,13,15,11,12,13,12,6,3,12,30,24,21,4,6,6,6,9,NA,NA,NA,NA,NA,NA
6,6,6,8,3,14,4,13,10,17,14,9,9,7,8,10,16,23,14,17,2,6,6,6,6,NA,NA,NA
6,6,6,12,9,13,7,18,22,10,5,8,12,12,13,12,3,12,4,20,5,5,2,6,6,6,NA,NA
2,9,2,7,5,4,18,10,26,8,20,2,7,5,12,6,13,10,16,4,16,9,5,2,6,6,8,2
8,2,2,19,7,31,15,25,3,5,3,2,4,7,9,10,14,9,9,8,6,6,6,6,6,6,6,6"
n       = matrix(NA,8,28)
#n[]    = as.matrix(read.table("data_n.txt",sep="\t"))
n[]     = as.matrix(read.table(file=textConnection(data_n),sep=","))

# 테스트 대상자가 테스트 수행동안 지속간격을 분류하는 횟수
data_r = 
"0,0,0,0,0,0,2,4,4,3,1,2,3,0,3,17,9,13,4,6,4,6,5,8,2,2,8,NA
0,0,0,0,0,0,3,1,2,3,0,2,6,10,7,3,7,4,6,11,13,12,6,8,NA,NA,NA,NA
0,2,1,3,2,8,2,4,2,3,2,1,3,11,9,8,6,8,10,14,8,6,6,6,6,6,6,NA
0,0,0,0,0,4,5,2,2,4,3,4,1,4,22,17,19,4,6,6,6,6,NA,NA,NA,NA,NA,NA
0,0,0,0,1,4,1,1,3,5,4,3,3,4,4,8,10,19,9,16,2,6,6,6,6,NA,NA,NA
0,0,0,0,3,3,2,3,5,5,2,3,8,7,10,11,2,7,2,17,4,4,2,6,6,6,NA,NA
0,1,0,1,0,1,4,2,10,5,6,1,4,1,11,4,9,5,13,1,12,7,4,2,6,6,7,2
0,0,0,3,1,6,4,11,2,4,2,2,4,5,9,7,13,8,8,8,6,6,6,6,6,6,6,6"
r       = matrix(NA,8,28)
#r[]    = as.matrix(read.table("data_r.txt",sep="\t"))
r[]     = as.matrix(read.table(file=textConnection(data_r),sep=","))

# 장기 지속적 반응의 비율
data_rprop = 
"0,0,0,0,0,0,0.125,0.210526315789474,0.235294117647059,0.333333333333333,0.111111111111111,0.181818181818182,0.6,0,0.5,0.708333333333333,0.692307692307692,0.8125,1,0.666666666666667,0.666666666666667,0.75,0.833333333333333,0.727272727272727,1,1,1,NA
0,0,0,0,0,0,0.25,0.090909090909091,0.153846153846154,0.3,0,0.2,0.4,0.526315789473684,0.636363636363636,0.6,0.777777777777778,0.571428571428571,0.75,0.578947368421053,0.8125,1,0.857142857142857,1,NA,NA,NA,NA
0,0.133333333333333,0.25,0.333333333333333,0.166666666666667,0.4,0.285714285714286,0.363636363636364,0.4,0.428571428571429,0.666666666666667,0.333333333333333,0.3,0.6875,0.818181818181818,0.888888888888889,0.857142857142857,0.8,0.714285714285714,0.933333333333333,1,1,1,1,1,1,1,NA
0,0,0,0,0,0.307692307692308,0.333333333333333,0.181818181818182,0.166666666666667,0.307692307692308,0.25,0.666666666666667,0.333333333333333,0.333333333333333,0.733333333333333,0.708333333333333,0.904761904761905,1,1,1,1,0.666666666666667,NA,NA,NA,NA,NA,NA
0,0,0,0,0.333333333333333,0.285714285714286,0.25,0.0769230769230769,0.3,0.294117647058824,0.285714285714286,0.333333333333333,0.333333333333333,0.571428571428571,0.5,0.8,0.625,0.826086956521739,0.642857142857143,0.941176470588235,1,1,1,1,1,NA,NA,NA
0,0,0,0,0.333333333333333,0.230769230769231,0.285714285714286,0.166666666666667,0.227272727272727,0.5,0.4,0.375,0.666666666666667,0.583333333333333,0.769230769230769,0.916666666666667,0.666666666666667,0.583333333333333,0.5,0.85,0.8,0.8,1,1,1,1,NA,NA
0,0.111111111111111,0,0.142857142857143,0,0.25,0.222222222222222,0.2,0.384615384615385,0.625,0.3,0.5,0.571428571428571,0.2,0.916666666666667,0.666666666666667,0.692307692307692,0.5,0.8125,0.25,0.75,0.777777777777778,0.8,1,1,1,0.875,1
0,0,0,0.157894736842105,0.142857142857143,0.193548387096774,0.266666666666667,0.44,0.666666666666667,0.8,0.666666666666667,1,1,0.714285714285714,1,0.7,0.928571428571429,0.888888888888889,0.888888888888889,1,1,1,1,1,1,1,1,1"
rprop   <- matrix(NA,8,28)
#rprop[] <- as.matrix(read.table("data_rprop.txt",sep="\t"))
rprop[] <- as.matrix(read.table(file=textConnection(data_rprop),sep=","))

###################
library(ggplot2)
# 자극(소리)의 지속 시간(x)와  자극(소리)의 차이 식별(n)
for (i in 1:8) {
	xv = x[i,]
	nv = n[i,]
	x_n=data.frame(xv ,nv)
	print(ggplot(data=x_n,aes(x=xv,y=nv)) + geom_point())	
	par(ask=TRUE)
}
# 자극(소리)의 지속 시간(x)와 지속간격을 분류하는 횟수(r)
for (i in 1:8) {
	xv = x[i,]
	rv = n[i,]
	x_r=data.frame(xv ,rv)
	print(ggplot(data=x_r,aes(x=xv,y=rv)) + geom_point())	
	par(ask=TRUE)
}
# 자극(소리)의 지속 시간(x)와  장기 지속적 반응의 비율(rprop)
for (i in 1:8) {
	xv = x[i,]
	rpropv = rprop[i,]
	x_rprop=data.frame(xv ,rpropv)
	print(ggplot(data=x_rprop,aes(x=xv,y=rpropv)) + geom_point())
	par(ask=TRUE)	
}
for (i in 1:8) {
	xmean[i] = mean(x[i,],na.rm=T)
}
# xmean = c(318.888,311.0417,284.4444,301.5909,296.2000,305.7692,294.6429,280.3571)
nstim =c(27,24,27,22,25,26,28,28)
nsubjs = 8

data = list("x","xmean","n","r","nsubjs","nstim") # to be passed on to WinBUGS
# alpha의 경우 Uniform 분포(-2, 2)의 8개의 난수 발생 및 대음
# beta의 경우  Uniform 분포(0, 0.5)의 8개의 난수 발생 및 대음
myinits <-	list(
  list(alpha = runif(nsubjs,-2,2),beta = runif(nsubjs,0,.5)),  
  list(alpha = runif(nsubjs,-2,2),beta = runif(nsubjs,0,.5)),
  list(alpha = runif(nsubjs,-2,2),beta = runif(nsubjs,0,.5)))

################################################
# Part for model without contamination parameter
################################################

# 관측할 파라미터:	
parameters <- c("alpha","beta")

# 모델
model = "
# Logistic Psychophysical Function
model{
  for (i in 1:nsubjs){
    for (j in 1:nstim[i]){
      r[i,j] ~ dbin(thetalim[i,j],n[i,j])
      logit(thetalim[i,j]) <- lthetalim[i,j]
      lthetalim[i,j] <- min(999,max(-999,ltheta[i,j]))
      ltheta[i,j] <- alpha[i]+beta[i]*(x[i,j]-xmean[i])
    }
    beta[i] ~ dnorm(mub,lambdab)
    alpha[i] ~ dnorm(mua,lambdaa)
  }
  # Priors
  mub ~ dnorm(0,.001)    
  mua ~ dnorm(0,.001)
  sigmab ~ dunif(0,1000)
  sigmaa ~ dunif(0,1000)
  lambdab <- pow(sigmab,-2)
  lambdaa <- pow(sigmaa,-2)
}"

# The following command calls WinBUGS with specific options.
# For a detailed description see Sturtz,Ligges,& Gelman (2005).
samples <- bugs(data,inits=myinits,parameters,
	 			model.file ="Psychophysical_1.txt",
	 			n.chains=3,n.iter=10000,n.burnin=5000,n.thin=1,
	 			DIC=T,bugs.directory=bugsdir,
	 			codaPkg=F,debug=T)

# Now the values for the monitored parameters are in the "samples" object,
# ready for inspection.

# Extracting the parameters
alpha	    <- samples$sims.list$alpha
beta  	  <- samples$sims.list$beta
alphaMAP  <- c(rep(0,nsubjs))
betaMAP   <- c(rep(0,nsubjs))
alpha_sel <- matrix(NA,20,8) 
beta_sel  <- matrix(NA,20,8) 

# Constructing MAP-estimates and alpha/beta range
for (i in 1:nsubjs)
{
	alphaMAP[i]   <- density(alpha[,i])$x[which(density(alpha[,i])$y==max(density(alpha[,i])$y))]
	betaMAP[i]    <- density(beta[,i])$x[which(density(beta[,i])$y==max(density(beta[,i])$y))]
	alpha_sel[,i] <- sample(alpha[,i],20)
	beta_sel[,i]  <- sample(beta[,i],20)
}

############################## PSYCHOMETRIC FUNCTIONS ##############################

F1 <- function(X,s) # only the MAP estimate; use this to plot psychometric functions
{
  exp(alphaMAP[s] + betaMAP[s]*(X - xmean[s]))/(1+exp(alphaMAP[s] + betaMAP[s]*(X - xmean[s])))
}

F1inv <- function(Y,s)
{
  (log(-Y/(Y-1))-alphaMAP[s])/betaMAP[s]
}

F2 <- function(X,s) # function for all the posterior alpha/beta values; use this to calculate JND posterior
{
  exp(alpha[,s] + beta[,s]*(X - xmean[s]))/(1+exp(alpha[,s] + beta[,s]*(X - xmean[s])))
}
F2inv <- function(Y,s)
{
  (log(-Y/(Y-1))-alpha[,s])/beta[,s]
}

F3 <- function(X,s,g) # function for 20 grabbed posterior alpha/beta values; use this to plot overlapping sigmoids to visualize variance
{
  exp(alpha_sel[g,s] + beta_sel[g,s]*(X - xmean[s]))/(1+exp(alpha_sel[g,s] + beta_sel[g,s]*(X - xmean[s])))
}

##################################### JND/PSE calculation ########################################
JND 	 <- F2inv(0.84,c(1:nsubjs))-F2inv(0.5,c(1:nsubjs))
JNDmap <- F1inv(0.84,c(1:nsubjs))-F1inv(0.5,c(1:nsubjs))
								  				             
PSE 	 <- F2inv(0.5,c(1:nsubjs))+xmean
PSEmap <- F1inv(0.5,c(1:nsubjs))+xmean
		
################## PLOTS ####################

### Figure 12.2

dev.new(width=10,height=5)
layout(matrix(1:nsubjs,2,4,byrow=T))
par(mar=c(1,2,2,0),oma=c(5,5,1,1))
for (i in 1:nsubjs)
{
	scale <- seq(x[i,1],x[i,nstim[i]],by=.1)
	plot(x[i,],rprop[i,],main=paste("Subject",as.character(i)),xlab="",ylab="",pch=15,col="dark grey",ylim=c(0,1),yaxt="n",xaxt="n")
	lines(scale,F1(scale,i),type="l")
	segments(x0=x[i,1],x1=PSEmap[i]+JNDmap[i],y0=0.84,lty=2)
	segments(x0=x[i,1],x1=PSEmap[i],y0=0.5,lty=2)
	segments(y0=0,y1=0.84,x0=PSEmap[i]+JNDmap[i],lty=2)
	segments(y0=0,y1=0.5,x0=PSEmap[i],lty=2)
	if (i==1 | i==5) 
  {
		axis(2,las=1,yaxp=c(0,1,2))
		axis(2,at=0.84,las=1)
	}
	if (i>4) axis(1)
}
mtext("Proportion 'Long' Response",side=2,line=2,outer=T,cex=1.4)
mtext("Test Interval (ms)",side=1,outer=T,line=3,cex=1.4)

### WARNING: Do not close R window.

### NOTE:	 Answers to the exercises can be found in PsychometricFunction1_Answers.R
