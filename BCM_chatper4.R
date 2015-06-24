##############################
# 0. Cartoon guide to statistics
# http://urizen-geography.nsm.du.edu/~psutton/Sutton_Courses/Geog_2000_Intro_Geog_Stats/Geog2000_ppts/Gonick_Cartoon_Guide/
##############################
# [Def] Statistics deals with uncertainty 
# [Def] Statistics can quantify uncertainty
# Statisticsis supported by Data analysis, Probability, Statistical inference
# Creating a Graph
attach(mtcars)
plot(wt, mpg) 
abline(lm(mpg~wt)) # Add Straight Lines to a Plot & Fitting Linear Models
title("Regression of MPG on Weight")

? lm
# Simple Dotplot
par(mfrow=c(1,1))
dotchart(mtcars$mpg,labels=row.names(mtcars),cex=.7,
         main="Gas Milage for Car Models", 
         xlab="Miles Per Gallon")

# Line Chart
x <- c(1:5); y <- x # create some data 
par(pch=22, col="red") # plotting symbol and color 
par(mfrow=c(2,4)) # all plots on one page 
opts = c("p","l","o","b","c","s","S","h") 
for(i in 1:length(opts)){ 
  heading = paste("type=",opts[i]) 
  plot(x, y, type="n", main=heading) 
  lines(x, y, type=opts[i]) 
}

# Simple Pie Chart
par(mfrow=c(1,1))
slices <- c(10, 12,4, 16, 8)
lbls <- c("US", "UK", "Australia", "Germany", "France")
pie(slices, labels = lbls, main="Pie Chart of Countries")

# Colored Histogram with Different Number of Bins(12 breaks)
hist(mtcars$mpg)
hist(mtcars$mpg, breaks=12, col="red")

# Boxplot of MPG by Car Cylinders 
boxplot(mpg~cyl,data=mtcars, main="Car Milage Data", 
        xlab="Number of Cylinders", ylab="Miles Per Gallon")

# Simple Bar Plot 
counts <- table(mtcars$gear)
barplot(counts, main="Car Distribution", 
        xlab="Number of Gears", ylab="Count")

# Simple Horizontal Bar Plot with Added Labels 
counts <- table(mtcars$gear)
barplot(counts, main="Car Distribution", horiz=TRUE,
        names.arg=c("3 Gears", "4 Gears", "5 Gears"))

# Stacked Bar Plot with Colors and Legend(vs means V engine or a straight engine)
counts <- table(mtcars$vs, mtcars$gear)
barplot(counts, main="Car Distribution by Gears and VS",
        xlab="Number of Gears", col=c("darkblue","red"),
        legend = c("V Engine","Straight"))

# Grouped Bar Plot
counts <- table(mtcars$vs, mtcars$gear)
barplot(counts, main="Car Distribution by Gears and VS",
        xlab="Number of Gears", col=c("darkblue","red"),
        legend = c("V Engine","Straight"), beside=TRUE)
##############################
# 4 Inference with Gaussians
##############################
# http://www.r-tutor.com/elementary-statistics/probability-distributions
help("distribution")
# Uniform Distribution
# p(x) = choose(n, x) p^x (1-p)^(n-x)
# choose(n, x) = n!/x!(n-x)!
theta = rbinom(1000,1000,1/2)
#theta = dbinom(1:100, 100, 1/2)
summary(theta)
par(mfrow=c(1,2))
plot(theta, main='Uniform trace',xlab='iteration', col='blue', lwd=1)
plot(density(theta), main='Uniform distribution',xlab='theta', col='skyblue', lwd=2)

# Beta Distribution
#theta = dbeta(1:100,1,1)
theta = rbeta(1000,1,1)
summary(theta)
par(mfrow=c(1,2))
plot(theta, main='Beta trace',xlab='iteration', col='blue', lwd=1)
plot(density(theta), main='Beta distribution',xlab='theta', col='skyblue', lwd=2)

# Normal Distribution
# f(x) = 1/(√(2 π) σ) e^-((x - μ)^2/(2 σ^2)) 
theta = rnorm(n=1000, mean=0, sd=5)   # 데이터 생성
summary(theta)
par(mfrow=c(1,2))
plot(theta, main='Normal trace',xlab='iteration', col='blue', lwd=1)
plot(density(theta), main='Normal distribution',xlab='theta', col='skyblue', lwd=2)
quantile(theta, c(0.15,0.85))
# Normal Distribution, mean & standard deviation
pnorm(5, mean=0, sd=5) - pnorm(-5, mean=0, sd=5) # 68% of the values are within 1 standard deviation 
pnorm(10, mean=0, sd=5) - pnorm(-10, mean=0, sd=5) # 95% of the values are within 2 standard deviation 
pnorm(15, mean=0, sd=5) - pnorm(-15, mean=0, sd=5) # 99.7% of the values are within 3 standard deviation 

? pnorm

# 평균키가 175이고 표준편차가 5인 사람들에 대한 변수 생성하고 Graph 
height=rnorm(n=10000,mean = 175, sd=5)   # 데이터 생성
plot(density(height))
hist(height, breaks=20, probability=TRUE)    # 사진1
hist(height, breaks=100, probability=TRUE)   # 사진2
hist(height, breaks=300, probability=TRUE)   # 더욱 더 확률밀도 함수에 가까운 히스토그램

###########################################################################
# Reference
# Variational Bayesian Multinomial Probit Regression
# http://www.r-tutor.com/gpu-computing/gaussian-process/rvbm
# http://www.bioconductor.org/packages/release/bioc/vignettes/vbmp/inst/doc/vbmp.pdf
source("http://bioconductor.org/biocLite.R")
biocLite("vbmp")
