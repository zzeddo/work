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

# Probability
#Mutually exclusive ?? not overlap
#    P(E OR F) = P(E) + P(F)
#Overlap of elementary outcomes
#    P(E OR F) = P(E) + P(F) ?? P(E AND F)
#When P(NOT E) is easier to compute use subtraction rule.
#    P(E) = 1 ?? P(NOT E)
#Conditional Probability
#  The ?쐏robability of A, given C??
#     P(E|F) = P(E AND F)/P(F)
#  When E and F are mutually exclusive
#     P(E|F) = 0, once F has occurred E is impossible
# Rearranging the definition get multiplication rule
#     P(E AND F) = P(E|F)P(F)
# Independence
# Two events E and F are independent of each other if the occurrence of one had no influence on the probability of the other.
#   P(E AND F) = P(E)P(F)


# http://www.r-tutor.com/elementary-statistics/probability-distributions
help("distribution")
# Uniform Distribution
# f(x) = 1/(max-min)
theta = runif(1000, min=0, max=10)
summary(theta)
par(mfrow=c(1,2))
plot(theta, main='Uniform trace',xlab='iteration', col='blue', lwd=1)
plot(density(theta), main='Uniform distribution',xlab='theta', col='skyblue', lwd=2)


# Binomial Distribution
# p(x) = choose(n, x) p^x (1-p)^(n-x)
# choose(n, x) = n!/x!(n-x)!
theta = rbinom(1000,1000,1/2)
#theta = dbinom(1:100, 100, 1/2)
summary(theta)
par(mfrow=c(1,2))
plot(theta, main='Uniform trace',xlab='iteration', col='blue', lwd=1)
plot(density(theta), main='Uniform distribution',xlab='theta', col='skyblue', lwd=2)

# Normal Distribution
# f(x) = 1/(?닖(2 ?) ??) e^-((x - 關)^2/(2 ??^2)) 
theta = rnorm(n=1000, mean=0, sd=5)   # ?뜲?씠?꽣 ?깮?꽦
summary(theta)
par(mfrow=c(1,2))
plot(theta, main='Normal trace',xlab='iteration', col='blue', lwd=1)
plot(density(theta), main='Normal distribution',xlab='theta', col='skyblue', lwd=2)
quantile(theta, c(0.15,0.85))
# Normal Distribution, mean & standard deviation
pnorm(5, mean=0, sd=5) - pnorm(-5, mean=0, sd=5) # 68% of the values are within 1 standard deviation 
pnorm(10, mean=0, sd=5) - pnorm(-10, mean=0, sd=5) # 95% of the values are within 2 standard deviation 
pnorm(15, mean=0, sd=5) - pnorm(-15, mean=0, sd=5) # 99.7% of the values are within 3 standard deviation 


# ?룊洹좏궎媛 175?씠怨? ?몴以?렪李④? 5?씤 ?궗?엺?뱾?뿉 ???븳 蹂?닔 ?깮?꽦?븯怨? Graph 
height=rnorm(n=10000,mean = 175, sd=5)   # ?뜲?씠?꽣 ?깮?꽦
plot(density(height))
hist(height, breaks=20, probability=TRUE)    # ?궗吏?1
hist(height, breaks=100, probability=TRUE)   # ?궗吏?2
hist(height, breaks=300, probability=TRUE)   # ?뜑?슧 ?뜑 ?솗瑜좊??룄 ?븿?닔?뿉 媛源뚯슫 ?엳?뒪?넗洹몃옩


# Beta Distribution
#theta = dbeta(1:100,1,1)
theta = rbeta(1000,1,1)
summary(theta)
par(mfrow=c(1,2))
plot(theta, main='Beta trace',xlab='iteration', col='blue', lwd=1)
plot(density(theta), main='Beta distribution',xlab='theta', col='skyblue', lwd=2)

# The Poisson Distribution
# ?젙?빐吏? ?떆媛? ?븞?뿉 ?뼱?뼡 ?궗嫄댁씠 ?씪?뼱?궇 ?슏?닔?뿉 ???븳 湲곕뙎媛믪쓣 貫(lambda)?씪怨? ?뻽?쓣 ?븣, 洹? ?궗嫄댁씠 x?쉶 ?씪?뼱?궇 ?솗瑜?
# p(x) = 貫^x exp(-貫)/x!
# for x = 0, 1, 2, ?? . The mean and variance are E(X) = Var(X) = 貫. 
#theta = rpois(n=1000, lambda=10)
theta = dpois(0:10, lambda = 1)
theta[3]
summary(theta)
par(mfrow=c(1,2))
plot(theta, main='Poisson trace',xlab='iteration', col='blue', lwd=1)
plot(density(theta), main='Poisson distribution',xlab='theta', col='skyblue', lwd=2)

# Gamma Distribution
#  ?뿰?냽 ?솗瑜좊텇?룷濡?, ?몢 媛쒖쓽 留ㅺ컻蹂?닔瑜? 諛쏆쑝硫? ?뼇?쓽 ?떎?닔瑜? 媛吏? ?닔 ?엳?떎.
? dgamma
theta = rgamma(1000, shape=0, scale=0)
summary(theta)
par(mfrow=c(1,2))
plot(theta, main='Gamma trace',xlab='iteration', col='blue', lwd=1)
plot(density(theta), main='Gamma distribution',xlab='theta', col='skyblue', lwd=2)

##############################
# 4.1 Inference a mean and standard deviation
##############################
# clears workspace:  
rm(list=ls()) 

#JARS
library(R2jags)

# Prior and Posterior Prediction
# 1. Model
modelString = "
# Inferring the Mean and Standard Deviation of a Gaussian
model{
# Data Come From A Gaussian
for (i in 1:n){
x[i] ~ dnorm(mu,lambda)
}
# Priors
mu ~ dnorm(0,0.001)
sigma ~ dunif(0,10)
lambda <- 1/pow(sigma,2)
}"

# 2. Data
x = c(1.1, 1.9, 2.3, 1.8)
n = length(x)
mydata = list("x", "n")

# 3. Start Values
myinits = list(
  list(mu =0, sigma = 1)
)

# 4. Paramters to be monitored
parameters = c("mu", "sigma","x")

# The following command calls JAGS with specific options.
# For a detailed description see the R2jags documentation.
mcmc_samples <- jags(mydata, inits=myinits, parameters,
                     model.file =textConnection(modelString), n.chains=1, n.iter=1000, 
                     n.burnin=1, n.thin=1, DIC=T)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.
x <- mcmc_samples$BUGSoutput$sims.list$x
mu <- mcmc_samples$BUGSoutput$sims.list$mu
sigma <- mcmc_samples$BUGSoutput$sims.list$sigma
summary(mu)
summary(sigma)

# Graph
par(mfrow=c(2,2))
plot(mu, main='mean trace',xlab='iteration', col='blue', lwd=1)
plot(density(mu), main='mean distribution',xlab='mu', col='skyblue', lwd=2)
plot(sigma, main='standard deviation trace',xlab='iteration', col='blue', lwd=1)
plot(density(sigma), main='standard deviation distribution',xlab='sigma', col='skyblue', lwd=2)

##############################
# 4.1 Inference a mean and standard deviation
# Exercise 4.1.1  mu ~ dnorm(0,0.001) -> mu ~ dnorm(0,1)
##############################
# clears workspace:  
rm(list=ls()) 

#JARS
library(R2jags)

# Prior and Posterior Prediction
# 1. Model
modelString = "
# Inferring the Mean and Standard Deviation of a Gaussian
model{
# Data Come From A Gaussian
for (i in 1:n){
x[i] ~ dnorm(mu,lambda)
}
# Priors
mu ~ dnorm(0,1)
sigma ~ dunif(0,10)
lambda <- 1/pow(sigma,2)
}"

# 2. Data
x = c(1.1, 1.9, 2.3, 1.8)
n = length(x)
mydata = list("x", "n")

# 3. Start Values
myinits = list(
  list(mu =0, sigma = 1)
)

# 4. Paramters to be monitored
parameters = c("mu", "sigma","x")

# The following command calls JAGS with specific options.
# For a detailed description see the R2jags documentation.
mcmc_samples <- jags(mydata, inits=myinits, parameters,
                     model.file =textConnection(modelString), n.chains=1, n.iter=1000, 
                     n.burnin=1, n.thin=1, DIC=T)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.
x <- mcmc_samples$BUGSoutput$sims.list$x
mu <- mcmc_samples$BUGSoutput$sims.list$mu
sigma <- mcmc_samples$BUGSoutput$sims.list$sigma
summary(mu)
summary(sigma)

# Graph
par(mfrow=c(2,2))
plot(mu, main='mean trace',xlab='iteration', col='blue', lwd=1)
plot(density(mu), main='mean distribution',xlab='mu', col='skyblue', lwd=2)
plot(sigma, main='standard deviation trace',xlab='iteration', col='blue', lwd=1)
plot(density(sigma), main='standard deviation distribution',xlab='sigma', col='skyblue', lwd=2)

##############################
# 4.1 Inference a mean and standard deviation
# Exercise 4.1.3  lambda <- 1/pow(sigma,2) -> lambda <- 1
##############################
# clears workspace:  
rm(list=ls()) 

#JARS
library(R2jags)

# Prior and Posterior Prediction
# 1. Model
modelString = "
# Inferring the Mean and Standard Deviation of a Gaussian
model{
# Data Come From A Gaussian
for (i in 1:n){
x[i] ~ dnorm(mu,lambda)
}
# Priors
mu ~ dnorm(0,1)
sigma ~ dunif(0,10)
lambda <- 1
}"

# 2. Data
x = c(1.1, 1.9, 2.3, 1.8)
n = length(x)
mydata = list("x", "n")

# 3. Start Values
myinits = list(
  list(mu =0, sigma = 1)
)

# 4. Paramters to be monitored
parameters = c("mu", "sigma","x")

# The following command calls JAGS with specific options.
# For a detailed description see the R2jags documentation.
mcmc_samples <- jags(mydata, inits=myinits, parameters,
                     model.file =textConnection(modelString), n.chains=1, n.iter=1000, 
                     n.burnin=1, n.thin=1, DIC=T)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.
x <- mcmc_samples$BUGSoutput$sims.list$x
mu <- mcmc_samples$BUGSoutput$sims.list$mu
sigma <- mcmc_samples$BUGSoutput$sims.list$sigma
summary(mu)
summary(sigma)

# Graph
par(mfrow=c(2,2))
plot(mu, main='mean trace',xlab='iteration', col='blue', lwd=1)
plot(density(mu), main='mean distribution',xlab='mu', col='skyblue', lwd=2)
plot(sigma, main='standard deviation trace',xlab='iteration', col='blue', lwd=1)
plot(density(sigma), main='standard deviation distribution',xlab='sigma', col='skyblue', lwd=2)

##############################
# 4.1 Inference a mean and standard deviation
# Exercise 4.1.4  dnorm(mu,lambda) -> dnorm(1, lambda)
##############################
# clears workspace:  
rm(list=ls()) 

#JARS
library(R2jags)

# Prior and Posterior Prediction
# 1. Model
modelString = "
# Inferring the Mean and Standard Deviation of a Gaussian
model{
# Data Come From A Gaussian
for (i in 1:n){
x[i] ~ dnorm(1,lambda)
}
# Priors
mu ~ dnorm(0,1)
sigma ~ dunif(0,10)
lambda <- 1/pow(sigma, 2)
}"

# 2. Data
x = c(1.1, 1.9, 2.3, 1.8)
n = length(x)
mydata = list("x", "n")

# 3. Start Values
myinits = list(
  list(mu =0, sigma = 1)
)

# 4. Paramters to be monitored
parameters = c("mu", "sigma","x")

# The following command calls JAGS with specific options.
# For a detailed description see the R2jags documentation.
mcmc_samples <- jags(mydata, inits=myinits, parameters,
                     model.file =textConnection(modelString), n.chains=1, n.iter=1000, 
                     n.burnin=1, n.thin=1, DIC=T)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.
x <- mcmc_samples$BUGSoutput$sims.list$x
mu <- mcmc_samples$BUGSoutput$sims.list$mu
sigma <- mcmc_samples$BUGSoutput$sims.list$sigma
summary(mu)
summary(sigma)

# Graph
par(mfrow=c(2,2))
plot(mu, main='mean trace',xlab='iteration', col='blue', lwd=1)
plot(density(mu), main='mean distribution',xlab='mu', col='skyblue', lwd=2)
plot(sigma, main='standard deviation trace',xlab='iteration', col='blue', lwd=1)
plot(density(sigma), main='standard deviation distribution',xlab='sigma', col='skyblue', lwd=2)

##############################
# 4.2 The seven scientiests
# lambda[i] ~ dgamma(.001,.001)
##############################
# clears workspace:  
rm(list=ls()) 

#JARS
library(R2jags)

# Prior and Posterior Prediction
# 1. Model
modelString = "
# The Seven Scientists
model{
  # Data Come From Gaussians With Common Mean But Different Precisions
  for (i in 1:n){
    x[i] ~ dnorm(mu,lambda[i])
  }
  # Priors
  mu ~ dnorm(0,0.001)
  for (i in 1:n){
    lambda[i] ~ dgamma(.001,.001)
    sigma[i] <- 1/sqrt(lambda[i])  
  }     
}"

# 2. Data
x <- c(-27.020,3.570,8.191,9.898,9.603,9.945,10.056)
n <- length(x)

mydata = list("x", "n")

# 3. Start Values
myinits = list(
  list(mu =0, lambda = rep(1,n))
)

# 4. Paramters to be monitored
parameters = c("mu", "sigma","x")

# The following command calls JAGS with specific options.
# For a detailed description see the R2jags documentation.
mcmc_samples <- jags(mydata, inits=myinits, parameters,
                     model.file =textConnection(modelString), n.chains=1, n.iter=1000, 
                     n.burnin=1, n.thin=1, DIC=T)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.
x <- mcmc_samples$BUGSoutput$sims.list$x
mu <- mcmc_samples$BUGSoutput$sims.list$mu
sigma <- mcmc_samples$BUGSoutput$sims.list$sigma
summary(mu)
summary(sigma)

# Graph
par(mfrow=c(3,2))
plot(mu, main='mean trace',xlab='iteration', col='blue', lwd=1)
plot(density(mu), main='mean distribution',xlab='mu', col='skyblue', lwd=2)
plot(sigma[,1], main='standard deviation 1 trace',xlab='iteration', col='blue', lwd=1)
plot(sigma[,2], main='standard deviation 2 trace',xlab='iteration', col='blue', lwd=1)
plot(sigma[,3], main='standard deviation 3 trace',xlab='iteration', col='blue', lwd=1)
plot(density(sigma), main='standard deviation distribution',xlab='sigma', col='skyblue', lwd=2)

##############################
# 4.2 The seven scientiests
# lambda[i] ~ dgamma(1, 1)
##############################
# clears workspace:  
rm(list=ls()) 
sqrt(1/0.001)
#JARS
library(R2jags)

# Prior and Posterior Prediction
# 1. Model
modelString = "
# The Seven Scientists
model{
# Data Come From Gaussians With Common Mean But Different Precisions
for (i in 1:n){
x[i] ~ dnorm(mu,lambda[i])
}
# Priors
mu ~ dnorm(0,0.001)
for (i in 1:n){
lambda[i] ~ dgamma(1, 1)
sigma[i] <- 1/sqrt(lambda[i])  
}     
}"

# 2. Data
x <- c(-27.020,3.570,8.191,9.898,9.603,9.945,10.056)
n <- length(x)

mydata = list("x", "n")

# 3. Start Values
myinits = list(
  list(mu =0, lambda = rep(1,n))
)

# 4. Paramters to be monitored
parameters = c("mu", "sigma","x")

# The following command calls JAGS with specific options.
# For a detailed description see the R2jags documentation.
mcmc_samples <- jags(mydata, inits=myinits, parameters,
                     model.file =textConnection(modelString), n.chains=1, n.iter=1000, 
                     n.burnin=1, n.thin=1, DIC=T)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.
x <- mcmc_samples$BUGSoutput$sims.list$x
mu <- mcmc_samples$BUGSoutput$sims.list$mu
sigma <- mcmc_samples$BUGSoutput$sims.list$sigma
summary(mu)
summary(sigma)

# Graph
par(mfrow=c(3,2))
plot(mu, main='mean trace',xlab='iteration', col='blue', lwd=1)
plot(density(mu), main='mean distribution',xlab='mu', col='skyblue', lwd=2)
plot(sigma[,1], main='standard deviation 1 trace',xlab='iteration', col='blue', lwd=1)
plot(sigma[,2], main='standard deviation 2 trace',xlab='iteration', col='blue', lwd=1)
plot(sigma[,3], main='standard deviation 3 trace',xlab='iteration', col='blue', lwd=1)
plot(density(sigma), main='standard deviation distribution',xlab='sigma', col='skyblue', lwd=2)

##############################
# 4.1 The seven scientiests
# 4.2.1 Draw posterior samples
##############################
# clears workspace:  
rm(list=ls()) 

#JARS
library(R2jags)

# Prior and Posterior Prediction
# 1. Model
modelString = "
# The Seven Scientists
model{
# Data Come From Gaussians With Common Mean But Different Precisions
for (i in 1:n){
x[i] ~ dnorm(mu,lambda[i])
}
# Priors
mu ~ dnorm(0,0.001)
for (i in 1:n){
lambda[i] ~ dgamma(.001,.001)
sigma[i] <- 1/sqrt(lambda[i])  
}     
}"

# 2. Data
x <- c(-27.020,3.570,8.191,9.898,9.603,9.945,10.056)
n <- length(x)

mydata = list("x", "n")

# 3. Start Values
myinits = list(
  list(mu =0, lambda = rep(1,n))
)

# 4. Paramters to be monitored
parameters = c("mu", "sigma","lambda", "x")

# The following command calls JAGS with specific options.
# For a detailed description see the R2jags documentation.
mcmc_samples <- jags(mydata, inits=myinits, parameters,
                     model.file =textConnection(modelString), n.chains=1, n.iter=1000, 
                     n.burnin=1, n.thin=1, DIC=T)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.
x <- mcmc_samples$BUGSoutput$sims.list$x
mu <- mcmc_samples$BUGSoutput$sims.list$mu
sigma <- mcmc_samples$BUGSoutput$sims.list$sigma
lambda <- mcmc_samples$BUGSoutput$sims.list$lambda
summary(mu)
summary(sigma)

# Graph
par(mfrow=c(3,2))
plot(sigma[,1], main='sd 1 trace',xlab='iteration', col='blue', lwd=1)
plot(density(sigma[,1]), main='sd 1 distribution',xlab='sigma', col='skyblue', lwd=2)
plot(sigma[,2], main='sd 2 trace',xlab='iteration', col='blue', lwd=1)
plot(density(sigma[,2]), main='sd 2 distribution',xlab='sigma', col='skyblue', lwd=2)
plot(sigma[,3], main='sd 3 trace',xlab='iteration', col='blue', lwd=1)
plot(density(sigma[,3]), main='sd 3 distribution',xlab='sigma', col='skyblue', lwd=2)
par(mfrow=c(3,2))
plot(sigma[,4], main='sd 4 trace',xlab='iteration', col='blue', lwd=1)
plot(density(sigma[,4]), main='sd 4 distribution',xlab='sigma', col='skyblue', lwd=2)
plot(sigma[,5], main='sd 5 trace',xlab='iteration', col='blue', lwd=1)
plot(density(sigma[,5]), main='sd 5 distribution',xlab='sigma', col='skyblue', lwd=2)
plot(sigma[,6], main='sd 6 trace',xlab='iteration', col='blue', lwd=1)
plot(density(sigma[,6]), main='sd 6 distribution',xlab='sigma', col='skyblue', lwd=2)
par(mfrow=c(3,2))
plot(sigma[,7], main='sd 7 trace',xlab='iteration', col='blue', lwd=1)
plot(density(sigma[,7]), main='sd 7 distribution',xlab='sigma', col='skyblue', lwd=2)

##############################
# 4.1 The seven scientiests
# 4.2.2 Change gamma to uniform(prior over the standard deviation)
# sigma ~ Uniform(0,10)
##############################
# clears workspace:  
rm(list=ls()) 

#JARS
# install.packages("R2jags")
library(R2jags)
# Prior and Posterior Prediction
# 1. Model
modelString = "
# The Seven Scientists
model{
  # Data Come From Gaussians With Common Mean But Different Precisions
  for (i in 1:n){
    x[i] ~ dnorm(mu,lambda[i])
  }
  # Priors
  mu ~ dnorm(0,.001)
  for (i in 1:n){
    sigma[i] ~ dunif(0, 10)  
    lambda[i] <- 1/pow(sigma[i],2)
  }     
}"

# 2. Data
x <- c(-27.020,3.570,8.191,9.898,9.603,9.945,10.056)
n <- length(x)

mydata = list("x", "n")

# 3. Start Values
myinits = list(
  list(mu =0, sigma = rep(1,n))
)

# 4. Paramters to be monitored
parameters = c("mu", "sigma","lambda", "x")

# The following command calls JAGS with specific options.
# For a detailed description see the R2jags documentation.
mcmc_samples <- jags(mydata, inits=myinits, parameters,
                     model.file =textConnection(modelString), n.chains=1, n.iter=1000, 
                     n.burnin=1, n.thin=1, DIC=T)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.
x <- mcmc_samples$BUGSoutput$sims.list$x
mu <- mcmc_samples$BUGSoutput$sims.list$mu
sigma <- mcmc_samples$BUGSoutput$sims.list$sigma
lambda <- mcmc_samples$BUGSoutput$sims.list$lambda
summary(mu)
summary(sigma)
summary(lambda)

# Graph
par(mfrow=c(3,2))
plot(mu, main='mean trace',xlab='iteration', col='blue', lwd=1)
plot(density(mu), main='mean distribution',xlab='mu', col='skyblue', lwd=2)
plot(sigma[,1], main='standard deviation 1 trace',xlab='iteration', col='blue', lwd=1)
plot(sigma[,2], main='standard deviation 2 trace',xlab='iteration', col='blue', lwd=1)
plot(sigma[,3], main='standard deviation 3 trace',xlab='iteration', col='blue', lwd=1)
plot(density(sigma), main='standard deviation distribution',xlab='sigma', col='skyblue', lwd=2)

##############################
# 4.1 The seven scientiests
# 4.2.2 Change gamma to uniform(prior over the standard deviation)
# sigma ~ Uniform(0,100)
##############################
# clears workspace:  
rm(list=ls()) 

#JARS
# install.packages("R2jags")
library(R2jags)
# Prior and Posterior Prediction
# 1. Model
modelString = "
# The Seven Scientists
model{
# Data Come From Gaussians With Common Mean But Different Precisions
for (i in 1:n){
x[i] ~ dnorm(mu,lambda[i])
}
# Priors
mu ~ dnorm(0,.001)
for (i in 1:n){
sigma[i] ~ dunif(0, 100)  
lambda[i] <- 1/pow(sigma[i],2)
}     
}"

# 2. Data
x <- c(-27.020,3.570,8.191,9.898,9.603,9.945,10.056)
n <- length(x)

mydata = list("x", "n")

# 3. Start Values
myinits = list(
  list(mu =0, sigma = rep(1,n))
)

# 4. Paramters to be monitored
parameters = c("mu", "sigma","lambda", "x")

# The following command calls JAGS with specific options.
# For a detailed description see the R2jags documentation.
mcmc_samples <- jags(mydata, inits=myinits, parameters,
                     model.file =textConnection(modelString), n.chains=1, n.iter=1000, 
                     n.burnin=1, n.thin=1, DIC=T)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.
x <- mcmc_samples$BUGSoutput$sims.list$x
mu <- mcmc_samples$BUGSoutput$sims.list$mu
sigma <- mcmc_samples$BUGSoutput$sims.list$sigma
lambda <- mcmc_samples$BUGSoutput$sims.list$lambda
summary(mu)
summary(sigma)
summary(lambda)

# Graph
par(mfrow=c(3,2))
plot(mu, main='mean trace',xlab='iteration', col='blue', lwd=1)
plot(density(mu), main='mean distribution',xlab='mu', col='skyblue', lwd=2)
plot(sigma[,1], main='standard deviation 1 trace',xlab='iteration', col='blue', lwd=1)
plot(sigma[,2], main='standard deviation 2 trace',xlab='iteration', col='blue', lwd=1)
plot(sigma[,3], main='standard deviation 3 trace',xlab='iteration', col='blue', lwd=1)
plot(density(sigma), main='standard deviation distribution',xlab='sigma', col='skyblue', lwd=2)

##############################
# 4.3 Repeated mesurement of IQ
#   sigma ~ dunif(0,100)
#   mu[i] ~ dunif(0,300)
##############################
# clears workspace:  
rm(list=ls()) 

#JARS
# install.packages("R2jags")
library(R2jags)
# Prior and Posterior Prediction
# 1. Model
modelString = "
# Repeated Measures of IQ
model{
  # Data Come From Gaussians With Different Means But Common Precision
  for (i in 1:n){ #Person
    for (j in 1:m){ # IQ Test
      x[i,j] ~ dnorm(mu[i],lambda)
    }
  }
  # Priors
  sigma ~ dunif(0,100)
  lambda <- 1/pow(sigma,2)     
  for (i in 1:n){
    mu[i] ~ dunif(0,300)
  }
}"

# 2. Data
x <- matrix(c(90,95,100,105,110,115,150,155,160),nrow=3,ncol=3,byrow=T) 
n <- nrow(x) # number of people
m <- ncol(x) # number of repeated measurements

mydata = list("x", "n", "m")

# 3. Start Values
myinits = list(
  list(mu =rep(100,n), sigma = 1)
)

# 4. Paramters to be monitored
parameters = c("mu", "sigma","lambda")

# The following command calls JAGS with specific options.
# For a detailed description see the R2jags documentation.
mcmc_samples <- jags(mydata, inits=myinits, parameters,
                     model.file =textConnection(modelString), n.chains=1, n.iter=1000, 
                     n.burnin=1, n.thin=1, DIC=T)
# Now the values for the monitored parameters are in the "samples" object, 

# ready for inspection.
mu <- mcmc_samples$BUGSoutput$sims.list$mu
sigma <- mcmc_samples$BUGSoutput$sims.list$sigma
lambda <- mcmc_samples$BUGSoutput$sims.list$lambda
summary(mu)
summary(sigma)
summary(lambda)

# Graph
par(mfrow=c(3,2))
plot(mu[,1], main='mean 1 trace',xlab='iteration', col='blue', lwd=1)
plot(density(mu[,1]), main='mean 1 distribution',xlab='mu 1', col='skyblue', lwd=2)
plot(mu[,2], main='mean 2 trace',xlab='iteration', col='blue', lwd=1)
plot(density(mu[,2]), main='mean 2 distribution',xlab='mu 2', col='skyblue', lwd=2)
plot(mu[,3], main='mean 3 trace',xlab='iteration', col='blue', lwd=1)
plot(density(mu[,3]), main='mean 3 distribution',xlab='mu 3', col='skyblue', lwd=2)

par(mfrow=c(1,2))
plot(sigma, main='standard deviation trace',xlab='sigma', col='skyblue', lwd=2)
plot(density(sigma), main='standard deviation distribution',xlab='sigma', col='skyblue', lwd=2)

##############################
# 4.3 Repeated mesurement of IQ
# Exercise 4.3.2
#    sigma ~ dunif(0,100)
#    mu[i] ~ dnorm(100, 0.0044)
##############################
# clears workspace:  
rm(list=ls()) 

#JARS
# install.packages("R2jags")
library(R2jags)
# Prior and Posterior Prediction
# 1. Model
modelString = "
# Repeated Measures of IQ
model{
  # Data Come From Gaussians With Different Means But Common Precision
  for (i in 1:n){ #Person
    for (j in 1:m){ # IQ Test
      x[i,j] ~ dnorm(mu[i],lambda)
  }
  }
  # Priors
  sigma ~ dunif(0,100)
  lambda <- 1/pow(sigma,2)     
  for (i in 1:n){
    mu[i] ~ dnorm(100,0.0044)
  }
}"

# 2. Data
x <- matrix(c(90,95,100,105,110,115,150,155,160),nrow=3,ncol=3,byrow=T) 
n <- nrow(x) # number of people
m <- ncol(x) # number of repeated measurements

mydata = list("x", "n", "m")

# 3. Start Values
myinits = list(
  list(mu =rep(100,n), sigma = 1)
)

# 4. Paramters to be monitored
parameters = c("mu", "sigma","lambda")

# The following command calls JAGS with specific options.
# For a detailed description see the R2jags documentation.
mcmc_samples <- jags(mydata, inits=myinits, parameters,
                     model.file =textConnection(modelString), n.chains=1, n.iter=1000, 
                     n.burnin=1, n.thin=1, DIC=T)
# Now the values for the monitored parameters are in the "samples" object, 

# ready for inspection.
mu <- mcmc_samples$BUGSoutput$sims.list$mu
sigma <- mcmc_samples$BUGSoutput$sims.list$sigma
lambda <- mcmc_samples$BUGSoutput$sims.list$lambda
summary(mu)
summary(sigma)
summary(lambda)

# Graph
par(mfrow=c(3,2))
plot(mu[,1], main='mean 1 trace',xlab='iteration', col='blue', lwd=1)
plot(density(mu[,1]), main='mean 1 distribution',xlab='mu 1', col='skyblue', lwd=2)
plot(mu[,2], main='mean 2 trace',xlab='iteration', col='blue', lwd=1)
plot(density(mu[,2]), main='mean 2 distribution',xlab='mu 2', col='skyblue', lwd=2)
plot(mu[,3], main='mean 3 trace',xlab='iteration', col='blue', lwd=1)
plot(density(mu[,3]), main='mean 3 distribution',xlab='mu 3', col='skyblue', lwd=2)

par(mfrow=c(1,2))
plot(sigma, main='standard deviation trace',xlab='sigma', col='skyblue', lwd=2)
plot(density(sigma), main='standard deviation distribution',xlab='sigma', col='skyblue', lwd=2)

##############################
# 4.3 Repeated mesurement of IQ
# Exercise 4.3.3
#    sigma ~ dunif(0, 100)
#    mu[i] ~ dunif(0, 100)
##############################
# clears workspace:  
rm(list=ls()) 

#JARS
# install.packages("R2jags")
library(R2jags)
# Prior and Posterior Prediction
# 1. Model
modelString = "
  # Repeated Measures of IQ
  model{
  # Data Come From Gaussians With Different Means But Common Precision
  for (i in 1:n){ #Person
  for (j in 1:m){ # IQ Test
  x[i,j] ~ dnorm(mu[i],lambda)
  }
  }
  # Priors
  sigma ~ dunif(0,100)
  lambda <- 1/pow(sigma,2)     
  for (i in 1:n){
  mu[i] ~ dunif(0,100)
  }
}"

# 2. Data
x <- matrix(c(94,95,96,109,110,111,154,155,156),nrow=3,ncol=3,byrow=T) 
n <- nrow(x) # number of people
m <- ncol(x) # number of repeated measurements

mydata = list("x", "n", "m")

# 3. Start Values
myinits = list(
  list(mu =rep(100,n), sigma = 1)
)

# 4. Paramters to be monitored
parameters = c("mu", "sigma","lambda")

# The following command calls JAGS with specific options.
# For a detailed description see the R2jags documentation.
mcmc_samples <- jags(mydata, inits=myinits, parameters,
                     model.file =textConnection(modelString), n.chains=1, n.iter=1000, 
                     n.burnin=1, n.thin=1, DIC=T)
# Now the values for the monitored parameters are in the "samples" object, 

# ready for inspection.
mu <- mcmc_samples$BUGSoutput$sims.list$mu
sigma <- mcmc_samples$BUGSoutput$sims.list$sigma
lambda <- mcmc_samples$BUGSoutput$sims.list$lambda
summary(mu)
summary(sigma)
summary(lambda)

# Graph
par(mfrow=c(3,2))
plot(mu[,1], main='mean 1 trace',xlab='iteration', col='blue', lwd=1)
plot(density(mu[,1]), main='mean 1 distribution',xlab='mu 1', col='skyblue', lwd=2)
plot(mu[,2], main='mean 2 trace',xlab='iteration', col='blue', lwd=1)
plot(density(mu[,2]), main='mean 2 distribution',xlab='mu 2', col='skyblue', lwd=2)
plot(mu[,3], main='mean 3 trace',xlab='iteration', col='blue', lwd=1)
plot(density(mu[,3]), main='mean 3 distribution',xlab='mu 3', col='skyblue', lwd=2)

par(mfrow=c(1,2))
plot(sigma, main='standard deviation trace',xlab='sigma', col='skyblue', lwd=2)
plot(density(sigma), main='standard deviation distribution',xlab='sigma', col='skyblue', lwd=2)
##############################
# 4.3 Repeated mesurement of IQ
# Exercise 4.3.3
#    sigma ~ dunif(0, 100)
#    mu[i] ~ dnorm(100, 0.0044)
##############################
# clears workspace:  
rm(list=ls()) 

#JARS
# install.packages("R2jags")
library(R2jags)
# Prior and Posterior Prediction
# 1. Model
modelString = "
# Repeated Measures of IQ
model{
# Data Come From Gaussians With Different Means But Common Precision
for (i in 1:n){ #Person
for (j in 1:m){ # IQ Test
x[i,j] ~ dnorm(mu[i],lambda)
}
}
# Priors
sigma ~ dunif(0,100)
lambda <- 1/pow(sigma,2)     
for (i in 1:n){
mu[i] ~ dnorm(100, 0.0044)
}
}"

# 2. Data
x <- matrix(c(94,95,96,109,110,111,154,155,156),nrow=3,ncol=3,byrow=T) 
n <- nrow(x) # number of people
m <- ncol(x) # number of repeated measurements

mydata = list("x", "n", "m")

# 3. Start Values
myinits = list(
  list(mu =rep(100,n), sigma = 1)
)

# 4. Paramters to be monitored
parameters = c("mu", "sigma","lambda")

# The following command calls JAGS with specific options.
# For a detailed description see the R2jags documentation.
mcmc_samples <- jags(mydata, inits=myinits, parameters,
                     model.file =textConnection(modelString), n.chains=1, n.iter=1000, 
                     n.burnin=1, n.thin=1, DIC=T)
# Now the values for the monitored parameters are in the "samples" object, 

# ready for inspection.
mu <- mcmc_samples$BUGSoutput$sims.list$mu
sigma <- mcmc_samples$BUGSoutput$sims.list$sigma
lambda <- mcmc_samples$BUGSoutput$sims.list$lambda
summary(mu)
summary(sigma)
summary(lambda)

# Graph
par(mfrow=c(3,2))
plot(mu[,1], main='mean 1 trace',xlab='iteration', col='blue', lwd=1)
plot(density(mu[,1]), main='mean 1 distribution',xlab='mu 1', col='skyblue', lwd=2)
plot(mu[,2], main='mean 2 trace',xlab='iteration', col='blue', lwd=1)
plot(density(mu[,2]), main='mean 2 distribution',xlab='mu 2', col='skyblue', lwd=2)
plot(mu[,3], main='mean 3 trace',xlab='iteration', col='blue', lwd=1)
plot(density(mu[,3]), main='mean 3 distribution',xlab='mu 3', col='skyblue', lwd=2)

par(mfrow=c(1,2))
plot(sigma, main='standard deviation trace',xlab='sigma', col='skyblue', lwd=2)
plot(density(sigma), main='standard deviation distribution',xlab='sigma', col='skyblue', lwd=2)

###########################################################################
# Reference
# Variational Bayesian Multinomial Probit Regression
# http://www.r-tutor.com/gpu-computing/gaussian-process/rvbm
# http://www.bioconductor.org/packages/release/bioc/vignettes/vbmp/inst/doc/vbmp.pdf
source("http://bioconductor.org/biocLite.R")
biocLite("vbmp")
