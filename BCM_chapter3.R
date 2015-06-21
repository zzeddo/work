###########################
# Chapter 1
###########################
x=c(1:10)
plot(x, dbinom(x,10,0.9))
?seq

x=seq(from=1, to=10, by=0.01)
plot(x, dbinom(x,10,0.9), type="p" )
plot(x, dbinom(x,10,0.5) )
plot(x, dbinom(x,10,0.3) )


?plot
plot(0:10, dbinom(0:10,10,0.1), type="l")
?dbinom


###########################
# Chapter 2
###########################
#install.packages("R2WinBUGS")
# WINGBUGS
setwd("C:/Work")
library(R2WinBUGS)
bugsdir = "C:/Program Files/WinBUGS14"
k = 5
n = 10
data = list("k", "n")
myinits = list(
  list(theta = 0.1),
  list(theta = 0.9)
)
parameters = c("theta")
samples = bugs(data, inits=myinits, parameters,
               model.file="Rate_1.txt",
               n.chain = 2, n.iter = 20000, n.burnin = 1, n.thin = 1,
               DIC=T, bugs.directory=bugsdir,
               codaPkg=F, debug=F)

samples$sims.array
samples$sims.list
plot(samples)



#JARS
#install.packages("rjags")
library(rjags)
# 1. Model
modelString = "
model {
  theta ~ dbeta(1, 1)
  k ~ dbin(theta, n)
}"
# 2. Data
mydata = list(k=5, n=10)
# 3. Start Values
myinits = list(
  theta = 0.1
)
jags = jags.model(
  textConnection(modelString), 
  data = mydata, inits=myinits, n.chains = 2, 
  n.adapt= 20000)
update(jags, 20000)
mcmc_samples = coda.samples(jags, 
                            variable.names = c("theta", "k") , 
                            n.iter=20000)
plot(mcmc_samples)
##############################
library(rjags)
simpleModelString = ""
model {
  # Likelihood
  for (i in 1:N) {
    Y[i] ~ dbern(theta)
  }
  # Prior
  theta ~ dbeta(3, 3)
}
m = jags.model(
  textConnection(simpleModelString),
  data = list(
    Y = c(1, 0, 0, 1, 1),
    N = 5
  )
)
plot(coda.samples(jags.model(
  textConnection("
model {
theta ~ dbeta(1, 1)
k ~ dbin(theta, n)
}
"),
  data = list(n=2)
), c("theta", "k"), 500))
###########################
# Binomials Distribution (BUGS)
###########################
# clears workspace:  
rm(list=ls()) 

# sets working directories:
setwd("C:/Work")

library(R2WinBUGS)
bugsdir <- "C:/Program Files/WinBUGS14"

k <- 5
n <- 10

data <- list("k", "n") # to be passed on to WinBUGS

myinits <-  list(
  list(theta = 0.1), #chain 1 starting value
  list(theta = 0.9)) #chain 2 starting value

# parameters to be monitored:	
parameters <- c("theta")

# The following command calls WinBUGS with specific options.
# For a detailed description see Sturtz, Ligges, & Gelman (2005).
samples <- bugs(data, inits=myinits, parameters,
                model.file ="Rate_1.txt",
                n.chains=2, n.iter=20000, n.burnin=1, n.thin=1,
                DIC=T, bugs.directory=bugsdir,
                codaPkg=F, debug=F)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

# The commands below are useful for a quick overview:
print(samples)  # a rough summary
plot(samples)   # a visual representation
names(samples)  # summarizes the variables
samples$summary # more detailed summary
samples$sims.array[1:15,,2]# array: sample, chain, parameter 

# Collect posterior samples across all chains:
theta <- samples$sims.list$theta 

# Now let's plot a histogram for theta. 
# NB. Some the plots will not look good in RStudio.
# First, some options to make the plot look better:
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
    font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
Nbreaks <- 80
y       <- hist(theta, Nbreaks, plot=F)
plot(c(y$breaks, max(y$breaks)), c(0,y$density,0), type="S", lwd=2, lty=1,
     xlim=c(0,1), xlab="Rate", ylab="Posterior Density") 

###########################
# Binomials Distribution (JAGS)
###########################
#JARS
#install.packages("rjags")
library(rjags)
# 1. Model
modelString = "
model {
  theta ~ dbeta(1, 1)
  k ~ dbin(theta, n)
}"
# 2. Data
mydata = list(k=5, n=10)
# 3. Start Values
myinits = list(
  list(theta = 0.1),
  list(theta = 0.9)
)
jags = jags.model(
  textConnection(modelString),
  data = mydata, inits=myinits, 
  n.chains = 2, n.adapt= 10000)
update(jags, 10000)
mcmc_samples = coda.samples(jags, 
                            variable.names = c("theta", "k") , 
                            n.iter=10000)
print(mcmc_samples)
summary(mcmc_samples)

plot(mcmc_samples)
traceplot(mcmc_samples)

###########################
# Binomials Distribution (R2jags)
###########################

# clears workspace:  
rm(list=ls()) 

# sets working directories:
setwd("C:/Work")

library(R2jags)
#install.packages("R2jags")
k <- 5
n <- 10

data <- list("k", "n") # to be passed on to JAGS

myinits <-  list(
  list(theta = 0.1), #chain 1 starting value
  list(theta = 0.9)) #chain 2 starting value

# parameters to be monitored:	
parameters <- c("theta")

# The following command calls JAGS with specific options.
# For a detailed description see the R2jags documentation.
samples <- jags(data, inits=myinits, parameters,
                model.file ="Rate_1.txt", n.chains=2, n.iter=20000, 
                n.burnin=1, n.thin=1, DIC=T)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

# The commands below are useful for a quick overview:
print(samples)  # a rough summary
plot(samples)   # a visual representation
traceplot(samples) # traceplot (press <enter> repeatedly to see the chains)

#more info on what is returned:
summary(samples)
summary(samples$BUGSoutput)

samples$BUGSoutput$sims.array[1:15,,2]# array: sample, chain, parameter 

# Collect posterior samples across all chains:
theta <- samples$BUGSoutput$sims.list$theta

# Now let's plot a histogram for theta. 
# NB. Some the plots will not look good in RStudio.
# First, some options to make the plot look better:
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
    font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
Nbreaks <- 80
y       <- hist(theta, Nbreaks, plot=F)
plot(c(y$breaks, max(y$breaks)), c(0,y$density,0), type="S", lwd=2, lty=1,
     xlim=c(0,1), xlab="Rate", ylab="Posterior Density") 

# Additional option: use some plots in coda
# first use as.mcmmc to convert rjags object into mcmc.list:
samples.mcmc <- as.mcmc(samples)
# then use the plotting methods from coda:
plot(samples.mcmc)
densityplot(samples.mcmc)
###########################
# 3.1 Inferring a rate (R2jags)
###########################
# clears workspace:  
rm(list=ls()) 

#JARS
#install.packages("R2jags")
library(R2jags)

# 1. Model
modelString = "
model {
# Observed Counts
theta ~ dbeta(1, 1)
# Prior on Rates
k ~ dbin(theta, n)
}"
# 2. Data
mydata = list(k=5, n=10)

# 3. Start Values
myinits = list(
  list(theta = 0.1), #chain 1 starting value
  list(theta = 0.9) #chain 2 starting value
)

# 4. Paramters to be monitored
parameters = c("theta")

# The following command calls JAGS with specific options.
# For a detailed description see the R2jags documentation.
mcmc_samples <- jags(mydata, inits=myinits, parameters,
                     model.file =textConnection(modelString), n.chains=2, n.iter=10000, 
                     n.burnin=1, n.thin=1, DIC=T)

# Collect posterior samples:
names(mcmc_samples$BUGSoutput$sims.list)
theta <- mcmc_samples$BUGSoutput$sims.list$theta

summary(theta)
# mean of delta:
mean(theta)
# median of delta:
median(theta)
# mode of delta, estimated from the "density" smoother:
density(theta)$x[which(density(theta)$y==max(density(theta)$y))]
#95% credibility probability
quantile(theta, c(.025,.975))

# Graph
par(mfrow=c(1,2))
plot(theta, main='Posterior trace',xlab='iteration', col='blue', lwd=1)
plot(density(theta), main='Posterior distribution',xlab='theta', col='skyblue', lwd=2)

# Now let's plot a histogram for theta 
# First, some options to make the plot look better:
par(mfrow=c(1,1))
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
    font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
Nbreaks <- 80
y       <- hist(theta, Nbreaks, plot=F)
plot(c(y$breaks, max(y$breaks)), c(0,y$density,0), type="S", lwd=2, lty=1,
     xlim=c(0,1), ylim=c(0,3), xlab="Difference in Rates", ylab="Posterior Density") 

###########################
# 3.2 Difference between two rates (JAGS)
###########################
# clears workspace:  
rm(list=ls()) 

#JARS
library(rjags)
# 1. Model
modelString = "
model {
# Observed Counts
theta1 ~ dbeta(1, 1)
theta2 ~ dbeta(1, 1)
# Prior on Rates
k1 ~ dbin(theta1, n1)
k2 ~ dbin(theta2, n2)
# Difference Between Rates 
delta <- theta1-theta2
}"
# 2. Data
mydata = list(k1=5, n1=10, k2=7, n2=10)

# 3. Start Values
myinits = list(
  list(theta1 = 0.1, theta2 = 0.9)
)

# 4. Paramters to be monitored
parameter = c("delta", "theta1", "theta2")

jags = jags.model(
  textConnection(modelString),
  data = mydata, inits = myinits, 
  n.chains = 1, n.adapt= 10000)
update(jags, 10000)
mcmc_samples = coda.samples(jags, 
                            variable.names = parameter , 
                            n.iter=10000)

summary(mcmc_samples)$stat
summary(mcmc_samples)$quant
summary(mcmc_samples)
plot(mcmc_samples, font.main=4)

# Collect posterior samples across all chains:

delta <- mcmc_samples[,1]
summary(delta)
plot(delta)
theta1 <- mcmc_samples[,2]
summary(theta1)
plot(theta1)
theta2 <- mcmc_samples[,3]
summary(theta2)
plot(theta2)

###########################
# 3.2 Difference between two rates (R2jags)
###########################
# clears workspace:  
rm(list=ls()) 

#JARS
#install.packages("R2jags")
library(R2jags)

# 1. Model
modelString = "
model {
# Observed Counts
theta1 ~ dbeta(1, 1)
theta2 ~ dbeta(1, 1)
# Prior on Rates
k1 ~ dbin(theta1, n1)
k2 ~ dbin(theta2, n2)
# Difference Between Rates 
delta <- theta1-theta2
}"
# 2. Data
mydata = list(k1=5, n1=10, k2=7, n2=10)

# 3. Start Values
myinits = list(
  list(theta1 = 0.1, theta2 = 0.9)
)

# 4. Paramters to be monitored
parameter = c("delta", "theta1", "theta2")

# The following command calls JAGS with specific options.
# For a detailed description see the R2jags documentation.
mcmc_samples <- jags(mydata, inits=myinits, parameter,
                     model.file =textConnection(modelString), n.chains=1, n.iter=10000, 
                     n.burnin=1, n.thin=1, DIC=T)

#95% credibility probability
mcmc_samples$summary
#mymat <- mcmc_samples$BUGSoutput$sims.array[,,"delta"]

# Collect posterior samples:
names(mcmc_samples$BUGSoutput$sims.list)
delta <- mcmc_samples$BUGSoutput$sims.list$delta
theta1 <- mcmc_samples$BUGSoutput$sims.list$theta1
theta2 <- mcmc_samples$BUGSoutput$sims.list$theta2

cbind(summary(theta1), summary(theta2), summary(delta)
# mean of delta:
mean(delta)
# median of delta:
median(delta)
# mode of delta, estimated from the "density" smoother:
density(delta)$x[which(density(delta)$y==max(density(delta)$y))]
# 95% credible interval for delta:
quantile(delta, c(.025,.975))

# Graph
par(mfrow=c(3,2))
plot(delta, main='Posterior trace',xlab='iteration', col='blue', lwd=1)
plot(density(delta), main='Posterior distribution',xlab='delta', col='skyblue', lwd=2)
plot(theta1, main='Posterior trace',xlab='iteration', col='blue', lwd=1)
plot(density(theta1), main='Posterior distribution',xlab='theta1', col='skyblue', lwd=2)
plot(theta2, main='Posterior trace',xlab='iteration', col='blue', lwd=1)
plot(density(theta2), main='Posterior distribution',xlab='theta2', col='skyblue', lwd=2)

#Graph Mapping
par(mfrow=c(1,1))
plot(density(delta), main='Posterior distribution', xlab ='Distributoin', xlim=c(-1,1),ylim=c(0,5), lty=1, col='1')
par(new=TRUE)
plot(density(theta1),xlim=c(-1,1),ylim=c(0,5), main="", xlab ="", ylab="",axes=FALSE,lty=2, col='2')
par(new=TRUE)
plot(density(theta2),xlim=c(-1,1),ylim=c(0,5), main="", xlab ="", ylab="",axes=FALSE, lty=3, col='3')
legend(0.5,5,legend=c("delta","theta1", "theta2"), lty=1:3, col=1:3)

# Now let's plot a histogram for delta. 
# First, some options to make the plot look better:
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
    font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
Nbreaks <- 80
y       <- hist(delta, Nbreaks, plot=F)
plot(c(y$breaks, max(y$breaks)), c(0,y$density,0), type="S", lwd=2, lty=1,
     xlim=c(-1,1), ylim=c(0,4), xlab="Difference in Rates", ylab="Posterior Density") 

###########################
# 3.3 Inferring a common rate (R2jags)
###########################
# clears workspace:  
rm(list=ls()) 

#JARS
library(R2jags)

# 1. Model
modelString = "
model {
# Observed Counts
k1 ~ dbin(theta, n1)
k2 ~ dbin(theta, n2)
# Prior on Single Rate Theta
theta ~ dbeta(1, 1)
}"
# 2. Data
mydata = list(k1=5, n1=10, k2=7, n2=10)

# 3. Start Values
myinits = list(
  list(theta1 =0.5)
)

# 4. Paramters to be monitored
parameters = c("theta")

# The following command calls JAGS with specific options.
# For a detailed description see the R2jags documentation.
mcmc_samples <- jags(mydata, inits=myinits, parameters,
                     model.file =textConnection(modelString), n.chains=1, n.iter=1000, 
                     n.burnin=1, n.thin=1, DIC=T)
theta <- mcmc_samples$BUGSoutput$sims.list$theta

#95% credibility probability

summary(theta)
# mean of delta:
mean(theta)
# median of delta:
median(theta)
# mode of delta, estimated from the "density" smoother:
density(theta)$x[which(density(theta)$y==max(density(theta)$y))]
#95% credibility probability
quantile(theta, c(.025,.975))

# Graph
par(mfrow=c(1,2))
plot(theta, main='Posterior trace',xlab='iteration', col='blue', lwd=1)
plot(density(theta), main='Posterior distribution',xlab='theta', col='skyblue', lwd=2)

###########################
# 3.3 Inferring a common rate (R2jags)
# Exercise 3.3.1 : k1=14, n1=20, k2=16, n2=20
###########################
# clears workspace:  
rm(list=ls()) 

#JARS
library(R2jags)

# 1. Model
modelString = "
model {
# Observed Counts
k1 ~ dbin(theta, n1)
k2 ~ dbin(theta, n2)
# Prior on Single Rate Theta
theta ~ dbeta(1, 1)
}"
# 2. Data
mydata = list(k1=14, n1=20, k2=16, n2=20)

# 3. Start Values
myinits = list(
  list(theta1 =0.5)
)

# 4. Paramters to be monitored
parameters = c("theta")

# The following command calls JAGS with specific options.
# For a detailed description see the R2jags documentation.
mcmc_samples <- jags(mydata, inits=myinits, parameters,
                     model.file =textConnection(modelString), n.chains=1, n.iter=1000, 
                     n.burnin=1, n.thin=1, DIC=T)
theta <- mcmc_samples$BUGSoutput$sims.list$theta

#95% credibility probability

summary(theta)
# mean of delta:
mean(theta)
# median of delta:
median(theta)
# mode of delta, estimated from the "density" smoother:
density(theta)$x[which(density(theta)$y==max(density(theta)$y))]
#95% credibility probability
quantile(theta, c(.025,.975))

# Graph
par(mfrow=c(1,2))
plot(theta, main='Posterior trace',xlab='iteration', col='blue', lwd=1)
plot(density(theta), main='Posterior distribution',xlab='theta', col='skyblue', lwd=2)

###########################
# 3.3 Inferring a common rate (R2jags)
# Exercise 3.3.2 : k1=0, n1=10, k2=10, n2=10
###########################
# clears workspace:  
rm(list=ls()) 

#JARS
library(R2jags)

# 1. Model
modelString = "
model {
# Observed Counts
k1 ~ dbin(theta, n1)
k2 ~ dbin(theta, n2)
# Prior on Single Rate Theta
theta ~ dbeta(1, 1)
}"
# 2. Data
mydata = list(k1=0, n1=10, k2=10, n2=10)

# 3. Start Values
myinits = list(
  list(theta1 =0.5)
)

# 4. Paramters to be monitored
parameters = c("theta")

# The following command calls JAGS with specific options.
# For a detailed description see the R2jags documentation.
mcmc_samples <- jags(mydata, inits=myinits, parameters,
                     model.file =textConnection(modelString), n.chains=1, n.iter=1000, 
                     n.burnin=1, n.thin=1, DIC=T)
theta <- mcmc_samples$BUGSoutput$sims.list$theta

#95% credibility probability

summary(theta)
# mean of delta:
mean(theta)
# median of delta:
median(theta)
# mode of delta, estimated from the "density" smoother:
density(theta)$x[which(density(theta)$y==max(density(theta)$y))]
#95% credibility probability
quantile(theta, c(.025,.975))

# Graph
par(mfrow=c(1,2))
plot(theta, main='Posterior trace',xlab='iteration', col='blue', lwd=1)
plot(density(theta), main='Posterior distribution',xlab='theta', col='skyblue', lwd=2)

###########################
# 3.3 Inferring a common rate (R2jags)
# Exercise 3.3.3 : Compare the CASE1, CASE2
# CASE 1 : k1=7, n1=10, k2=3, n2=10
# CASE 2 : k1=5, n1=10, k2=5, n2=10
###########################
# clears workspace:  
rm(list=ls()) 

#JARS
library(R2jags)

# 1. Model
modelString = "
model {
# Observed Counts
k1 ~ dbin(theta, n1)
k2 ~ dbin(theta, n2)
# Prior on Single Rate Theta
theta ~ dbeta(1, 1)
}"
# 2. Data
mydata1 = list(k1=3, n1=10, k2=7, n2=10)
mydata2 = list(k1=5, n1=10, k2=5, n2=10)

# 3. Start Values
myinits = list(
  list(theta =0.5)
)

# 4. Paramters to be monitored
parameters = c("theta")

# The following command calls JAGS with specific options.
# For a detailed description see the R2jags documentation.
mcmc_samples1 <- jags(mydata1, inits=myinits, parameters,
                      model.file =textConnection(modelString), n.chains=1, n.iter=1000, 
                      n.burnin=1, n.thin=1, DIC=T)

mcmc_samples2 <- jags(mydata2, inits=myinits, parameters,
                      model.file =textConnection(modelString), n.chains=1, n.iter=1000, 
                      n.burnin=1, n.thin=1, DIC=T)
theta1 <- mcmc_samples1$BUGSoutput$sims.list$theta
theta2 <- mcmc_samples2$BUGSoutput$sims.list$theta

# Summary
cbind(summary(theta1), summary(theta2))
#95% credibility probability
cbind(quantile(theta1, c(.025,.975)), quantile(theta2, c(.025,.975)))

# Graph
par(mfrow=c(2,2))
plot(theta1, main='Posterior trace',xlab='iteration', col='blue', lwd=1)
plot(density(theta1), main='Posterior distribution',xlab='theta1', col='skyblue', lwd=2)
plot(theta2, main='Posterior trace',xlab='iteration', col='blue', lwd=1)
plot(density(theta1), main='Posterior distribution',xlab='theta2', col='skyblue', lwd=2)

#Graph Mapping
par(mfrow=c(1,1))
plot(density(theta1), main='Posterior distribution', xlab ='Distributoin', xlim=c(0,1),ylim=c(0,4), lty=1, col='1')
par(new=TRUE)
plot(density(theta2),xlim=c(0,1),ylim=c(0,4), main="", xlab ="", ylab="",axes=FALSE,lty=2, col='2')
legend(0,3,legend=c("theta1", "theta2"), lty=1:2, col=1:2)

###########################
# 3.4 Prior and posterior prediction (R2jags)
###########################
# clears workspace:  
rm(list=ls()) 

#JARS
library(R2jags)

# Prior and Posterior Prediction
# 1. Model
modelString = "
model{
  # Observed Data
  k ~ dbin(theta,n)
  # Prior on Rate Theta
  theta ~ dbeta(1,1)
  # Posterior Predictive
  postpredk ~ dbin(theta,n)
  # Prior Predictive
  thetaprior ~ dbeta(1,1)
  priorpredk ~ dbin(thetaprior,n)
}"
# 2. Data
k=1
n=15
mydata = list("k", "n")

# 3. Start Values
myinits = list(
  list(theta =0.5, thetaprior = 0.5)
)

# 4. Paramters to be monitored
parameters = c("theta", "thetaprior", "k", "postpredk", "priorpredk")

# The following command calls JAGS with specific options.
# For a detailed description see the R2jags documentation.
mcmc_samples <- jags(mydata, inits=myinits, parameters,
                     model.file =textConnection(modelString), n.chains=1, n.iter=5000, 
                     n.burnin=1, n.thin=1, DIC=T)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.
theta <- mcmc_samples$BUGSoutput$sims.list$theta
thetaprior <- mcmc_samples$BUGSoutput$sims.list$thetaprior
k <- mcmc_samples$BUGSoutput$sims.list$postpredk
postpredk <- mcmc_samples$BUGSoutput$sims.list$postpredk
priorpredk <- mcmc_samples$BUGSoutput$sims.list$priorpredk

layout(matrix(c(1,2),2,1))
layout.show(2)
#Prior and posterior of theta
plot(density(theta, from=0, to=1), zero.line=F, axes=F, main="", xlab="", ylab="", xlim=c(0,1), ylim=c(0,6))
axis(1, at=c(0,0.2,0.4,0.6,0.8,1), lab=c("0","0.2","0.4","0.6","0.8","1"),cex.axis=0.8)
mtext("Rate", side=1, line=2.25, cex=1.2)
axis(2, at=c(0,2,4,6),cex.axis=0.8)
mtext("Density", side=2, line=2.25, cex=1.2)
lines(density(thetaprior, from=0, to=1), lty=3, col="gray")
#legend(0.6,5.75, legend=c("Prior", "Posterior"), lty=c(3,1), col=c ("grey", "black"))

#Prior and posterior predictive
mybreaks <- seq(from=-.5,to=n+1,by=1)
my.at    <- seq(from=0,to=n,by=1)
hist(postpredk,breaks=mybreaks,freq=F, right=F, ylab="", xlab="", ylim=c(0,0.3),main="", axes=F )
axis(1, at=my.at,lab=my.at,cex.axis=0.8)
mtext("Success Count", side=1, line=2.25, cex=1.2)
axis(2,at=c(0,0.1,0.2,0.3),lab=c("0","0.1","0.2","0.3"),cex.axis=0.8)
mtext("Mass", side=2, line=2.25, cex=1.2)
hist(priorpredk, breaks=mybreaks,freq=F,right=F,add=T, lty=3,border="grey")
#legend(8,0.3, legend=c("Prior", "Posterior"), lty=c(3,1),col=c("grey", "black"))

summary(theta)
summary(thetaprior)

# Graph
par(mfrow=c(3,2))
plot(theta, main='Posterior trace',xlab='iteration', col='blue', lwd=1)
plot(density(theta), main='Posterior distribution',xlab='theta', col='skyblue', lwd=2)
plot(thetaprior, main='Prior Prediction trace',xlab='iteration', col='blue', lwd=1)
plot(density(thetaprior), main='Prior Prediction distribution',xlab='thetaprior', col='skyblue', lwd=2)
plot(theta, main='Posterior Prediction trace',xlab='iteration', col='blue', lwd=1)
plot(density(theta), main='Posterior Prediction distribution',xlab='theta', col='skyblue', lwd=2)

par(mfrow=c(3,2))
plot(k, main='Posterior trace',xlab='iteration', col='blue', lwd=1)
plot(density(k), main='Posterior distribution',xlab='k', col='skyblue', lwd=2)
plot(priorpredk, main='Prior Prediction trace',xlab='iteration', col='blue', lwd=1)
plot(density(priorpredk), main='Prior Prediction distribution',xlab='priorpredk', col='skyblue', lwd=2)
plot(postpredk, main='Posterior Prediction trace',xlab='iteration', col='blue', lwd=1)
plot(density(postpredk), main='Posterior Prediction distribution',xlab='postpredk', col='skyblue', lwd=2)

###########################
# 3.4 Prior and posterior prediction (R2jags)
# Exercise 3.4.1 k=24, n=121
###########################
# clears workspace:  
rm(list=ls()) 

#JARS
library(R2jags)

# Prior and Posterior Prediction
# 1. Model
modelString = "
model{
# Observed Data
k ~ dbin(theta,n)
# Prior on Rate Theta
theta ~ dbeta(1,1)
# Posterior Predictive
postpredk ~ dbin(theta,n)
# Prior Predictive
thetaprior ~ dbeta(1,1)
priorpredk ~ dbin(thetaprior,n)
}"
# 2. Data
k=24
n=121
mydata = list("k", "n")

# 3. Start Values
myinits = list(
  list(theta =0.5, thetaprior = 0.5)
)

# 4. Paramters to be monitored
parameters = c("theta", "thetaprior", "k", "postpredk", "priorpredk")

# The following command calls JAGS with specific options.
# For a detailed description see the R2jags documentation.
mcmc_samples <- jags(mydata, inits=myinits, parameters,
                     model.file =textConnection(modelString), n.chains=1, n.iter=5000, 
                     n.burnin=1, n.thin=1, DIC=T)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.
theta <- mcmc_samples$BUGSoutput$sims.list$theta
thetaprior <- mcmc_samples$BUGSoutput$sims.list$thetaprior
k <- mcmc_samples$BUGSoutput$sims.list$postpredk
postpredk <- mcmc_samples$BUGSoutput$sims.list$postpredk
priorpredk <- mcmc_samples$BUGSoutput$sims.list$priorpredk

layout(matrix(c(1,2),2,1))
layout.show(2)
#Prior and posterior of theta
plot(density(theta, from=0, to=1), zero.line=F, axes=F, main="", xlab="", ylab="", xlim=c(0,1), ylim=c(0,15))
axis(1, at=c(0,0.2,0.4,0.6,0.8,1), lab=c("0","0.2","0.4","0.6","0.8","1"),cex.axis=0.8)
mtext("Rate", side=1, line=2.25, cex=1.2)
axis(2, at=c(0,2,4,6,8,10,12,14,16),cex.axis=0.8)
mtext("Density", side=2, line=2.25, cex=1.2)
lines(density(thetaprior, from=0, to=1), lty=3, col="gray")
#legend(0.6,5.75, legend=c("Prior", "Posterior"), lty=c(3,1), col=c ("grey", "black"))

#Prior and posterior predictive
mybreaks <- seq(from=-.5,to=n+1,by=1)
my.at    <- seq(from=0,to=n,by=1)
hist(postpredk,breaks=mybreaks,freq=F, right=F, ylab="", xlab="", ylim=c(0,0.1),main="", axes=F )
axis(1, at=my.at,lab=my.at,cex.axis=0.8)
mtext("Success Count", side=1, line=2.25, cex=1.2)
axis(2,at=c(0,0.1),lab=c("0","0.1"),cex.axis=0.8)
mtext("Mass", side=2, line=2.25, cex=1.2)
hist(priorpredk, breaks=mybreaks,freq=F,right=F,add=T, lty=3,border="grey")
#legend(8,0.3, legend=c("Prior", "Posterior"), lty=c(3,1),col=c("grey", "black"))

summary(theta)
summary(thetaprior)

# Graph
par(mfrow=c(3,2))
plot(theta, main='Posterior trace',xlab='iteration', col='blue', lwd=1)
plot(density(theta), main='Posterior distribution',xlab='theta', col='skyblue', lwd=2)
plot(thetaprior, main='Prior Prediction trace',xlab='iteration', col='blue', lwd=1)
plot(density(thetaprior), main='Prior Prediction distribution',xlab='thetaprior', col='skyblue', lwd=2)
plot(theta, main='Posterior Prediction trace',xlab='iteration', col='blue', lwd=1)
plot(density(theta), main='Posterior Prediction distribution',xlab='theta', col='skyblue', lwd=2)

par(mfrow=c(3,2))
plot(k, main='Posterior trace',xlab='iteration', col='blue', lwd=1)
plot(density(k), main='Posterior distribution',xlab='k', col='skyblue', lwd=2)
plot(priorpredk, main='Prior Prediction trace',xlab='iteration', col='blue', lwd=1)
plot(density(priorpredk), main='Prior Prediction distribution',xlab='priorpredk', col='skyblue', lwd=2)
plot(postpredk, main='Posterior Prediction trace',xlab='iteration', col='blue', lwd=1)
plot(density(postpredk), main='Posterior Prediction distribution',xlab='postpredk', col='skyblue', lwd=2)

###########################
# 3.4 Prior and posterior prediction (R2jags)
# Exercise 3.4.2 theta ~ dbeta(3, 1)
###########################
# clears workspace:  
rm(list=ls()) 

#JARS
library(R2jags)

# Prior and Posterior Prediction
# 1. Model
modelString = "
model{
# Observed Data
k ~ dbin(theta,n)
# Prior on Rate Theta
theta ~ dbeta(3,1)
# Posterior Predictive
postpredk ~ dbin(theta,n)
# Prior Predictive
thetaprior ~ dbeta(1,1)
priorpredk ~ dbin(thetaprior,n)
}"
# 2. Data
k=5
n=10
mydata = list("k", "n")

# 3. Start Values
myinits = list(
  list(theta =0.5, thetaprior = 0.5)
)

# 4. Paramters to be monitored
parameters = c("theta", "thetaprior", "k", "postpredk", "priorpredk")

# The following command calls JAGS with specific options.
# For a detailed description see the R2jags documentation.
mcmc_samples <- jags(mydata, inits=myinits, parameters,
                     model.file =textConnection(modelString), n.chains=1, n.iter=5000, 
                     n.burnin=1, n.thin=1, DIC=T)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.
theta <- mcmc_samples$BUGSoutput$sims.list$theta
thetaprior <- mcmc_samples$BUGSoutput$sims.list$thetaprior
k <- mcmc_samples$BUGSoutput$sims.list$postpredk
postpredk <- mcmc_samples$BUGSoutput$sims.list$postpredk
priorpredk <- mcmc_samples$BUGSoutput$sims.list$priorpredk

layout(matrix(c(1,2),2,1))
layout.show(2)
#Prior and posterior of theta
plot(density(theta, from=0, to=1), zero.line=F, axes=F, main="", xlab="", ylab="", xlim=c(0,1), ylim=c(0,6))
axis(1, at=c(0,0.2,0.4,0.6,0.8,1), lab=c("0","0.2","0.4","0.6","0.8","1"),cex.axis=0.8)
mtext("Rate", side=1, line=2.25, cex=1.2)
axis(2, at=c(0,2,4,6),cex.axis=0.8)
mtext("Density", side=2, line=2.25, cex=1.2)
lines(density(thetaprior, from=0, to=1), lty=3, col="gray")
#legend(0.6,5.75, legend=c("Prior", "Posterior"), lty=c(3,1), col=c ("grey", "black"))

#Prior and posterior predictive
mybreaks <- seq(from=-.5,to=n+1,by=1)
my.at    <- seq(from=0,to=n,by=1)
hist(postpredk,breaks=mybreaks,freq=F, right=F, ylab="", xlab="", ylim=c(0,0.3),main="", axes=F )
axis(1, at=my.at,lab=my.at,cex.axis=0.8)
mtext("Success Count", side=1, line=2.25, cex=1.2)
axis(2,at=c(0,0.1,0.2,0.3),lab=c("0","0.1","0.2","0.3"),cex.axis=0.8)
mtext("Mass", side=2, line=2.25, cex=1.2)
hist(priorpredk, breaks=mybreaks,freq=F,right=F,add=T, lty=3,border="grey")
#legend(8,0.3, legend=c("Prior", "Posterior"), lty=c(3,1),col=c("grey", "black"))

summary(theta)
summary(thetaprior)

# Graph
par(mfrow=c(3,2))
plot(theta, main='Posterior trace',xlab='iteration', col='blue', lwd=1)
plot(density(theta), main='Posterior distribution',xlab='theta', col='skyblue', lwd=2)
plot(thetaprior, main='Prior Prediction trace',xlab='iteration', col='blue', lwd=1)
plot(density(thetaprior), main='Prior Prediction distribution',xlab='thetaprior', col='skyblue', lwd=2)
plot(theta, main='Posterior Prediction trace',xlab='iteration', col='blue', lwd=1)
plot(density(theta), main='Posterior Prediction distribution',xlab='theta', col='skyblue', lwd=2)

par(mfrow=c(3,2))
plot(k, main='Posterior trace',xlab='iteration', col='blue', lwd=1)
plot(density(k), main='Posterior distribution',xlab='k', col='skyblue', lwd=2)
plot(priorpredk, main='Prior Prediction trace',xlab='iteration', col='blue', lwd=1)
plot(density(priorpredk), main='Prior Prediction distribution',xlab='priorpredk', col='skyblue', lwd=2)
plot(postpredk, main='Posterior Prediction trace',xlab='iteration', col='blue', lwd=1)
plot(density(postpredk), main='Posterior Prediction distribution',xlab='postpredk', col='skyblue', lwd=2)

###########################
# 3.4 Prior and posterior prediction (R2jags)
# Exercise 3.4.3 n' <> n observed data is not same. 
###########################
# clears workspace:  
rm(list=ls()) 

#JARS
library(R2jags)

# Prior and Posterior Prediction
# 1. Model
modelString = "
model{
# Observed Data
k ~ dbin(theta,n)
# Prior on Rate Theta
theta ~ dbeta(1,1)
# Posterior Predictive
postpredk ~ dbin(theta,npred)
# Prior Predictive
thetaprior ~ dbeta(1,1)
priorpredk ~ dbin(thetaprior,npred)
}"
# 2. Data
k=5
n=10
npred=15
mydata = list("k", "n", "npred")

# 3. Start Values
myinits = list(
  list(theta =0.5, thetaprior = 0.5)
)

# 4. Paramters to be monitored
parameters = c("theta", "thetaprior", "k", "postpredk", "priorpredk")

# The following command calls JAGS with specific options.
# For a detailed description see the R2jags documentation.
mcmc_samples <- jags(mydata, inits=myinits, parameters,
                     model.file =textConnection(modelString), n.chains=1, n.iter=5000, 
                     n.burnin=1, n.thin=1, DIC=T)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.
theta <- mcmc_samples$BUGSoutput$sims.list$theta
thetaprior <- mcmc_samples$BUGSoutput$sims.list$thetaprior
k <- mcmc_samples$BUGSoutput$sims.list$postpredk
postpredk <- mcmc_samples$BUGSoutput$sims.list$postpredk
priorpredk <- mcmc_samples$BUGSoutput$sims.list$priorpredk

layout(matrix(c(1,2),2,1))
layout.show(2)
#Prior and posterior of theta
plot(density(theta, from=0, to=1), zero.line=F, axes=F, main="", xlab="", ylab="", xlim=c(0,1), ylim=c(0,6))
axis(1, at=c(0,0.2,0.4,0.6,0.8,1), lab=c("0","0.2","0.4","0.6","0.8","1"),cex.axis=0.8)
mtext("Rate", side=1, line=2.25, cex=1.2)
axis(2, at=c(0,2,4,6),cex.axis=0.8)
mtext("Density", side=2, line=2.25, cex=1.2)
lines(density(thetaprior, from=0, to=1), lty=3, col="gray")
#legend(0.6,5.75, legend=c("Prior", "Posterior"), lty=c(3,1), col=c ("grey", "black"))

#Prior and posterior predictive
mybreaks <- seq(from=-.5,to=npred+1,by=1)
my.at    <- seq(from=0,to=npred,by=1)
hist(postpredk,breaks=mybreaks,freq=F, right=F, ylab="", xlab="", ylim=c(0,0.3),main="", axes=F )
axis(1, at=my.at,lab=my.at,cex.axis=0.8)
mtext("Success Count", side=1, line=2.25, cex=1.2)
axis(2,at=c(0,0.1,0.2,0.3),lab=c("0","0.1","0.2","0.3"),cex.axis=0.8)
mtext("Mass", side=2, line=2.25, cex=1.2)
hist(priorpredk, breaks=mybreaks,freq=F,right=F,add=T, lty=3,border="grey")
#legend(8,0.3, legend=c("Prior", "Posterior"), lty=c(3,1),col=c("grey", "black"))

summary(theta)
summary(thetaprior)

# Graph
par(mfrow=c(3,2))
plot(theta, main='Posterior trace',xlab='iteration', col='blue', lwd=1)
plot(density(theta), main='Posterior distribution',xlab='theta', col='skyblue', lwd=2)
plot(thetaprior, main='Prior Prediction trace',xlab='iteration', col='blue', lwd=1)
plot(density(thetaprior), main='Prior Prediction distribution',xlab='thetaprior', col='skyblue', lwd=2)
plot(theta, main='Posterior Prediction trace',xlab='iteration', col='blue', lwd=1)
plot(density(theta), main='Posterior Prediction distribution',xlab='theta', col='skyblue', lwd=2)

par(mfrow=c(3,2))
plot(k, main='Posterior trace',xlab='iteration', col='blue', lwd=1)
plot(density(k), main='Posterior distribution',xlab='k', col='skyblue', lwd=2)
plot(priorpredk, main='Prior Prediction trace',xlab='iteration', col='blue', lwd=1)
plot(density(priorpredk), main='Prior Prediction distribution',xlab='priorpredk', col='skyblue', lwd=2)
plot(postpredk, main='Posterior Prediction trace',xlab='iteration', col='blue', lwd=1)
plot(density(postpredk), main='Posterior Prediction distribution',xlab='postpredk', col='skyblue', lwd=2)

###########################
# 3.4 Prior and posterior prediction (R2jags)
# Exercise 3.4.4 bullied elderly. 
###########################
par(mfrow=c(1,1))
theta = dbinom(1:121, 121, 24/121)
#95% credibility probability
Low=0
i=1
Temp=0
while(Temp <= 0.025){
  i=i+1
  Low=i
  Temp=sum(dbinom(1:i, 121, 24/121))
}
High=0
i=1
Temp=0
while(Temp <= 0.975){
  i=i+1
  High=i
  Temp=sum(dbinom(1:i, 121, 24/121))
}
Low
High

###########################
# 3.5 Posterior prediction (R2jags) 
###########################
# clears workspace:  
rm(list=ls()) 

#JARS
library(R2jags)

# Prior and Posterior Prediction
# 1. Model
modelString = "
model{  
   # Observed Counts
   k1 ~ dbin(theta,n1)
   k2 ~ dbin(theta,n2)
   # Prior on Single Rate Theta
   theta ~ dbeta(1,1)
   # Posterior Predictive
   postpredk1 ~ dbin(theta,n1)
   postpredk2 ~ dbin(theta,n2)   
}"
# 2. Data
k1=0
k2=10
n1=10
n2=10
mydata = list("k1", "n1", "k2", "n2")

# 3. Start Values
myinits = list(
  list(theta =0.5)
)

# 4. Paramters to be monitored
parameters = c("theta", "k1", "k2", "postpredk1", "postpredk2")

# The following command calls JAGS with specific options.
# For a detailed description see the R2jags documentation.
mcmc_samples <- jags(mydata, inits=myinits, parameters,
                     model.file =textConnection(modelString), n.chains=1, n.iter=5000, 
                     n.burnin=1, n.thin=1, DIC=T)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.
theta <- mcmc_samples$BUGSoutput$sims.list$theta
k1 <- mcmc_samples$BUGSoutput$sims.list$k1
k2 <- mcmc_samples$BUGSoutput$sims.list$k2
postpredk1 <- mcmc_samples$BUGSoutput$sims.list$postpredk1
postpredk2 <- mcmc_samples$BUGSoutput$sims.list$postpredk2

# Two-panel plot. 
layout(matrix(c(1,2),1,2))
layout.show(2)
# First, a histogram for theta.
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
    font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
Nbreaks <- 80
y       <- hist(theta, Nbreaks, plot=F)
plot(c(y$breaks, max(y$breaks)), c(0,y$density,0), type="S", lwd=2, lty=1,
     xlim=c(0,1), ylim=c(0,4), xlab="Theta", ylab="Density") 
# let's plot a density estimate over this:
lines(density(theta), col="red", lwd=2)

# Second plot, the data space (predictives)
plot(k1,k2,type="p", pch=4, cex=2, lwd=2, xlab="Success Count 1", ylab="Success Count 2",
     xlim=c(-1, n1+1), ylim=c(-1,n2+1))        
nsamples <- length(theta)
sc <- 10
for (i in 0:n1)
{
  for (j in 0:n2) 
  {
    match.preds <- sum(postpredk1==i & postpredk2==j)/nsamples
    if (match.preds > 0)
    {
      points(i,j, pch=0, cex=sc*sqrt(match.preds)) 
    }
  }
}

summary(theta)

# Graph
par(mfrow=c(1,2))
plot(theta, main='Posterior trace',xlab='iteration', col='blue', lwd=1)
plot(density(theta), main='Posterior distribution',xlab='theta', col='skyblue', lwd=2)

par(mfrow=c(2,2))
plot(k1, main='Posterior trace',xlab='iteration', col='blue', lwd=1)
plot(density(k1), main='Posterior distribution',xlab='k1', col='skyblue', lwd=2)
plot(k2, main='Posterior trace',xlab='iteration', col='blue', lwd=1)
plot(density(k2), main='Posterior distribution',xlab='k2', col='skyblue', lwd=2)

par(mfrow=c(2,2))
plot(postpredk1, main='Posterior Prediction trace',xlab='iteration', col='blue', lwd=1)
plot(density(postpredk1), main='Posterior Prediction distribution',xlab='postpredk1', col='skyblue', lwd=2)
plot(postpredk2, main='Posterior Prediction trace',xlab='iteration', col='blue', lwd=1)
plot(density(postpredk2), main='Posterior Prediction distribution',xlab='postpredk2', col='skyblue', lwd=2)
###########################
# 3.6 Joint Distribution (R2jags) 
###########################
# clears workspace:  
rm(list=ls()) 

#JARS
library(R2jags)

# Inferring Return Rate and Number of Surveys from Observed Returns
# 1. Model
modelString = "
model{
  # Observed Returns
  for (i in 1:m){
     k[i] ~ dbin(theta,n)
  }   
  # Priors on Rate Theta and Number n
  theta ~ dbeta(1,1)
  n ~ dcat(p[])
  for (i in 1:nmax){
    p[i] <- 1/nmax
  }
}"
# 2. Data
nmax <- 500
k    <- c(16,18,22,25,27)
m    <- length(k)
data <- list("nmax", "k", "m") # to be passed on to JAGS

# 3. Start Values
myinits = list(
  list(theta =0.5, n = nmax)
)

# 4. Paramters to be monitored
parameters = c("theta", "n")

# The following command calls JAGS with specific options.
# For a detailed description see the R2jags documentation.
samples <- jags(data, inits=myinits, parameters,
                model.file =textConnection(modelString), n.chains=1, n.iter=5000, 
                n.burnin=1, n.thin=1, DIC=T)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

theta      <- samples$BUGSoutput$sims.list$theta
n          <- samples$BUGSoutput$sims.list$n 

par(mfrow=c(1,1))
plot(density(theta))
?dcat
## First calculate MLE:
cc <- -Inf
ind <- 0

for (i in 1:length(n))
{
  logL <- 0
  for(j in 1:m)
  {		
    logL <- logL+lgamma(n[i]+1)-lgamma(k[j]+1)-lgamma(n[i]-k[j]+1)
    logL <- logL+k[j]*log(theta[i])+(n[i]-k[j])*log(1-theta[i])
  }
  if (logL>cc) 
  {
    ind <- i
    cc <- logL
  }
}
# end MLE

######################Plots########################################################################
layout(matrix(c(2,0,1,3),2,2,byrow=T), width=c(2/3, 1/3), heights=c(1/3,2/3))
xhist <- hist(n,plot=F)
yhist <- hist(theta,plot=F)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,nmax)
yrange <- c(0,1)

par(mar=c(5,5,1,1))
plot(n,theta,xlim=xrange, ylim=yrange,ylab="", xlab="")
axis(1)
mtext("Number of Surveys", side=1,line=2.25, cex=1.2)
axis(2, cex=1.2)
las=0
mtext("Rate of Return", side=2 ,line=2.25, cex=1.2)
las=1
points(mean(n),mean(theta), col="red", lwd=3, pch=4) #expectation
points(n[ind],theta[ind], col="green", lwd=3, pch=10) #Maximum Likelihood

par(mar=c(0,4,1,1))
barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0,col="lightblue")

par(mar=c(4,0,1,3))
barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE,col="lightblue")
##################
