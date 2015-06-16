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

# Collect posterior samples:
delta <- mcmc_samples$BUGSoutput$sims.list$delta
#plot(delta)
summary(delta)
# mean of delta:
mean(delta)
# median of delta:
median(delta)
# mode of delta, estimated from the "density" smoother:
density(delta)$x[which(density(delta)$y==max(density(delta)$y))]
# 95% credible interval for delta:
quantile(delta, c(.025,.975))

# Now let's plot a histogram for delta. 
# First, some options to make the plot look better:
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
    font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
Nbreaks <- 80
y       <- hist(delta, Nbreaks, plot=F)
plot(c(y$breaks, max(y$breaks)), c(0,y$density,0), type="S", lwd=2, lty=1,
     xlim=c(-1,1), ylim=c(0,4), xlab="Difference in Rates", ylab="Posterior Density") 
