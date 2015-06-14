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
               codaPkg = F, debug=F)

samples$sims.array
plot(samples)
samples$sims.list



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
  data = mydata, inits=myinits, n.chains = 2, n.adapt= 20000)
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
