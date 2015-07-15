# clears workspace:  
rm(list=ls()) 

#JARS
library(R2jags)

# Prior and Posterior Prediction
# 1. Model
modelString = "
# Pearson Correlation
model{
  # Data
  for (i in 1:n){
    x[i,1:2] ~ dmnorm(mu[],TI[,])
  }
  # Priors
  mu[1] ~ dnorm(0,.001)
  mu[2] ~ dnorm(0,.001)
  lambda[1] ~ dgamma(.001,.001)
  lambda[2] ~ dgamma(.001,.001)
  r ~ dunif(-1,1)
  # Reparameterization
  sigma[1] <- 1/sqrt(lambda[1])
  sigma[2] <- 1/sqrt(lambda[2])
  T[1,1] <- 1/lambda[1]
  T[1,2] <- r*sigma[1]*sigma[2]
  T[2,1] <- r*sigma[1]*sigma[2]
  T[2,2] <- 1/lambda[2]
  TI[1:2,1:2] <- inverse(T[1:2,1:2])
}"

# 2. Data(Interface Counts, Duration(Delay) Time(Minutes))
x = matrix(c(0.8, 102,
             1, 98,
             5, 100,
             0.9, 105,
             0.7,103,
             0.4, 100,
             1.2, 99,
             1.4, 87,
             0.6, 113,
             1.1, 89,
             1.3, 93), nrow=11, ncol=2, byrow=T)
cor(x[,1], x[,2]) # corelation is -0.13

n = nrow(x)
mydata = list("x", "n") # to be passed on to JAGS

# 3. Start Values
myinits = list(
  list(r = 0, mu = c(0,0), lambda = c(1,1))
)

# 4. Paramters to be monitored
parameters = c("r", "mu","sigma")

# The following command calls JAGS with specific options.
# For a detailed description see the R2jags documentation.
mcmc_samples <- jags(mydata, inits=myinits, parameters,
                     model.file =textConnection(modelString), n.chains=1, n.iter=5000, 
                     n.burnin=1, n.thin=1, DIC=T)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.
r <- mcmc_samples$BUGSoutput$sims.list$r
#Frequentist point-estimate of r:
freq.r = cor(x[,1],x[,2])
mu <- mcmc_samples$BUGSoutput$sims.list$mu
sigma <- mcmc_samples$BUGSoutput$sims.list$sigma
summary(mu)
summary(sigma)
#make the two panel plot:
par(mfrow=c(1,2))
#some plotting options to make things look better:
#par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
#     cex.axis = 1.3, bty = "n", las=1)
# data panel:    
plot(x[,1],x[,2], type="p", pch=19, cex=1)
# correlation panel:
plot(density(r, from=-1,to=1), main="", ylab="Posterior Density", xlab="Correlation", lwd=1)
lines(c(freq.r, freq.r), c(0,100), lwd=1, lty=2)
#-----------------------------------------
# Graph
par(mfrow=c(1,2))
plot(r, main='correlation trace',xlab='iteration', col='blue', lwd=1)
plot(density(r), main='correlation distribution',xlab='r', col='skyblue', lwd=2)
lines(c(freq.r, freq.r), c(0,100), lwd=1, lty=2)
#-----------------------------------------
par(mfrow=c(2,2))
plot(mu[,1], main='mean1 trace',xlab='iteration', col='blue', lwd=1)
plot(density(mu[,1]), main='mean distribution',xlab='mu[1]', col='skyblue', lwd=2)
plot(sigma[,1], main='standard deviation trace',xlab='iteration', col='blue', lwd=1)
plot(density(sigma[,1]), main='standard deviation distribution',xlab='sigma1', col='skyblue', lwd=2)

#-----------------------------------------
par(mfrow=c(2,2))
plot(mu[,2], main='mean2 trace',xlab='iteration', col='blue', lwd=1)
plot(density(mu[,2]), main='mean distribution',xlab='mu[2]', col='skyblue', lwd=2)
plot(sigma[,2], main='standard deviation trace',xlab='iteration', col='blue', lwd=1)
plot(density(sigma[,2]), main='standard deviation distribution',xlab='sigma2', col='skyblue', lwd=2)
