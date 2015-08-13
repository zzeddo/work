# Reference : http://ezstat.snu.ac.kr/textbook_sources/chapter_03.pdf
# t-test example
# 1. Group Selection(2 Group of Monthly Salary)
Winter = c(-0.05,0.41,0.17,-0.13,0.00,-0.05,0.00,0.17,0.29,0.04,0.21,0.08,0.37,0.17,0.08,-0.04,-0.04,0.04,-0.13,-0.12,0.04,0.21,0.17,
            0.17,0.17,0.33,0.04,0.04,0.04,0.00,0.21,0.13,0.25,-0.05,0.29,0.42,-0.05,0.12,0.04,0.25,0.12)

Summer = c(0.00,0.38,-0.12,0.12,0.25,0.12,0.13,0.37,0.00,0.50,0.00,0.00,-0.13,-0.37,-0.25,-0.12,0.50,0.25,0.13,0.25,0.25,0.38,0.25,0.12,
            0.00,0.00,0.00,0.00,0.25,0.13,-0.25,-0.38,-0.13,-0.25,0.00,0.00,-0.12,0.25,0.00,0.50,0.00)
summary(Winter)
summary(Summer)
sd(Winter) # Standard Deviation Winter
sd(Summer) # Standard Deviation Summer

# 2. Density
plot(density(Winter))
lines(density(Summer), lty=2)

# 3. hypothesis
# H0 : u2-u1 vs H1 : u2-u1 <> 0

# 4. t-test
t.test(Winter,Summer, var.equal=T, alt='two.sided')

# Cauchy Distributon
x=seq(from=-2, to=2, by=0.001)
y=1/(pi*(1+x^2))
plot(x,y)


##############################
#  8.1 One Sample comparison
##############################
# clears workspace:  
rm(list=ls()) 

#JARS
library(R2jags)

# Read data Dr. Smith
Winter = c(-0.05,0.41,0.17,-0.13,0.00,-0.05,0.00,0.17,0.29,0.04,0.21,0.08,0.37,0.17,0.08,-0.04,-0.04,0.04,-0.13,-0.12,0.04,0.21,0.17,
            0.17,0.17,0.33,0.04,0.04,0.04,0.00,0.21,0.13,0.25,-0.05,0.29,0.42,-0.05,0.12,0.04,0.25,0.12)

Summer = c(0.00,0.38,-0.12,0.12,0.25,0.12,0.13,0.37,0.00,0.50,0.00,0.00,-0.13,-0.37,-0.25,-0.12,0.50,0.25,0.13,0.25,0.25,0.38,0.25,0.12,
            0.00,0.00,0.00,0.00,0.25,0.13,-0.25,-0.38,-0.13,-0.25,0.00,0.00,-0.12,0.25,0.00,0.50,0.00)

x = Winter-Summer # allowed because it is a within-subjects design
x = x/sd(x)       # standardize
ndata = length(Winter) # number of subjects
data  = list("x", "ndata") # to be passed on to JAGS

# inital value
myinits <- list(
  list(delta = rnorm(1,0,3), sigmatmp = rnorm(1,0,1)),
  list(delta = rnorm(1,0,3), sigmatmp = rnorm(1,0,1)),
  list(delta = rnorm(1,0,3), sigmatmp = rnorm(1,0,1)))

# Parameters to be monitored
parameters <- c("delta")

# 1. Model
# One-Sample Comparison of Means
modelString = "
model{
  # Data
  for (i in 1:ndata){
    x[i] ~ dnorm(mu,lambda)
  } 
  mu <- delta*sigma
  lambda <- pow(sigma,-2)
  # delta and sigma Come From (Half) Cauchy Distributions
  lambdadelta ~ dchisqr(1)
  delta ~ dnorm(0,lambdadelta)
  lambdasigma ~ dchisqr(1)
  sigmatmp ~ dnorm(0,lambdasigma)
  sigma <- abs(sigmatmp)
  # Sampling from Prior Distribution for Delta
  deltaprior ~ dnorm(0,lambdadeltaprior)
  lambdadeltaprior ~ dchisqr(1)
}"
