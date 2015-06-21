##############################
# 4 Inference with Gaussians
##############################
# Uniform Distribution
# p(x) = choose(n, x) p^x (1-p)^(n-x)
# choose(n, x) = n!/x!(n-x)!
theta = dbinom(1:100,100,1/2)
summary(theta)
par(mfrow=c(1,2))
plot(theta, main='Uniform trace',xlab='iteration', col='blue', lwd=1)
plot(density(theta), main='Uniform distribution',xlab='theta', col='skyblue', lwd=2)

?dbinom

# Beta Distribution
theta = dbeta(1:100,1,1,1/2)
summary(theta)
par(mfrow=c(1,2))
plot(theta, main='Uniform trace',xlab='iteration', col='blue', lwd=1)
plot(density(theta), main='Uniform distribution',xlab='theta', col='skyblue', lwd=2)

exp(1)
?exp
