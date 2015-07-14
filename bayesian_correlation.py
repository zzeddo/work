%matplotlib inline
%reset -f
from pymc import *
import numpy as np
from numpy.linalg import inv
import pandas as pd
import matplotlib.pyplot as plt

# 1. Data
x=np.array([[87030, 2549],
     [79571, 2244],
     [67385, 2555],
     [72402, 2110],
     [70816, 1406],
     [62985, 942],
     [131463, 2504]])
n = len(x)
plt.plot([x[i][0] for i in range(n)], [x[i][1] for i in range(n)], 'ro')
mean = x.mean(1)
# 2. Model
# Pearson Correlation
#Priors
mu1 = Normal('mu1', mu=0, tau=0.001, value=mean[0], observed=True)
mu2 = Normal('mu2', mu=0, tau=0.001, value=mean[1], observed=True)
lambda1 = Gamma('lambda1', alpha=0.001, beta=0.001)
lambda2 = Gamma('lambda2', alpha=0.001, beta=0.001)
r = Uniform('r', lower=-1, upper = 1, value=0)

@pymc.deterministic
def mean(mu1=mu1, mu2=mu2):
       return np.array([mu1, mu2])

@pymc.deterministic
def precision(lambda1=lambda1, lambda2=lambda2, rho=r):
    sigma1 = 1/sqrt(lambda1)
    T11 = 1/lambda1
    T12 = r*sigma1*sigma2
    T21 = r*sigma1*sigma2
    T22 = 1/lambda2
    #return np.power(np.mat([[ss1, rss], [rss, ss2]]), -1)
    return np.mat(inv([[T11, T12], [T21, T22]]))

# Observed Counts
xy = MvNormal('xy', mu=mean, tau=precision, value=x.T, observed=True)
    
#3. MCMC sampling
S = pymc.MCMC(locals())
S.sample(iter = 1000, burn = 1, thin = 1)
S.db.close()
pymc.Matplot.plot(S)
