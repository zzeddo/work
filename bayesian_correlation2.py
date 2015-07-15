%matplotlib inline
%reset -f
from pymc import *
import numpy as np
from numpy.linalg import inv
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
# 1. Data
x = np.array([87030, 79571, 67385, 72402, 70816, 62985, 131463])
y = np.array([2549, 2244, 2555, 2110, 1406, 942, 2504])
data = np.array([x,y])
n = len(x)
mean = [x.mean(), y.mean()]
plt.plot(x, y, 'ro');
pearsonr(x, y)

# 2. Model
# Pearson Correlation
#Priors
mu1 = Normal('mu1', mu=0, tau=0.001, value=mean[0])
mu2 = Normal('mu2', mu=0, tau=0.001, value=mean[1])
lambda1 = Gamma('lambda1', alpha=0.001, beta=0.001)
lambda2 = Gamma('lambda2', alpha=0.001, beta=0.001)
r = Uniform('r', lower=-1, upper = 1, value=0)

@pymc.deterministic
def mean(mu1=mu1, mu2=mu2):
       return np.array([mu1, mu2])

@pymc.deterministic
def precision(lambda1=lambda1, lambda2=lambda2, rho=r):
    sigma1 = 1/sqrt(lambda1)
    sigma2 = 1/sqrt(lambda2)
    T11 = 1/lambda1
    T12 = r*sigma1*sigma2
    T21 = r*sigma1*sigma2
    T22 = 1/lambda2
    return np.mat(inv([[T11, T12], [T21, T22]]))

# Observed Counts
xy = MvNormal('xy', mu=mean, tau=precision, value=data.T, observed=True)
    
#3. MCMC sampling
S = pymc.MCMC(locals())
S.sample(iter = 10000, burn = 1, thin = 1)
