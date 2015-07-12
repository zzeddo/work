
# coding: utf-8

# In[18]:

get_ipython().magic('matplotlib inline')


# In[19]:

import numpy as np
import pandas as pd
from scipy import stats, integrate
import matplotlib.pyplot as plt


# In[20]:

import seaborn as sns
sns.set(color_codes=True)


# In[21]:

# setting the seed number(defult mean : 0, standard devariance : 1)
np.random.seed(sum(map(ord, "distributions")))


# In[88]:

# normal distribution(standard normal distribution)
mu = 0
sigma = 1
x = np.random.normal(loc=mu, scale=sigma, size=100)
sns.distplot(x);


# In[57]:

# Uniform distribution
x=np.random.uniform(low=0, high=1, size=1000)
sns.distplot(x);


# In[68]:

# binomial distribution
# n is number of trial
# p is possibility
x=np.random.binomial(n=10, p=9/10, size=1000)
sns.distplot(x)


# In[81]:

# Sine, Consine Graph
x= np.linspace(-np.pi, np.pi, 256, endpoint=True)
sns.set()
plt.plot(x, np.sin(x), x, np.cos(x));


# In[86]:

data = np.random.multivariate_normal([0, 0], [[5, 2], [2, 2]], size=2000)
data = pd.DataFrame(data, columns=['x', 'y'])
print('xy')
for col in 'xy':
    plt.hist(data[col], normed=True, alpha=0.5)


# In[ ]:



