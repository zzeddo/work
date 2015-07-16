
# coding: utf-8

# In[1]:

from math import *
a=[pi, e, -1.0,2.0]
print(a)


# In[2]:

b=list(enumerate(a))
print(b)
{i:j for (i,j) in b}


# In[13]:

S="The rain in Spain falls mainly on the plain"
L=S.split()
wordcount={}
for word in L :
    if word not in wordcount :
        wordcount[word] = 1
    else :
        wordcount[word] += 1
print(wordcount)


# In[17]:

get_ipython().magic(u'matplotlib inline')
import pandas as pd
import matplotlib.pyplot as plt

L=[[2,2], [3,2], [1.75,1], [2,1], [2.25,1], [2.5, 1], [2.75, 1], [3,1], [3.25,1]]
limit = 4
n = len(L)
plt.plot([L[i][0] for i in range(n)], [L[i][1] for i in range(n)], 'ro')
plt.xlim(0,limit)
plt.ylim(0,limit)
plt.xlabel('x')	
plt.ylabel('y')
plt.grid(True)


# In[8]:

a=[1,2]
b=[3,4]
print(a+b)
def add2(v, w) :
    return [v[0]+w[0], v[1]+w[1]]
print(add2(a,b))


# In[13]:

L=[[2,2], [3,2], [1.75,1], [2,1], [2.25,1], [2.5, 1], [2.75, 1], [3,1], [3.25,1]]
LP=[[i+1, j+2] for (i,j) in L]
print(LP)
def add2(v, w) :
    return [v[0]+w[0], v[1]+w[1]]
LP = [add2(v, [1,2]) for v in L]
limit = 5
n = len(LP)
plt.plot([LP[i][0] for i in range(n)], [LP[i][1] for i in range(n)], 'ro')
plt.xlim(0,limit)
plt.ylim(0,limit)
plt.xlabel('x')	
plt.ylabel('y')
plt.grid(True)


# In[17]:

def addn(v, w):
    return [v[i]+w[i] for i in range(len(v)) ]
a=[1,2]
b=[3,4]
print(addn(a,b))


# In[52]:

def add2(v, w) :
    return [v[0]+w[0], v[1]+w[1]]
limit = 4
#plt.plot([LP[i][0] for i in range(n)], [LP[i][1] for i in range(n)], 'ro')
plt.figure()
P = [3,1.5]
PP=add2(P, [-1, -1])
plt.quiver(0,0, 3, 1.5, scale_units='xy', angles='xy', scale=1)
plt.quiver(-1, -1,2, 0.5, scale_units='xy', angles='xy', scale=1)

plt.xlim(0,limit)
plt.ylim(0,limit)
plt.xlabel('x')	
plt.ylabel('y')
plt.grid(True)


# In[6]:

get_ipython().magic(u'matplotlib inline')
import pandas as pd
import matplotlib.pyplot as plt
L=[[2,2], [3,2], [1.75,1], [2,1], [2.25,1], [2.5, 1], [2.75, 1], [3,1], [3.25,1]]
def prd2(v, w) :
    return [v[0]*w[0], v[1]*w[1]]
LP = [prd2(v, [-1,-1]) for v in L]
limit = 5
n = len(LP)
plt.plot([L[i][0] for i in range(n)], [L[i][1] for i in range(n)], 'ro')
plt.plot([LP[i][0] for i in range(n)], [LP[i][1] for i in range(n)], 'bo')
plt.xlim(-limit, limit)
plt.ylim(-limit, limit)
plt.xlabel('x')	
plt.ylabel('y')
plt.grid(True)


# In[1]:

x=[1,2,3,4,5,6,7]
print(x[:3])
print(x[5:])


# In[2]:

[1,2]*2


# In[9]:

def scalar_vector_mult(alpha, v):
    return [v[i]*alpha for i in range(len(v))]

v=[1,2,3,4,5]
scalar_vector_mult(2, v)


# In[12]:

get_ipython().magic(u'matplotlib inline')
import pandas as pd
import matplotlib.pyplot as plt
L=[[2,2], [3,2], [1.75,1], [2,1], [2.25,1], [2.5, 1], [2.75, 1], [3,1], [3.25,1]]
def scalar_vector_mult(alpha, v):
    return [v[i]*alpha for i in range(len(v))]
LP1 = [scalar_vector_mult(0.5, v) for v in L]
LP2 = [scalar_vector_mult(-0.5, v) for v in L]
limit = 4
n = len(LP1)
plt.plot([L[i][0] for i in range(n)], [L[i][1] for i in range(n)], 'ro')
plt.plot([LP1[i][0] for i in range(n)], [LP1[i][1] for i in range(n)], 'bo')
plt.plot([LP2[i][0] for i in range(n)], [LP2[i][1] for i in range(n)], 'bo')
plt.xlim(-limit, limit)
plt.ylim(-limit, limit)
plt.xlabel('x')	
plt.ylabel('y')
plt.grid(True)


# In[16]:

get_ipython().magic(u'matplotlib inline')
import pandas as pd
import matplotlib.pyplot as plt
def scalar_vector_mult(alpha, v):
    return [v[i]*alpha for i in range(len(v))]
limit = 4
P = [3.0,1.5]
PS = [scalar_vector_mult(i/10.0, P) for i in range(11)]
print(PS)
plt.quiver(0,0, 3, 1.5, scale_units='xy', angles='xy', scale=1)
plt.plot([PS[i][0] for i in range(n)], [PS[i][1] for i in range(n)], 'ro')

plt.xlim(0,limit)
plt.ylim(0,limit)
plt.xlabel('x')	
plt.ylabel('y')
plt.grid(True)


# In[ ]:



