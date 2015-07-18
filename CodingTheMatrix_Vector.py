
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




# coding: utf-8

# In[11]:

class Vec:
    def __init__(self, labels, function):
        self.D=labels
        self.f=function
v=Vec(['A', 'B', 'C'], {'A':1})

for d in v.D:
    if d in v.f:
        print(v.f[d])

def zero_vec(D):
    return Vec(D, {d:0 for d in D})

d={1,2,3,4,5,6,7,8}
dp = zero_vec(d)
print(dp)


# In[34]:

class Vec:
    def __init__(self, labels, function):
        self.D=labels
        self.f=function
v=Vec(['A', 'B', 'C'], {'A':1})

def setitem(v, d, val) :
    if d in v.D :
        v.f[d] = val
    else :
        print('out of range')

def getitem(v,d) :
    return v.f[d] if d in v.f else 0

setitem(v,'B',2)
setitem(v,'C',5)

print(v.D)
print(getitem(v, 'B'))


# In[33]:

def scalar_mul(v, alpha) :
    return ([v[i]*alpha for i in range(len(v))])
v=[1,2,3,4,5,6]
print(scalar_mul(v, 3))


# In[55]:

def setitem(v, d, val) :
    if d in v.D :
        v.f[d] = val
    else :
        print('out of range')

def getitem(v,d) :
    return v.f[d] if d in v.f else 0

def scalar_mul(v, alpha) :
    return Vec(v.D, {d:alpha*getitem(v, d) for d in v.D})
def scalar_mul2(v, alpha) :
    return Vec(v.D, {d:alpha*value for d, value in v.f.items() })
v=Vec(['A', 'B','C'], {'A':1, 'B':2, 'C':3})
vp=scalar_mul2(v,3)
getitem(vp, 'B')

def add(u, v) :
    return Vec(u.D, {getitem(u,d) + getitem(v, d) for d in u.D})
add(v,vp)
add(v,vp).f
add(v,vp).D


# In[65]:

def dot_prod(v,w) :
    return (sum(v[i]*w[i] for i in range(len(v))))
def list_dot(u, v):
    return sum(i*j for (i,j) in zip(u,v))
v=[1.,2.,3.,4.,5.]
w=[1/5, 1/5, 1/5, 1/5,1/5]
print(dot_prod(v,w))
print(list_dot(v,w))


# In[ ]:



