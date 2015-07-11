%matplotlib inline
import pylab as pl
import numpy as np

X = np.linspace(-np.pi, np.pi, 256, endpoint=True)
C, S = np.cos(X), np.sin(X)

pl.plot(X, C)
pl.plot(X, S)

pl.show()
#-----------------------------------------------------------
%matplotlib inline
import numpy as np
import matplotlib.pyplot as plt
limit=4
s=[2+2j, 3+2j, 1.75+1j, 2+1j, 2.25+1j, 2.5+1j, 2.75+1j, 3+1j, 3.25+1j]
for x in range(len(s)) :
	plt.plot(s[x].real, s[x].imag, 'ro-', label='complex')
	##  plt.plot([0,s[x].real], [0, s[x].imag])
plt.xlim(-limit,limit)
plt.ylim(-limit,limit)
plt.ylabel('Imaginary')
plt.xlabel('Real')	
plt.grid(True)
plt.show()
#-----------------------------------------------------------
# CSV file reader
import csv

cars=[]
f = open('auto-mpg.csv', 'rb')
csvReader  = csv.reader(f)
for row in csvReader :
   cars.append(row)
#print(cars[0])
print(cars[len(cars)-1])
f.close()
#-----------------------------------------------------------
%matplotlib inline
import matplotlib.pyplot as plt
import numpy as np
import pylab as plt
import csv as csv
import pandas as pd

readdata = csv.reader(open('C:/Work/auto-mpg.csv'))
data=[]
for row in readdata :
    data.append(row)
Header=data[0]
data.pop(0)

cars = pd.DataFrame(data, columns=Header)
plt.plot(cars.mpg, cars.displ, 'ro', label='cars')
plt.ylabel('displ')
plt.xlabel('mpg')

plt.show()
#------------------------------------------
# Task 1.6.6
def makeInverseIndex(strlist) :
    return {x:y for (x,y) in enumerate(list(strlist.split()))}
makeInverseIndex("Ask not what you can do for your country")
