#pyparsing-2.0.3.tar.gz
#python-dateutil-2.4.2.tar.gz
# six-1.9.0.tar.gz
#numpy-1.9.2-win32-superpack-python3.4.exe
#matplotlib-1.4.3.win32-py3.4.exe
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
