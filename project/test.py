import numpy as np
from numpy.fft import ifft
import matplotlib.pyplot as plt



def f(x, A = 0.07, k = 1):
    f = A*( np.cos(k*x) - A*k*(1 - np.cos(2*k*x) ) 
    return f

x = np.linspace(-2*np.pi,2*np.pi,1000)
A = 0.15
k = 1
plt.cla()
plt.plot(x,f(x,A=A,k=k),label='CWM')
plt.plot(x,A*np.cos(k*x))
plt.legend()
plt.show()
