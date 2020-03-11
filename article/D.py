import numpy as np
import matplotlib.pyplot as plt 
A = 0.08
pi = np.pi
k0 = 8*pi
def D(k):
    return  2*A/pi*(k0*k0)/(k**2-k0**2)*np.cos(pi*k/(2*k0))
k = np.linspace(-100*k0,100*k0,10000)
S = A*pi*(D(k-k0*np.ones(k.size)) + D(k+k0*np.ones(k.size)))
plt.plot(k,S)
plt.show()
