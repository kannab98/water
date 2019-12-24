import numpy as np
import matplotlib.pyplot as plt
from numpy import fft
end = 3.2
step = 100000
x = np.arange(-3,4,1)
k0 = 1
a = 0.08

def Sg(x):
    Sg = np.zeros(len(x))
    for i in range(x.size):
        if x[i]%k0 < 1e-4:
            Sg[i] = np.pi * np.sinc(np.pi*x[i]/(2*k0))
    return Sg
def S(x):
    S = -a**2/2*(k0) * (Sg(x-2*k0) + Sg(x+2*k0) + 2*Sg(x) )
    return S

S = S(x)
fix,ax =plt.subplots()
plt.stem(x,S)
# plt.stem([0],[a])
plt.show()
