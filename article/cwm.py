import numpy as np
import matplotlib.pyplot as plt
from numpy import fft
end = 3
step = 1000
x = np.linspace(0,end,step)
k = 2*np.pi
a = 0.08
y = a * np.cos(k*x) 
N = y.size


def D(x,a=1):
    D = np.zeros(x.size)
    for i in range(x.size):
        D[i] = a*np.sin(k*x[i])
    return D

def parD(x,a=1):
    D = np.zeros(x.size)
    for i in range(x.size):
        D[i] = a*np.cos(k*x[i])
    return D


fig, ax = plt.subplots(nrows = 1, ncols = 1)
# Surface
D = D(x,a = a)
parD = parD(x, a=a)
for i in range(D.size):
    if D[i] < 0:   
        D[i] = 0

for i in range(D.size):
    if parD[i] < 0:   
        parD[i] = 0
y = a * np.cos(k*(x)) 
ax.plot(x-D, y,'-',label='CWM')
ax.plot(x,y,'-',label='Default')
ax.legend()
ax.set_xlabel('k/k0')
ax.set_ylabel('z(x), м')
print(np.trapz(y,x = x))
print(np.trapz(y,x = x - D))
# plt.savefig('cwm_surface.pdf')
fig, ax = plt.subplots(nrows = 1, ncols = 1)
# Spectrum
freq = np.fft.rfftfreq(n = N, d = end/step)

# S = fft.rfft(y)
# S[0]*=0.5
# ax.plot(freq,2*np.abs(S)/x.size,'.',label='Standart', )
S = fft.rfft(y*(1-parD))
S[0]*=0.5
ax.stem(freq,2*np.abs(S)/x.size,markerfmt = '.',label='Численно', )

freq = [0,1,2]
S = [a**2/2, a, a**2/2] 
ax.stem(freq,S,label='Аналитически',markerfmt ='or')
ax.legend()
ax.set_xlabel('k/k0')
ax.set_ylabel('z(x)/pi, м')
ax.set_xlim(-0.2,2.2)
plt.legend()
# plt.savefig('cwm_spectrum.pdf')
plt.show()
