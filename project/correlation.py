import numpy as np
import matplotlib.pyplot as plt
from numpy import pi
from scipy import interpolate,integrate
from tqdm import tqdm
import water
N = 128
M = 1

t = 0
x0 = np.linspace(0,400,400)
y0 = 0


surface = water.Surface(N=N, M=M, whitening=1,KT=[0.05,2000])
k = surface.k
k0 = surface.k0
S = surface.spectrum(k)


def angles(rho):
    k= np.logspace(np.log10(surface.KT[0]), np.log10(surface.KT[-1]), 5*10**4)
    integral=np.zeros(len(rho))
    y=lambda k: k**2*surface.spectrum(k)
    for i in range(len(rho)):
        integral[i]=np.trapz(y(k)*np.cos(k*rho[i]),x=k)
    return integral


def angles_sum(k,rho):
    f=0
    A=surface.amplitude(k)
    k = k[:-1]
    f=np.zeros(len(rho))
    for j in range(len(rho)):
            f[j]=sum( k**2*A**2/2*np.cos(k*rho[j]) )
    return f

def height(rho,k, fourier = 'real'):
    if fourier == 'real':
        S=surface.spectrum(k)
        integral=np.zeros(len(rho))
        for i in range(len(rho)):
            integral[i]=np.trapz(S*np.cos(k*rho[i]),x=k)
        return integral

    if fourier == 'complex':
        S=surface.spectrum(k)
        integral=np.fft.fft(S,n=k.size)
        return integral
        
def height_sum(k,rho):
    f=0
    f=np.zeros(len(rho))
    A=surface.amplitude(k)
    k = k[:-1]
    for j in range(len(rho)):
            f[j]=sum( A**2/2*np.cos(k*rho[j]) )
    return f

from numpy.linalg import norm
fix,ax = plt.subplots(1,2)
rho = np.linspace(0,100,1000)
klog= np.logspace(np.log10(surface.KT[0]), np.log10(surface.KT[-1]), N)
heights = height(rho,k0)
ax[0].plot(rho,heights)
ax[0].plot(rho,height_sum(k,rho))
ax[0].plot(rho,height_sum(klog,rho))
print(np.std(heights)-np.std(height_sum(k,rho)), 
np.std(heights)-np.std(height_sum(klog,rho)))
slopes = angles(rho)
ax[1].plot(rho,slopes)
ax[1].plot(rho,angles_sum(k,rho))
ax[1].plot(rho,angles_sum(klog,rho))
print(np.std(slopes)-np.std(angles_sum(k,rho)),
np.std(slopes)-np.std(angles_sum(k,rho))) 

plt.show()