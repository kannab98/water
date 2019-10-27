import numpy as np
import matplotlib.pyplot as plt
from numpy import pi
from scipy import interpolate,integrate
from tqdm import tqdm
import water
N = 2000
M = 1

surface = water.Surface(N=N,M=M)
x0 = np.linspace(0,200,10**5)
y0 = 0
t = 0
x, y = np.meshgrid(x0, y0)
z_real = surface.model([x,y],t)[0]
fig,ax = plt.subplots(nrows = 1, ncols = 2)

a = surface.amplitudes
k = surface.k
print(k,a)

def fft(t,x):
    N = x.size
    m = 1
    X = np.fft.rfft(x, n = m*N)
    # 2/(surface.KT[0]*x.size)
    freq = np.fft.rfftfreq(n= m*N, d = 1 )
    return freq[1:],X[1:]


ax[0].plot(x0, z_real)
freq,S = fft(x0,z_real)

S = 2*np.abs(S)/x0.size
q = np.argmax(S)
ax[1].loglog(freq,S)

plt.show()