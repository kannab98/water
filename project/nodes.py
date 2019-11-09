import numpy as np
import matplotlib.pyplot as plt
from numpy import pi
from scipy import interpolate,integrate
from tqdm import tqdm
import water
N = 39
M = 1

t = 0
x0 = np.linspace(0,400,400)
y0 = 0

surface = water.Surface(N=N, M=M, whitening=1)
k = surface.k
k0 = surface.k0
fix,ax = plt.subplots(1,2)
S = surface.spectrum(k)
ax[0].loglog(k,S,'.')
ax[0].stem(surface.k_heights, surface.spectrum(surface.k_heights), 
            use_line_collection=True, markerfmt = ' ',linefmt='C0-', label='Высоты',
            bottom=min(S))

ax[0].stem(surface.k_slopes, surface.spectrum(surface.k_slopes), 
            use_line_collection=True, markerfmt = ' ', linefmt='C1-', label='Наклоны',
            bottom=min(S))
ax[1].loglog(k,k**2*surface.spectrum(k),'.')
ax[1].stem(surface.k_heights, surface.k_heights**2*surface.spectrum(surface.k_heights), 
            use_line_collection=True, markerfmt = ' ',linefmt='C0-', label='Высоты',
            bottom=min(S))

ax[1].stem(surface.k_slopes, surface.k_slopes**2*surface.spectrum(surface.k_slopes), 
            use_line_collection=True, markerfmt = ' ', linefmt='C1-', label='Наклоны',
            bottom=min(S))
ax[0].legend()
ax[1].legend()
plt.show()

