import numpy as np
import matplotlib.pyplot as plt
from numpy import pi
from scipy import interpolate,integrate
from tqdm import tqdm
import water
from matplotlib import rc
plt.rcParams['axes.labelsize'] = 20
rc('text',usetex=True)
rc('text.latex',preamble=[r'\usepackage[russian]{babel}',r'\usepackage{amsmath}'])
rc('font',family = 'serif')
N = 64
M = 1

t = 0
x0 = np.linspace(0,400,400)
y0 = 0

surface = water.Surface(N=N, M=M, whitening='hs')
k = surface.k
k0 = surface.k0
S = surface.spectrum(k)

plt.figure(figsize=[8,6])
plt.loglog(k0,surface.spectrum(k0),'-',color='black')
plt.stem(surface.k_slopes, surface.spectrum(surface.k_slopes),
            use_line_collection=True, markerfmt = ' ', linefmt='r', label='Наклоны',
            bottom=0)

plt.stem(surface.k_heights[0:-2], surface.spectrum(surface.k_heights[0:-2]),
            use_line_collection=True, markerfmt = ' ', linefmt='darkblue',
            label='Высоты',
            bottom=0)

plt.xlabel(r'$k, \text{ рад}\cdot\text{м}^{-1}$')
plt.ylabel(r'S,~\text{a.u.}')
plt.grid(which='major', linestyle='-')
plt.grid(which='minor', linestyle=':')
plt.minorticks_on()
plt.legend()
plt.savefig('heights_slopes.png',dpi=1000)
plt.savefig('heights_slopes.pdf')

plt.show()

