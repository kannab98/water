import numpy as np
import matplotlib.pyplot as plt
from numpy import pi
from scipy import interpolate,integrate
from tqdm import tqdm
import water

N = 1
M = 1
for i in range(2):
    surface = water.Surface(N=N,M=M)
    x0 = np.linspace(0,200,10**5)
    y0 = 0
    t = 0
    x, y = np.meshgrid(x0, y0)
    z = surface.model([x,y],t)[0]

    plt.plot(x0,z)
plt.show()
# ax[0].plot(t,heights)
# ax[1].plot(freq,np.abs(S)/t.size)
