import numpy as np
import matplotlib.pyplot as plt
from numpy import pi
from scipy import interpolate,integrate
from tqdm import tqdm
import water
N = 1000
M = 1
surface = water.Surface(N=N,M=M)
print(surface.k_m)
print(surface.sigma_sqr)
# x0 = np.linspace(0,400,10**5)
# y0 = 0
# t = 0
# x, y = np.meshgrid(x0, y0)
# z_real = surface.model([x,y],t)[0]
# fig,ax = plt.subplots(nrows = 1, ncols = 1)

# Dx,Dy = surface.D([x,y])
# ax.plot(x0+Dx[0],z_real)
# ax.plot(x0,z_real)
# # # ax.plot(ans)
# plt.show()

