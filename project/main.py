import numpy as np
import matplotlib.pyplot as plt
from numpy import pi
from scipy import interpolate,integrate
from tqdm import tqdm
import water
N = 2000
M = 100

t = 0


def plot_surface(x , y, t):
    fig,ax = plt.subplots(nrows = 1, ncols = 1)
    surface = water.Surface(N=N,M=M)
    x, y = np.meshgrid(x, y)
    z = surface.model([x,y],t)
    print(z.shape)
    from matplotlib.cm import winter
    plt.contourf(z,100,cmap=winter)
    # plt.colorbar()
    # # plt.ylabel(r'Y, \text{м}',fontsize=16)
    # plt.xlabel(r'X, \text{м}',fontsize=16)

x0 = np.linspace(0,200,200)
y0 = x0
plot_surface(x0,y0,t)
plt.show()
