import numpy as np
import matplotlib.pyplot as plt
from numpy import pi
from scipy import interpolate,integrate
from tqdm import tqdm
import water
N = 256
M = 256

t = 0


def plot_surface(x , y, t):
    fig,ax = plt.subplots(nrows = 1, ncols = 1)
    surface = water.Surface(N=N,M=M,U10=10)
    x, y = np.meshgrid(x, y)
    z = surface.model([x,y],t)
    print(z.shape)
    from matplotlib.cm import winter
    plt.contourf(z,100,cmap=winter)
    plt.colorbar()
    plt.ylabel('Y, м',fontsize=16)
    plt.xlabel('X, м',fontsize=16)
    plt.savefig('water5.png', transparent=True, pdi=10**6,)

x0 = np.linspace(0,400,400)
y0 = x0
plot_surface(x0,y0,t)
plt.show()
