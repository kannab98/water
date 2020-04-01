import matplotlib.pyplot as plt
from numpy.linalg import norm
import numpy as np
from numpy import pi
from scipy import interpolate,integrate
from tqdm import tqdm
import water

xrad = yrad = 0
H0 = 1e3;
def ionsphere():
    pass

def troposphere_dry():
    pass
def troposphere_wet():
    pass

def spherical_earth():
    pass


def U10(sigma0):
    if sigma0 < 10.12:
        A =  0.080074
        B = -0.12465
    elif sigma0 < 10.9:
        A = 0.03989
        B = -0.03199
    else:
        A = 0.01595
        B = 0.017215

    return (np.exp(10**(0.21 + 0.1*sigma0)) - B)/A
    


M = 256
N = 256
x0 = np.linspace(0,200,50)
y0 = x0
t = 0


fig,ax = plt.subplots(nrows = 1, ncols = 1)
surface = water.Surface(N=N,M=M,U10=5,wind= np.pi/6)
x, y = np.meshgrid(x0, y0)
z = surface.model([x,y],t)
print(z.shape)
from matplotlib.cm import winter
plt.contourf(x,y,z,levels=100,cmap=winter)
plt.colorbar()

plt.ylabel('Y, м',fontsize=16)
plt.xlabel('X, м',fontsize=16)
plt.show()