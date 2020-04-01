import numpy as np
from numpy import pi

from matplotlib.cm import winter
import matplotlib.pyplot as plt


# Вот тут должны быть нужные тебе массивы. (Если они одномерные)
x0 = np.linspace(0,200,200)
y0 = x0

# Тут эти массивы превращаются в сетку
x, y = np.meshgrid(x0, y0)

# Вот тут твой двумерный массив
z = np.sin(x)**2



fig,ax = plt.subplots(nrows = 1, ncols = 1)
plt.contourf(x,y,z,levels=100,cmap=winter)
plt.colorbar()
plt.ylabel(r'$Y$, м',fontsize=16)
plt.xlabel(r'$X$, м',fontsize=16)

plt.show()
