import matplotlib.pyplot as plt
from numpy.linalg import norm
import numpy as np
from numpy import pi
from scipy import interpolate,integrate
from tqdm import tqdm
import water
from matplotlib import rc
plt.rcParams['axes.labelsize'] = 20
rc('text',usetex=False)
rc('text.latex',preamble=[r'\usepackage[russian]{babel}',r'\usepackage{amsmath}'])
rc('font',family = 'serif')



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
NN = [i for i in range(64,512,2)]
H=[]
S0=[]
Hlog=[]
Slog=[]
for i in NN:
    N = i
    M = 1

    t = 0
    x0 = np.linspace(0,400,400)
    y0 = 0
    surface = water.Surface(N=N, M=M, whitening='s', KT=[0.05,2000])
    k = surface.k
    k0 = surface.k0
    S = surface.spectrum(k)
    rho = np.linspace(0,100,1000)
    klog= np.logspace(np.log10(surface.KT[0]), np.log10(surface.KT[-1]), N)
    heights = height(rho,k0)


    # plt.plot(rho,height_sum(k,rho),label='a')
    # plt.plot(rho,height_sum(klog,rho),label='b')
    # plt.plot(rho,heights,label='c')
    # plt.legend()
    # plt.savefig('/home/kannab/documents/water/poster/fig/corr1.pdf')
    H.append(np.linalg.norm(heights-height_sum(k,rho)))
    Hlog.append(np.linalg.norm(heights-height_sum(klog,rho)))

    # plt.figure('2')
    slopes = angles(rho)
    # plt.plot(rho,angles_sum(k,rho),label='a')
    # plt.plot(rho,angles_sum(klog,rho),label='b')
    # plt.plot(rho,slopes,label='c')
    S0.append(np.linalg.norm(slopes-angles_sum(k,rho)))
    Slog.append(np.linalg.norm(slopes-angles_sum(klog,rho)))
    # plt.legend()
    # plt.savefig('/home/kannab/documents/water/poster/fig/corr2.pdf')
    plt.show()
plt.figure('1')
plt.plot(H,label='white')
plt.plot(Hlog,label='log')
plt.legend()

plt.figure('2')
plt.plot(S0,label='white')
plt.plot(Slog,label='log')
plt.legend()
plt.show()
