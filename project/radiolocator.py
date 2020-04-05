import matplotlib.pyplot as plt
from numpy.linalg import norm
import numpy as np
from numpy import pi
from scipy import interpolate,integrate
from scipy.optimize import curve_fit
from scipy.special import erf
from scipy.optimize import anderson
from scipy.optimize import root

from pandas import read_csv
xrad = yrad = 0

def U10(sigma0):
    if sigma0 < 10.12:
        A =  0.080074
        B = -0.12465
    elif sigma0 < 19.9:
        A = 0.03989
        B = -0.03199
    else:
        A = 0.01595
        B = 0.017215

    return (np.exp(10**-(0.21 + 0.1*sigma0)) - B)/A

class Radiolocator():
    def __init__(self, h=1e6, xi=0.0, theta=1.5, c=3e8, sigma=1, 
                    angles_in='degrees', pulse = np.array([None]), t = None):
        self.R = 6370e3 # Радиус Земли в метрах
        self.c = c # Скорость света в м/с
        if angles_in=='degrees':
            self.xi = np.deg2rad(xi) # Отклонение антенны в радианах
            self.theta = np.deg2rad(theta) # Ширина диаграммы направленности в радианах
            self.Gamma = self.gamma(self.theta)

        if pulse.any() != None:
            params = radiolocator.calc(t,pulse)
            print('xi = {},\nh = {},\nsigma_s = {}'.format(params[0],params[1],params[2]))

        else:
            self.h = h # Высота орбиты в метрах
            self.sigma_s = sigma # Дисперсия наклонов


        self.T = 3e-9 # Временное разрешение локатора


    def H(self,h):
        return h*( 1+ h/self.R )
    
    def A(self,gamma,xi,A0=1):
        return A0*np.exp(-4/gamma * np.sin(xi)**2 )

    def u(self,t,alpha,sigma_c):
        return (t - alpha * sigma_c**2) / (np.sqrt(2) * sigma_c)

    def v(self,t,alpha,sigma_c):
        return alpha*(t - alpha/2 * sigma_c**2)

    def alpha(self,beta,delta):
        return delta - beta**2/4

    def delta(self,gamma,xi,h):
        return 4/gamma * self.c/self.H(h) * np.cos(2 * xi)
    
    def gamma(self,theta):
        return 2*np.sin(theta/2)**2/np.log(2)

    def beta(self,gamma,xi,h):
        return 4/gamma * np.sqrt(self.c/self.H(h)) * np.sin(2*xi)


    def sigma_c(self,sigma_s):
        sigma_p = 0.425 * self.T 
        return np.sqrt(sigma_p**2 + (2*sigma_s/self.c)**2 )

    def pulse(self,t, dim = 1):

        self.dim = dim
        gamma = self.Gamma
        delta = self.delta(gamma,self.xi,self.h)
        beta  = self.beta(gamma,self.xi,self.h)

        if dim == 1:
            alpha = self.alpha(beta,delta)
        else:
            alpha = self.alpha(beta/np.sqrt(2),delta)

        sigma_c = self.sigma_c(self.sigma_s)

        u = self.u(t, alpha, sigma_c)
        v = self.v(t, alpha, sigma_c)

        A = self.A(gamma,self.xi)
        pulse = A*np.exp(-v)*( 1 + erf(u) )
        
        if self.dim == 2:
            alpha = gamma
            u = self.u(t, alpha, sigma_c)
            v = self.v(t, alpha, sigma_c)
            pulse -= A/2*np.exp(-v)*( 1 + erf(u) )

        return pulse
    def pulse_min(self):
        gamma = self.Gamma
        delta = self.delta(gamma,self.xi,self.h)
        beta  = self.beta(gamma,self.xi,self.h)

        alpha = self.alpha(beta,delta)

        sigma_c = self.sigma_c(self.sigma_s)

        return np.sqrt(2*sigma_c**2 * np.log(1/(np.sqrt(2*np.pi*sigma_c**2 *alpha**2))))  + alpha*sigma_c**2
        


    def calc(self,t,pulse):
        self.brown = lambda t,A,alpha,sigma:  A*np.exp(-alpha*(t - alpha/2 * sigma**2))*(1 + erf( (t - alpha*sigma**2)/(2**0.5*sigma))) 
        
        theta = np.deg2rad(1.5)
        gamma = 2*np.sin(theta/2)**2/np.log(2)

        popt,pcov = curve_fit(self.brown, 
                            xdata=t,
                            ydata=pulse,
                            p0=[1,2e6,0],
                            bounds=( (0,1e5,0), (2,4e6,np.inf) )
                        ) 
        print(popt)
        c = self.c
        R = self.R

        tmp = np.log(popt[0])
        if tmp > 0 :
            xi = 0
        else:
            xi = np.arcsin( np.sqrt( -gamma/4 * tmp ))


        h0 = 4*c/gamma/popt[1]  * (np.cos(2*xi) - np.sin(2*xi)**2/gamma)
        h  = -R/2 + np.sqrt(R**2 +4*h0*R)/2

        sigma_s = (popt[2]**2 - (self.T*0.425) **2)*c**2/4
        return np.rad2deg(xi), h, sigma_s, P,pcov


df = read_csv('imp_5.dat',sep='\s+',header=None).T

t = df.iloc[0].values*1e-9
pulse = df.iloc[1].values
pulse=2*pulse/max(pulse)
tmax = t[np.argmax(pulse)]

plt.plot(t- tmax,pulse,label='практика')


radiolocator = Radiolocator(xi = 0, sigma=0.0)
# t = np.linspace(-5e-8, 2e-7, 25600) 

# pulse = radiolocator.pulse(t)     

# tmax1 = radiolocator.pulse_min()
# A = max(pulse) - min(pulse)
# err = 0.01

# plt.plot(t - tmax1,pulse,label='теория')
# pulse +=  np.random.uniform(-err*A,err*A,size=t.size)
# plt.plot(t,pulse,':',label='теория + шум')

params = radiolocator.calc(t=t,pulse=pulse)


# radiolocator = Radiolocator(xi = params[0], h = params[1], sigma = params[2])
# pulse = radiolocator.pulse(t)     
# plt.plot(t,pulse,label='аппроксимация')
# plt.legend()
# # plt.title( 'xi = {:.2},\nh = {:.1},\nsigma_s = {:.3}, sigma0 = {:.2}'.format(params[0],params[1],params[2],params[3]))
plt.show()
# # print(U10(params[3]))

