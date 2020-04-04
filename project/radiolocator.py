import matplotlib.pyplot as plt
from numpy.linalg import norm
import numpy as np
from numpy import pi
from scipy import interpolate,integrate
from scipy.optimize import curve_fit
from scipy.special import erf

xrad = yrad = 0
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
    elif sigma0 < 19.9:
        A = 0.03989
        B = -0.03199
    else:
        A = 0.01595
        B = 0.017215

    return (np.exp(10**(0.21 + 0.1*sigma0)) - B)/A

class Radiolocator():
    def __init__(self, h=1e6, xi=0, theta=1.5, c=3e8, sigma=1, angles_in='degrees'):
        if angles_in=='degrees':
            self.xi = np.deg2rad(xi) # Отклонение антенны в радианах
            self.theta = np.deg2rad(theta) # Ширина диаграммы направленности в радианах
        
        self.h = h # Высота орбиты в метрах
        self.c = c # Скорость света в м/с
        self.sigma_s = sigma # Дисперсия наклонов
        self.R = 6370e3 # Радиус Земли в метрах
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
        gamma = self.gamma(self.theta)
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

    def pulse1(self, v, dim = 1):

        self.dim = dim
        gamma = self.gamma(self.theta)
        delta = self.delta(gamma,self.xi,self.h)
        beta  = self.beta(gamma,self.xi,self.h)

        alpha = self.alpha(beta,delta)

        sigma_c = self.sigma_c(self.sigma_s)

        u = np.sqrt(2)*v/(alpha*sigma_c) - alpha*sigma_c/np.sqrt(2)

        A = self.A(gamma,self.xi)
        print(A, sigma_c)
        pulse = A*np.exp(-v)*( 1 + erf(u) )

        return pulse

    def pulse2(self,t, theta, xi, h, sigma_s):


        xi = np.deg2rad(xi) # Отклонение антенны в радианах
        theta = np.deg2rad(theta) # Ширина диаграммы направленности в радианах
        gamma = self.gamma(theta)
        delta = self.delta(gamma,xi,h)
        beta  = self.beta(gamma,xi,h)

        alpha = self.alpha(beta,delta)

        sigma_c = self.sigma_c(sigma_s)

        u = self.u(t, alpha, sigma_c)
        v = self.v(t, alpha, sigma_c)

        A = self.A(gamma,xi)
        pulse = A*np.exp(-v)*( 1 + erf(u) )

        return pulse



alpha = 2097955
sigma = 6.7874e-9
brown = lambda t,A,alpha,sigma,T:  A*np.exp(-alpha*(t - alpha/2 * sigma**2))*(1 + erf( (t - alpha*sigma**2)/(2**0.5*sigma))) + T
ice = lambda t,A,tau,sigma,S:  A*np.exp(S*(t - tau))*(1 + erf( (t - tau)/(sigma)))
t = np.linspace(-5e-8, 2e-7, 256)
v = np.linspace(-0.1,0.5,1000)
radiolocator = Radiolocator()
plt.figure(1)
pulse = radiolocator.pulse(t,dim=1) + np.random.uniform(-0.3,0.3,size=t.size)


p, pcov = curve_fit(brown, xdata=t,ydata=pulse,p0=[0,alpha*0.01,0,0])
# p1, pcov = curve_fit(ice, xdata=t,ydata=pulse,p0=[1,(alpha*sigma)**2,0,alpha])
# print(p1)
# print(popt)



# plt.plot(t,brown(t,p[0],p[1],p[2]))
plt.plot(t,pulse,label='теория + белый шум')
t = np.linspace(-5e-8, 2e-7, 512)
pulse0 = brown(t,p[0],p[1],p[2],p[3])
sigma0 = max(pulse0) - min(pulse0)
print(sigma0)
plt.plot(t,pulse0,label='аппроксимация')
# plt.plot(t,ice(t,p1[0],p1[1],p1[2],p1[3]))
plt.plot(t,brown(t,1,alpha,sigma,0),label='теория')
plt.legend()
plt.show()



# M = 256


# N = 256
# x0 = np.linspace(0,200,50)
# y0 = x0
# t = 0


# fig,ax = plt.subplots(nrows = 1, ncols = 1)
# surface = water.Surface(N=N,M=M,U10=5,wind= np.pi/6)
# x, y = np.meshgrid(x0, y0)
# z = surface.model([x,y],t)
# print(z.shape)
# from matplotlib.cm import winter
# plt.contourf(x,y,z,levels=100,cmap=winter)
# plt.colorbar()

# plt.ylabel('Y, м',fontsize=16)
# plt.xlabel('X, м',fontsize=16)
# plt.show()