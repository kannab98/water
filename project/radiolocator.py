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
import pandas as pd
# import pandas.DataFrame


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
    def __init__(self, h=1e6, xi=0.0, theta=1.5, c=299792458, sigma=1, 
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
    
    def A(self,gamma,xi,A0=1.):
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

    def pulse_v(self, v, dim = 1):

        self.dim = dim
        gamma = self.gamma(self.theta)

        delta = self.delta(gamma,self.xi,self.h)
        beta  = self.beta(gamma,self.xi,self.h)

        alpha = self.alpha(beta,delta)

        sigma_c = self.sigma_c(self.sigma_s)

        u = np.sqrt(2)*v/(alpha*sigma_c) - alpha*sigma_c/np.sqrt(2)

        A = self.A(gamma,self.xi)
        pulse = A*np.exp(-v)*( 1 + erf(u) )

        return pulse

    def calc(self,t,pulse):
        self.brown = lambda t,A,alpha,sigma:  A*np.exp(-alpha*(t - alpha/2 * sigma**2))*(1 + erf( (t - alpha*sigma**2)/(2**0.5*sigma))) 
        
        theta = np.deg2rad(1.5)
        gamma = 2*np.sin(theta/2)**2/np.log(2)

        popt = curve_fit(self.brown, 
                            xdata=t,
                            ydata=pulse,
                            p0=[1,2e6,0],
                            bounds=( (0,1e5,0), (2,4e6,np.inf) )
                        ) [0]
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
        P = max(pulse) - min(pulse)
        return np.rad2deg(xi), h, sigma_s, 10*np.log10(P)



radiolocator = Radiolocator(xi = 0, sigma = 1)
t = np.linspace(-2e-7, 2e-7, 256)
pulse = radiolocator.pulse(t)     
data = pd.DataFrame({"col1":t,"col2":pulse})
data.to_csv(r'test.csv',index=False,header=None)
plt.plot(t,pulse)


# df = read_csv('imp04_10.dat', sep ='\s+', header = None).T
# t0 = df.iloc[0].values*1e-9
# pulse0 = df.iloc[1].values
# pulse0 = 2*pulse0/max(pulse0)
# # pulse = pulse0
# # t = t0



# T  =  t[np.argmax(pulse)]
# # T0 = t0[np.argmax(pulse0)]
# # print(T,T0)

# # plt.plot(t,pulse,label='теория')
# # t00 = t0 -T0 + T


# # params = radiolocator.calc(t=t00,pulse=pulse0)


# # radiolocator = Radiolocator(xi = params[0], h = params[1], sigma = params[2])
# # pulse = radiolocator.pulse(t)     








# # plt.plot(t,pulse,label='аппроксимация')
# # plt.plot(t,pulse,label='аппроксимация')
# # plt.plot(t0 - T0 + T , pulse0, label='практика')
# # plt.plot(t0 - T0 + T , pulse0, label='практика')

# # print('xi = {},\nh = {},\nsigma_s = {}'.format(params[0],params[1],params[2]))
# # plt.legend()

# ############3 Ищу коэффициент наклона заднегот фронта
# Nmax = np.argmax(pulse)
# Pulse = np.log(pulse[Nmax:])
# T = t[Nmax:]
# line = lambda t,alpha,b: -alpha*t + b   

# popt = curve_fit(line, 
#                     xdata=T,
#                     ydata=Pulse,
#                     p0=[2e6,0],
#                 )[0]


# print(popt)

# theta = np.deg2rad(1.5)
# gamma = 2*np.sin(theta/2)**2/np.log(2)
# plt.figure(1)
# plt.plot(T,Pulse)
# plt.plot(T, line(T,popt[0],popt[1]))
# print('наклон заднего фронта= ', popt[0])
# alpha = popt[0]
# #### Оценка амплитуды
# A = ( line(T[0],popt[0],popt[1]) - line(T[-1],popt[0],popt[1]) )/2 * np.exp(+popt[0] *(T[-1]-T[0]))
# print('оценка амплитуды=', A)
# ###### Оценка xi
# # квадратное уравнение
# h = radiolocator.H(radiolocator.h)
# # h = 1e6
# a = 1
# b = - gamma
# c = alpha * radiolocator.Gamma**2 * h/4/radiolocator.c - 1
# xi =  (-b + np.sqrt(b**2 - 4*a*c))/2
# xi = np.abs(xi)
# print(xi)
# print('оценка xi=', np.rad2deg(np.arccos(xi)/2))
# plt.figure(2)
# ############# Время искать ширину переднего фронта.

# args = (np.argwhere(pulse < 0.01*np.max(pulse)))
# args = np.max(args)
# norm = []
# left_edge = args
# for i in range(2,args):
#     line = lambda t,T: T   
#     popt = curve_fit(line, 
#                         xdata=t[0:i],
#                         ydata=pulse[0:i],
#                         p0=[0],
#                     )[0]
#     # if abs(pulse[i] - line(t[i],popt[0],popt[1])) < 0.05:
#     #     count = i
#     if np.linalg.norm(pulse[0:i] - line(t[0:i],popt[0])) > 0.01:
#         left_edge = i
#         break
#     # print(popt[0]*1e-9)

# # count = np.argmin(norm)
# count = left_edge

# print(count)
# sigmal = t[Nmax] - t[left_edge]
# print('ширина переднего фронта= ', sigmal)
# print('оценка амплитуды', (max(pulse) - min(pulse))/2 ) 
# print('эпоха=', alpha/2*sigmal**2)

# rad = radiolocator
# h =  
# alpha_new =  4/gamma * rad.c/h * (cos )
# brown = lambda t,,alpha,sigma:  A*np.exp(-alpha*(t - alpha/2 * sigma**2))*(1 + erf( (t - alpha*sigma**2)/(2**0.5*sigma))) 


# theta = np.deg2rad(1.5)
# gamma = 2*np.sin(theta/2)**2/np.log(2)

# popt = curve_fit(brown, 
#                     xdata=t,
#                     ydata=pulse,
#                     p0=[1,2e6,0],
#                     bounds=( (0,1e5,0), (2,4e6,np.inf) )
#                 ) [0]

# plt.plot(t[0:count+1],pulse[0:count+1],label='нули')
# plt.plot(t[count:Nmax+1],pulse[count:Nmax+1],label='передний фронт')
# plt.plot(t[Nmax:],pulse[Nmax:],label='задний фронт')
# plt.legend()
# plt.show()