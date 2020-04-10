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

df = read_csv('imp_15a.dat', sep ='\s+', header = None).T
# df = read_csv('test.csv', sep =',', header = None).T

t0 = df.iloc[0].values*1e-9
pulse0 = df.iloc[1].values
pulse = pulse0
t = t0

Nmax = np.argmax(pulse)
Pulse = np.log(pulse[Nmax:])
T = t[Nmax:]
line = lambda t,alpha,b: -alpha*t + b   

popt = curve_fit(line, 
                    xdata=T,
                    ydata=Pulse,
                    p0=[2e6,0],
                )[0]

plt.plot(T,Pulse)
plt.plot(T,line(T,popt[0],popt[1]))
alpha = popt[0] # 2096503.7683083578
b = popt[1]

c = 299792458
R = 6370e3
theta = np.deg2rad(1.5)
gamma = 2*np.sin(theta/2)**2/np.log(2)


print(alpha) # 2095703 Почти идеально
ice = lambda t, A, tau,sigma_l:   A * (1 + erf( (t-tau)/sigma_l ))
popt = curve_fit(ice, 
                    xdata=t[0:Nmax],
                    ydata=pulse[0:Nmax]*2,
                    p0=[(max(pulse) - min(pulse))/2,5e-11,6e-9],
                )[0]
plt.figure()
A = popt[0]
tau = popt[1]
sigma_l = popt[2]
# plt.plot(t,ice(t,popt[0],popt[1],popt[2]))
# plt.plot(t,pulse*2)
ice0 = A * np.exp(-alpha*(t-tau/2))* (1 + erf( (t-tau)/sigma_l ))
# plt.plot(t, ice(t))
# plt.show()

N0max = np.argmax(ice0)

# # # #############################################################
# arg = np.max(np.argwhere(pulse < 0.01*np.max(pulse)))
# left_edge = arg
# for i in range(2,arg):
#     line = lambda t,T: T   
#     popt = curve_fit(line, 
#                         xdata=t[0:i],
#                         ydata=pulse[0:i],
#                         p0=[0],
#                     )[0]
#     if np.linalg.norm(pulse[0:i] - line(t[0:i],popt[0])) > 0.01:
#         left_edge = i
#         break

# sigma_l = (t[Nmax] - t[left_edge])/3.5
# print(sigma_l)




# # ice = lambda t, A, tau,sigma_l:   A * np.exp(-alpha*(t-tau/2))* (1 + erf( (t-tau)/sigma_l ))
# # popt = curve_fit(ice, 
# #                     xdata=t[0:Nmax],
# #                     ydata=pulse[0:Nmax]/max(pulse)*2,
# #                     p0=[1,0,1e-8],
# #                 )[0]
# # plt.figure()
# # print(popt)
# # plt.plot(t,ice(t,popt[0],popt[1],popt[2]))
# # plt.plot(t[:Nmax],pulse[:Nmax]/max(pulse)*2)
# # plt.show()

# # ########################## оптимизация по амплитуде
# A = 1
# sigma_c = sigma_l / np.sqrt(2)
# tau = alpha * sigma_c**2
# print('tau= ', tau)
# print('sigma_l= ', sigma_l)

brown = lambda t,A,alpha,sigma:  A*np.exp(-alpha*(t - alpha/2 * sigma**2))*(1 + erf( (t - alpha*sigma**2)/(2**0.5*sigma))) 
ice = lambda t, A,alpha,tau,sigma_l:  A * np.exp(-alpha*(t-tau/2)) * (1 + erf((t-tau)/sigma_l))


t_true = t - t[Nmax] + t[N0max]
popt = curve_fit(ice, 
                    xdata=t_true,
                    ydata=pulse,
                    p0=[A,alpha,tau,sigma_l],
                )[0]





plt.figure(3)
plt.plot(t,pulse)
# y = brown(t, 1, alpha,sigma_c)
# plt.plot(t,y)
y = ice(t_true,popt[0],popt[1],popt[2],popt[3])
plt.plot(t, y)
plt.show()



