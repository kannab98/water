
import numpy as np
from numpy import pi
from scipy import interpolate,integrate
from tqdm import tqdm
# from water.spectrum import Spectrum
from spectrum import Spectrum

class Surface(Spectrum):
    def __init__(self,  N=256, M=100, space='log',
        random_phases = 1, 
        whitening = None, wind = 0,cwm = None, **kwargs):
        Spectrum.__init__(self,**kwargs)
        self.get_spectrum()
        self.N = N
        self.M = M
        KT = self.KT
        self.wind = wind # Направление ветра
        if space=='log':
            self.k = np.logspace(np.log10(KT[0]), np.log10(KT[-1]),self.N + 1)
        else:
            self.k = np.linspace(KT[0], KT[-1],self.N + 1)


        if whitening != None:
            if 'h' in whitening:
                interspace = self.interspace(self.k, N, power=0)
                self.k_heights = self.nodes(interspace,power=0)
                self.k = self.k_heights

            if 's' in whitening:
                interspace = self.interspace(self.k, N, power=2)
                self.k_slopes = self.nodes(interspace,power=2)
                self.k = self.k_slopes

            if 'h' in whitening and 's' in whitening:
                self.k = np.hstack((self.k_heights,self.k_slopes))
                self.k = np.sort(self.k)

        self.phi = np.linspace(0,2*np.pi,self.M + 1)
        self.phi_c = self.phi + self.wind

        # случайные фазы
        if random_phases == 0:
            self.psi = np.array([
                    [0 for m in range(self.M) ] for n in range(self.N) ])
        elif random_phases == 1:
            self.psi = np.array([
                [ np.random.uniform(0,2*pi) for m in range(self.M)]
                            for n in range(self.N) ])


        # массив с амплитудами i-ой гармоники
        # self.A = self.amplitude1(self.k)
        # угловое распределение
        # self.F = self.angle1(self.k,self.phi)

    def B(self,k):
          def b(k):
              b=(
                  -0.28+0.65*np.exp(-0.75*np.log(k/self.k_m))
                  +0.01*np.exp(-0.2+0.7*np.log10(k/self.k_m))
                )
              return b
          B=10**b(k)
          return B

    def Phi(self,k,phi):
        # Функция углового распределения
        phi = phi -self.wind
        normalization = lambda B: B/np.arctan(np.sinh(2* (pi)*B))
        B0 = self.B(k)
        A0 = normalization(B0)
        Phi = A0/np.cosh(2*B0*(phi) )
        return Phi


    def angle(self,k,phi,method='h'):
        M = self.M
        # N = self.N
        N = self.N
        # print(k.size)
        if method =='h':
            Phi = lambda phi,k: self.Phi(k,phi)
        elif method == 'xx':
            Phi = lambda phi,k: self.Phi(k,phi)*np.cos(phi)**2
        elif method == 'yy':
            Phi = lambda phi,k: self.Phi(k,phi)*np.sin(phi)**2
        else:
            Phi = lambda phi,k: self.Phi(k,phi)

        integral = np.zeros((N,M))
        # print(integral.shape)
        for i in range(N):
            for j in range(M):
                # integral[i][j] = integrate.quad( Phi, phi[j], phi[j+1], args=(k[i],) )[0]
                integral[i][j] = np.trapz( Phi( phi[j:j+2],k[i] ), phi[j:j+2])
        amplitude = np.sqrt(2 *integral )
        return amplitude

    def amplitude(self, k,method='h'):
        N = len(k)
        if method == 'h':
            S = self.spectrum
        else:
            S = lambda k: self.spectrum(k) * k**2

        integral = np.zeros(k.size-1)
        # progress_bar = tqdm( total = N-1)
        for i in range(1,N):
            integral[i-1] = integrate.quad(S,k[i-1],k[i])[0] 
        #     progress_bar.update(1)
        # progress_bar.clear()
        # progress_bar.close()
        # integral = np.array([ integrate.quad(S,k[i-1],k[i])[0] for i in range(1,N) ])
        amplitude = np.sqrt(2 *integral )
        return np.array(amplitude)


    
    def calc_amplitude(self,k,phi):
        N = k.size
        M = phi.size 
        # print(N,M)
        integral = np.zeros((N-1,M-1))
        progress_bar = tqdm(total=(N-1)*M)
        spectrum_k = lambda k: self.spectrum(k)
        spectrum_phi = lambda k,phi: self.Phi(k,phi)
        func = lambda k,phi: spectrum_k(k)*spectrum_phi(k,phi)
        for i in range(1,N):
            for j in range(1,M):
                integral[i-1][j-1] = integrate.dblquad(func,phi[j-1],phi[j],lambda x: k[i-1],lambda x: k[i] )[0]
            progress_bar.update(1)
        progress_bar.clear()
        progress_bar.close()
        return np.sqrt(2*integral)


    def model(self,r,t,method='h'):
        N = self.N
        M = self.M
        # self.k = self.k[:N]
        k = self.k
        phi = self.phi
        A = self.amplitude(k,method=method)
        F = self.angle(k,phi,method=method)
        psi = self.psi
        self.surface = 0
        # self.amplitudes = np.array([ A[i]*sum(F[i])  for i in range(N)])
        progress_bar = tqdm( total = N*M,  leave = False )
        for n in range(N):
            for m in range(M):
                self.surface += A[n] * \
                np.cos(
                    +k[n]*(r[0]*np.cos(phi[m])+r[1]*np.sin(phi[m]))
                    +psi[n][m]
                    +self.omega_k(k[n])*t) \
                    * F[n][m]
                progress_bar.update(1)
        progress_bar.close()
        progress_bar.clear()
        print()
        return self.surface

    # def slopesxx(self,r,t):
    #     N = self.N
    #     M = self.M
    #     spectrum = self.spectrum(self.k) * (k*cos(phi))**2
    #     angles_distrib = self.Phi(self.k,phi)


    #     self.k = self.k[:N]
    #     k = self.k
    #     phi = self.phi
    #     A = self.A
    #     F = self.F
    #     psi = self.psi
    #     self.surface = 0
    #     self.amplitudes = np.array([ A[i]*sum(F[i])  for i in range(N)])
    #     progress_bar = tqdm( total = N*M,  leave = False )
    #     for n in range(N):
    #         for m in range(M):
    #             self.surface += A[n] * \
    #             np.cos(
    #                 +k[n]*(r[0]*np.cos(phi[m])+r[1]*np.sin(phi[m]))
    #                 +psi[n][m]
    #                 +self.omega_k(k[n])*t) \
    #                 * F[n][m]
    #             progress_bar.update(1)
    #     progress_bar.close()
    #     progress_bar.clear()
    #     print()
    #     return self.surface

    def D(self,r,t):
        N = self.N
        M = self.M
        k = self.k
        phi = self.phi
        A = self.A
        F = self.F
        psi = self.psi
        Dx = 0
        Dy = 0
        for n in range(N):
            for m in range(M):
                Dx += A[n]  * np.cos(phi[m]) *\
                np.cos(
                        +k[n]*(r[0]*np.cos(phi[m])+r[1]*np.sin(phi[m]))
                        +psi[n][m] + np.pi/2+self.omega_k(k[n])*t)  \
                        * F[n][m]

                Dy += A[n]  * np.sin(phi[m]) *\
                np.cos(
                        +k[n]*(r[0]*np.cos(phi[m])+r[1]*np.sin(phi[m]))
                        +psi[n][m] + np.pi/2 + self.omega_k(k[n])*t)  \
                        * F[n][m]


        return np.array(Dx),np.array(Dy)

# Экспериментальный способ. Метод <<отбеливания>> спектра
    def interspace(self, k, N, power=0):
        #     Функция разбивает заданный k-интервал на N участков, в каждом
        # из которых интеграл по спектру принимает одно и тоже значение.
        full_spectrum=self.spectrum
        y = lambda k: full_spectrum(k)*k**power
        sigma = self.sigma_sqr
        b0 = sigma/(N)
        k_new = k[0]
        k_end = k[-1]
        k = np.zeros(N+1)
        k[0] = k_new
#        err = np.zeros(N)
        epsabs = 1.49e-8
#        epsrel = 1.49e-8
        progress_bar = tqdm(total = N-1 )
        for i in range(N-1):
            integral = 0
            n = 1
            m = 1
            while integral < b0:
                k_new*= 1 + 10**(-m)
                integral,error = integrate.quad(y, k[i], k_new)
                if integral> b0 and abs(integral-b0) > epsabs:
                    k_new*= 1 - 10**(-m)
                    m+= 1
                    integral, error = integrate.quad(y, k[i], k_new)
                n+= 1
            progress_bar.update(1)
            k[i+1] = k_new
        k = k[0:N+1]
        k[-1] = k_end
        # Возвращает k размерности N+1, т.о. имеем как раз N интервалов.
        # err -- ошибка вычисения интегралов
        progress_bar.close()
        return k

    def nodes(self,ki, power=0):
        # Функции interspace и nodes используются последовательно.
        # Первая возвращает интервалы, а вторая вычисляет в возвращенных
        # интервалах координаты узлов
        sigma = self.sigma_sqr
        b0 = sigma/(len(ki)+1)
        progress_bar = tqdm(total = len(ki)-1)
        full_spectrum=self.spectrum

        y=lambda k: k**(2+power)*full_spectrum(k)

        nodes=np.zeros(len(ki))
        A=np.sqrt(2*b0)

        for i in range(1,len(ki)):
            integral,error=integrate.quad(y,ki[i-1],ki[i])
            B=(np.sqrt(integral/A**2))
            nodes[i-1]=B
            progress_bar.update(1)
        nodes = nodes[:-1]
        if nodes[-1]>self.KT[-1]:
            nodes[-1]=self.KT[-1]

        return nodes
    
    def model1(self,r,t):
        N = self.N
        M = self.M
        print(N,M)
        self.k = self.k[:N]
        k = self.k
        phi = self.phi
        A = self.calc_amplitude(k,phi)
        psi = self.psi
        self.surface = 0
        progress_bar = tqdm( total = N*M,  leave = False )
        for n in range(N-1):
            for m in range(M-1):
                self.surface += A[n][m] * \
                np.cos(
                    +k[n]*(r[0]*np.cos(phi[m])+r[1]*np.sin(phi[m]))
                    +psi[n][m]
                    +self.omega_k(k[n])*t) 
                progress_bar.update(1)
        progress_bar.close()
        progress_bar.clear()
        print()
        return self.surface

import matplotlib.pyplot as plt

surface = Surface(N=2048,M=1,U10=10)
x0 = np.linspace(0,100,10000)
y0 = x0
t=0
# x, y = np.meshgrid(x0, y0)
from matplotlib.cm import winter
# fig,ax = plt.subplots(nrows = 1, ncols = 3)
methods=['h','xx']
# titles=['высоты', 'наклоны xx', 'наклоны yy']
# for i in range(len(methods)):
fig, ax1 = plt.subplots()
z = surface.model([x0,0],t,method='h')
# slopes=np.diff(z[0])
slopes = surface.model([x0,0],t,method='x')

ax2 = ax1.twinx()

# slopes = np.diff(z)/np.diff(x0)
# ax2.plot(x0[1:],np.rad2deg(np.arctan(slopes)))
# slopes = surface.model([x0,0],t,method='x')
ax2.plot(x0,np.rad2deg(np.arctan(slopes)))

ax2.set_ylabel('Наклоны, градусы')
ax1.plot(x0,z,'r')
ax1.set_ylabel('Высоты, м')
fig.tight_layout()
# plt.savefig(titles[i]+'.png',dpi=600)
plt.show()
# print("\n")
# A = surface.calc_amplitude(surface.k,surface.phi)
# print(A.shape)
# A0 = surface.amplitude(surface.k)*surface.angle(surface.k,surface.phi)
# print(A0.shape)
