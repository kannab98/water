import multiprocessing as mp
import matplotlib.pyplot as plt
import water
import numpy as np
import time
import os

def mpCalc(nCPU = mp.cpu_count() ):
    x = np.linspace(0,100,10)
    y = np.linspace(0,100,10)
    x,y = np.meshgrid(x,y)
    x1 = x[0:50]
    x2 = x[50:-1]
    jobs = [[x1,x1], [x2,x2] ]
    p = mp.Pool()

    t = 0
    self = water.Surface()
    N = self.N
    M = self.M
    self.k = self.k[:N]
    k = self.k
    phi = self.phi
    A = self.A
    F = self.F
    psi = self.psi

    p.map(model,jobs)
    p.close()

def model(r):
    self.surface = 0
    self.amplitudes = np.array([ A[i]*sum(F[i])  for i in range(N)])
    for n in range(N):
        for m in range(M):
            self.surface += A[n] * \
            np.cos(
                +k[n]*(r[0]*np.cos(phi[m])+r[1]*np.sin(phi[m]))
                +psi[n][m]
                +self.omega_k(k[n])*t) \
                * F[n][m]
    print('Done')
    return self.surface

mpCalc()
