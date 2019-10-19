import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import tqdm




df = pd.read_csv('signal.dat', sep="\s+", header = 0 ).T
t = np.array(df.iloc[0])
# heights = np.array(df.iloc[1])
# T = 1/2000
# N = 2*np.pi/T
t = np.linspace(0,1,100)
# heights = sum([np.sin(i*2*np.pi*t) for i in range(1,2)])
# slopes = np.array(df.iloc[2])
heights = np.sin(2*np.pi*t)
def fft(t, x, fm = None, ftype='real',short = True):

    N = x.size
    # freq = np.linspace(-N/(2*t[-1]), N/(2*t[-1]), 3*N+1)
    freq = np.linspace(-2,2,3*N+1)

    # freq = np.linspace(-1.5,+1.5,10000)
    X = np.zeros(freq.size, dtype='complex')
    progress_bar = tqdm.tqdm(total = freq.size)

    for j in range(freq.size):
        for i in range(t.size):
            X[j] += x[i]*np.exp(2j*np.pi*freq[j]*t[i])  
        progress_bar.update(1)
    progress_bar.close()

    # X = np.fft.fft(x,n=N)
    # freq = np.linspace(-1/(2*T),+1/(2*T),N)
    if ftype == 'real':
        N = freq.size
        freq = freq[N//2:]
        S = 2*np.absolute(X[N//2:])/t.size
    return freq,S

fig,ax = plt.subplots(nrows = 1, ncols = 2)
ax[0].plot(t,heights)
freq,S = fft(t,heights,ftype='real',fm=5,short=False)
ax[1].plot(freq,S)
# ax[1][0].plot(t,slopes)
# freq,S = fft(t,slopes,'real')
# ax[1][1].plot(freq,S)
plt.show()