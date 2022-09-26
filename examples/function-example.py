import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt

#fs = 10e3
#N = 1e5

fs = 100
N = 1e3

amp = 2*np.sqrt(2)
#freq = 1270.0
freq = 100.0
noise_power = 0.001 * fs / 2
time = np.arange(N) / fs
x = amp*np.sin(2*np.pi*freq*time)
#x += np.random.normal(scale=np.sqrt(noise_power), size=time.shape)

plt.plot(x)
plt.show()

