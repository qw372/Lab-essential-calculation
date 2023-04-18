import numpy as np
import matplotlib.pyplot as plt

sigma = 20 # um
U1 = -107.4
U2 = -309
x = np.linspace(-100, 100, 1000)
y1 = U1*np.exp(-x**2/2/sigma**2)
y2 = U2*np.exp(-x**2/2/sigma**2)

plt.plot(x, y1, linewidth=5)
plt.plot(x, y2, linewidth=5)

T = 8
limit = np.sqrt(T/np.abs(U1))*sigma
x = np.linspace(-limit*4, limit*4, 200)
y1 = np.abs(U1)/2*np.exp(-(np.abs(U1)/T)*x**2/2/sigma**2)-np.abs(U1)
plt.fill_between(x, y1, U1, color='green', alpha=0.45, linewidth=0)
limit = np.sqrt(T/np.abs(U2))*sigma
x = np.linspace(-limit*4, limit*4, 200)
y2 = np.abs(U1)/2*np.exp(-(np.abs(U2)/T)*x**2/2/sigma**2)-np.abs(U2)
plt.fill_between(x, y2, U2, color='green', alpha=0.45, linewidth=0)


plt.show()