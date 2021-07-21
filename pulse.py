from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import scqubits as scq
from matplotlib import pyplot, animation
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import cmath
import parameters
import scipy as sci

A = 9
s = 0.05

def amps(t):
    return A*np.exp(-((t-0.5)/s)**2)

tlist = np.linspace(0,1,1000)
energy = sci.integrate.quad(amps,0,2)
print(energy[0])

fig, ax = plt.subplots(figsize=(12,6))
plt.rcParams.update({'font.size': 22})
ax.plot(tlist, amps(tlist), 'b')
ax.set_xlabel('Time',fontsize=24)
ax.set_ylabel('Amplitude',fontsize=24)
plt.title("Qubit Drive",fontsize=34)

plt.show()