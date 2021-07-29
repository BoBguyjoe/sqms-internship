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

Ej = parameters.Ej
Ec = parameters.Ec
ng = parameters.ng
frequency = parameters.frequency
anharmonicity = parameters.anharmonicity
g = parameters.g
A = 18
s = 0.05
kappa = np.power(1.0e1,-1)
options = Options(atol=1e-12, rtol=1e-12,max_step=1e-4)

N = 2
wc = frequency*2*np.pi
wq = anharmonicity*2*np.pi
g = g*2*np.pi
a = tensor(qeye(2),destroy(N))
sm = tensor(destroy(2),qeye(N))
sx = tensor(sigmax(),qeye(N))
sy = tensor(sigmay(),qeye(N))
sz = tensor(sigmaz(),qeye(N))
down = tensor(basis(N,1),basis(N,0))
up = tensor(basis(N,0),basis(N,1))

def pulse(t, args):
    return args['A'] * np.exp(-((t-.5) / args['s']) ** 2)

A = 20
s = 0.05
width = 8
delay = 0
wd = 1

def square(t, args):
    return args['A']*(t > args['delay']) - args['A']*(t>(args['width']+args['delay']))

H0 = -0.5*wq*sz # qubit Hamiltonian
Hc = wc*(a.dag()*a) # cavity Hamiltonian
Hi = g*(a.dag()*sm + a*sm.dag()) # interaction Hamiltonian
Hdrive1 = sm + sm.dag()

def Hdrive2(t, args):
    return square(t, args)*np.exp(-1j*wd*t)*(a.dag() + a)


H = [H0,[Hdrive1,square]]
tlist = np.linspace(0,50,500)
c_ops = []
#c_ops.append(np.sqrt(kappa)*a)
e_ops = [down*down.dag(),up*up.dag(),(down+up).unit()*(down+up).unit().dag(),(down-up).unit()*(down-up).unit().dag()]
result = mesolve(H,down,tlist,c_ops,e_ops,args={'A': A,'s': s, 'width': width, 'delay': delay})

fig, ax = plt.subplots(1, 1, figsize=(7,5))
plt.rcParams.update({'font.size': 22})
ax.plot(tlist, pulse(tlist,args={'A': A,'s': s}), 'b')
ax.set_xlabel('Time',fontsize=24)
ax.set_ylabel('Amplitude',fontsize=24)
plt.title("Qubit Drive",fontsize=34)
ax.plot(tlist, [square(t,args={'A': A,'s': s, 'width': width, 'delay': delay}) for t in tlist], 'r')

fig, ax = plt.subplots(figsize=(12,6))
plt.rcParams.update({'font.size': 22})
ax.plot(tlist, np.real(result.expect[0]), 'b')
ax.plot(tlist, np.real(result.expect[1]), 'r')
ax.plot(tlist, np.real(result.expect[2]), 'm+')
ax.plot(tlist, np.real(result.expect[3]), 'm--')
ax.legend(("|0]", "|1]", "|0]+|1]", "|0]-|1]"))
ax.set_xlabel('Time',fontsize=24)
ax.set_ylabel('Probability',fontsize=24)
plt.ylim([-.1,1.1])

plt.show()