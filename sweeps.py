from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import scqubits as scq
from matplotlib import pyplot, animation
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import cmath
import parameters

# --Set Parameters--
Ej = parameters.Ej
Ec = parameters.Ec
ng = parameters.ng
frequency = parameters.frequency
anharmonicity = parameters.anharmonicity
g = parameters.g
# --/Set Parameters--

# --Take User Inputs
mode_p = int(input("Vary the {1} amplitude or the {2} length?"))
range = int(input("Enter the range of values to sweep over: "))
params = np.arange(0,range,range/500.0)
if (mode_p == 1):
    width = float(input("Enter the pulse width (this will be constant): "))
if (mode_p == 2):
    A = float(input("Enter the pulse amplitude (this will be constant): "))

# Create the qubit
qubit = scq.Transmon(EJ=Ej, EC=Ec, ng=ng, ncut=150) # ncut seems to be just a number that needs to be big
N = 2
wc = frequency*2*np.pi
wq = anharmonicity*2*np.pi
a = tensor(qeye(2),destroy(N))
sm = tensor(destroy(2),qeye(N))
sz = tensor(sigmaz(),qeye(N))
down = tensor(basis(N,1),basis(N,0))
up = tensor(basis(N,0),basis(N,1))

# Defining the drive pulse
def pulse(t, args):
    return args['A']*(t > args['delay']) - args['A']*(t>(args['width']+args['delay']))

# Hamiltonians
H0 = -0.5*wq*sz # qubit Hamiltonian
Hc = wc*(a.dag()*a) # cavity Hamiltonian
Hi = g*(a.dag()*sm + a*sm.dag()) # interaction Hamiltonian
Hdrive = sm+sm.dag() # drive Hamiltonian
H = [H0,[Hdrive,pulse]]

results = np.zeros(len(params))
c_ops = []
e_ops = [down*down.dag()]

for i in np.arange(0,len(params)):
    if (mode_p == 1):
        tlist = np.linspace(0, width + 5, 50)
        result = mesolve(H, down, tlist, c_ops, e_ops, args={'A': params[i], 'width': width, 'delay': 0}, options=Options(nsteps=5000))
    if (mode_p == 2):
        tlist = np.linspace(0, params[i] + 5, 50)
        result = mesolve(H, down, tlist, c_ops, e_ops, args={'A': A, 'width': params[i], 'delay': 0}, options=Options(nsteps=5000))
    results[i] = result.expect[0][len(tlist)-1]
    print(result.expect[0][len(tlist)-1])

fig, ax = plt.subplots(figsize=(12,6))
plt.rcParams.update({'font.size': 22})
ax.plot(params, results, 'b')
ax.set_xlabel('Pulse Parameter',fontsize=24)
ax.set_ylabel('Probability',fontsize=24)
plt.ylim([-.1,1.1])
plt.show()

fig = pyplot.figure()
ax = Axes3D(fig, azim=-40, elev=30)
sphere = qutip.Bloch(axes=ax)

# Convert expectation values to spherical coordinates
theta = [i * np.pi for i in results]
phi = [0 for i in results]

def animate(i):
    sphere.clear()
    sphere.add_vectors([np.sin(theta[i]) * np.cos(phi[i]), np.sin(theta[i]) * np.sin(phi[i]), np.cos(theta[i])])
    sphere.make_sphere()
    return ax

def init():
    sphere.vector_color = ['r']
    return ax

ani = animation.FuncAnimation(fig, animate, np.arange(len(results)), init_func=init, blit=False, repeat=False)
ani.save("rabi.gif", fps=50)