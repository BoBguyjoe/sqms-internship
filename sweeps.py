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

sample = 50
T1 = 10.0

# --Take User Inputs
mode_p = int(input("Vary the {1} amplitude, {2} length, or {3} both?"))
if (mode_p == 1 or mode_p == 2):
    range = int(input("Enter the range of values to sweep over: "))
    if (mode_p == 1):
        params = np.linspace(-range, range, sample)
        width = float(input("Enter the pulse width (this will be constant): "))
    if (mode_p == 2):
        params = np.linspace(0, range, sample)
        A = float(input("Enter the pulse amplitude (this will be constant): "))
    print()
    makeSphere = int(input("Make a Bloch sphere animation? {0} for no, {1} for yes"))
if (mode_p == 3):
    rangeA = int(input("Enter the range of amplitudes to sweep over: "))
    amplitudes = np.linspace(-rangeA,rangeA, sample)
    rangeW = int(input("Enter the range of lengths to sweep over: "))
    widths = np.linspace(0,rangeW, sample)
    makeSphere = 0 # Making sure that this is declared

# Create the qubit
qubit = scq.Transmon(EJ=Ej, EC=Ec, ng=ng, ncut=150) # ncut seems to be just a number that needs to be big
N = 2
wc = frequency*2*np.pi
wq = anharmonicity*2*np.pi
a = tensor(qeye(2),destroy(N))
sm = tensor(destroy(2),qeye(N))
sx = tensor(sigmax(),qeye(N))
sy = tensor(sigmay(),qeye(N))
sz = tensor(sigmaz(),qeye(N))
down = tensor(basis(N,1),basis(N,0))
up = tensor(basis(N,0),basis(N,1))

# Defining the drive pulses
def square(t, args):
    return args['A']*(t > args['delay']) - args['A']*(t>(args['width']+args['delay']))

def gauss(t, args):
    return args['A'] * np.exp(-((t-.5) / args['width']) ** 2)

# Hamiltonians
H0 = -0.5*wq*sz # qubit Hamiltonian
Hc = wc*(a.dag()*a) # cavity Hamiltonian
Hi = g*(a.dag()*sm + a*sm.dag()) # interaction Hamiltonian
Hdrive = a+a.dag() # drive Hamiltonian
H = [H0,[Hdrive,gauss]]

c_ops = []
kappa_di = np.power(T1,-1)
#c_ops.append(np.sqrt(kappa_di)*sm)
e_ops = [down*down.dag(), sx, sy, sz]

if (mode_p == 1 or mode_p == 2):
    results = np.zeros(len(params))
    for i in np.arange(0,len(params)):
        if (mode_p == 1):
            tlist = np.linspace(0, width + 5, 50)
            result = mesolve(H, down, tlist, c_ops, e_ops, args={'A': params[i], 'width': width, 'delay': 0}, options=Options(nsteps=5000))
        if (mode_p == 2):
            tlist = np.linspace(0, params[i] + 5, 50)
            result = mesolve(H, down, tlist, c_ops, e_ops, args={'A': A, 'width': params[i], 'delay': 0}, options=Options(nsteps=5000))
        results[i] = result.expect[0][len(tlist)-1]

    fig, ax = plt.subplots(figsize=(12, 6))
    plt.rcParams.update({'font.size': 22})
    ax.plot(params, results, 'b')
    ax.set_xlabel('Pulse Parameter', fontsize=20)
    ax.set_ylabel('Probability', fontsize=20)
    plt.ylim([-.1, 1.1])
    plt.show()

if (mode_p == 3):
    tlist = np.linspace(0, widths[sample - 1] + 5, 50)
    results = [[0 for x in np.arange(0,len(amplitudes))] for y in np.arange(0,len(widths))]
    for i in np.arange(0,len(amplitudes)):
        for j in np.arange(0,len(widths)):
            result = mesolve(H, down, tlist, c_ops, e_ops, args={'A': amplitudes[i], 'width': widths[j], 'delay': 0},options=Options(nsteps=5000))
            results[j][i] = result.expect[0][len(tlist) - 1]
            print(str(i) + ", " + str(j) + ": " + str(results[i][j]))

    # rotating the axes for the 2D plot
    Results = [[0 for x in np.arange(0,len(amplitudes))] for y in np.arange(0,len(widths))]
    for i in np.arange(0,len(results)):
        for j in np.arange(0,len(results[0])):
            Results[i][j] = 1-results[sample-1-i][j]

    fig, ax = plt.subplots(figsize=(12,6))
    img = ax.imshow(Results,cmap = 'plasma', interpolation='none')
    plt.title('Rabi Sweeps')
    ax.set_xlabel('Amplitude')
    ax.set_ylabel('Length')
    plt.xticks([sample/4, sample/2, sample*3/4],[-rangeA/2, 0, rangeA/2])
    plt.yticks([sample/4, sample/2, sample*3/4],[rangeW*3/4, rangeW/2, rangeW/4])
    fig.colorbar(img)
    plt.savefig("rabi sweeps.png")
    plt.show()

if (makeSphere == 1):
    fig = pyplot.figure()
    ax = Axes3D(fig, azim=-40, elev=30)
    sphere = qutip.Bloch(axes=ax)

    # Convert expectation values to spherical coordinates
    theta = [i * np.pi for i in results]
    phi = [0 * np.pi for i in results]

    def animate(i):
        sphere.clear()
        sphere.add_vectors([np.sin(theta[i]) * np.cos(phi[i]), np.sin(theta[i]) * np.sin(phi[i]), np.cos(theta[i])])
        sphere.make_sphere()
        return ax

    def init():
        sphere.vector_color = ['r']
        return ax

    ani = animation.FuncAnimation(fig, animate, np.arange(len(x)), init_func=init, blit=False, repeat=False)
    ani.save("rabi.gif", fps=10)