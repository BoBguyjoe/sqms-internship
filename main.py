from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import scqubits as scq
from matplotlib import pyplot, animation
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import cmath
import parameters

scq.settings.T1_DEFAULT_WARNING=False
scq.set_units("GHz") # All units in GHZ (or 1/GHz = ns)

# --Set Parameters--
Ej = parameters.Ej
Ec = parameters.Ec
ng = parameters.ng
frequency = parameters.frequency
anharmonicity = parameters.anharmonicity
g = parameters.g
A = 0 # declaring pulse variables to prevent from running into a "variable not defined" error
width = 0
delay = 0
# --/Set Parameters--

# --Take User's Specifications--
mode_dr = int(input("Drive the qubit? {0} for no, {1} for yes"))
if (mode_dr != 1 and mode_dr != 0):
    print("Yo, get your act together. " + str(mode_dr) + " wasn't an option.")
    exit()
elif (mode_dr == 1):
    print("The qubit will be driven with a square pulse.")
    A = float(input("Amplitude: "))
    width = float(input("Width: "))
    delay = float(input("Delay: "))
mode_di = int(input("Include Dissipation? {0} for no, {1} for yes"))
if (mode_di != 1 and mode_di != 0):
    print("Yo, get your act together. " + str(mode_di) + " wasn't an option.")
    exit()
#mode_de = int(input("Include Dephasing? {0} for no, {1} for yes"))
#if (mode_de != 1 and mode_de != 0):
#    print("Yo, get your act together. " + str(mode_de) + " wasn't an option.")
#    exit()

print()
print("Decoherence Times")
print("{1} Calculate from parameters")
print("{2} Input manually")
manualInput = int(input()) - 1
if (manualInput != 0 and manualInput != 1):
    print("Yo, get your act together. " + str(manualInput + 1) + " wasn't an option.")
    exit()
if (manualInput == 0):
    print()
#    if (mode_de == 1):
#        mode_ng = int(input("Source dephasing noise from {1} Critical Current, or {2} Offset Charge")) - 1
#        if (mode_ng != 0 and mode_ng != 1):
#            print("Yo, get your act together. " + str(mode_ng) + " wasn't an option.")
#            exit()
#        else:
#            mode_cc = 1 - mode_ng

print()
print("Initialize the qubit at:")
print("1: |1]")
print("2: |0]")
print("3: 1/sqrt(2)(|0] + |1])")
mode_in = int(input())
if (mode_in != 1 and mode_in != 2 and mode_in != 3 and mode_in != 4):
    print("Yo, get your act together. " + str(mode_in) + " wasn't an option.")
    exit()

# Specify Output(s)
print()
printExpect = int(input("Print expectation values to console? {0} for no, {1} for yes"))
if (printExpect != 1 and printExpect != 0):
    print("Yo, get your act together. " + str(printExpect) + " wasn't an option.")
    exit()
makePlot = int(input("Create probability plot? {0} for no, {1} for yes"))
if (makePlot != 1 and makePlot != 0):
    print("Yo, get your act together. " + str(makePlot) + " wasn't an option.")
    exit()
elif (makePlot == 1):
    print("The plot will save to a .png file. What would you like to name it?")
    plotName = str(input())
makeSphere = int(input("Create Bloch sphere animation? {0} for no, {1} for yes"))
if (makeSphere != 1 and makeSphere != 0):
    print("Yo, get your act together. " + str(makeSphere) + " wasn't an option.")
    exit()
elif (makeSphere == 1):
    print("The sphere animation will save to a .gif file. What would you like to name it?")
    sphereName = str(input())
# --/Take User's Specifications--

# Create the Qubit
qubit = scq.Transmon(EJ=Ej, EC=Ec, ng=ng, ncut=150) # ncut seems to be just a number that needs to be big
N = 2
wc = frequency*2*np.pi # converting inputted frequencies to angular frequency
wq = anharmonicity*2*np.pi
g = g*2*np.pi
a = tensor(qeye(2),destroy(N))
sm = tensor(destroy(2),qeye(N))
sx = tensor(sigmax(),qeye(N))
sy = tensor(sigmay(),qeye(N))
sz = tensor(sigmaz(),qeye(N))
ground = tensor(basis(N,1),basis(N,0))
excite = tensor(basis(N,0),basis(N,1))

H0 = -0.5*wq*sz # qubit Hamiltonian
Hc = wc*(a.dag()*a) # cavity Hamiltonian
Hi = g*(a.dag()*sm + a*sm.dag()) # interaction Hamiltonian
Hdrive = a+a.dag() # drive Hamiltonian

if (mode_in == 1):
    psi0 = excite
elif (mode_in == 2):
    psi0 = ground
elif (mode_in == 3):
    psi0 = (excite+ground).unit()

# Defining the drive pulse
def pulse(t, args):
    return args['A']*(t > args['delay']) - args['A']*(t>(args['width']+args['delay']))

# Calculating Decoherence Times
if (manualInput == 0):
    if (mode_di == 1):
        T1 = qubit.t1_effective() / (2*np.pi)
#    if (mode_de == 1):
#        Tphi = (mode_ng*qubit.tphi_1_over_f_ng() + mode_cc*qubit.tphi_1_over_f_cc()) / (2*np.pi)
elif (manualInput == 1):
    print()
    if (mode_di == 1):
        T1 = float(input("Enter T1 value: "))
        if (T1 < 0):
            print("Yo, get your act together. This can't be a negative number")
            exit()
#    if (mode_de == 1):
#        Tphi = float(input("Enter Tphi value: "))
#        if (Tphi < 0):
#            print("Yo, get your act together. This can't be a negative number")
#            exit()

print()
if (mode_di == 1):
    print("T1 = " + str(T1) + " microseconds")
#if (mode_de == 1):
#    print("Tphi = " + str(Tphi) + " microseconds")

print()
range = float(input("Enter which t value this should calculate to: "))
if (range <= 0):
    print("Yo, get your act together. This has to be a positive number")
    exit()

# --Doing the thing--
c_ops = []
tlist = np.linspace(0,range,200)
# Add dissipation to collapse operators
if (mode_di == 1):
    kappa_di = np.power(T1,-1)
    c_ops.append(np.sqrt(kappa_di)*a)

# Add dephasing to collapse operators
#if (mode_de == 1):
#    kappa_de = np.power(Tphi/2.0,-1)
#    c_ops.append(np.sqrt(kappa_de)*a*a.dag())

# Set expectation value output to: [0] excited state, [1] ground state, and [2] phase
e_ops = [excite*excite.dag(),ground*ground.dag(),(ground+excite).unit()*(ground+excite).unit().dag()]

# Set the hamiltonian
if (mode_dr == 1):
    H = [H0,[Hdrive,pulse]] # This adds on the driving term. The driving term has a coefficient that is the pulse
else:
    H = H0

result = mesolve(H, psi0, tlist, c_ops, e_ops, args={'A': A, 'width': width, 'delay': delay}, options = Options(nsteps=5000))
# --/Doing the thing--

# --Print expectation values--
if (printExpect == 1):
    for i in np.arange(0,len(result.expect)):
        print(result.expect[i])
# --/Print expectation values--

# --Create Bloch sphere animation--
if (makeSphere == 1):
    fig = pyplot.figure()
    ax = Axes3D(fig, azim=-40, elev=30)
    sphere = qutip.Bloch(axes=ax)

    # Convert expectation values to spherical coordinates
    if (mode_dr == 0):
        theta = [i * np.pi for i in result.expect[0]]
        phi = [i * np.pi for i in result.expect[2]]
    if (mode_dr == 1 and mode_di == 0):
        theta = [(1-i) * np.pi for i in result.expect[1]]
        phi = [(1-i) * np.pi for i in result.expect[2]]
    if (mode_dr == 1 and mode_di == 1):
        theta = [i * np.pi for i in result.expect[0]]
        phi = [(1-i) * np.pi for i in result.expect[2]]

    def animate(i):
        sphere.clear()
        sphere.add_vectors([np.sin(theta[i]) * np.cos(phi[i]), np.sin(theta[i]) * np.sin(phi[i]), np.cos(theta[i])])
        sphere.make_sphere()
        return ax

    def init():
        sphere.vector_color = ['r']
        return ax

    ani = animation.FuncAnimation(fig, animate, np.arange(len(result.expect[1])),
                                  init_func=init, blit=False, repeat=False)
    ani.save(sphereName + ".gif", fps=50)
    plt.close()
    plt.close()
    print("Animation saved")
# --/Create Bloch sphere animation--

# --Create probability plot--
if (makePlot == 1):
    plt.close()
    fig, ax = plt.subplots(figsize=(12, 6))
    plt.rcParams.update({'font.size': 22})
    if (mode_dr == 0):
        ax.plot(tlist, np.real(result.expect[0]), 'r')
        ax.plot(tlist, np.real(result.expect[2]), 'm--')
    if (mode_dr == 1 and mode_di == 0):
        ax.plot(tlist, np.real(1-result.expect[1]), 'r')
        ax.plot(tlist, np.real(1-result.expect[2]), 'm--')
    if (mode_dr == 1 and mode_di == 1):
        ax.plot(tlist, np.real(result.expect[0]), 'r')
        ax.plot(tlist, np.real(result.expect[2]), 'm--')
    ax.legend(("|1]", "|0]+|1]"))
    ax.set_xlabel('Time', fontsize=20)
    ax.set_ylabel('Probability', fontsize=20)
    plt.ylim([-.1, 1.1])
    plt.savefig(plotName + ".png")

    # Plot pulse (if applicable)
    if (mode_dr == 1):
        fig, ax = plt.subplots(1, 1, figsize=(7, 5))
        plt.rcParams.update({'font.size': 22})
        ax.set_xlabel('Time', fontsize=24)
        ax.set_ylabel('Amplitude', fontsize=24)
        plt.title("Qubit Drive", fontsize=34)
        ax.plot(tlist, [pulse(t, args={'A': A, 'width': width, 'delay': delay}) for t in tlist], 'r')
    plt.show()
# --/Create probability plot--