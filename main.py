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
# --/Set Parameters--

# --Take User's Specifications--
mode_dr = int(input("Drive the qubit? {0} for no, {1} for yes"))
if (mode_dr != 1 and mode_dr != 0):
    print("Yo, get your act together. " + str(mode_dr) + " wasn't an option.")
    exit()
mode_di = int(input("Include Dissipation? {0} for no, {1} for yes"))
if (mode_di != 1 and mode_di != 0):
    print("Yo, get your act together. " + str(mode_di) + " wasn't an option.")
    exit()
mode_de = int(input("Include Dephasing? {0} for no, {1} for yes"))
if (mode_de != 1 and mode_de != 0):
    print("Yo, get your act together. " + str(mode_de) + " wasn't an option.")
    exit()

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
    if (mode_de == 1):
        mode_ng = int(input("Source dephasing noise from {1} Critical Current, or {2} Offset Charge")) - 1
        if (mode_ng != 0 and mode_ng != 1):
            print("Yo, get your act together. " + str(mode_ng) + " wasn't an option.")
            exit()
        else:
            mode_cc = 1 - mode_ng

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
makeSphere = int(input("Create Bloch sphere animation? {0} for no, {1} for yes"))
if (makeSphere != 1 and makeSphere != 0):
    print("Yo, get your act together. " + str(makeSphere) + " wasn't an option.")
    exit()
elif (makeSphere == 1):
    print("The sphere animation will save to a .gif file. What would you like to name it?")
    name = str(input())
# --/Take User's Specifications--

# Create the Qubit
qubit = scq.Transmon(EJ=Ej, EC=Ec, ng=ng, ncut=150) # ncut seems to be just a number that needs to be big
N = 2
wc = frequency*2*np.pi
wa = anharmonicity*2*np.pi
g = g*2*np.pi
a = tensor(qeye(2),destroy(N))
sm = tensor(destroy(2),qeye(N))
down = tensor(basis(N,1),basis(N,0))
up = tensor(basis(N,0),basis(N,1))

#H = wc*((a.dag()*a)+0.5) + 0.5*wa*sm.dag()*sm + g*(a.dag()+a)*(sm+sm.dag()) #JC Hamiltonian
#H = wc*a.dag()*a + wa*sm.dag()*sm + g*(a.dag()*sm + a*sm.dag()) # RWA Hamiltonian
#H = 0.5*wa*a.dag()*a.dag()*a*a # some other Hamiltonian I found in a qutip example code
H = wc*a.dag()*a + g*(a.dag()*sm + a*sm.dag())

if (mode_in == 1):
    psi0 = up
elif (mode_in == 2):
    psi0 = down
elif (mode_in == 3):
    psi0 = (up+down).unit()

# Calculating Decoherence Times
if (manualInput == 0):
    if (mode_di == 1):
        T1 = qubit.t1_effective() / (2*np.pi)
    if (mode_de == 1):
        Tphi = (mode_ng*qubit.tphi_1_over_f_ng() + mode_cc*qubit.tphi_1_over_f_cc()) / (2*np.pi)
elif (manualInput == 1):
    print()
    if (mode_di == 1):
        T1 = float(input("Enter T1 value: "))
        if (T1 < 0):
            print("Yo, get your act together. This can't be a negative number")
            exit()
    if (mode_de == 1):
        Tphi = float(input("Enter Tphi value: "))
        if (Tphi < 0):
            print("Yo, get your act together. This can't be a negative number")
            exit()

print()
if (mode_di == 1):
    print("T1 = " + str(T1) + " microseconds")
if (mode_de == 1):
    print("Tphi = " + str(Tphi) + " microseconds")

print()
range = float(input("Enter which t value this should calculate to: "))
if (range <= 0):
    print("Yo, get your act together. This has to be a positive number")
    exit()

# Doing the thing
c_ops = []
tlist = np.linspace(0,range,200)
if (mode_di == 1):
    kappa_di = np.power(T1,-1)
    c_ops.append(np.sqrt(kappa_di)*a)
if (mode_de == 1):
    kappa_de = np.power(Tphi/2.0,-1)
    c_ops.append(np.sqrt(kappa_de)*a*a.dag())
e_ops = [down*down.dag(),up*up.dag(),(down+up).unit()*(down+up).unit().dag(),(down-up).unit()*(down-up).unit().dag()]
if (mode_dr == 1):
    H = H + wa*sm.dag()*sm

result = mesolve(H, psi0, tlist, c_ops, e_ops, args={'A': 14.0,'s':0.05})

# Print expectation values
if (printExpect == 1):
    i = 0
    while (i < 4):
        print(result.expect[i])
        i=i+1

# Create Bloch sphere animation
if (makeSphere == 1):
    fig = pyplot.figure()
    ax = Axes3D(fig, azim=-40, elev=30)
    sphere = qutip.Bloch(axes=ax)

    theta = [i * np.pi for i in result.expect[1]]
    phi = [i * np.pi for i in result.expect[3]]

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
    ani.save(name + ".gif", fps=50)
    plt.close()
    plt.close()
    print("Animation saved")

# Create probability plot
if (makePlot == 1):
    plt.close()
    fig, ax = plt.subplots(figsize=(12, 6))
    plt.rcParams.update({'font.size': 22})
    ax.plot(tlist, np.real(result.expect[0]), 'b')
    ax.plot(tlist, np.real(result.expect[1]), 'r')
    ax.plot(tlist, np.real(result.expect[2]), 'm+')
    ax.plot(tlist, np.real(result.expect[3]), 'm--')
    ax.legend(("|0]", "|1]", "|0]+|1]", "|0]-|1]"))
    ax.set_xlabel('Time', fontsize=20)
    ax.set_ylabel('Probability', fontsize=20)
    plt.ylim([-.1, 1.1])
    plt.show()