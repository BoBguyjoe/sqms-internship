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
scq.set_units("MHz") # All units in MHZ (or 1/MHz = us)

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
print("4: Initialize yourself")
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
wc = frequency*2*np.pi
wa = anharmonicity*2*np.pi
g = g*2*np.pi
ac = qutip.create(2)
ad = qutip.destroy(2)
sz = qutip.sigmaz()
sp = qutip.sigmap()
sm = qutip.sigmam()
N = 2
down = basis(N,0)
up = basis(N,1)

cavity = scq.Oscillator(E_osc=wc, truncated_dim=2)
space = scq.HilbertSpace([qubit, cavity])
space.add_interaction(g_strength=g, op1=qubit.n_operator, op2=cavity.creation_operator, add_hc=True)

H = wc*ad*ac + 0.5*wa*sz + g*(ad*sm + ac*sp) # JC Hamiltonian
#H = 0.5*wa*ac*ac*ad*ad # some other Hamiltonian I found in a qutip example code

if (mode_in == 1):
    psi0 = up
elif (mode_in == 2):
    psi0 = down
elif (mode_in == 3):
    psi0 = (up+down).unit()
elif (mode_in == 4):
    psi0 = basis(2, 0)
    i = True
    while (i == True):
        sphere = qutip.Bloch()
        sphere.vector_color = ['r']
        print("The qubit's state is: " + str(psi0))
        print("Which operator would you like to apply?")
        print("1: Rotate around x axis")
        print("2: Rotate around y axis")
        print("3: Rotate around z axis")
        print("4: Reset to (|1], |0])")
        print("5: Finish and exit")
        action = int(input())
        if (action != 1 and action != 2 and action != 3 and action != 4 and action != 5):
            print("Yo, get your act together. " + str(action) + " wasn't an option.")
            exit()

        if (action == 1 or action == 2 or action == 3):
            rotate = float(input("How much rotation: {1} for half-pi, {2} for pi")) / 2.0
        if (action == 1):
            psi0 = rotate * sigmax() * psi0
            theta = np.pi * psi0[1][0][0].real
            phi = np.pi * psi0[1][0][0].imag
            sphere.add_vectors([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)])
            sphere.make_sphere()
            plt.show()
        if (action == 2):
            psi0 = rotate * sigmay() * psi0
            theta = np.pi * psi0[1][0][0].real
            phi = np.pi * psi0[1][0][0].imag
            sphere.add_vectors([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)])
            sphere.make_sphere()
            plt.show()
        if (action == 3):
            psi0 = rotate * sigmaz() * psi0
            theta = np.pi * psi0[1][0][0].real
            phi = np.pi * psi0[1][0][0].imag
            sphere.add_vectors([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)])
            sphere.make_sphere()
            plt.show()
        if (action == 4):
            psi0 = basis(2,0)
            theta = np.pi * psi0[1][0][0].real
            phi = np.pi * psi0[1][0][0].imag
            sphere.add_vectors([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)])
            sphere.make_sphere()
            plt.show()
        if (action == 5):
            i = False

# Calculating Decoherence Times
if (manualInput == 0):
    if (mode_di == 1):
        T1 = qubit.t1_effective() / (2*np.pi)
    if (mode_de == 0):
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
if (range <=0):
    print("Yo, get your act together. This has to be a positive number")
    exit()

# Doing the thing
c_ops = []
tlist = np.linspace(0,range,200)
if (mode_di == 1):
    kappa_di = np.power(T1,-1)
    c_ops.append(np.sqrt(kappa_di)*ad)
if (mode_de == 1):
    kappa_de = Tphi/2.0
    c_ops.append(np.sqrt(kappa_de)*ad*ac)
e_ops = [down*down.dag(),up*up.dag(),(down+up).unit()*(down+up).unit().dag(),(down-up).unit()*(down-up).unit().dag()]
if (mode_dr == 1):
    H = H + 1*(ad+ac)

result = mesolve(H, psi0, tlist, c_ops, e_ops, options=Options(nsteps=5000))

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