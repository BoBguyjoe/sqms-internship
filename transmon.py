from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import scqubits as scq

scq.settings.T1_DEFAULT_WARNING=False
scq.get_units()
scq.set_units("MHz")
options=Options(nsteps=5000) # a (failed) attempt to get the transmon Hamiltonian to work with driving

# --Set Parameters--
Ej = 34.0e3 # Transmon Josephson energy
Ec = 264.0 # Transmon charging energy
ng = 0.3 # cutoff charge (something from scqubits, I don't know)
frequency = 6.02e3 # Resonance frequency (in MHz)
anharmonicity = 3e2 # Anharmonicity (in MHz)
c = 130.0 # coupling strength
# --/Set Parameters--

# Create the qubit
qubit = scq.Transmon(EJ=Ec, EC=Ec, ng=ng, ncut=150) # ncut seems to be just a number that needs to be big

# --Driving, Decoherence, and Noise Types--
print("DRIVING")
print("[0] for no driving")
print("[1] for position driving")
print("[2] for momentum driving")
mode_dr = int(input())

print("DECOHERENCE & Noise")
mode_di = int(input("Include Photon Dissipation? [0] for no, [1] for yes"))
if (mode_di == 1):
    mode_cap = int(input("Include Capacitive Noise? [0] for no, [1] for yes"))
    mode_ci = int(input("Include Charge Impedance Noise? [0] for no, [1] for yes"))
elif (mode_di == 0):
    mode_cap = 0
    mode_ci = 0

mode_de = int(input("Include Dephasing? [0] for no, [1] for yes"))
if (mode_de == 1):
    mode_ng = int(input("Dephasing Noise from [1] Critical Current, or [2] Offset Charge")) - 1
    mode_cc = 1 - mode_ng
elif (mode_de == 0):
    mode_cc = 0
    mode_ng = 0
# --/Driving and Noise Types--

# Calculate T1 and Tphi from selected noise types
if (mode_cap == 1 and mode_ci == 0):
    T = qubit.t1_effective(noise_channels=["t1_capacitive"])/(2*np.pi)
elif (mode_cap == 0 and mode_ci == 1):
    T = qubit.t1_effective(noise_channels=["t1_charge_impedance"])/(2*np.pi)
elif (mode_cap == 1 and mode_ci == 1):
    T = qubit.t1_effective(noise_channels=["t1_capacitive", "t1_charge_impedance"])/(2*np.pi)
elif (mode_di == 0):
    T = 1.0

if (mode_ng == 1):
    dephasing = qubit.tphi_1_over_f_ng(get_rate=False)/(2*np.pi)
elif (mode_cc == 1):
    dephasing = qubit.tphi_1_over_f_cc(get_rate=False)/(2*np.pi)
elif (mode_de == 0):
    dephasing = 1.0

print("T1: " + str(T) + " microseconds")
print("Tphi: " + str(dephasing) + " microseconds")

# x-axis range
scalingRange = True
if (scalingRange == False):
    range = 10
else:
    range = (T/5.0)/(dephasing/30.0)
sample = 1000

# Hamiltonians (use only one, the JC one doesn't work, the transmon one only works if you're not driving)
wc = 2*np.pi*frequency #Resonance Frequency
wa = 2*np.pi*anharmonicity #Transmon Anharmonicity
N = 2 #Levels in the system
g = basis(N,0)
e = basis(N,1)
A = destroy(N)
a = tensor(destroy(N), qeye(2))
sm = tensor(qeye(N), destroy(2))

#H = wc*a.dag()*a + wa*sm.dag()*sm/2 + c*(a.dag()*sm + a*sm.dag()) # JC Hamiltonian
#H = wc*A*A.dag()+0.5*wa*A*A*A.dag()*A.dag() #transmon Hamiltonian
H = -.5*wa*A.dag()*A.dag()*A*A #some other Hamiltonian that the example uses

# --Coherency Plots--
# Perfect Transmon
if ((mode_dr == 0) and (mode_di == 0) and (mode_de == 0)):
    tlist = np.linspace(0, 10, 1000)
    c_ops = []  # Collapse operators
    e_ops = [g * g.dag(), e * e.dag(), (g + e), (g - e)]  # Expectation values

    result = mesolve(H, (g + e).unit(), tlist, c_ops, e_ops)

    fig, ax = plt.subplots(figsize=(12,6))
    plt.rcParams.update({'font.size': 22})
    ax.plot(tlist, np.real(result.expect[0]), 'b')
    ax.plot(tlist, np.real(result.expect[1]), 'r')
    ax.plot(tlist, np.real(result.expect[2]), 'm')
    ax.plot(tlist, np.real(result.expect[3]), 'k')
    ax.legend(("|0]", "|1]", "|0]+|1]", "|0]-|1]"))
    ax.set_xlabel('Time',fontsize=24)
    ax.set_ylabel('Probability',fontsize=24)
    plt.ylim([-.1,1.1])
    plt.title("Perfect Coherence",fontsize=34)

# Dissipation
if ((mode_dr == 0) and (mode_di == 1) and (mode_de == 0)):
    tlist = np.linspace(0,range,sample)
    c_ops_dis = []
    kappa_dis = np.power(T,-1)
    c_ops_dis.append(np.sqrt(kappa_dis)*A)
    e_ops_dis = [g*g.dag(),e*e.dag(),(g+e).unit()*(g+e).unit().dag(),(g-e).unit()*(g-e).unit().dag()]

    result_dis = mesolve(H,e,tlist,c_ops_dis,e_ops_dis)

    fig, ax = plt.subplots(figsize=(12,6))
    plt.rcParams.update({'font.size': 22})
    ax.plot(tlist, np.real(result_dis.expect[0]), 'b')
    ax.plot(tlist, np.real(result_dis.expect[1]), 'r')
    ax.plot(tlist, np.real(result_dis.expect[2]), 'm+')
    ax.plot(tlist, np.real(result_dis.expect[3]), 'm--')
    ax.legend(("|0]", "|1]", "|0]+|1]", "|0]-|1]"))
    ax.set_xlabel('Time',fontsize=24)
    ax.set_ylabel('Probability',fontsize=24)
    plt.ylim([-.1,1.1])
    plt.title("Photon Dissipation",fontsize=34)

#Dephasing
if ((mode_dr == 0) and (mode_di == 0) and (mode_de == 1)):
    tlist = np.linspace(0,range,sample)
    c_ops_dep = []
    kappa_dep = dephasing/2.0 # Half the dephasing time, in microseconds
    c_ops_dep.append(np.sqrt(kappa_dep)*A.dag()*A)
    e_ops_dep = [g*g.dag(),e*e.dag(),(g+e).unit()*(g+e).unit().dag(),(g-e).unit()*(g-e).unit().dag()]

    result_dep = mesolve(H,(g+e).unit(),tlist,c_ops_dep,e_ops_dep)

    fig, ax = plt.subplots(figsize=(12,6))
    plt.rcParams.update({'font.size': 22})
    ax.plot(tlist, np.real(result_dep.expect[0]), 'b')
    ax.plot(tlist, np.real(result_dep.expect[1]), 'r')
    ax.plot(tlist, np.real(result_dep.expect[2]), 'm+')
    ax.plot(tlist, np.real(result_dep.expect[3]), 'm--')
    ax.legend(("|0]", "|1]", "|0]+|1]", "|0]-|1]"))
    ax.set_xlabel('Time',fontsize=24)
    ax.set_ylabel('Probability',fontsize=24)
    plt.ylim([-.1,1.1])
    plt.title("Dephasing",fontsize=34)

#Dissipation and Dephasing
if ((mode_dr == 0) and (mode_di == 1) and (mode_de == 1)):
    tlist = np.linspace(0,range,sample)
    c_ops = []
    kappa_dep = dephasing/2.0 # Half the dephasing time (in microseconds)
    c_ops.append(np.sqrt(kappa_dep)*A.dag()*A)
    kappa_dis = T
    c_ops.append(np.sqrt(kappa_dis)*A)
    e_ops = [g*g.dag(),e*e.dag(),(g+e).unit()*(g+e).unit().dag(),(g-e).unit()*(g-e).unit().dag()]

    result = mesolve(H,(g+e).unit(),tlist,c_ops,e_ops)

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
    plt.title("Dissapation and Dephasing",fontsize=34)

# Position Driving, No Noise
if ((mode_dr == 1) and (mode_di == 0) and (mode_de == 0)):
    tlist = np.linspace(0,10,1000)
    c_ops = []
    e_ops_re = [g * g.dag(), e * e.dag(), (g + e).unit() * (g + e).unit().dag(),
                (g - e).unit() * (g - e).unit().dag(), (g + 1j * e).unit() * (g + 1j * e).unit().dag(),
                (g - 1j * e).unit() * (g - 1j * e).unit().dag()]
    H_dre = 1 * (A + A.dag())
    H_re = H + H_dre

    result_dre = mesolve(H_re, g, tlist, c_ops, e_ops_re)

    fig, ax = plt.subplots(figsize=(12, 6))
    plt.rcParams.update({'font.size': 22})
    ax.plot(tlist, np.real(result_dre.expect[0]), 'b')
    ax.plot(tlist, np.real(result_dre.expect[1]), 'r')
    ax.plot(tlist, np.real(result_dre.expect[2]), 'm+')
    ax.plot(tlist, np.real(result_dre.expect[3]), 'm--')
    ax.plot(tlist, np.real(result_dre.expect[4]), 'g+')
    ax.plot(tlist, np.real(result_dre.expect[5]), 'g--')
    ax.legend(("|0]", "|1]", "|0]+|1]", "|0]-|1]","|0]+i|1]","|0]-i|1]"))
    ax.set_xlabel('Time', fontsize=24)
    ax.set_ylabel('Probability', fontsize=24)
    plt.ylim([-.1, 1.1])
    plt.title("Position Drive with No Noise",fontsize=34)

# Position Driving, Dissipation
if ((mode_dr == 1) and (mode_di == 1) and (mode_de == 0)):
    tlist = np.linspace(0,range,sample)
    c_ops = []
    kappa_dis = np.power(T, -1)  # First number is T1 in microseconds
    c_ops.append(np.sqrt(kappa_dis) * A)
    e_ops_re = [g * g.dag(), e * e.dag(), (g + e).unit() * (g + e).unit().dag(),
                (g - e).unit() * (g - e).unit().dag(), (g + 1j * e).unit() * (g + 1j * e).unit().dag(),
                (g - 1j * e).unit() * (g - 1j * e).unit().dag()]
    H_dre = 1 * (A + A.dag())
    H_re = H + H_dre

    result_dre = mesolve(H_re, g, tlist, c_ops, e_ops_re)

    fig, ax = plt.subplots(figsize=(12, 6))
    plt.rcParams.update({'font.size': 22})
    ax.plot(tlist, np.real(result_dre.expect[0]), 'b')
    ax.plot(tlist, np.real(result_dre.expect[1]), 'r')
    ax.plot(tlist, np.real(result_dre.expect[2]), 'm+')
    ax.plot(tlist, np.real(result_dre.expect[3]), 'm--')
    ax.plot(tlist, np.real(result_dre.expect[4]), 'g+')
    ax.plot(tlist, np.real(result_dre.expect[5]), 'g--')
    ax.legend(("|0]", "|1]", "|0]+|1]", "|0]-|1]", "|0]+i|1]", "|0]-i|1]"))
    ax.set_xlabel('Time', fontsize=24)
    ax.set_ylabel('Probability', fontsize=24)
    plt.ylim([-.1, 1.1])
    plt.title("Position Drive with Dissipation", fontsize=34)

# Position Driving, Dephasing
if ((mode_dr == 1) and (mode_di == 0) and (mode_de == 1)):
    tlist = np.linspace(0,range,sample)
    c_ops = []
    kappa_dep = dephasing / 2  # Half the dephasing time (in microseconds)
    c_ops.append(np.sqrt(kappa_dep) * A.dag() * A)
    kappa_dis = T  # T1 in microseconds
    c_ops.append(np.sqrt(kappa_dis) * A)
    e_ops_re = [g * g.dag(), e * e.dag(), (g + e).unit() * (g + e).unit().dag(),
                (g - e).unit() * (g - e).unit().dag(), (g + 1j * e).unit() * (g + 1j * e).unit().dag(),
                (g - 1j * e).unit() * (g - 1j * e).unit().dag()]
    H_dre = 1 * (A + A.dag())
    H_re = H + H_dre

    result_dre = mesolve(H_re, g, tlist, c_ops, e_ops_re)

    fig, ax = plt.subplots(figsize=(12, 6))
    plt.rcParams.update({'font.size': 22})
    ax.plot(tlist, np.real(result_dre.expect[0]), 'b')
    ax.plot(tlist, np.real(result_dre.expect[1]), 'r')
    ax.plot(tlist, np.real(result_dre.expect[2]), 'm+')
    ax.plot(tlist, np.real(result_dre.expect[3]), 'm--')
    ax.plot(tlist, np.real(result_dre.expect[4]), 'g+')
    ax.plot(tlist, np.real(result_dre.expect[5]), 'g--')
    ax.legend(("|0]", "|1]", "|0]+|1]", "|0]-|1]", "|0]+i|1]", "|0]-i|1]"))
    ax.set_xlabel('Time', fontsize=24)
    ax.set_ylabel('Probability', fontsize=24)
    plt.ylim([-.1, 1.1])
    plt.title("Position Drive with Dephasing", fontsize=34)

# Position Driving, Dissipation and Dephasing
if ((mode_dr == 1) and (mode_di == 1) and (mode_de == 1)):
    tlist = np.linspace(0,range,sample)
    c_ops = []
    kappa_dep = dephasing / 2  # Half the dephasing time (in microseconds)
    c_ops.append(np.sqrt(kappa_dep) * A.dag() * A)
    kappa_dis = T  # T1 in microseconds
    c_ops.append(np.sqrt(kappa_dis) * A)
    e_ops_re = [g * g.dag(), e * e.dag(), (g + e).unit() * (g + e).unit().dag(),
                (g - e).unit() * (g - e).unit().dag(), (g + 1j * e).unit() * (g + 1j * e).unit().dag(),
                (g - 1j * e).unit() * (g - 1j * e).unit().dag()]
    H_dre = 1 * (A + A.dag())
    H_re = H + H_dre

    result_dre = mesolve(H_re, g, tlist, c_ops, e_ops_re)

    fig, ax = plt.subplots(figsize=(12, 6))
    plt.rcParams.update({'font.size': 22})
    ax.plot(tlist, np.real(result_dre.expect[0]), 'b')
    ax.plot(tlist, np.real(result_dre.expect[1]), 'r')
    ax.plot(tlist, np.real(result_dre.expect[2]), 'm+')
    ax.plot(tlist, np.real(result_dre.expect[3]), 'm--')
    ax.plot(tlist, np.real(result_dre.expect[4]), 'g+')
    ax.plot(tlist, np.real(result_dre.expect[5]), 'g--')
    ax.legend(("|0]", "|1]", "|0]+|1]", "|0]-|1]", "|0]+i|1]", "|0]-i|1]"))
    ax.set_xlabel('Time', fontsize=24)
    ax.set_ylabel('Probability', fontsize=24)
    plt.ylim([-.1, 1.1])
    plt.title("Position Drive with Dissipation and Dephasing", fontsize=34)

# Momentum Driving, No Noise
if ((mode_dr == 2) and (mode_di == 0) and (mode_de == 0)):
    tlist = np.linspace(0,10,1000)
    c_ops = []
    e_ops_re = [g * g.dag(), e * e.dag(), (g + e).unit() * (g + e).unit().dag(),
                (g - e).unit() * (g - e).unit().dag(), (g + 1j * e).unit() * (g + 1j * e).unit().dag(),
                (g - 1j * e).unit() * (g - 1j * e).unit().dag()]
    H_dim = 10j * (A - A.dag())
    H_im = H + H_dim

    result_dim = mesolve(H_im, g, tlist, c_ops, e_ops_re)

    fig, ax = plt.subplots(figsize=(12, 6))
    plt.rcParams.update({'font.size': 22})
    ax.plot(tlist, np.real(result_dim.expect[0]), 'b')
    ax.plot(tlist, np.real(result_dim.expect[1]), 'r')
    ax.plot(tlist, np.real(result_dim.expect[2]), 'm+')
    ax.plot(tlist, np.real(result_dim.expect[3]), 'm--')
    ax.plot(tlist, np.real(result_dim.expect[4]), 'g+')
    ax.plot(tlist, np.real(result_dim.expect[5]), 'g--')
    ax.legend(("|0]", "|1]", "|0]+|1]", "|0]-|1]", "|0]+i|1]", "|0]-i|1]"))
    ax.set_xlabel('Time', fontsize=24)
    ax.set_ylabel('Probability', fontsize=24)
    plt.ylim([-.1, 1.1])
    plt.title("Momentum Drive with No Noise",fontsize=34)

# Momentum Driving, Dissipation
if ((mode_dr == 2) and (mode_di == 1) and (mode_de == 0)):
    tlist = np.linspace(0,range,sample)
    c_ops = []
    kappa_dis = np.power(T, -1)  # First number is T1 in microseconds
    c_ops.append(np.sqrt(kappa_dis) * A)
    e_ops_re = [g * g.dag(), e * e.dag(), (g + e).unit() * (g + e).unit().dag(),
                (g - e).unit() * (g - e).unit().dag(), (g + 1j * e).unit() * (g + 1j * e).unit().dag(),
                (g - 1j * e).unit() * (g - 1j * e).unit().dag()]
    H_dim = 10j * (A - A.dag())
    H_im = H + H_dim

    result_dim = mesolve(H_im, g, tlist, c_ops, e_ops_re)

    fig, ax = plt.subplots(figsize=(12, 6))
    plt.rcParams.update({'font.size': 22})
    ax.plot(tlist, np.real(result_dim.expect[0]), 'b')
    ax.plot(tlist, np.real(result_dim.expect[1]), 'r')
    ax.plot(tlist, np.real(result_dim.expect[2]), 'm+')
    ax.plot(tlist, np.real(result_dim.expect[3]), 'm--')
    ax.plot(tlist, np.real(result_dim.expect[4]), 'g+')
    ax.plot(tlist, np.real(result_dim.expect[5]), 'g--')
    ax.legend(("|0]", "|1]", "|0]+|1]", "|0]-|1]", "|0]+i|1]", "|0]-i|1]"))
    ax.set_xlabel('Time', fontsize=24)
    ax.set_ylabel('Probability', fontsize=24)
    plt.ylim([-.1, 1.1])
    plt.title("Momentum Drive with Dissipation", fontsize=34)

# Momentum Driving, Dephasing
if ((mode_dr == 2) and (mode_di == 0) and (mode_de == 1)):
    tlist = np.linspace(0,range,sample)
    c_ops = []
    kappa_dep = dephasing / 2  # Half the dephasing time (in microseconds)
    c_ops.append(np.sqrt(kappa_dep) * A.dag() * A)
    kappa_dis = T  # T1 in microseconds
    c_ops.append(np.sqrt(kappa_dis) * A)
    e_ops_re = [g * g.dag(), e * e.dag(), (g + e).unit() * (g + e).unit().dag(),
                (g - e).unit() * (g - e).unit().dag(), (g + 1j * e).unit() * (g + 1j * e).unit().dag(),
                (g - 1j * e).unit() * (g - 1j * e).unit().dag()]
    H_dim = 10j * (A - A.dag())
    H_im = H + H_dim

    result_dim = mesolve(H_im, g, tlist, c_ops, e_ops_re)

    fig, ax = plt.subplots(figsize=(12, 6))
    plt.rcParams.update({'font.size': 22})
    ax.plot(tlist, np.real(result_dim.expect[0]), 'b')
    ax.plot(tlist, np.real(result_dim.expect[1]), 'r')
    ax.plot(tlist, np.real(result_dim.expect[2]), 'm+')
    ax.plot(tlist, np.real(result_dim.expect[3]), 'm--')
    ax.plot(tlist, np.real(result_dim.expect[4]), 'g+')
    ax.plot(tlist, np.real(result_dim.expect[5]), 'g--')
    ax.legend(("|0]", "|1]", "|0]+|1]", "|0]-|1]", "|0]+i|1]", "|0]-i|1]"))
    ax.set_xlabel('Time', fontsize=24)
    ax.set_ylabel('Probability', fontsize=24)
    plt.ylim([-.1, 1.1])
    plt.title("Momentum Drive with Dephasing", fontsize=34)

# Momentum Driving, Dissipation and Dephasing
if ((mode_dr == 2) and (mode_di == 1) and (mode_de == 1)):
    tlist = np.linspace(0,range,sample)
    c_ops = []
    kappa_dep = dephasing / 2  # Half the dephasing time (in microseconds)
    c_ops.append(np.sqrt(kappa_dep) * A.dag() * A)
    kappa_dis = T  # T1 in microseconds
    c_ops.append(np.sqrt(kappa_dis) * A)
    e_ops_re = [g * g.dag(), e * e.dag(), (g + e).unit() * (g + e).unit().dag(),
                (g - e).unit() * (g - e).unit().dag(), (g + 1j * e).unit() * (g + 1j * e).unit().dag(),
                (g - 1j * e).unit() * (g - 1j * e).unit().dag()]
    H_dim = 10j * (A - A.dag())
    H_im = H + H_dim

    result_dim = mesolve(H_im, g, tlist, c_ops, e_ops_re)

    fig, ax = plt.subplots(figsize=(12, 6))
    plt.rcParams.update({'font.size': 22})
    ax.plot(tlist, np.real(result_dim.expect[0]), 'b')
    ax.plot(tlist, np.real(result_dim.expect[1]), 'r')
    ax.plot(tlist, np.real(result_dim.expect[2]), 'm+')
    ax.plot(tlist, np.real(result_dim.expect[3]), 'm--')
    ax.plot(tlist, np.real(result_dim.expect[4]), 'g+')
    ax.plot(tlist, np.real(result_dim.expect[5]), 'g--')
    ax.legend(("|0]", "|1]", "|0]+|1]", "|0]-|1]", "|0]+i|1]", "|0]-i|1]"))
    ax.set_xlabel('Time', fontsize=24)
    ax.set_ylabel('Probability', fontsize=24)
    plt.ylim([-.1, 1.1])
    plt.title("Momentum Drive with Dissipation and Dephasing", fontsize=34)
# --/Coherency Plots--

plt.show()