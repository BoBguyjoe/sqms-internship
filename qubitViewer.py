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

wc = frequency*2*np.pi
wa = anharmonicity*2*np.pi
g = g*2*np.pi
qubit = scq.Transmon(EJ=Ej, EC=Ec, ng=ng, ncut=150)
cavity = scq.Oscillator(E_osc=wc, truncated_dim=2)
space = scq.HilbertSpace([qubit, cavity])
space.add_interaction(g_strength=g, op1=qubit.n_operator, op2=cavity.creation_operator, add_hc=True)

plotStates = int(input("Plot energy states? {0} for no, {1} for yes"))
if (plotStates != 1 and plotStates != 0):
    print("Yo, get your act together. " + str(plotStates) + " wasn't an option.")
    exit()
elif (plotStates == 1):
    print("Plot energy states against which parameter:")
    print("1: Offset Charge (ng)")
    print("2: Josephson Energy (Ej)")
    print("3: Charging Energy (Ec)")
    parameter1 = int(input())
    if (parameter1 != 1 and parameter1 != 2 and parameter1 != 3):
        print("Yo, get your act together. " + str(parameter1) + " wasn't an option.")
        exit()
    elif (parameter1 == 1):
        parameter = "ng"
    elif (parameter1 == 2):
        parameter = "EJ"
    elif (parameter1 == 3):
        parameter = "EC"
    param_max1 = float(input("Enter range to sweep over (energy states will be plotted against +/- the entered value for the parameter you specified): "))
    param_range1 = np.linspace(-1*param_max1,param_max1,100)

print()
plotT1 = int(input("Plot effective T1? {0} for no, {1} for yes"))
if (plotT1 != 1 and plotT1 != 0):
    print("Yo, get your act together. " + str(plotT1) + " wasn't an option.")
    exit()
elif (plotT1 == 1):
    print("Plot effective T1 against which parameter:")
    print("1: Offset Charge (ng)")
    print("2: Josephson Energy (Ej)")
    print("3: Charging Energy (Ec)")
    parameter2 = int(input())
    if (parameter2 != 1 and parameter2 != 2 and parameter2 != 3):
        print("Yo, get your act together. " + str(parameter2) + " wasn't an option.")
        exit()
    elif (parameter2 == 1):
        parameter2 = "ng"
    elif (parameter2 == 2):
        parameter2 = "EJ"
    elif (parameter2 == 3):
        parameter2 = "EC"
    param_max2 = float(input("Enter range to sweep over (effective T1 will be plotted against +/- the entered value for the parameter you specified): "))
    param_range2 = np.linspace(-1*param_max2,param_max2,100)

print()
plotTransitions = int(input("Plot transition energies? {0} for no, {1} for yes."))
if (plotTransitions != 1 and plotTransitions != 0):
    print("Yo, get your act together. " + str(plotTransitions) + " wasn't an option.")
    exit()
elif (plotTransitions == 1):
    print("Plot transition energies against which parameter:")
    print("1: Offset Charge (ng)")
    print("2: Josephson Energy (Ej)")
    print("3: Charging Energy (Ec)")
    parameter3 = int(input())
    if (parameter3 != 1 and parameter3 != 2 and parameter3 != 3):
        print("Yo, get your act together. " + str(parameter3) + " wasn't an option.")
        exit()
    elif (parameter3 == 1):
        parameter3 = "ng"
    elif (parameter3 == 2):
        parameter3 = "EJ"
    elif (parameter3 == 3):
        parameter3 = "EC"
    param_max3 = float(input(
        "Enter range to sweep over (transition energies will be plotted against the parameter from 0 to this number): "))
    param_range3 = np.linspace(0, param_max3, 100)


# Create energy states plot
if (plotStates == 1):
    qubit.plot_evals_vs_paramvals(parameter1, param_range1, evals_count=6, subtract_ground=False)

# Create effective T1 plot
if (plotT1 == 1):
    qubit.plot_t1_effective_vs_paramvals(parameter2, param_range2)

# Create transition energies plot
if (plotTransitions == 1):
    qubit.plot_dispersion_vs_paramvals('ng',parameter3,param_range3,transitions=(((0,1), (0,2))))

plt.show()