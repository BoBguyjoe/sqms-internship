PACKAGES
qutip
scqubits
numpy
matplotlib
ffmpeg
h5py
ipywidgets

HOW TO USE
The "Set Parameters" section of the script is where you can input the parameters of the qubit to simulate. The parameters.txt file has a set of such parameters used in an actual qubit.

When you run the script, it will ask you a series of questions that will allow you to specify how to initialize the qubit, what type(s) of decoherence and noise you want, and what form(s) you want the output to be in.

TO FIX
 - When the qubit is initialized at halfway between |0] and |1] and only dissipation is used, the state still dephases
 - Initializing the qubit halfway between |0] and |1] puts and moves the vector along the xz plane, but initializing it at |1] makes it move along the yz plane
 - The bloch sphere animation will throw a warning and an error. It still works but it's annoying. (UPDATE: the given way of suppressing the warning will make the figure blank)

TO DO
 - Add using pulses (pi and pi/2 and the like) to initialize the qubit instead of simply starting the qubit in the desired state
 - Add more parameter sets to the parameters.txt file
 - Replace the current Hamiltonian with the more complete JC hamiltonian (it currently throws a variety of errors when used)
 - Show the current time (as it passes) in Bloch sphere animations

CHANGELOG
7/15/2021:
 - User can now specify a name for the sphere animation file
 - Printing the expectation values now works
 - Script now indicates to the user when the sphere animation has been saved
 - Re-added driving (it simply sums a creation and destruction operator to the hamiltonian, as per a piece of qutip example code)
 - Added objects for the resonator and hilbert space (setup for parameter sweeps)
 - Can now plot eigenstates