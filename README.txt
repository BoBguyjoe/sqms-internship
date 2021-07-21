PACKAGES
qutip
scqubits
numpy
matplotlib
ffmpeg
h5py
ipywidgets
cmath

HOW TO USE
The 'parameters' file is where you enter the parameters of the qubit you want to simulate. The example parameters.txt file has a set of such parameters used in an actual qubit.

The main script asks you a series of questions to specify what type(s) of decoherence and noise you want, and what form(s) you want the output to be in.

The initializer script allows you to initialize a qubit manually, starting in the |0] state, and using sigma-x, sigma-y, and sigma-z operators

The qubitViewer script allows you to plot various characteristics of the qubit against ranges of certain parameters.

TO FIX
 - When the qubit is initialized at halfway between |0] and |1] and only dissipation is used, the state still dephases
 - Initializing the qubit halfway between |0] and |1] puts and moves the vector along the xz plane, but initializing it at |1] makes it move along the yz plane
 - The bloch sphere animation will throw a warning and an error. It still works but it's annoying.

TO DO
 - Add more sweeps to qubitViewer
 - Add more parameter sets to the parameters.txt file
 - Replace the current Hamiltonian with the more complete JC hamiltonian (it currently throws a variety of errors when used)
 - FINAL GOAL: Add in a mode that simulates a driving pulse leading into Rabi oscillations

CHANGELOG
7/20/2021:
 - Moved the initializer script into the simulation script
 - Fixed initializer only showing bloch sphere plot once
 - Added pulse script, which is start of the rabi oscillations simulation (not for use)
 - Added transition energies plot to qubitViewer (will hopefully be helpful for the rabi oscillator pulse)

7/19/2021:
 - Fixed the qubit setting its Ej value to the specified Ec value (whoops)
 - Moved the eigenstate plot to a separate python file 'qubitViewer' and added other plot options
 - Renamed the parameters txt file to 'example parameters' and made a new 'parameters' file that the other scripts read from
 - Added an 'initializer' script that allows user to set up a qubit state manually for the main script to use

7/15/2021:
 - User can now specify a name for the sphere animation file
 - Printing the expectation values now works
 - Script now indicates to the user when the sphere animation has been saved
 - Re-added driving (it simply sums a creation and destruction operator to the hamiltonian, as per a piece of qutip example code)
 - Added objects for the resonator and hilbert space (setup for parameter sweeps)
 - Can now plot eigenstates