PACKAGES
qutip
scqubits
numpy
matplotlib
ffmpeg
h5py
ipywidgets
cmath
scipy

HOW TO USE
The 'parameters' file is where you enter the parameters of the qubit you want to simulate. The example parameters.txt file has a set of such parameters used in an actual qubit.

main.py is the simulation. It initializes a qubit that you can apply dissipation and/or a drive pulse to to see what happens

sweeps.py sweeps over a range of amplitudes and/or lengths of pulses and outputs what are called Rabi oscillations

CHANGELOG
8/2/2021:
 - Added a 2D rabi sweep to sweeps.py, where it sweeps amplitude and length at the same time
 - Fixed up which expectation values main.py plots (having a drive and noise makes it behave unexpectedly)
 - Added capability for gaussian pulse

7/30/2021:
 - Folded the single drive pulse changing the state to the main simulation
 - Folded the pulse sweep "experiment" to qubitViewer
 - Added Bloch sphere animation to the rabi oscillations sweep
 - Renamed 'qubitViewer' to 'sweeps'

7/29/2021:
 - Added a 'pulse.py' file that collects a bunch of state measurements of a qubit when driven by a square pulse of varying length.

7/23/2021:
 - Switched to the JC hamiltonian in the rotating frame
 - Corrected the various operators to be tensor products, so that they actually act on the right parts of the system
 - The default hamiltonian is the first and third terms of this JC hamiltonian. Driving adds in the usual second term
 - Fixed bug that would prevent the simulation from using calculated Tphi value
 - The initializer has been removed for the time being (until I fix it or decide to cut it)

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