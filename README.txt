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
 - The bloch sphere animation will throw a warning and an error. It still works but it's annoying.

TO DO
 - The script was restructured recently, and the qubit driving hasn't been added back in yet
 - Allow user to specify the name of the file to save the bloch sphere animation as so that it doesn't keep overwriting the same one
 - Add using pulses (pi and pi/2 and the like) to initialize the qubit instead of simply starting the qubit in the desired state
 - Add more parameter sets to the parameters.txt file
