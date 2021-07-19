from qutip import *
import numpy as np
import cmath
import matplotlib.pyplot as plt

sphere = qutip.Bloch()
sphere.vector_color = ['r']
psi0 = basis(2,0)

i = True
while (i == True):
    print("The qubit's state is: " + str(psi0))
    print("Which operator would you like to apply?")
    print("1: Rotate around x axis")
    print("2: Rotate around y axis")
    print("3: Rotate around z axis")
    print("4: Finish and exit")
    action = int(input())
    if (action != 1 and action != 2 and action != 3 and action != 4):
        print("Yo, get your act together. " + str(action) + " wasn't an option.")
        exit()

    if (action == 1 or action == 2 or action == 3):
        rotate = float(input("How much rotation: {1} for half-pi, {2} for pi"))/2.0
    if (action == 1):
        psi0 = rotate*sigmax()*psi0
    if (action == 2):
        psi0 = rotate*sigmay()*psi0
    if (action == 3):
        psi0 = rotate*sigmaz()*psi0
    if (action == 4):
        i = False

    theta = np.pi * psi0[1][0][0].real
    phi = np.pi * psi0[1][0][0].imag
    sphere.add_vectors([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)])

    sphere.make_sphere()
    print("Here I am!")
    plt.show()
    print("There I am!")

print("The final state of the qubit is: " + str(psi0))