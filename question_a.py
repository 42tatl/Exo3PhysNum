import numpy as np
import subprocess
import matplotlib.pyplot as plt
import os
import functions as fct

executable = './Exe3'  # Remove .exe for Mac
repertoire = r"/Users/lilimouelle/Desktop/PHYSNUM/Exo3PhysNum"  # Modify for correct directory
os.chdir(repertoire)

input_filename = "configa.in"
output_template = "output_{paramstr}_{value}.out"

params = fct.read_in_file(input_filename)

tFin, m1, m2, x0, v0x, v0y, a, tol = fct.get_params(params)

outputs_nsteps = fct.run_param_sweep(executable, input_filename, "nsteps", [8000], {"adapt": 0, "tol": 30})
#outputs_tol = fct.run_param_sweep(executable, input_filename, "tol", [1,0.1,0.01], {"adapt": True, "nsteps": 1000})


def read_output_file(filename):
    data = np.loadtxt(filename)
    t = data[:, 0]
    x = data[:, 3]
    y = data[:, 4]
    vx = data[:, 1]
    vy = data[:, 2]
    energy = data[:, 5]
    return t, x, y, vx, vy, energy

def plot_trajectory(x, y, label=""):
    plt.plot(x, y, label=label)
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    #plt.axis("equal")
    plt.grid(True)

def plot_energy(t, energy, label=""):
    plt.plot(t, energy, label=label)
    plt.xlabel("t [s]")
    plt.ylabel("Energie per mass [J/kg]")
    plt.grid(True)


#Fixed time step
for output in outputs_nsteps:
    t, x, y, vx, vy, energy = read_output_file(output)
    label = output.split("/")[-1]
    plot_trajectory(x, y, label=label)
plt.legend()
plt.savefig('trajectoires.pdf')
plt.show()


for output in outputs_nsteps:
    t, x, y, vx, vy, energy = read_output_file(output)
    plot_energy(t, energy, label=output)
plt.legend()
plt.savefig('energie.pdf')
plt.show()



'''
#Adaptive time step
for output in outputs_tol:
    t, x, y, vx, vy, energy = read_output_file(output)
    label = output.split("/")[-1]
    plot_trajectory(x, y, label=label)
plt.legend()
plt.savefig('trajectoires_adapt.pdf')
plt.show()

for output in outputs_tol:
    t, x, y, vx, vy, energy = read_output_file(output)
    plot_energy(t, energy, label=output)
plt.legend()
plt.savefig('energie_adapt.pdf')
plt.show()
'''