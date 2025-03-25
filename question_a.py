import numpy as np
import subprocess
import matplotlib.pyplot as plt
import os
import functions as fct

executable = './Exe'  # Remove .exe for Mac
repertoire = r"/Users/lilimouelle/Desktop/PHYSNUM/Exo3PhysNum"  # Modify for correct directory
os.chdir(repertoire)

input_filename = "configa.in"
output_template = "output_{paramstr}_{value}.out"

params = fct.read_in_file(input_filename)

tFin, m1, m2, x0, v0x, v0y, a, tol = fct.get_params(params)

outputs_nsteps = fct.run_param_sweep(executable, input_filename, "nsteps", [10, 20, 50, 100], {"adapt": False, "tol": 30})
#outputs_tol = fct.run_param_sweep(executable, input_filename, "tol", [10, 20, 30], {"adapt": True, "nsteps": 100})


def read_output_file(filename):
    data = np.loadtxt(filename)
    t = data[:, 0]
    x = data[:, 1]
    y = data[:, 2]
    vx = data[:, 3]
    vy = data[:, 4]
    energy = data[:, 5]
    return t, x, y, vx, vy, energy

def plot_trajectory(x, y, label=""):
    plt.plot(x, y, label=label)
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    plt.axis("equal")
    plt.grid(True)


#Fixed time step
for output in outputs_nsteps:
    t, x, y, vx, vy, energy = read_output_file(output)
    label = output.split("/")[-1]
    plot_trajectory(x, y, label=label)
plt.title("Trajectoires avec pas de temps fixe")
plt.legend()
plt.show()



