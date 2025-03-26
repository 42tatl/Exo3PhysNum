import numpy as np
import subprocess
import matplotlib.pyplot as plt
import os
import functions as fct

repertoire = r"C:\Users\Avril\Desktop\Exo3PhysNum"  
executable = os.path.join(repertoire, "Exe42.exe")
os.chdir(repertoire)
input_filename = "configb.txt"
output_template = "output_{paramstr}_{value}.out"

params = fct.read_in_file(input_filename)

tFin, m1, m2, x0, v0x, v0y, a, tol = fct.get_params(params)

#outputs_nsteps = fct.run_param_sweep(executable, input_filename, "nsteps", [1000], {"adapt": 0, "tol": 30})
outputs_tol = fct.run_param_sweep(executable, input_filename, "tol", [1e-3], {"adapt": 1, "nsteps": 8000})

plt.figure()
plt.title("Trajectories")
for output in outputs_tol:
    t, vx, vy, x, y, E, nsteps = np.loadtxt(output, unpack=True)
    plt.plot(x, y, label=output)
    plt.legend()
    plt.grid()

plt.show()
