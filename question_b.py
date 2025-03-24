import numpy as np
import subprocess
import matplotlib.pyplot as plt
import os
import functions as fct

repertoire = r"C:\Users\Avril\Desktop\Exo3PhysNum"  
executable = os.path.join(repertoire, "Exe.exe")
os.chdir(repertoire)
input_filename = "config_a"
output_template = "output_{paramstr}_{value}.out"

params = fct.read_in_file(input_filename)

tFin, m1, m2, x0, v0x, v0y, a, tol = fct.get_params(params)

outputs_nsteps = fct.run_param_sweep(executable, input_filename, "nsteps", [10, 20, 50, 100], {"adapt": False, "tol": 30})
outputs_tol = fct.run_param_sweep(executable, input_filename, "tol", [10, 20, 30], {"adapt": True, "nsteps": 100})


