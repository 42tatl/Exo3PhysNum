import numpy as np
import matplotlib.pyplot as plt
import os
import functions as fct

repertoire = r"C:\Users\Avril\Desktop\Exo3PhysNum"  
executable = os.path.join(repertoire, "Exe.exe")
os.chdir(repertoire)
input_filename = "configb.txt"

params = fct.read_in_file(input_filename)

tFin, m1, m2, x0, v0x, v0y, a, tol = fct.get_params(params)

outputs_nsteps, params = fct.run_param_sweep(executable, input_filename, "nsteps", [1e2, 1e3, 1e4, 1e5, 1e6], {"adapt": 0, "tol": 30})

alpha = m2/(m1 + m2)
beta = m1/(m1 + m2)

x0 = 2*a - alpha*a
xS_prime = -alpha*a
xJ_prime = beta*a
rad_S = 6.96240e8
rad_J = 6.9911e7

# Plotting the trajectories
plt.figure()

for output, param in zip(outputs_nsteps, params):  # zip them properly now
    t, vx, vy, x, y, E, nsteps = np.loadtxt(output, unpack=True)
    plt.plot(x, y, label=rf"$nsteps$ = {int(param['nsteps'])}")
    
plt.grid()
plt.plot(x0, 0, "ro", label="Initial asteroid position")
plt.plot
plt.plot(xS_prime, 0, "bo", label="Sun")
plt.plot(xJ_prime, 0, "go", label="Jupiter")
#Zoom1
'''
plt.xlim(-1e11, 1.8e12)
plt.ylim(-9e11, 7e11)
'''

#Zoom2
plt.xlim(-5e12, 6.5e12)
plt.ylim(-5.5e12, 6.5e12)

plt.xlabel(r"$x'\ [\mathrm{m}]$")
plt.ylabel(r"$y'\ [\mathrm{m}]$")
plt.legend()
fct.save_figure("traj_b_Rprime_nsteps_zoom2.png")
plt.show()

