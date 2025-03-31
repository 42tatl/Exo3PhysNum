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

outputs_tol, params = fct.run_param_sweep(executable, input_filename, "tol", [1e5, 1e-1, 1e-5], {"adapt": 1})

alpha = m2/(m1 + m2)
beta = m1/(m1 + m2)

x0 = 2*a - alpha*a
xS_prime = -alpha*a
xJ_prime = beta*a
rad_S = 6.96240e8
rad_J = 6.9911e7

#TRAJECTORY

plt.figure()
for output, param in zip(outputs_tol, params):  
    t, vx, vy, x, y, E, nsteps = fct.read_output_file(output)
    plt.plot(x, y, label=rf"$tol$ = {int(param['tol'])}")

plt.grid()
plt.plot(x0, 0, "ro", label="Initial asteroid position")
plt.plot
plt.plot(xS_prime, 0, "bo", label="Sun")
plt.plot(xJ_prime, 0, "go", label="Jupiter")
plt.xlabel(r"$x'\ [\mathrm{m}]$")
plt.ylabel(r"$y'\ [\mathrm{m}]$")
plt.legend()
fct.save_figure("traj_b_Rprime_tol.png")
plt.show()

#ENERGY VS TIME

plt.figure()
for output, param in zip(outputs_tol, params):  
    t, vx, vy, x, y, E, nsteps = fct.read_output_file(output)
    plt.plot(t, E, label=rf"$tol$ = {int(param['tol'])}")

plt.grid()
plt.xlabel(r"$t'\ [\mathrm{s}]$")
plt.ylabel(r"$E'\ [\mathrm{J}\cdot K^{-1}]$")
plt.legend()
fct.save_figure("E_b_Rprime_tol.png")
plt.show()

#ENERGY CONSERVATION

plt.figure()
for output, param in zip(outputs_tol, params):  
    t, vx, vy, x, y, E, nsteps = fct.read_output_file(output)
    plt.plot(t, E, label=rf"$tol$ = {int(param['tol'])}")
    dE = max(E) - min(E)

plt.grid()
plt.xlabel(r"$nsteps (adaptive)'\ [\mathrm{s}]$")
plt.ylabel(r"$dE'\ [\mathrm{J}\cdot K^{-1}]$")
plt.legend()
fct.save_figure("dE_b_Rprime_tol.png")
plt.show()

