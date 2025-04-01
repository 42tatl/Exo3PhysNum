import numpy as np
import matplotlib.pyplot as plt
import os
import functions as fct

repertoire = r"C:\Users\Avril\Desktop\Exo3PhysNum"  
executable = os.path.join(repertoire, "Exe.exe")
os.chdir(repertoire)
input_filename = "configb.txt"

params = fct.read_in_file(input_filename)

#--------------------PARAMETERS--------------------

tFin, m1, m2, x0, v0x, v0y, a, tol = fct.get_params(params)

GM = 6.674e-11
mtot = m1 + m2
Omega = np.sqrt(GM*(mtot)/(a**3))
alpha = m2/(m1 + m2)
beta = m1/(m1 + m2)

x0 = 2*a - alpha*a
xS_prime = -alpha*a
xJ_prime = beta*a
rad_S = 6.96240e8
rad_J = 6.9911e7

outputs_tol, params1 = fct.run_param_sweep(executable, input_filename, "tol", [1e-1], {"adapt": 1, "nsteps" : 5000})

#--------------------R' REFERENCE FRAME--------------------

dE1 = []
dE2 = []
jsteps = []

for output, param in zip(outputs_tol, params):  
    t, x, y, vx, vy, E, nsteps = fct.read_output_file(output)
    dE1 = max(E) - min(E)
    jsteps.append(nsteps[-1])

outputs_nsteps, params2 = fct.run_param_sweep(executable, input_filename, "nsteps", jsteps, {"adapt": 0})

plt.figure()
plt.title("Non-inertial frame")
plt.plot(dE1, jsteps, label=r"$\Delta E_{tol}$", marker="o", markersize=5, linestyle="--")
plt.plot(dE1, jsteps, label=r"$\Delta E_{nsteps}$", marker="o", markersize=5, linestyle="--")
plt.grid()
plt.xlabel(r"$nsteps$")
plt.ylabel(r"$\Delta E$")
plt.legend()
fct.save_figure("dE_b_R_1.png")
plt.show()

plt.figure()
plt.title("Non-inertial frame")
plt.plot(dE1/dE2, jsteps, label=r"$\Delta E_{tol}/\Delta E_{nsteps}$", marker="o", markersize=5, linestyle="--")
plt.grid()
plt.xlabel(r"$nsteps$")
plt.ylabel(r"$\Delta E_{tol}/\Delta E_{nsteps}$")
plt.legend()
fct.save_figure("dE_b_R_1.png")
plt.show()

#---------------------R REFERENCE FRAME--------------------

xS_prime = -alpha * a
xJ_prime = beta * a

# dE vs nsteps in R

dE1 = []
dE2 = []
jsteps = []

for output, param in zip(outputs_tol, params1):
    t, x_prime, y_prime, vx_prime, vy_prime, E_rot, nsteps = fct.read_output_file(output)

    cosO = np.cos(Omega * t)
    sinO = np.sin(Omega * t)
    x = x_prime * cosO - y_prime * sinO
    y = x_prime * sinO + y_prime * cosO

    vx_eff = vx_prime - Omega * y_prime
    vy_eff = vy_prime + Omega * x_prime
    vx = vx_eff * cosO - vy_eff * sinO
    vy = vx_eff * sinO + vy_eff * cosO

    rS = np.sqrt((x + alpha * a)**2 + y**2)
    rJ = np.sqrt((x - beta * a)**2 + y**2)
    U = -GM * (m1 / rS + m2 / rJ)

    E_inertial = 0.5 * (vx**2 + vy**2) + U

    dE = max(E_inertial) - min(E_inertial)
    jsteps.append(nsteps[-1])
    dE1.append(dE)

for output, param in zip(outputs_nsteps, params2):
    t, x_prime, y_prime, vx_prime, vy_prime, E_rot, nsteps = fct.read_output_file(output)

    cosO = np.cos(Omega * t)
    sinO = np.sin(Omega * t)
    x = x_prime * cosO - y_prime * sinO
    y = x_prime * sinO + y_prime * cosO

    vx_eff = vx_prime - Omega * y_prime
    vy_eff = vy_prime + Omega * x_prime
    vx = vx_eff * cosO - vy_eff * sinO
    vy = vx_eff * sinO + vy_eff * cosO

    rS = np.sqrt((x + alpha * a)**2 + y**2)
    rJ = np.sqrt((x - beta * a)**2 + y**2)
    U = -GM * (m1 / rS + m2 / rJ)

    E_inertial = 0.5 * (vx**2 + vy**2) + U

    dE = max(E_inertial) - min(E_inertial)
    dE2.append(dE)

plt.figure()
plt.title("Inertial frame")
plt.plot(dE1, jsteps, label=r"$\Delta E_{tol}$", marker="o", markersize=5, linestyle="--")
plt.plot(dE1, jsteps, label=r"$\Delta E_{nsteps}$", marker="o", markersize=5, linestyle="--")
plt.grid()
plt.xlabel(r"$nsteps$")
plt.ylabel(r"$\Delta E$")
plt.legend()
fct.save_figure("dE_b_R_1.png")
plt.show()

plt.figure()
plt.title("Inertial frame")
plt.plot(dE1/dE2, jsteps, label=r"$\Delta E_{tol}/\Delta E_{nsteps}$", marker="o", markersize=5, linestyle="--")
plt.grid()
plt.xlabel(r"$nsteps$")
plt.ylabel(r"$\Delta E_{tol}/\Delta E_{nsteps}$")
plt.legend()
fct.save_figure("dE_b_R_1.png")
plt.show()

