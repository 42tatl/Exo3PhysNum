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

outputs_tol, params = fct.run_param_sweep(executable, input_filename, "tol", [1e-1], {"adapt": 1, "nsteps" : 5000})

#--------------------R' REFERENCE FRAME--------------------

#-----------TRAJECTORY---------------

plt.figure()
for output, param in zip(outputs_tol, params):  
    t, x, y, vx, vy, E, nsteps = fct.read_output_file(output)
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

#-----------CONVERGENCE FINAL POSITION VS JSTEPS----------
plt.figure()
for output, param in zip(outputs_tol, params):  
    t, x, y, vx, vy, E, nsteps = fct.read_output_file(output)
    norm_final = np.sqrt(pow(x[-1],2) + pow(y[-1],2))
    jsteps = nsteps[-1]
    plt.plot(jsteps, norm_final, label=rf"$tol$ = {int(param['tol'])}")

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

#---------ENERGY VS TIME----------------

plt.figure()
for output, param in zip(outputs_tol, params):  
    t, x, y, vx, vy, E, nsteps = fct.read_output_file(output)
    plt.plot(t, E, label=rf"$tol$ = {int(param['tol'])}")

plt.grid()
plt.xlabel(r"$t'\ [\mathrm{s}]$")
plt.ylabel(r"$E'\ [\mathrm{J}\cdot K^{-1}]$")
plt.legend()
fct.save_figure("E_b_Rprime_tol.png")
plt.show()


#---------------------R REFERENCE FRAME--------------------
# Transform R′ (rotating frame) → R inertial frame)
#The asteroid is at (x, y) in the inertial frame.

plt.figure()
# Time-independent Sun/Jupiter positions in R'
xS_prime = -alpha * a
xJ_prime = beta * a

# Compute their motion in R
for output, param in zip(outputs_tol, params):
    t, x_prime, y_prime, vx_prime, vy_prime, E_rot, nsteps = fct.read_output_file(output)

    # Transform position: rotate back
    cosO = np.cos(Omega * t)
    sinO = np.sin(Omega * t)
    x = x_prime * cosO - y_prime * sinO
    y = x_prime * sinO + y_prime * cosO

    # Sun and Jupiter trajectories in R
    x_S = xS_prime * cosO
    y_S = xS_prime * sinO
    x_J = xJ_prime * cosO
    y_J = xJ_prime * sinO

    # Velocity transformation: inertial velocity = rotated (v' + Omega × r')
    vx_eff = vx_prime - Omega * y_prime
    vy_eff = vy_prime + Omega * x_prime
    vx = vx_eff * cosO - vy_eff * sinO
    vy = vx_eff * sinO + vy_eff * cosO
    
    plt.plot(x0, 0, "ro", label="Initial asteroid position" if param is params[0] else "")
    plt.plot(x, y, label=rf"$nsteps$ = {int(param['nsteps'])}")
    plt.plot(x_S, y_S, "y--", label="Sun trajectory" if param is params[0] else "")
    plt.plot(x_J, y_J, "g--", label="Jupiter trajectory" if param is params[0] else "")

plt.grid()
plt.xlabel(r"$x\ [\mathrm{m}]$")
plt.ylabel(r"$y\ [\mathrm{m}]$")
plt.legend()
fct.save_figure("traj_b_R_tol.png")
plt.show()

# ENERGY in inertial frame
plt.figure()
plt.title("Energy in Inertial Frame R")

for output, param in zip(outputs_tol, params):
    t, x_prime, y_prime, vx_prime, vy_prime, E_rot, nsteps = fct.read_output_file(output)

    # Convert positions and velocities
    cosO = np.cos(Omega * t)
    sinO = np.sin(Omega * t)
    x = x_prime * cosO - y_prime * sinO
    y = x_prime * sinO + y_prime * cosO

    vx_eff = vx_prime - Omega * y_prime
    vy_eff = vy_prime + Omega * x_prime
    vx = vx_eff * cosO - vy_eff * sinO
    vy = vx_eff * sinO + vy_eff * cosO

    # Compute gravitational potential in inertial frame
    rS = np.sqrt((x + alpha * a)**2 + y**2)
    rJ = np.sqrt((x - beta * a)**2 + y**2)
    U = -6.674e-11 * (m1 / rS + m2 / rJ)

    E_inertial = 0.5 * (vx**2 + vy**2) + U

    plt.plot(t, E_inertial, label=rf"$nsteps$ = {int(param['nsteps'])}")

plt.grid()
plt.xlabel(r"$t\ [\mathrm{s}]$")
plt.ylabel(r"$E\ [\mathrm{J/kg}]$")
plt.legend()
fct.save_figure("E_b_R_nsteps.png")
plt.show()

