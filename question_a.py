import numpy as np
import subprocess
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.lines import Line2D
import os
import functions as fct

executable = './Exe3'  # Remove .exe for Mac
repertoire = r"/Users/lilimouelle/Desktop/PHYSNUM/Exo3PhysNum"  # Modify for correct directory
os.chdir(repertoire)

input_filename = "configa.in"
output_template = "output_{paramstr}_{value}.out"

params = fct.read_in_file(input_filename)

tFin, m1, m2, x0, v0x, v0y, a, tol = fct.get_params(params)

#outputs_nsteps = fct.run_param_sweep(executable, input_filename, "nsteps", [8000,10000,90000], {"adapt": 0, "tol": 30})
#outputs_tol = fct.run_param_sweep(executable, input_filename, "tol", [1,0.1,0.01], {"adapt": True, "nsteps": 1000})

outputs_nsteps = fct.run_param_sweep(executable, input_filename, "nsteps", [90000], {"adapt": 0, "tol": 30})
outputs_tol = fct.run_param_sweep(executable, input_filename, "tol", [0.01], {"adapt": True, "nsteps": 1000})


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

'''
#TRAJECTORY AND ENERGY PLOTS
#Fixed time step
for output in outputs_nsteps:
    t, x, y, vx, vy, energy = read_output_file(output)
    label = output.split("_")[-1].replace(".out", "")
    n_label = f"nsteps = {label}"   
    plot_trajectory(x, y, label=n_label)
    ax = plt.gca()  # get current axis
    soleil = Circle((0, 0), radius=0.05e12, color='yellow', zorder=10)
    pos_init = Circle((x[0], y[0]), radius=0.05e12, color='red', zorder=10)
    pos_final = Circle((x[-1], y[-1]), radius=0.05e12, color='blue', zorder=10)
    ax.add_patch(soleil)
    ax.add_patch(pos_init)
    ax.add_patch(pos_final)
    legend_elements = [
    Line2D([0], [0], marker='o', color='yellow', label='Sun', markersize=10, linestyle='None'),
    Line2D([0], [0], marker='o', color='red', label='Initial position', markersize=10, linestyle='None'),
    Line2D([0], [0], marker='o', color='blue', label='Final position', markersize=10, linestyle='None'),
    Line2D([], [], label=n_label)
]

    plt.legend(handles=legend_elements)
    plt.axis("equal")
    plt.savefig(f"traja_{label}.pdf")
    plt.show()


for output in outputs_nsteps:
    t, x, y, vx, vy, energy = read_output_file(output)
    label = output.split("_")[-1].replace(".out", "")
    n_label = f"nsteps = {label}"  
    plot_energy(t, energy, label=n_label)
    plt.legend()
    plt.savefig(f"energiea_{label}.pdf")
    plt.show()



#Adaptive time step
for output in outputs_tol:
    t, x, y, vx, vy, energy = read_output_file(output)
    label = output.split("_")[-1].replace(".out", "")
    eps_label = f"epsilon = {label}"   
    plot_trajectory(x, y, label=eps_label)
    ax = plt.gca()  # get current axis
    soleil = Circle((0, 0), radius=0.05e12, color='yellow', zorder=10)
    pos_init = Circle((x[0], y[0]), radius=0.05e12, color='red', zorder=10)
    pos_final = Circle((x[-1], y[-1]), radius=0.05e12, color='blue', zorder=10)
    ax.add_patch(soleil)
    ax.add_patch(pos_init)
    ax.add_patch(pos_final)
    legend_elements = [
    Line2D([0], [0], marker='o', color='yellow', label='Sun', markersize=10, linestyle='None'),
    Line2D([0], [0], marker='o', color='red', label='Initial position', markersize=10, linestyle='None'),
    Line2D([0], [0], marker='o', color='blue', label='Final position', markersize=10, linestyle='None'),
    Line2D([], [], label=eps_label)
]

    plt.legend(handles=legend_elements)
    plt.axis("equal")
    plt.savefig(f"traja_adapt_{label}.pdf")
    plt.show()

for output in outputs_tol:
    t, x, y, vx, vy, energy = read_output_file(output)
    label = output.split("_")[-1].replace(".out", "")
    eps_label = f"epsilon = {label}"   
    plot_energy(t, energy, label=eps_label)
    plt.legend()
    plt.savefig(f"energiea_adapt_{label}.pdf")
    plt.show()
'''




#COMPARISON OF NUMERICAL AND THEORETICAL VALUES FOR r_min,r_max,v_min,v_max
G = 6.67430e-11
E = 0.5*(v0x**2+v0y**2)-(G*m1)/x0
L = x0*v0y
e = np.sqrt(1 + 2*E*L**2/(G**2*m1**2))
a_orb = -G*m1/(2*E)
r_min_th = a_orb*(1-e)
r_max_th = a_orb*(1+e)
v_min_th = np.sqrt(2*(E+G*m1/r_max_th))
v_max_th = np.sqrt(2*(E+G*m1/r_min_th))

def compute_min_max(x,y,vx,vy):
    r = np.sqrt(x**2+y**2)
    v = np.sqrt(vx**2+vy**2)
    return np.min(r), np.max(r), np.min(v), np.max(v)

for outputs in outputs_nsteps:
    t, x, y, vx, vy, energy = read_output_file(outputs)
    r_min, r_max, v_min, v_max = compute_min_max(x,y,vx,vy)
    label = outputs.split("_")[-1].replace(".out", "")
    print(f"\n--- nsteps = {label} ---")
    print(f"  r_min num = {r_min:.3e} m   | r_min th = {r_min_th:.3e} m")
    print(f"  r_max num = {r_max:.3e} m   | r_max th = {r_max_th:.3e} m")
    print(f"  v_min num = {v_min:.2f} m/s | v_min th = {v_min_th:.2f} m/s")
    print(f"  v_max num = {v_max:.2f} m/s | v_max th = {v_max_th:.2f} m/s")

for outputs in outputs_tol:
    t, x, y, vx, vy, energy = read_output_file(outputs)
    r_min, r_max, v_min, v_max = compute_min_max(x,y,vx,vy)
    label = outputs.split("_")[-1].replace(".out", "")
    print(f"\n--- epsilon = {label} ---")
    print(f"  r_min num = {r_min:.3e} m   | r_min th = {r_min_th:.3e} m")
    print(f"  r_max num = {r_max:.3e} m   | r_max th = {r_max_th:.3e} m")
    print(f"  v_min num = {v_min:.2f} m/s | v_min th = {v_min_th:.2f} m/s")
    print(f"  v_max num = {v_max:.2f} m/s | v_max th = {v_max_th:.2f} m/s")



#VARIATION OF TIME STEP WITH ADAPTIVE TIME STEP
for output in outputs_tol:
    t, x, y, vx, vy, energy = read_output_file(output)
    label = output.split("_")[-1].replace(".out", "")
    eps_label = f"epsilon = {label}"

    dt = np.diff(t)
    plt.plot(t[:-1], dt, label=eps_label)
    plt.xlabel("t [s]")
    plt.ylabel("Δt [s]")
    plt.legend()
    plt.grid(True)
    plt.savefig(f"dt_adapt_{label}.pdf")
    plt.show()
