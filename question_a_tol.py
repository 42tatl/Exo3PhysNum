import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.lines import Line2D
import os
import functions as fct
from matplotlib.ticker import ScalarFormatter

executable = './Exe'  # Remove .exe for Mac
repertoire = r"/Users/lilimouelle/Desktop/PHYSNUM/Exo3PhysNum"  # Modify for correct directory
os.chdir(repertoire)

input_filename = "configa.in"
output_template = "output_{paramstr}_{value}.out"

params = fct.read_in_file(input_filename)

tFin, m1, m2, x0, v0x, v0y, a, tol = fct.get_params(params)

outputs_tol, parameters = fct.run_param_sweep(executable, input_filename, "tol", [1e-3,100,1e7], {"adapt": True, "nsteps": 1000})

def plot_trajectory(x, y, label=""):
    plt.plot(x, y, label=label)
    plt.xlabel(r'$x$ [m]',fontsize=14)
    plt.ylabel(r'$y$ [m]',fontsize=14)
    plt.xticks(fontsize=12)  
    plt.yticks(fontsize=12) 
    #plt.axis("equal")
    plt.grid(True)
    plt.tight_layout()

def plot_energy(t, energy, label=""):
    plt.plot(t, energy, label=label)
    plt.xlabel(r'$t$ [s]',fontsize=14)
    plt.xticks(fontsize=12)  
    plt.yticks(fontsize=12)
    plt.ylabel(r'$\frac{E}{m}$ [J $\cdot$ kg$^{-1}$]',fontsize=14)
    plt.grid(True)
    plt.tight_layout()


#TRAJECTORY
for output, param in zip(outputs_tol, parameters):
    t, x, y, vx, vy, energy, nsteps = fct.read_output_file(output)
    label = output.split("_")[-1].replace(".out", "")
    nsteps = int(nsteps[-1])
    eps_label = f"epsilon = {label}, nsteps = {nsteps}"
    plot_trajectory(x, y, label=eps_label)
    ax = plt.gca()  
    soleil = Circle((0, 0), radius=0.05e12, color='yellow', zorder=10)
    pos_init = Circle((x[0], y[0]), radius=0.05e12, color='red', zorder=10)
    pos_final = Circle((x[-1], y[-1]), radius=0.05e12, color='blue', zorder=10)
    print("Initial pos =", x[0], y[0])
    print("Final pos   =", x[-1], y[-1])
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
    fct.save_figure(f"traja_adapt_{label}.pdf")
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

for output, param in zip(outputs_tol, parameters):
    t, x, y, vx, vy, energy, nsteps = fct.read_output_file(output)
    r_min, r_max, v_min, v_max = compute_min_max(x,y,vx,vy)
    label = output.split("_")[-1].replace(".out", "")
    print(f"\n--- epsilon = {label} ---")
    print(f"  r_min num = {r_min:.3e} m   | r_min th = {r_min_th:.3e} m")
    print(f"  r_max num = {r_max:.3e} m   | r_max th = {r_max_th:.3e} m")
    print(f"  v_min num = {v_min:.2f} m/s | v_min th = {v_min_th:.2f} m/s")
    print(f"  v_max num = {v_max:.2f} m/s | v_max th = {v_max_th:.2f} m/s")



#VARIATION OF TIME STEP
for output, param in zip(outputs_tol, parameters):
    t, x, y, vx, vy, energy, nsteps = fct.read_output_file(output)
    label = output.split("_")[-1].replace(".out", "")
    eps_label = f"epsilon = {label}"
    dt = np.diff(t)
    plt.plot(t[:-1], dt, label=eps_label)
    #plt.yscale("log")
    plt.xlabel(r'$t$ [s]',fontsize=14)
    plt.ylabel(r'$\Delta t$ [s]',fontsize=14)
    plt.xticks(fontsize=12)  
    plt.yticks(fontsize=12) 
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    fct.save_figure(f"dt_adapt_{label}.pdf")
    plt.show()
'''


'''
#CONVERGENCE OF FINAL POSITION
outputs_tol, parameters = fct.run_param_sweep(executable, input_filename, "tol", np.logspace(-2,-1,num=10), {"adapt": True, "nsteps": 1000})
nsteps_tot = []
xfin = []
yfin = []
pos_final = []
delta_ts_four = []

for output, param in zip(outputs_tol, parameters):
    t, x, y, vx, vy, energy,nsteps = fct.read_output_file(output)
    nsteps_tot.append(nsteps[-1])
    xfin.append(x[-1])
    yfin.append(y[-1])
    pos_final.append(np.sqrt(x[-1]**2+y[-1]**2))
    delta_ts_four.append((tFin/nsteps[-1])**4)


plt.plot(delta_ts_four, pos_final, marker='o', linestyle='-', color='b', label="Adaptive")
plt.xlabel(r'$(\Delta t)^4$ [s$^4$]',fontsize=14)
plt.ylabel(r'$r_{fin}$ [m]',fontsize=14)
plt.xticks(fontsize=12)  
plt.yticks(fontsize=12)
plt.grid(True)
plt.legend()
plt.tight_layout()
fct.save_figure(f"conv_pos_tol_1.pdf")
plt.show()


plt.plot(delta_ts_four, xfin, marker='o', linestyle='-', color='b', label="Adaptive")
plt.xlabel(r'$(\Delta t)^4$ [s$^4$]',fontsize=14)
plt.ylabel(r'$x_{fin}$ [m]',fontsize=14)
plt.xticks(fontsize=12)  
plt.yticks(fontsize=12)
plt.grid(True)
plt.legend()
plt.tight_layout()
fct.save_figure(f"conv_pos_tol_2.pdf")
plt.show()


plt.plot(delta_ts_four, yfin, marker='o', linestyle='-', color='b', label="Adaptive")
plt.xlabel(r'$(\Delta t)^4$ [s$^4$]',fontsize=14)
plt.ylabel(r'$y_{fin}$ [m]',fontsize=14)
plt.xticks(fontsize=12)  
plt.yticks(fontsize=12)
plt.grid(True)
plt.legend()
plt.tight_layout()
fct.save_figure(f"conv_pos_tol_3.pdf")
plt.show()



#CONVERGENCE ENERGY
outputs_tol, parameters = fct.run_param_sweep(executable, input_filename, "tol", np.logspace(2,4,num=10), {"adapt": True, "nsteps": 1000})
nsteps_tot_e = []
diff_energy = []
Efinal = []
dts_four = []


for output, param in zip(outputs_tol, parameters):
    t, x, y, vx, vy, energy,nsteps = fct.read_output_file(output)
    nsteps_tot_e.append(nsteps[-1])
    Efinal.append(energy[-1])
    dts_four.append((tFin/nsteps[-1])**4)
    Emin = np.min(energy)
    Emax = np.max(energy)
    diff_energy.append(np.abs(Emax-Emin))

ref_n = np.array(nsteps_tot_e, dtype=float)
ref_y = ref_n**(-4)
ref_curve = 1e2 * (ref_n / ref_n[0])**(-4)


plt.loglog(nsteps_tot_e, diff_energy, marker='o', linestyle='-', color='b', label="Adaptive") 
plt.plot(ref_n, ref_curve, '--', color='r', label=r"$\propto \Delta n^{-4}$")
plt.xlabel(r'$n_{steps}$',fontsize=14)
plt.ylabel(r'$\frac{\left| E_{max} - E_{min} \right|}{m}$ [J $\cdot$ kg$^{-1}$]',fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid(True)
plt.legend()
plt.tight_layout()
fct.save_figure(f"conv_energy_tol_1.pdf")
plt.show()

plt.plot(dts_four, Efinal, marker='o', linestyle='-', color='b', label="Adaptive") 
plt.xlabel(r'$(\frac{1}{n_{steps}})^4$',fontsize=14)
plt.ylabel(r'$\frac{E_{fin}}{m}$ [J $\cdot$ kg$^{-1}$]',fontsize=14)
plt.xticks(fontsize=12)  
plt.yticks(fontsize=12)
plt.grid(True)
plt.legend()
plt.tight_layout()
fct.save_figure(f"conv_energy_tol_2.pdf")
plt.show()
'''

