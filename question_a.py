import numpy as np
import subprocess
import matplotlib.pyplot as plt
import os
import functions as fct

repertoire = r"C:\Users\Avril\Desktop\Exo2PhysNum-SLAY"  
executable = os.path.join(repertoire, "Exe.exe")

os.chdir(repertoire)

input_filename = "config.a"

params = fct.read_in_file(input_filename)

theta0, thetadot0, mu, B0, m, L, w0 = fct.get_params(params)

A = theta0
B = thetadot0 / w0

#analytical solution functions
def theta_a(t):
    return A * np.cos(w0 * t) + B * np.sin(w0 * t)

def theta_dot(t):
    return w0 * (B * np.cos(w0 * t) - A * np.sin(w0 * t))

nsteps = np.array([10, 20, 50, 100], dtype=int)
nsimul = len(nsteps)
paramstr = "nsteps"
param = nsteps

outputs = []
convergence_list = []

#creates the output files for different nsteps
for i in range(nsimul):
    output_file = f"{paramstr}={param[i]}.out"
    outputs.append(output_file)

    cmd = f'"{executable}" {input_filename} {paramstr}={param[i]:.15g} output={output_file}'

    print(f"\n Running command: {cmd}")

    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    print(" Command Output:", result.stdout)
    print(" Command Error:", result.stderr)

    if os.path.exists(output_file):
        print(f" SUCCESS: Output file '{output_file}' was created!")
    else:
        print(f" ERROR: The output file '{output_file}' was NOT created!")

print("\n🔍 Checking if output files exist before reading...\n")

plt.figure(figsize=(8, 5))

#plotting the simulated data
for i in range(nsimul):
    if os.path.exists(outputs[i]): 
        print(f"Reading file: {outputs[i]}")
        data = np.loadtxt(outputs[i])
        t = data[:, 0]
        theta = data[:, 1]
        plt.plot(t, theta, label=rf"$\theta(t)$ for {param[i]} steps", linewidth=1)

t_a = np.linspace(0, max(t), 1000)
plt.plot(t_a, theta_a(t_a), label=r"$\theta(t)$", color='black', linestyle='--', linewidth=2) #plotting the analytical solution
plt.xlabel("t [s]", fontsize=16)
plt.ylabel(r"$\theta$ [rad]", fontsize=16)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.grid(True)
plt.legend(prop={'size': 14}, loc="upper left")
plt.savefig("questiona1.png", dpi=300)
plt.show()
