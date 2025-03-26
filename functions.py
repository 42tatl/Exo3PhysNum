import numpy as np
import subprocess
import matplotlib.pyplot as plt
import os

def read_in_file(filename): 
    '''Reads in a file and returns the data as a list of floats'''
    variables = {}
    with open(filename, "r") as file:
        for line in file:
            line = line.strip()
            if line and not line.startswith("#"):  
                key, value = line.split("=")
                key = key.strip()
                value = value.split("//")[0].strip()  
                
                try:
                    if "." in value:
                        variables[key] = float(value)
                    else:
                        variables[key] = int(value)
                except ValueError:
                    print(f"Warning: Could not convert '{value}' to number. Check {filename}.")
    return variables

def get_params(params):
    '''Extracts the parameters from the dictionary'''
    tFin = params.get("tFin", 0.0)
    m1 = params.get("m1", 0.0)
    m2 = params.get("m2", 0.0)
    x0 = params.get("x0", 0.0)
    v0x = params.get("v0x", 0.0)
    v0y = params.get("v0y", 0.0)
    a = params.get("a", 0.0)
    tol = params.get("tol", 0.0)
    return tFin, m1, m2, x0, v0x, v0y, a, tol

def run_simulation(executable, input_filename, output_template, **params):
    '''Runs the simulation with the given parameters'''

    output_filename = output_template.format(**params)
    
    param_str = " ".join(f"{key}={value:.15g}" for key, value in params.items())
    cmd = f'"{executable}" {input_filename} {param_str} output={output_filename}'

    print(f"\n Running command: {cmd}")

    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    print(" Command Output:", result.stdout)
    print(" Command Error:", result.stderr)

    if os.path.exists(output_filename):
        print(f" SUCCESS: Output file '{output_filename}' was created!")
    else:
        print(f" ERROR: The output file '{output_filename}' was NOT created!")

    return output_filename, result

def run_param_sweep(executable, input_filename, param_name, values, fixed_params):
    outputs = []
    for val in values:
        params = fixed_params.copy()
        params[param_name] = val
        output_template = f"output_{param_name}_{{{param_name}}}.out"
        outname, result = run_simulation(executable, input_filename, output_template, **params)
        outputs.append(outname)
<<<<<<< HEAD

    return outputs
=======
    return outputs
>>>>>>> 3a4e8e688516d8466f6f22cf8394ff9f82dbfc7f
