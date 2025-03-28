# Project Overview

## Python Scripts

### `functions.py`
- **Purpose**: Contains reusable utility functions for numerical computations and data processing.
- **Details**:
  - Includes helper functions for mathematical operations, data formatting, and file handling.
  - Used as a shared library across other Python scripts to avoid code duplication.

### `question_a.py`
- **Purpose**: Solves "Question A" of the project.
- **Details**:
  - Implements a specific numerical simulation or data analysis task.
  - Uses the functions from `functions.py` to perform calculations.
  - Outputs results to a file or visualizes them using plots.

### `question_b_nsteps.py`
- **Purpose**: Solves "Question B" with a focus on varying the number of steps (`nsteps`).
- **Details**:
  - Explores how the number of time steps affects the accuracy and performance of the simulation.
  - Outputs results for comparison with other configurations.

### `question_b_tol.py`
- **Purpose**: Solves "Question B" with a focus on varying the tolerance (`tol`).
- **Details**:
  - Investigates the impact of different tolerance values on the simulation's accuracy and stability.
  - Outputs results for analysis and comparison.

---

## C++ Files

### `Ex3.cpp`
- **Purpose**: Implements the main C++ program for numerical simulations.
- **Details**:
  - Contains the core logic for solving a specific problem (e.g., "Exercise 3").
  - Reads input parameters from configuration files and performs computations.
  - Outputs results to files for further analysis or visualization.

### `ConfigFile.h` and `ConfigFile.hpp`
- **Purpose**: Provide functionality for handling configuration files in the C++ program.
- **Details**:
  - Define classes and methods for reading, parsing, and managing configuration settings.
  - Allow the program to dynamically adjust parameters like time steps, tolerances, and initial conditions without modifying the source code.

---

## Configuration Files

### `configa.in`
- **Purpose**: Provides input parameters for the C++ program.
- **Details**:
  - Defines initial conditions, physical constants, and simulation settings for a specific scenario.
  - Example parameters: `tFin` (final time), `dt` (time step), `m1` and `m2` (masses), etc.

### `configb.txt`
- **Purpose**: Another configuration file for the C++ program.
- **Details**:
  - Similar to `configa.in`, but defines a different set of parameters for another experiment or simulation.

---

## Output Files

### `output_tol_*.out`
- **Purpose**: Store the results of simulations or computations for different tolerance values.
- **Details**:
  - Each file corresponds to a specific tolerance setting (e.g., `0.001`, `0.01`, `1e-06`).
  - Contains numerical data such as positions, velocities, energies, or errors over time.

---

## How to Run

1. **Python Scripts**:
   - Run the Python scripts using the following command:
     ```sh
     python <script_name>.py
     ```
   - Example:
     ```sh
     python question_a.py
     ```

2. **C++ Program**:
   - Compile the C++ program using a C++ compiler (e.g., `g++`):
     ```sh
     g++ Ex3.cpp -o Exe
     ```
   - Run the compiled executable with a configuration file:
     ```sh
     ./Exe configa.in
     ```

---

## Notes

- Ensure all required dependencies are installed before running the scripts or programs.
- Modify the configuration files (`configa.in`, [configb.txt](http://_vscodecontentref_/2)) to adjust simulation parameters as needed.
- Use the [functions.py](http://_vscodecontentref_/3) script to add or modify utility functions for custom computations.

---