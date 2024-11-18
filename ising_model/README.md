# Ising Model Package

## Overview

This package provides tools for simulating and analyzing the 2-dimensional Ising model. It includes functionalities for performing simulations using the Monte Carlo Metropolis algorithm on a cluster of spins. Additionally, it provides functions for analyzing results and visualizing the data.

## File Structure

```plaintext
ising_model/
├── dist/                     # Distribution files for package
├── src/                      # Source code
│   └── ising2d/             # Main package directory
│       ├── __init__.py      # Package initialization
│       ├── model/           # Core Ising model implementation
│       │   ├── __init__.py
│       │   └── ising.py     # Main Ising model class
│       ├── analysis/        # Analysis tools
│       │   ├── __init__.py
│       │   └── observables.py # Observable calculations
│       ├── visualization/   # Visualization tools
│       │   ├── __init__.py
│       │   └── plotting.py  # Plotting functions
│       ├── config/         # Configuration
│       │   ├── __init__.py
│       │   └── parameters.py # Simulation parameters
│       └── cli/            # Command Line Interface
│           ├── __init__.py
│           └── run_simulation.py  # Main simulation script
├── LICENSE.txt            # License file
├── README.md             # Project documentation
├── pyproject.toml        # Project configuration

Installation
To run this package, follow the steps below:

1 Create a Virtual Environment
   Open your terminal or PowerShell and navigate to the project directory. Create a new virtual environment using:
      python -m venv venv_ising
2 Activate the Virtual Environment
   For windowns:
      .\venv_ising\Scripts\Activate
   On macOS/Linux:
      source venv_ising/bin/activate
3 Install the package
   Run the following command to install the package:
      pip install dist/ising_model-0.1.0-py3-none-any.whl
4 Run the simulation
   Once installed, run the simulation with the following command:
      run-ising

Change of Parameters
To change the parameters, such as magnetic field, lattice size, temperature range, number of steps, and intervals, please go to the ising2d/config/parameters.py file.





