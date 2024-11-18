##Overview

This package provides tools for simulating and analyzing 2 dimensional ising model. It includes functionalities for performing simulations using Monte Carlo metropolis algorithm on a cluster of spins. It also includes functions for analyzing reults and visualising the data.

##FILE STRUCTURE

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
│           └── run_simulation.py  # Main simul
├── LICENSE.txt            # License file
├── README.md             # Project documentation
├── pyproject.toml        # Project configuration

##Installation
To run this package you will follow these steps:

1. **Create a Virtual Environment**
   Open your terminal or PowerShell and navigate to the project directory. Create a new virtual environment using:

   
bash
   python -m venv venv_ising
   
2. Activate the Virtual Environment Activate the virtual environment with the following command:

    On Windows: .\venv_ising\Scripts\Activate
    On macOS/Linux: source venv_ising/bin/activate
   
3. Installing the package:
bash
   pip install dist/ising_model-0.1.0-py3-none-any.whl
4. run the following command:
   run-ising
