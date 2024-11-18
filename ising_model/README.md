<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Ising Model Package Documentation</title>
</head>
<body>
    <h1>Overview</h1>
    <p>
        This package provides tools for simulating and analyzing the 2D Ising model. It includes functionalities for performing simulations using the Monte Carlo Metropolis algorithm on a cluster of spins. Additionally, it provides functions for analyzing results and visualizing the data.
    </p>

    <h2>File Structure</h2>
    <pre>
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
    </pre>

    <h2>Installation</h2>
    <p>To run this package, follow the steps below:</p>

    <ol>
        <li>
            <strong>Create a Virtual Environment</strong><br>
            Open your terminal or PowerShell and navigate to the project directory. Create a new virtual environment using:
            <pre><code>python -m venv venv_ising</code></pre>
        </li>
        <li>
            <strong>Activate the Virtual Environment</strong><br>
            Activate the virtual environment with the following command:
            <ul>
                <li>On Windows: <code>.\\venv_ising\\Scripts\\Activate</code></li>
                <li>On macOS/Linux: <code>source venv_ising/bin/activate</code></li>
            </ul>
        </li>
        <li>
            <strong>Installing the package</strong><br>
            Run the following command to install the package:
            <pre><code>pip install dist/ising_model-0.1.0-py3-none-any.whl</code></pre>
        </li>
        <li>
            <strong>Run the simulation</strong><br>
            Once installed, run the simulation with the following command:
            <pre><code>run-ising</code></pre>
        </li>
    </ol>
</body>
</html>
