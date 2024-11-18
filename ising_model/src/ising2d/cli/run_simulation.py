"""Main script to run the Ising model simulation."""

import numpy as np
from ising2d.model import IsingModel2D
from ising2d.analysis import analyze_system
from ising2d.analysis import track_thermalization
from ising2d.visualization import (
    plot_magnetization,
    plot_specific_heat,
    plot_susceptibility,
    plot_correlation_length,
    plot_thermalization,
    plot_spin_configuration
)
from ising2d.config.parameters import (
    LATTICE_SIZE,
    TEMPERATURE_RANGE,
    MAGNETIC_FIELD,
    SIMULATION_STEPS
)

def main():
    """Main function to run the Ising model simulation."""
    print(f"Running Ising model analysis for {LATTICE_SIZE}x{LATTICE_SIZE} lattice...")
    print(f"Temperature range: {TEMPERATURE_RANGE[0]:.2f} to {TEMPERATURE_RANGE[-1]:.2f}")
    print(f"External magnetic field: {MAGNETIC_FIELD}")
    print(f"Number of steps: {SIMULATION_STEPS}")
    
    results = analyze_system()
    
    
    # Plots for Observables 
    plot_magnetization(results)
    plot_specific_heat(results)
    plot_susceptibility(results)
    plot_correlation_length(results)
    
    
    model = IsingModel2D(temperature=2.0, h=MAGNETIC_FIELD)
    model.run_simulation(steps=SIMULATION_STEPS)
    max_thermalization_steps = 1000
    energies, magnetizations = track_thermalization(model, max_steps=max_thermalization_steps)

    #Burn-in period plots
    plot_thermalization(energies, magnetizations)
    plot_spin_configuration(model)

if __name__ == "__main__":
    main()
