"""Analysis functions for the Ising model simulation."""

import numpy as np
from tqdm import tqdm
from scipy.optimize import curve_fit

from ..model.ising import IsingModel2D
from ..config.parameters import (
    LATTICE_SIZE, TEMPERATURE_RANGE, MAGNETIC_FIELD, 
    SIMULATION_STEPS
)

def analyze_system(size=LATTICE_SIZE, T_range=TEMPERATURE_RANGE, h=MAGNETIC_FIELD, steps=SIMULATION_STEPS):
    #Analyze system properties across temperatures with fixed magnetic field
    results = {
        'T': T_range,
        'h': h,
        'M': np.zeros_like(T_range),
        'E': np.zeros_like(T_range),
        'C': np.zeros_like(T_range),
        'X': np.zeros_like(T_range),
        'correlation_length': np.zeros_like(T_range)
    }
    
    for i, T in enumerate(tqdm(T_range, desc=f"Analyzing {size}x{size} lattice with h={h}")):
        model = IsingModel2D(size=size, temperature=T, h=h)
        energies, mags = model.run_simulation(steps=steps)
        
        results['M'][i] = np.mean(mags)
        results['E'][i] = np.mean(energies)
        results['C'][i] = (np.var(energies) * size**2) / T**2
        results['X'][i] = (np.var(mags) * size**2) / T
        
        correlations = model.calculate_correlation()
        try:
            def exp_decay(x, xi, a):
                return a * np.exp(-x/xi)
            popt, _ = curve_fit(exp_decay, np.arange(len(correlations)), 
                              correlations, p0=[1.0, 1.0])
            results['correlation_length'][i] = popt[0]
        except:
            results['correlation_length'][i] = np.nan
            
    return results

def track_thermalization(model, max_steps=1000):
    #Track energy and magnetization during thermalization
    energies = []
    magnetizations = []
    
    for step in range(max_steps):
        model.cluster_flip()
        energies.append(model.energy / (model.size ** 2))  # Normalize energy
        magnetizations.append(np.mean(model.spins))       # Magnetization
    
    return np.array(energies), np.array(magnetizations)
