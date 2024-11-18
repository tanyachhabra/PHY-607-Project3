"""Visualization functions for the Ising model simulation."""

import numpy as np
import matplotlib.pyplot as plt
from ..config.parameters import LATTICE_SIZE, CRITICAL_TEMP

def plot_thermalization(energies, magnetizations):
    # Plot for burn-in 
    steps = np.arange(len(energies))
    
    plt.figure(figsize=(12, 6))
    
    # Energy plot
    plt.subplot(1, 2, 1)
    plt.plot(steps, energies, label='Energy per Spin', color='blue')
    plt.xlabel('Steps')
    plt.ylabel('Energy per Spin')
    plt.title('Thermalization: Energy vs Steps')
    plt.grid(True)
    plt.legend()
    
    # Magnetization plot
    plt.subplot(1, 2, 2)
    plt.plot(steps, magnetizations, label='Magnetization', color='red')
    plt.xlabel('Steps')
    plt.ylabel('Magnetization')
    plt.title('Thermalization: Magnetization vs Steps')
    plt.grid(True)
    plt.legend()
    
    plt.tight_layout()
    plt.show()


def plot_magnetization(results_dict):
    # Plot for magnetization 
    plt.figure(figsize=(8, 6))
    plt.plot(results_dict['T'], results_dict['M'], 'o-')
    plt.axvline(x=CRITICAL_TEMP, color='r', linestyle='--', alpha=0.5, 
                label=f'Critical Temperature (T={CRITICAL_TEMP})')
    plt.axhline(y=0, color='k', linestyle=':', alpha=0.5)
    plt.xlabel('Temperature (T)')
    plt.ylabel('Magnetization (M)')
    plt.grid(True)
    plt.legend()
    plt.title(f'Magnetization vs Temperature\n{LATTICE_SIZE}x{LATTICE_SIZE} lattice, h={results_dict["h"]}')
    plt.show()


def plot_specific_heat(results_dict):
    #Plot specific heat
    plt.figure(figsize=(8, 6))
    plt.plot(results_dict['T'], results_dict['C'], 'o-')
    plt.axvline(x=CRITICAL_TEMP, color='r', linestyle='--', alpha=0.5,
                label=f'Critical Temperature (T={CRITICAL_TEMP})')
    plt.xlabel('Temperature (T)')
    plt.ylabel('Specific Heat (C)')
    plt.grid(True)
    plt.legend()
    plt.title(f'Specific Heat vs Temperature\n{LATTICE_SIZE}x{LATTICE_SIZE} lattice, h={results_dict["h"]}')
    plt.show()

def plot_susceptibility(results_dict):
    #Plot susceptibility
    plt.figure(figsize=(8, 6))
    plt.plot(results_dict['T'], results_dict['X'], 'o-')
    plt.axvline(x=CRITICAL_TEMP, color='r', linestyle='--', alpha=0.5,
                label=f'Critical Temperature (T={CRITICAL_TEMP})')
    plt.xlabel('Temperature (T)')
    plt.ylabel('Susceptibility (χ)')
    plt.grid(True)
    plt.legend()
    plt.title(f'Susceptibility vs Temperature\n{LATTICE_SIZE}x{LATTICE_SIZE} lattice, h={results_dict["h"]}')
    plt.show()

def plot_correlation_length(results_dict):
    #Plot correlation length
    plt.figure(figsize=(8, 6))
    plt.plot(results_dict['T'], results_dict['correlation_length'], 'o-')
    plt.axvline(x=CRITICAL_TEMP, color='r', linestyle='--', alpha=0.5,
                label=f'Critical Temperature (T={CRITICAL_TEMP})')
    plt.xlabel('Temperature (T)')
    plt.ylabel('Correlation Length (ξ)')
    plt.grid(True)
    plt.legend()
    plt.title(f'Correlation Length vs Temperature\n{LATTICE_SIZE}x{LATTICE_SIZE} lattice, h={results_dict["h"]}')
    plt.show()
def plot_spin_configuration(model):
    # plot for spin configuration 
    plt.figure(figsize=(6, 6))
    plt.imshow(model.spins, cmap='coolwarm')
    plt.colorbar(label='Spin')
    plt.title(f'Spin Configuration\nT={model.T:.2f}, h={model.h:.3f}')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()
