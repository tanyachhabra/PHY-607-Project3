#importing libraries 
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.optimize import curve_fit

# Define all parameters upfront
LATTICE_SIZE = 20
TEMPERATURE_RANGE = np.linspace(1.5, 6,30 )
SIMULATION_STEPS = 500
MEASUREMENT_INTERVAL = 10
THERMALIZATION_FRACTION = 0.2
CRITICAL_TEMP = 2.27
J_COUPLING = 1.0
MAGNETIC_FIELD = 0.1  # External magnetic field strength

class IsingModel2D:
    def __init__(self, size=LATTICE_SIZE, temperature=2.0, J=J_COUPLING, h=MAGNETIC_FIELD):
        self.size = size
        self.T = temperature
        self.J = J
        self.h = h  # External magnetic field
        self.spins = np.random.choice([-1, 1], size=(size, size))
        self.energy = self.calculate_total_energy()
        
    def calculate_local_energy(self, i, j):
        """Calculate energy of a single spin with its neighbors and magnetic field"""
        # Exchange interaction term
        neighbors_sum = (self.spins[i, (j-1)%self.size] + 
                       self.spins[i, (j+1)%self.size] +
                       self.spins[(i-1)%self.size, j] + 
                       self.spins[(i+1)%self.size, j])
        exchange_energy = -self.J * self.spins[i,j] * neighbors_sum
        
        # Magnetic field term
        field_energy = -self.h * self.spins[i,j]
        
        return exchange_energy + field_energy
    
    def calculate_total_energy(self):
        """Calculate total energy of the system including magnetic field term"""
        exchange_energy = np.sum([self.calculate_local_energy(i, j) 
                                for i in range(self.size) 
                                for j in range(self.size)]) / 2
        
        # Add total magnetic field contribution
        field_energy = -self.h * np.sum(self.spins)
        
        return exchange_energy + field_energy
    
    def metropolis_step(self):
        """Perform one Metropolis step"""
        i, j = np.random.randint(0, self.size, 2)
        # Delta E includes both exchange and field terms
        delta_E = -2 * self.calculate_local_energy(i, j)
        
        if delta_E <= 0 or np.random.random() < np.exp(-delta_E/self.T):
            self.spins[i,j] *= -1
            self.energy += delta_E
    
    # Rest of the methods remain the same
    def calculate_correlation(self, max_distance=None):
        if max_distance is None:
            max_distance = self.size // 2
            
        correlations = np.zeros(max_distance)
        for r in range(max_distance):
            corr = np.mean([self.spins[i,j] * self.spins[i,(j+r)%self.size]
                          for i in range(self.size)
                          for j in range(self.size)])
            correlations[r] = corr
        return correlations
            
    def run_simulation(self, steps=SIMULATION_STEPS, measure_every=MEASUREMENT_INTERVAL):
        """Run simulation and collect measurements"""
        energies = []
        magnetizations = []
        
        # Thermalization
        thermalization_steps = int(steps * THERMALIZATION_FRACTION)
        for _ in range(thermalization_steps):
            for _ in range(self.size * self.size):
                self.metropolis_step()
        
        # Measurement phase
        for step in range(steps):
            for _ in range(self.size * self.size):
                self.metropolis_step()
            
            if step % measure_every == 0:
                energies.append(self.energy/(self.size**2))
                magnetizations.append(np.mean(self.spins))  # Note: removed abs() to show field direction
        
        return np.array(energies), np.array(magnetizations)

def analyze_system(size=LATTICE_SIZE, T_range=TEMPERATURE_RANGE, h=MAGNETIC_FIELD, steps=SIMULATION_STEPS):
    """Analyze system properties across temperatures with fixed magnetic field"""
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
        
        results['M'][i] = np.mean(mags)  # Note: can be negative now
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

def plot_magnetization(results_dict):
    """Plot magnetization"""
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


def plot_energy(results_dict):
    plt.figure(figsize=(8, 6))
    plt.plot(results_dict['T'], results_dict['E'], 'o-')
    plt.axvline(x=CRITICAL_TEMP, color='r', linestyle='--', alpha=0.5,
                label=f'Critical Temperature (T={CRITICAL_TEMP})')
    plt.xlabel('Temperature (T)')
    plt.ylabel('Energy per spin (E)')
    plt.grid(True)
    plt.legend()
    plt.title(f'Energy vs Temperature\n{LATTICE_SIZE}x{LATTICE_SIZE} lattice, h={results_dict["h"]}')
    plt.show()


def plot_magnetization(results_dict):
    """Plot magnetization"""
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
    """Plot specific heat"""
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
    """Plot susceptibility"""
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
    """Plot correlation length"""
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
    """Plot current spin configuration"""
    plt.figure(figsize=(6, 6))
    plt.imshow(model.spins, cmap='coolwarm')
    plt.colorbar(label='Spin')
    plt.title(f'Spin Configuration\nT={model.T:.2f}, h={model.h:.3f}')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()

if __name__ == "__main__":
    print(f"Running Ising model analysis for {LATTICE_SIZE}x{LATTICE_SIZE} lattice...")
    print(f"Temperature range: {TEMPERATURE_RANGE[0]:.2f} to {TEMPERATURE_RANGE[-1]:.2f}")
    print(f"External magnetic field: {MAGNETIC_FIELD}")
    print(f"Number of steps: {SIMULATION_STEPS}")
    
    results = analyze_system()
    
    
    plot_magnetization(results)
    plot_specific_heat(results)
    plot_susceptibility(results)
    plot_correlation_length(results)
    
    
    model = IsingModel2D(temperature=2.0, h=MAGNETIC_FIELD)
    model.run_simulation(steps=SIMULATION_STEPS)
    plot_spin_configuration(model)

