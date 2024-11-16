import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.optimize import curve_fit

class IsingModel2D:
    def __init__(self, size=20, temperature=2.0, J=1.0):
        self.size = size
        self.T = temperature
        self.J = J
        self.spins = np.random.choice([-1, 1], size=(size, size))
        self.energy = self.calculate_total_energy()
        
    def calculate_local_energy(self, i, j):
        """Calculate energy of a single spin with its neighbors"""
        neighbors_sum = (self.spins[i, (j-1)%self.size] + 
                       self.spins[i, (j+1)%self.size] +
                       self.spins[(i-1)%self.size, j] + 
                       self.spins[(i+1)%self.size, j])
        return -self.J * self.spins[i,j] * neighbors_sum
    
    def calculate_total_energy(self):
        """Calculate total energy of the system"""
        return np.sum([self.calculate_local_energy(i, j) 
                      for i in range(self.size) 
                      for j in range(self.size)]) / 2
    
    def calculate_correlation(self, max_distance=None):
        """Calculate spin-spin correlation function"""
        if max_distance is None:
            max_distance = self.size // 2
            
        correlations = np.zeros(max_distance)
        for r in range(max_distance):
            # Calculate correlation along rows
            corr = np.mean([self.spins[i,j] * self.spins[i,(j+r)%self.size]
                          for i in range(self.size)
                          for j in range(self.size)])
            correlations[r] = corr
        return correlations
    
    def metropolis_step(self):
        """Perform one Metropolis step"""
        i, j = np.random.randint(0, self.size, 2)
        delta_E = -2 * self.calculate_local_energy(i, j)
        
        if delta_E <= 0 or np.random.random() < np.exp(-delta_E/self.T):
            self.spins[i,j] *= -1
            self.energy += delta_E
            
    def run_simulation(self, steps, measure_every=10):
        """Run simulation and collect measurements"""
        energies = []
        magnetizations = []
        
        # Thermalization
        for _ in range(steps//5):  # 20% of steps for thermalization
            for _ in range(self.size * self.size):
                self.metropolis_step()
        
        # Measurement phase
        for step in range(steps):
            for _ in range(self.size * self.size):
                self.metropolis_step()
            
            if step % measure_every == 0:
                energies.append(self.energy/(self.size**2))
                magnetizations.append(np.abs(np.mean(self.spins)))
        
        return np.array(energies), np.array(magnetizations)

def analyze_system(size=20, T_range=None, steps=500):
    """Analyze system properties across temperatures"""
    if T_range is None:
        T_range = np.linspace(1.5, 3.0, 15)
    
    results = {
        'T': T_range,
        'M': np.zeros_like(T_range),  # Magnetization
        'E': np.zeros_like(T_range),  # Energy
        'C': np.zeros_like(T_range),  # Specific Heat
        'X': np.zeros_like(T_range),  # Susceptibility
        'correlation_length': np.zeros_like(T_range)
    }
    
    for i, T in enumerate(tqdm(T_range, desc=f"Analyzing {size}x{size} lattice")):
        model = IsingModel2D(size=size, temperature=T)
        energies, mags = model.run_simulation(steps=steps)
        
        # Calculate observables
        results['M'][i] = np.mean(mags)
        results['E'][i] = np.mean(energies)
        results['C'][i] = (np.var(energies) * size**2) / T**2  # Specific heat
        results['X'][i] = (np.var(mags) * size**2) / T  # Susceptibility
        
        # Calculate correlation length
        correlations = model.calculate_correlation()
        # Fit exponential decay to get correlation length
        try:
            def exp_decay(x, xi, a):
                return a * np.exp(-x/xi)
            popt, _ = curve_fit(exp_decay, np.arange(len(correlations)), 
                              correlations, p0=[1.0, 1.0])
            results['correlation_length'][i] = popt[0]
        except:
            results['correlation_length'][i] = np.nan
            
    return results

def plot_results(results_dict):
    """Plot all observables"""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Magnetization and Energy
    axes[0,0].plot(results_dict['T'], results_dict['M'], 'o-')
    axes[0,0].set_xlabel('Temperature (T)')
    axes[0,0].set_ylabel('Magnetization |M|')
    axes[0,0].axvline(x=2.27, color='r', linestyle='--', alpha=0.5)
    axes[0,0].grid(True)
    
    # Specific Heat
    axes[0,1].plot(results_dict['T'], results_dict['C'], 'o-')
    axes[0,1].set_xlabel('Temperature (T)')
    axes[0,1].set_ylabel('Specific Heat C')
    axes[0,1].axvline(x=2.27, color='r', linestyle='--', alpha=0.5)
    axes[0,1].grid(True)
    
    # Susceptibility
    axes[1,0].plot(results_dict['T'], results_dict['X'], 'o-')
    axes[1,0].set_xlabel('Temperature (T)')
    axes[1,0].set_ylabel('Susceptibility χ')
    axes[1,0].axvline(x=2.27, color='r', linestyle='--', alpha=0.5)
    axes[1,0].grid(True)
    
    # Correlation Length
    axes[1,1].plot(results_dict['T'], results_dict['correlation_length'], 'o-')
    axes[1,1].set_xlabel('Temperature (T)')
    axes[1,1].set_ylabel('Correlation Length ξ')
    axes[1,1].axvline(x=2.27, color='r', linestyle='--', alpha=0.5)
    axes[1,1].grid(True)
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    # Optimized parameters for faster runtime
    lattice_size = 20  # Reduced from 50
    T_range = np.linspace(1.5, 3.0, 15)  # Fewer temperature points
    steps = 500  # Reduced number of steps
    
    print("Running complete Ising model analysis...")
    results = analyze_system(size=lattice_size, T_range=T_range, steps=steps)
    plot_results(results)