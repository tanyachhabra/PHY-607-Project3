import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.optimize import curve_fit

# Physical Parameters
LATTICE_SIZE = 20
TEMPERATURE_RANGE = np.linspace(1.5, 6, 30)
SIMULATION_STEPS = 500
MEASUREMENT_INTERVAL = 10
THERMALIZATION_FRACTION = 0.2
CRITICAL_TEMP = 2.27
J_COUPLING = 1.0
MAGNETIC_FIELD = 0.0

class IsingModel2D:
    def __init__(self, size=LATTICE_SIZE, temperature=2.0, J=J_COUPLING, h=MAGNETIC_FIELD):
        self.size = size
        self.T = temperature
        self.J = J
        self.h = h
        self.beta = 1.0 / temperature
        self.spins = np.random.choice([-1, 1], size=(size, size))
        self.energy = self.calculate_total_energy()
        
    def calculate_local_energy(self, i, j): #local energy of a spin
        neighbors_sum = (self.spins[i, (j-1)%self.size] + 
                       self.spins[i, (j+1)%self.size] +
                       self.spins[(i-1)%self.size, j] + 
                       self.spins[(i+1)%self.size, j])
        exchange_energy = -self.J * self.spins[i,j] * neighbors_sum
        field_energy = -self.h * self.spins[i,j]
        return exchange_energy + field_energy
    
    def calculate_total_energy(self): # total energy of the spin configration
        #Calculate total energy of the system including magnetic field term
        exchange_energy = np.sum([self.calculate_local_energy(i, j) 
                                for i in range(self.size) 
                                for j in range(self.size)]) / 2
        field_energy = -self.h * np.sum(self.spins)
        return exchange_energy + field_energy
    
    def identify_clusters(self):
        #Identify clusters using Swendsen-Wang algorithm
        # Bond probability
        p = 1 - np.exp(-2 * self.beta * self.J)
        
        # Initialize bonds arrays
        horizontal_bonds = np.random.random((self.size, self.size)) < p
        vertical_bonds = np.random.random((self.size, self.size)) < p
        
        # Only keep bonds between parallel spins
        for i in range(self.size):
            for j in range(self.size):
                if self.spins[i,j] != self.spins[i,(j+1)%self.size]:
                    horizontal_bonds[i,j] = False
                if self.spins[i,j] != self.spins[(i+1)%self.size,j]:
                    vertical_bonds[i,j] = False
        
        # Label clusters using Hoshen-Kopelman algorithm
        labels = np.arange(self.size * self.size).reshape(self.size, self.size)
        parent = np.arange(self.size * self.size)
        
        def find_root(x):
            #Find root label using path compression
            if parent[x] != x:
                parent[x] = find_root(parent[x])
            return parent[x]
        
        def union(x, y):
            #Union of two clusters
            root_x = find_root(x)
            root_y = find_root(y)
            if root_x != root_y:
                parent[root_y] = root_x
        
        # Connect spins through bonds
        for i in range(self.size):
            for j in range(self.size):
                if horizontal_bonds[i,j]:
                    union(labels[i,j], labels[i,(j+1)%self.size])
                if vertical_bonds[i,j]:
                    union(labels[i,j], labels[(i+1)%self.size,j])
        
        # Relabel clusters
        cluster_map = {}
        clusters = np.zeros((self.size, self.size), dtype=int)
        current_label = 0
        
        for i in range(self.size):
            for j in range(self.size):
                root = find_root(labels[i,j])
                if root not in cluster_map:
                    cluster_map[root] = current_label
                    current_label += 1
                clusters[i,j] = cluster_map[root]
                
        return clusters, current_label
    
    def cluster_flip(self):
        #Perform one Swendsen-Wang cluster flip step
        # Identify clusters
        clusters, num_clusters = self.identify_clusters()
        
        # Randomly flip clusters
        flip_decisions = np.random.choice([-1, 1], size=num_clusters)
        
        # Calculate magnetic field contribution
        if self.h != 0:
            for cluster_id in range(num_clusters):
                cluster_mask = (clusters == cluster_id)
                cluster_size = np.sum(cluster_mask)
                delta_E = 2 * self.h * self.spins[cluster_mask].sum() * flip_decisions[cluster_id]
                if delta_E > 0 and np.random.random() > np.exp(-self.beta * delta_E):
                    flip_decisions[cluster_id] = 1  # Don't flip
        
        # Apply flips
        for cluster_id in range(num_clusters):
            cluster_mask = (clusters == cluster_id)
            if flip_decisions[cluster_id] == -1:
                self.spins[cluster_mask] *= -1
        
        # Update energy
        self.energy = self.calculate_total_energy()
    
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
        #Run simulation using cluster updates and collect measurements
        energies = []
        magnetizations = []
        
        # Thermalization
        thermalization_steps = int(steps * THERMALIZATION_FRACTION)
        for _ in range(thermalization_steps):
            self.cluster_flip()
        
        # Measurement phase
        for step in range(steps):
            self.cluster_flip()
            
            if step % measure_every == 0:
                energies.append(self.energy/(self.size**2))
                magnetizations.append(np.mean(self.spins))
        
        return np.array(energies), np.array(magnetizations)



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

if __name__ == "__main__":
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