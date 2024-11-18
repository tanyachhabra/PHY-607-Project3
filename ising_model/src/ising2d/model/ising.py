"""Implementation of the 2D Ising Model."""

import numpy as np
from ..config.parameters import (
    LATTICE_SIZE, J_COUPLING, MAGNETIC_FIELD, 
    SIMULATION_STEPS, MEASUREMENT_INTERVAL, THERMALIZATION_FRACTION
)

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
