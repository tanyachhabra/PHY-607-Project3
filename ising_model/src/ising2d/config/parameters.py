"""Configuration parameters for the Ising model simulation."""

import numpy as np

# Physical Parameters
LATTICE_SIZE = 20
TEMPERATURE_RANGE = np.linspace(1.5, 6, 30)
SIMULATION_STEPS = 500
MEASUREMENT_INTERVAL = 10
THERMALIZATION_FRACTION = 0.2
CRITICAL_TEMP = 2.27
J_COUPLING = 1.0
MAGNETIC_FIELD = 0.0
