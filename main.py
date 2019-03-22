# Script for tuning and running simulation

import numpy as np

import verlet_int as vi

# Constants for book keeping
# m = 6.6e-26  # kg
# sigma = 3.405  # Å
# k_B = 1.38e-23  # J/K
# epsilon = 119.8 * k_B  # J

# Simulation parameters
time_step = 1e-3
t_max = 8e0
# temp = 1.0
# rho = 0.8

if __name__ == "__main__":
    # Gas
    # vi.run_sim(time_step, t_max, temp=3.0, rho=0.3, n=5, plot=False, write_to_file=True)
    # Liquid
    # vi.run_sim(time_step, t_max, temp=1.0, rho=0.8, n=5, plot=False, write_to_file=True)
    # Solid
    vi.run_sim(time_step, t_max, temp=0.5, rho=1.2, n=4, plot=True, write_to_file=True)
