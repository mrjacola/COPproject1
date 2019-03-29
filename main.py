# Script for tuning and running simulation

import numpy as np

import verlet_int as vi

# Constants for book keeping
# m = 6.6e-26  # kg
# sigma = 3.405  # Å
# k_B = 1.38e-23  # J/K
# epsilon = 119.8 K * k_B  # J

# Simulation parameters
time_step = 5e-3
t_max = 1e1

if __name__ == "__main__":
    # Run simulations with input parameters rho and N = 4*n^3

    # For corr and energy plots
    # Gas
    vi.run_sim(time_step, t_max, temp=3.0, rho=0.3, n=5, plot=False, write_to_file=True)
    # Liquid
    vi.run_sim(time_step, t_max, temp=1.0, rho=0.8, n=5, plot=False, write_to_file=True)
    # Solid
    vi.run_sim(time_step, t_max, temp=0.5, rho=1.2, n=5, plot=False, write_to_file=True)

    # For performance plot
    ns = np.array([2, 3, 4, 6])
    for n in ns:

        vi.run_sim(time_step, t_max, temp=3.0, rho=0.3, n=n, plot=False, write_to_file=True)