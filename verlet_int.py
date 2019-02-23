# Modules
import numpy as np
import matplotlib.pyplot as plt
import time

# Home made module
import remove_diag as rd
import init_lattice as il

# Constants for book keeping
m = 6.6e-26  # kg
sigma = 3.405  # Å
k_B = 1.38e-23  # J/K
epsilon = 119.8 * k_B  # J

# Sim parameters
L = 10      # Box size in units of sigma
N = 36      # Nr of particles (use a square for now...)
h = 1e-3    # Time step in units of sqrt(m*sigma^2/epsilon) approx 2e-12 s
d = 2       # Dimensions
v_range = 1.5  # v_max - v_min
t_max = 5     # Max sim time
nr_steps = int(t_max/h)

# Allocation
x = np.zeros([N, d, nr_steps])
v = np.zeros([N, d, nr_steps])
dx = np.zeros([N, N, d])
E_tot = np.zeros(nr_steps)

# Init
x[:,:,0] = il.init_sc(N, d, L)
v[:,:,0] = v_range*(np.random.rand(N, d) - 0.5)

# x[:, :, 0] = [[-0.6, 0], [0.6, 0]]
# v[:, :, 0] = [[2, 0], [-2, 0]]

# Calculate initial distances, forces etc
for i in range(N):
    for j in range(N):
        dx[i, j, :] = x[i, :, 0] - x[j, :, 0]

# Implement boundary nearest image convention
dx[dx > L / 2] = dx[dx > L / 2] - L
dx[dx < -L / 2] = dx[dx < -L / 2] + L

# Remove zeros on diagonal
dx_no_diag = np.zeros([N, N - 1, d])

for i in range(d):
    dx_no_diag[:, :, i] = rd.remove_diag(dx[:, :, i])

# Calc distance between
r = np.sqrt(np.sum(np.square(dx_no_diag), 2))

x_vec = dx_no_diag

# Calc force vector on each particle
F = -(np.sum((24 * (1 - 2 * np.power(r, -6)) * np.power(r, -8)).T * x_vec.T, 1)).T

# Run simulation in time
for t in range(nr_steps - 1):
    # Update positions
    x[:, :, t + 1] = x[:, :, t] + v[:, :, t] * h + F * np.square(h) / 2

    # Implement boundary conditions

    x[:, :, t + 1][x[:, :, t + 1] > L / 2] = x[:, :, t + 1][x[:, :, t + 1] > L / 2] - L
    x[:, :, t + 1][x[:, :, t + 1] < -L / 2] = x[:, :, t + 1][x[:, :, t + 1] < -L / 2] + L

    # Calculate force at new positions

    # Calc all coordinate wise distances between particles
    for i in range(N):
        for j in range(N):
            dx[i, j, :] = x[i, :, t + 1] - x[j, :, t + 1]

    # Implement boundary nearest image convention
    dx[dx > L / 2] = dx[dx > L / 2] - L
    dx[dx < -L / 2] = dx[dx < -L / 2] + L

    # Remove zeros on diagonal
    dx_no_diag = np.zeros([N, N - 1, d])

    for i in range(d):
        dx_no_diag[:, :, i] = rd.remove_diag(dx[:, :, i])

    # Calc distance between
    r = np.sqrt(np.sum(np.square(dx_no_diag), 2))

    x_vec = dx_no_diag

    # Calc force vector on each particle
    F_new = -(np.sum((24 * (1 - 2 * np.power(r, -6)) * np.power(r, -8)).T * x_vec.T, 1)).T

    # Update state variables
    # v[:,:,t+1] = v[:,:,t] + F*h
    # x[:,:,t+1] = x[:,:,t] + v[:,:,t]*h

    v[:, :, t + 1] = v[:, :, t] + (F + F_new) * h / 2
    F = F_new

    # Calculate energy
    E_kin = np.sum(np.square(v[:, :, t])) / 2
    r_6 = np.power(r, -6)
    E_pot = np.sum(4 * (r_6 - 1) * r_6) / 2
    E_tot[t] = E_kin + E_pot

fig = plt.figure(num=1, figsize=(10,5))
ax1 = fig.add_subplot(1, 2, 1)
ax2 = fig.add_subplot(1, 2, 2)


for t in range(0, nr_steps - 1, 20):

    ax1.clear()
    ax1.scatter(x[:, 0, t], x[:, 1, t], s=4, color="black")
    ax1.set_aspect("equal", "box")
    ax1.set_xlim(- L / 2, L / 2)
    ax1.set_ylim(- L / 2, L / 2)
    ax1.set_title("Progress: " + str(round(100*t*h/t_max)) + "% " + " t_max = " + str(t_max))

    ax2.clear()
    ax2.plot(np.linspace(0.0, (t + 1) * h, num=t + 1), E_tot[0: t + 1])
    ax2.set_ylim(np.min(E_tot) - 1, np.max(E_tot) + 1)
    ax2.set_xlim(0.0, t_max)
    ax2.set_title("Total energy")
    ax2.set_xlabel("Time")
    plt.pause(0.001)

plt.show()
