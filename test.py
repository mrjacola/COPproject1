# Modules
import numpy as np
import matplotlib.pyplot as plt
import time

# Home made module
import remove_diag as rd

# Constants for book keeping

m = 6.6e-26             # kg
sigma = 3.405           # Å
k_B = 1.38e-23          # J/K
epsilon = 119.8*k_B     # J

L = 10                  # Box size in units of sigma
N = 10                  # Nr of particles
h = 1e-2                # Time step in units of sqrt(m*sigma^2/epsilon) approx 2e-12 s
d = 2                   # Dimensions
v_range = 3             # v_max - v_min
t_max = int(100*h)      # Max sim time
nr_steps = int(t_max/h)

# Allocation
x = np.zeros([N, d, nr_steps])
v = np.zeros([N, d, nr_steps])
dx = np.zeros([N, N, d])

# Init
x[:,:,0] = L*(np.random.rand(N, d) - 0.5)
v[:,:,0] = v_range*(np.random.rand(N, d) - 0.5)


# Run simulation in time
for t in range(nr_steps - 1):

    # Calc all coordinate wise distances between particles
    for i in range(N):
        for j in range(N):

            dx[i, j, :] = x[i, :, t] - x[j, :, t]

    # Implement boundary nearest image convention
    dx[dx > L/2] = dx[dx > L/2] - L
    dx[dx < -L/2] = dx[dx < -L/2] + L

    # Remove zeros on diagonal
    dx_no_diag = np.zeros([N, N - 1, d])

    for i in range(d):
        dx_no_diag[:,:,i] = rd.remove_diag(dx[:,:,i])

    # Calc distance between
    r = np.sqrt(np.sum(np.square(dx_no_diag), 2))

    x_vec = dx_no_diag

    # Calc force vector on each particle
    F = -(np.sum((24*(1 - 2*np.power(r, -6))*np.power(r, -8)).T*x_vec.T, 1)).T

    # Update state variables
    v[:,:,t+1] = v[:,:,t] + F*h
    x[:,:,t+1] = x[:,:,t] + v[:,:,t]*h

    # Implement boundary conditions
    for dim in range(d):
        x[:,dim,t + 1][x[:,dim,t + 1] > L/2] = x[:,dim,t + 1][x[:,dim,t + 1] > L/2] - L
        x[:,dim,t + 1][x[:,dim,t + 1] < -L/2] = x[:,dim,t + 1][x[:,dim,t + 1] < -L/2] + L

print(np.max(x))
fig = plt.figure()
ax = fig.add_subplot(111)

for t in range(nr_steps - 1):
    ax.clear()
    ax.scatter(x[:,0,t], x[:,1,t])
    plt.pause(1)
