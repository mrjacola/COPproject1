import numpy as np
import matplotlib.pyplot as plt

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

# Init
x[:,:,0] = L*(np.random.rand(N, d) - 0.5)
v[:,:,0] = v_range*(np.random.rand(N, d) - 0.5)

r = np.zeros([N, N, d])

for i in range(N):
    for j in range(N):

        r[i, j, :] = x[i, :, 0] - x[j, :, 0]

r[r > L/2] = r[r > L/2] - L
r[r < -L/2] = r[r < -L/2] + L

print(r[r[:]!=0])
#F = np.sum(24*(1 - 2/np.power(r[r!=0], 6))/np.power(r[r!=0], 7), 0)

#print(F)



#plt.scatter(x[:,0,0], x[:,1,0])
#plt.show()

#print(type(x[0,0,0]))
#print(v_0)