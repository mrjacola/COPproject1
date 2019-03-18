# Modules
import numpy as np
import matplotlib.pyplot as plt
from itertools import product
import av_corr as ac
import time
%matplotlib notebook

# Home made module
import remove_diag as rd
import init_lattice as il

# Constants for book keeping
# m = 6.6e-26  # kg
# sigma = 3.405  # Å
# k_B = 1.38e-23  # J/K
# epsilon = 119.8 * k_B  # J

# Sim parameters
L = 10               # Box size in units of sigma
N = 108             # Nr of particles
h = 1e-3            # Time step in units of sqrt(m*sigma^2/epsilon) approx 2e-12 s
d = 3               # Dimensions


v_range = 1.5       # v_max - v_min
beta = 10           # eps/k_B*T
t_max = 5           # Max sim time
nr_steps = int(t_max/h)
rescale_vel = True
lam = 100
 

# Calculate lattice
t0 = time.time()

# x[:,:,0] = il.init_sc(N, d, L, rand_disp=False)
x_0, act_N = il.init_3dfcc(N, L, rand_disp=False)

t1 = time.time()

rho = act_N / L**d               # Number density
n_bins=int(act_N*(act_N-1)/20)   # Number of bins that will 
                                 #be used to calculate the correlation function.           #HERE NUMBER OF BINS CALCULATION (GO TO 62)
                                 # With this definition we will get 
                                 #10 pairs per bin on average

print("Time to init lattice: " + str(t1 - t0))
print("Simulation parameters")
print("Desired N: " + str(N))
print("Actual N: " + str(act_N))
print("Number density: " + str(rho))
print("Box size: " + str(L))

# Allocation
x = np.zeros([act_N, d, nr_steps])
v = np.zeros([act_N, d, nr_steps])
dx = np.zeros([act_N, act_N, d])
E_tot = np.zeros(nr_steps)
E_pot = np.zeros(nr_steps)
E_kin = np.zeros(nr_steps)
corr = np.zeros([n_bins, nr_steps])                                                                #DECLARATION HERE (GO TO 105)

# Initial lattice positions
x[:,:,0] = x_0

# Maxwell distribution (Gaussian in components)
v[:,:,0] = np.random.normal(loc=0.0, scale=1/np.sqrt(beta), size=(act_N, d))

# x[:, :, 0] = [[-0.6, 0], [0.6, 0]]
# v[:, :, 0] = [[2, 0], [-2, 0]]

# Calculate initial distances, forces etc
for i, j in product(range(act_N), repeat=2):
    dx[i, j, :] = x[i, :, 0] - x[j, :, 0]

# Implement boundary nearest image convention
dx[dx > L / 2] = dx[dx > L / 2] - L
dx[dx < -L / 2] = dx[dx < -L / 2] + L

# Remove zeros on diagonal
dx_no_diag = np.zeros([act_N, act_N - 1, d])

for i in range(d):
    dx_no_diag[:, :, i] = rd.remove_diag(dx[:, :, i])

# Calc distance between pairs of atoms
r = np.sqrt(np.sum(np.square(dx_no_diag), 2))

x_vec = dx_no_diag
r_6 = np.power(r, -6)

# Calc force vector on each particle
F = -(np.sum((24 * (1 - 2 * np.power(r, -6)) * np.power(r, -8)).T * x_vec.T, 1)).T


# Run simulation in time
for t in range(nr_steps - 1):
    if t % 100 == 0:
        print("Simulation progress: " + str(100 * t * h / t_max) + "%")

    # Calc kinetic energy
    E_kin[t] = np.sum(np.square(v[:, :, t])) / 2
    E_pot[t] = np.sum(4 * (r_6 - 1) * r_6) / 2
    E_tot[t] = E_kin[t] + E_pot[t]
    #Calc correlation histogram
    corr[:,t], bin_edges=np.histogram(r,n_bins,(0,1/2*L)) #np.sqrt(3)/2*L                               #HERE MEASURMENT (GO TO 162)
    
    

    # Rescale velocities for better temperature settings
    if rescale_vel:
        lam = np.sqrt(((act_N - 1) * (3 / 2)) / (E_kin[t] * beta))
        v[:, :, t] = lam * v[:, :, t]
        print("lambda(" + str(t) + "): " + str(lam))

    if abs(lam - 1) < 1e-3:
        rescale_vel = False

    # Update positions
    x[:, :, t + 1] = x[:, :, t] + v[:, :, t] * h + F * np.square(h) / 2

    # Implement boundary conditions

    x[:, :, t + 1][x[:, :, t + 1] > L / 2] = x[:, :, t + 1][x[:, :, t + 1] > L / 2] - L
    x[:, :, t + 1][x[:, :, t + 1] < -L / 2] = x[:, :, t + 1][x[:, :, t + 1] < -L / 2] + L

    # Calculate force at new positions

    # Calc all coordinate wise distances between particles
    for i, j in product(range(act_N), repeat=2):

        dx[i, j, :] = x[i, :, t + 1] - x[j, :, t + 1]

    # Implement boundary nearest image convention
    dx[dx > L / 2] = dx[dx > L / 2] - L
    dx[dx < -L / 2] = dx[dx < -L / 2] + L

    # Remove zeros on diagonal
    dx_no_diag = np.zeros([act_N, act_N - 1, d])

    for i in range(d):
        dx_no_diag[:, :, i] = rd.remove_diag(dx[:, :, i])

    # Calc distance between
    r = np.sqrt(np.sum(np.square(dx_no_diag), 2))

    x_vec = dx_no_diag

    # Calc force vector on each particle
    F_new = -(np.sum((24 * (1 - 2 * np.power(r, -6)) * np.power(r, -8)).T * x_vec.T, 1)).T

    # Calculate energy
    r_6 = np.power(r, -6)
    

    # Update velocities
    v[:, :, t + 1] = v[:, :, t] + (F + F_new) * h / 2
    F = F_new


    
#Calc observables
av_corr,R=ac.g(corr,bin_edges,L,act_N)                                                            #HERE CALCULATION
    
    
fig = plt.figure(num=1, figsize=(10,5))
if d == 2:
    ax1 = fig.add_subplot(1, 2, 1)
elif d == 3:
    ax1 = fig.add_subplot(1, 2, 1, projection="3d")

ax2 = fig.add_subplot(1, 2, 2)

for t in range(0, nr_steps - 1, 20):

    ax1.clear()
    if d == 2:
        ax1.scatter(x[:, 0, t], x[:, 1, t], s=4, color="black")
        ax1.set_aspect("equal", "box")
        ax1.set_xlim(- L / 2, L / 2)
        ax1.set_ylim(- L / 2, L / 2)
    elif d == 3:
        ax1.scatter(x[:, 0, t], x[:, 1, t], x[:, 2, t], s=4, color="black")
        ax1.set_aspect("equal", "box")
        ax1.set_xlim(- L / 2, L / 2)
        ax1.set_ylim(- L / 2, L / 2)
        ax1.set_zlim(- L / 2, L / 2)
    ax1.set_title("Progress: " + str(round(100*t*h/t_max)) + "% " + " t_max = " + str(t_max))

    ax2.clear()
    ax2.plot(np.linspace(0.0, (t + 1) * h, num=t + 1), E_tot[0: t + 1])
    ax2.plot(np.linspace(0.0, (t + 1) * h, num=t + 1), E_pot[0: t + 1])
    ax2.plot(np.linspace(0.0, (t + 1) * h, num=t + 1), E_kin[0: t + 1])
    ax2.set_xlim(0.0, t_max)
    ax2.set_title("Total energy")
    ax2.set_xlabel("Time")
    plt.pause(0.001)

# plt.hist(v[:,0,0], color='b')
# plt.hist(v[:,1,0], color='r')
    plt.show()
