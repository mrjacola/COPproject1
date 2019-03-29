# Modules
import numpy as np
import matplotlib.pyplot as plt
from itertools import product
import time

# Home made modules
import remove_diag as rd
import init_lattice as il
import create_results_file as crf
import observables as obs


def run_sim(h, t_max, temp, rho, d=3, n=5, plot=True, write_to_file=True):
    # h time step
    # t_max simulation time
    # temp temperature
    # rho number density
    # d dimensions

    nr_steps = int(t_max/h)         # Number of time steps in simulation
    N = 4 * n ** d                  # Number of particles consistent with FCC
    L = np.power(N / rho, 1 / d)    # Simulation box size
    n_bins = int(N * (N - 1) / 100)  # Number of bins that will
                                    # be used to calculate the correlation function.
                                    # With this definition we will get
                                    # 10 pairs per bin on average

    # Calculate lattice
    t0 = time.time()
    x_0 = il.init_3dfcc(N, L, rand_disp=False)
    t1 = time.time()

    print("Time to init lattice: " + str(t1 - t0))
    print("Simulation parameters:")
    print("Particle number N: " + str(N))
    print("Number density: " + str(rho))
    print("Box size: " + str(L))
    print("Temperature: " + str(temp))
    print("Simulation time: " + str(t_max))
    print("Time step: " + str(h))

    # Allocation
    x = np.zeros([N, d, nr_steps])
    v = np.zeros([N, d, nr_steps])
    dx = np.zeros([N, N, d])
    E_tot = np.zeros(nr_steps)
    E_pot = np.zeros(nr_steps)
    E_kin = np.zeros(nr_steps)
    corr = np.zeros([n_bins, nr_steps])
    press_minus_rho_T = np.zeros(nr_steps)

    # Initial lattice positions
    x[:,:,0] = x_0

    # Initial Maxwell distribution (Gaussian in components)
    v[:,:,0] = np.random.normal(loc=0.0, scale=np.sqrt(temp), size=(N, d))

    # Calculate initial distances
    for i, j in product(range(N), repeat=2):
        dx[i, j, :] = x[i, :, 0] - x[j, :, 0]

    # Implement boundary nearest image convention
    dx[dx > L / 2] = dx[dx > L / 2] - L
    dx[dx < -L / 2] = dx[dx < -L / 2] + L

    # Remove zeros on diagonal
    dx_no_diag = np.zeros([N, N - 1, d])

    for i in range(d):
        dx_no_diag[:, :, i] = rd.remove_diag(dx[:, :, i])

    # Calc initial distance between
    r = np.sqrt(np.sum(np.square(dx_no_diag), 2))
    r_6 = np.power(r, -6)

    x_vec = dx_no_diag

    # Calc initial force vector on each particle
    f = - 24 * (1 - 2 * np.power(r, -6)) * np.power(r, -8)
    F = (np.sum(f.T * x_vec.T, 1)).T

    vel_rescalings = 0
    last_scaling = 0
    eq_time = 0.8 / h
    lam_tol = 3e-2

    # Run simulation in time
    t0 = time.time()
    for t in range(nr_steps - 1):
        if t % int(nr_steps / 100) == 0:
            print("Simulation progress: " + str(np.around(100 * t * h / t_max, decimals=1)) + "%")

        # Calc kinetic energy
        E_kin[t] = np.sum(np.square(v[:, :, t])) / 2
        # Potential
        E_pot[t] = np.sum(4 * (r_6 - 1) * r_6) / 2
        # Total
        E_tot[t] = E_kin[t] + E_pot[t]
        # Correlation
        corr[:, t], bin_edges = np.histogram(r, n_bins, (0, 1 / 2 * L))  # np.sqrt(3)/2*L
        # Pressure
        press_minus_rho_T[t] = (1 / 6)/(L**3) * np.sum(f * r)

        # Rescale velocities for better temperature equilibration
        lam = np.sqrt(((N - 1) * (3 / 2) * temp) / (E_kin[t]))
        if t != 0 and t % eq_time == 0 and abs(lam - 1) > lam_tol:
            v[:, :, t] = lam * v[:, :, t]
            print("lambda(" + str(t) + "): " + str(lam))
            vel_rescalings += 1
            last_scaling = t

        # Update positions
        x[:, :, t + 1] = x[:, :, t] + v[:, :, t] * h + F * np.square(h) / 2
        # Implement boundary conditions
        x[:, :, t + 1][x[:, :, t + 1] > L / 2] = x[:, :, t + 1][x[:, :, t + 1] > L / 2] - L
        x[:, :, t + 1][x[:, :, t + 1] < -L / 2] = x[:, :, t + 1][x[:, :, t + 1] < -L / 2] + L

        # Calc all coordinate wise distances between particles
        for i, j in product(range(N), repeat=2):

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

        # For calculation of energy
        r_6 = np.power(r, -6)

        x_vec = dx_no_diag

        # Calc force vector on each particle
        f = - 24 * (1 - 2 * r_6) * np.power(r, -8)
        F_new = (np.sum(f.T * x_vec.T, 1)).T

        # Update velocities and forces
        v[:, :, t + 1] = v[:, :, t] + (F + F_new) * h / 2
        F = F_new

    t1 = time.time()
    real_sim_time = t1 - t0
    print("Simulation time: "
          + str(int(real_sim_time / 60)) + " min " + str(np.around(real_sim_time % 60, decimals=0)) + " s")

    t0 = time.time()
    # Index for start of averaging, should be after equilibration
    avg_init_cut = int(last_scaling + eq_time)
    # Calc correlation fcn as average
    av_corr, R, sigma_corr = obs.g(corr[:, avg_init_cut:nr_steps - 2], bin_edges, nr_steps, L, N)
    # Calc specific heat
    specific_heat, sigma_specific_heat, E_kin_mean = obs.c(E_kin[avg_init_cut:nr_steps - 2], nr_steps, N)
    temp_true = (2 / 3) / (N - 1) * E_kin_mean
    # Pressure
    av_press, sigma_press = obs.p(press_minus_rho_T[avg_init_cut:nr_steps - 2], temp_true, rho, nr_steps)

    t1 = time.time()
    t_obs = t1 - t0
    print("Time to calculate observables: " + str(t_obs))

    print("C_V = " + str(specific_heat))
    print("std C_V = " + str(sigma_specific_heat))
    print("p = " + str(av_press))
    print("std p = " + str(sigma_press))
    print("Last scaling: " + str(last_scaling))
    print("true_temp = " + str(temp_true))

    if write_to_file:
        # Write results to file
        crf.create_results_file("sim_data/T_" + str(temp) + "_rho_" + str(rho) + "_N_" + str(N) + ".txt",
                                rho,
                                N,
                                L,
                                temp,
                                h,
                                t_max,
                                real_sim_time,
                                t_obs,
                                vel_rescalings,
                                avg_init_cut,
                                spec_heat=specific_heat,
                                sigma_C_V=sigma_specific_heat,
                                pressure=av_press,
                                sigma_P=sigma_press,
                                corr_func=av_corr,
                                sigma_corr=sigma_corr,
                                r_vals=R,
                                E_kin=E_kin,
                                E_pot=E_pot,
                                T_true=temp_true)

    if plot:
        fig = plt.figure(num=1, figsize=(10,5))

        if d == 2:
            ax1 = fig.add_subplot(1, 2, 1)
        elif d == 3:
            ax1 = fig.add_subplot(1, 2, 1, projection="3d")

        ax2 = fig.add_subplot(1, 2, 2)

        atom_color = (243 / 255, 8 / 255, 71 / 255)

        for t in range(0, nr_steps - 1, 20):

            ax1.clear()
            if d == 2:
                ax1.scatter(x[:, 0, t], x[:, 1, t], s=10, color=atom_color)
                ax1.set_aspect("equal", "box")
                ax1.set_xlim(- L / 2, L / 2)
                ax1.set_ylim(- L / 2, L / 2)
            elif d == 3:
                ax1.scatter(x[:, 0, t], x[:, 1, t], x[:, 2, t], s=10, color=atom_color)
                ax1.set_aspect("equal", "box")
                ax1.set_xlim(- L / 2, L / 2)
                ax1.set_ylim(- L / 2, L / 2)
                ax1.set_zlim(- L / 2, L / 2)
            ax1.set_title("Progress: " + str(round(100 * t * h / t_max)) + "% " + " t_max = " + str(t_max))

            ax2.clear()
            ax2.plot(np.linspace(0.0, (t + 1) * h, num=t + 1), E_tot[0: t + 1] / N, color="black", label="$E_{tot}$")
            ax2.plot(np.linspace(0.0, (t + 1) * h, num=t + 1), E_pot[0: t + 1] / N, color="red", label="$E_{pot}$")
            ax2.plot(np.linspace(0.0, (t + 1) * h, num=t + 1), E_kin[0: t + 1] / N, color="blue", label="$E_{kin}$")
            ax2.set_xlim(0.0, t_max)
            ax2.set_title("Energy evolution of the simulation")
            ax2.set_ylabel("E / N (natural unit)")
            ax2.set_xlabel("Time (natural unit)")
            ax2.legend()
            plt.pause(0.001)

        plt.show()
