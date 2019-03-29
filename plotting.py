import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from scipy.stats import linregress
from matplotlib import rc


def plot_corr(filenames, nr_bins):

    corr_skip_header = 26
    corr_skip_footer = 4
    sigma_corr_skip_header = 28
    sigma_corr_skip_footer = 2
    r_skip_header = 30
    r_skip_footer = 0
    corr = np.zeros([len(filenames), nr_bins])
    sigma_corr = np.zeros([len(filenames), nr_bins])
    r = np.zeros([len(filenames), nr_bins])
    i = 0

    for filename in filenames:

        corr[i, :] = np.genfromtxt(filename, skip_header=corr_skip_header, skip_footer=corr_skip_footer)
        # print(corr[i, :].shape)
        sigma_corr[i, :] = np.genfromtxt(filename, skip_header=sigma_corr_skip_header, skip_footer=sigma_corr_skip_footer)
        r[i, :] = np.genfromtxt(filename, skip_header=r_skip_header, skip_footer=r_skip_footer)
        # print(r[i, :].shape)
        i += 1

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    fs = 16
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(2, 1, 1)

    ax.plot(r[0, :], corr[0, :], color="black", label="Solid")
    ax.plot(r[1, :], corr[1, :], color="blue", label="Liquid")
    ax.plot(r[2, :], corr[2, :], color="red", label="Gas")

    ax.set_title("Pair correlation function for different phases of Argon", fontsize=fs)
    # ax.set_xlabel("Distance $r$", fontsize=fs)
    ax.set_ylabel("Correlation $g(r)$", fontsize=fs)
    ax.set_xlim(0.0, max(r[2, :]))
    ax.set_xticks(np.arange(0, max(r[2, :]), 0.5))
    ax.legend(fontsize=fs)

    ax_sigma = fig.add_subplot(2, 1, 2)
    ax_sigma.plot(r[0, :], sigma_corr[0, :], color="black", label="Solid")
    ax_sigma.plot(r[1, :], sigma_corr[1, :], color="blue", label="Liquid")
    ax_sigma.plot(r[2, :], sigma_corr[2, :], color="red", label="Gas")
    # ax_sigma.set_title("STD", fontsize=fs)
    ax_sigma.set_xlabel("Distance $r$", fontsize=fs)
    ax_sigma.set_ylabel(r"Standard deviation $\sigma_{g}(r)$", fontsize=fs)
    ax_sigma.set_xlim(0.0, max(r[2, :]))
    ax_sigma.set_xticks(np.arange(0, max(r[2, :]), 0.5))
    # ax_sigma.legend(fontsize=fs)

    plt.show()


def plot_energy(filename, t_max, N):

    E_kin_skip_header = 22
    E_kin_skip_footer = 8
    E_pot_skip_header = 24
    E_pot_skip_footer = 6

    E_kin = np.genfromtxt(filename, skip_header=E_kin_skip_header, skip_footer=E_kin_skip_footer)
    E_pot = np.genfromtxt(filename, skip_header=E_pot_skip_header, skip_footer=E_pot_skip_footer)
    E_tot = E_kin + E_pot

    nr_steps = len(E_kin)
    time = np.linspace(0, t_max, num=nr_steps - 1)

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    fig = plt.figure(figsize=(6, 5))
    ax = fig.add_subplot(111)
    lw = 1.0
    fs = 16
    ax.plot(time, E_kin[0:nr_steps-1] / N, linestyle=(0, (5, 1)), color="blue", linewidth=lw, label=r'$E_{\textrm{kin}}$')
    ax.plot(time, E_pot[0:nr_steps-1] / N, linestyle=(0, (1, 1)), color="red", linewidth=lw, label=r"$E_{\textrm{pot}}$")
    ax.plot(time, E_tot[0:nr_steps-1] / N, linestyle=(0, ()), color="black", linewidth=lw, label=r"$E_{\textrm{tot}}$")
    ax.legend(fontsize=fs, loc=[0.72, 0.52])
    ax.set_title("Energy evolution of a simulation", fontsize=fs)
    ax.set_xlabel("Time", fontsize=fs)
    ax.set_ylabel("Average energy per particle", fontsize=fs)
    ax.set_xlim(0, t_max)
    axins = inset_axes(ax, width=1.5, height=.5, loc="center")
    axins.set_title(r"Fluctuations in $E_{\textrm{tot}}$")
    axins.plot(time, E_tot[0:nr_steps-1] / N, color="black", linewidth=lw)
    axins.set_xlim(4.5, 5.5)
    axins.set_ylim(-4.06, -4.055)
    plt.show()


def plot_sim_time(filenames):

    real_sim_time_skip_header = 10
    real_sim_time_skip_footer = 20
    N_skip_header = 4
    N_skip_footer = 26
    N = np.zeros(len(filenames))
    real_sim_time = np.zeros(len(filenames))

    i = 0
    for filename in filenames:
        N[i] = int(np.genfromtxt(filename, skip_header=N_skip_header, skip_footer=N_skip_footer))
        real_sim_time[i] = float(np.genfromtxt(filename, skip_header=real_sim_time_skip_header, skip_footer=real_sim_time_skip_footer))
        i += 1

    log_N = np.log(N)
    log_time = np.log(real_sim_time / 60)
    reg_params = linregress(log_N, log_time)
    slope = reg_params[0]
    intercept = reg_params[1]
    fit = slope * log_N + intercept

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    fig = plt.figure(figsize=(6, 5))
    ax = fig.add_subplot(111)
    ax.plot(log_N, log_time, "+", markeredgecolor="black", markersize=10, label="Simulation data")
    ax.plot(log_N, fit, "--", color="grey", label="Linear fit, slope = " + str(np.around(slope, 3)))
    fs = 16
    ax.legend(fontsize=fs)
    ax.set_title("Simulation times for different number of particles", fontsize=fs)
    ax.set_xlabel(r"$\log(N)$", fontsize=fs)
    ax.set_ylabel(r"$\log(t_{\textrm{sim}}/1\; \textrm{min})$", fontsize=fs)
    plt.show()


if __name__ == "__main__":

    N = 500
    n_bins = int(N * (N - 1) / 100)
    plot_corr(["sim_data/T_0.5_rho_1.2_N_500.txt", "sim_data/T_1.0_rho_0.8_N_500.txt", "sim_data/T_3.0_rho_0.3_N_500.txt"],
             nr_bins=n_bins)

    # plot_energy("sim_data/T_1.0_rho_0.8_N_500.txt", t_max=10.0, N=N)

    # plot_sim_time(["sim_data/T_3.0_rho_0.3_N_32.txt",
    #               "sim_data/T_3.0_rho_0.3_N_108.txt",
    #               "sim_data/T_3.0_rho_0.3_N_256.txt",
    #               "sim_data/T_3.0_rho_0.3_N_500.txt",
    #               "sim_data/T_3.0_rho_0.3_N_864.txt"])

