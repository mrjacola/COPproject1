import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":

    corr_1 = np.genfromtxt("sim_data/T_0.5_rho_1.2.txt", skip_header=15, skip_footer=2)
    r_1 = np.genfromtxt("sim_data/T_0.5_rho_1.2.txt", skip_header=18, skip_footer=0)

    corr_2 = np.genfromtxt("sim_data/T_1.0_rho_0.8.txt", skip_header=15, skip_footer=2)
    r_2 = np.genfromtxt("sim_data/T_1.0_rho_0.8.txt", skip_header=18, skip_footer=0)

    corr_3 = np.genfromtxt("sim_data/T_3.0_rho_0.3.txt", skip_header=15, skip_footer=2)
    r_3 = np.genfromtxt("sim_data/T_3.0_rho_0.3.txt", skip_header=18, skip_footer=0)

    print(r_1.shape)
    print(corr_1.shape)
    print(r_2.shape)
    print(corr_2.shape)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(r_1, corr_1, color="black")
    ax.plot(r_2, corr_2, color="blue")
    ax.plot(r_3, corr_3, color="red")
    plt.show()
