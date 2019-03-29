

def create_results_file(filename,
                        rho,
                        N,
                        L,
                        temp,
                        h,
                        t_max,
                        real_sim_time,
                        obs_time,
                        vel_rescalings,
                        avg_cut,
                        spec_heat,
                        sigma_C_V,
                        pressure,
                        sigma_P,
                        corr_func,
                        sigma_corr,
                        r_vals,
                        E_kin,
                        E_pot,
                        T_true):

    f = open(filename, "w+")
    f.write("Parameters\n")
    f.write("------------\n")
    f.write("rho = %.4f\n" % rho)
    f.write("N =\n")
    f.write("%d\n" % N)
    f.write("L = %.4f\n" % L)
    f.write("T = %.4f\n" % temp)
    f.write("time_step h = %.4f\n" % h)
    f.write("sim_time = %.4f\n" % t_max)
    f.write("real_sim_time (s) =\n")
    f.write("%.4f\n" % real_sim_time)
    f.write("obs_time = %.4f s\n" % obs_time)
    f.write("vel_rescalings = %d\n" % vel_rescalings)
    f.write("avg_cut = %d\n" % avg_cut)

    f.write("Sim results:\n")
    f.write("------------\n")
    f.write("C_V = %.4f\n" % spec_heat)
    f.write("sigma_C_V = %.4f\n" % sigma_C_V)
    f.write("P = %.4f\n" % pressure)
    f.write("sigma_P = %.4f\n" % sigma_P)
    f.write("T_true = %.4f\n" % T_true)

    f.write("E_kin: \n")
    for E_kin_val in E_kin:
        f.write("%.4f " % E_kin_val)

    f.write("\n")

    f.write("E_pot: \n")
    for E_pot_val in E_pot:
        f.write("%.4f " % E_pot_val)

    f.write("\n")

    f.write("corr_func:\n")
    for corr_val in corr_func:
        f.write("%.4f " % corr_val)

    f.write("\n")

    f.write("sigma_corr:\n")
    for sigma_corr_val in sigma_corr:
        f.write("%.4f " % sigma_corr_val)

    f.write("\n")

    f.write("r_vals:\n")
    for r_val in r_vals:
        f.write("%.4f " % r_val)

    f.close()
