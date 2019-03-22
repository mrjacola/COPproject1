

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
                        pressure,
                        corr_func,
                        r_vals,
                        E_kin,
                        E_pot):

    f = open(filename, "w+")
    f.write("Parameters\n")
    f.write("------------\n")
    f.write("rho = %.4f\n" % rho)
    f.write("N = %d\n" % N)
    f.write("L = %.4f\n" % L)
    f.write("T = %.4f\n" % temp)
    f.write("time_step h = %.4f\n" % h)
    f.write("sim_time = %.4f\n" % t_max)
    f.write("real_sim_time = %.4f s\n" % real_sim_time)
    f.write("obs_time = %.4f s\n" % obs_time)
    f.write("vel_rescalings = %d\n" % vel_rescalings)
    f.write("avg_cut = %d\n" % avg_cut)
    f.write("\n")
    f.write("Sim results:\n")
    f.write("------------\n")
    f.write("C_V = %.4f\n" % spec_heat)
    f.write("P = %.4f\n" % pressure)

    f.write("E_kin: \n")
    for E_kin_val in E_kin:
        f.write("%.4f " % E_kin_val)

    f.write("\n\n")

    f.write("E_pot: \n")
    for E_pot_val in E_pot:
        f.write("%.4f " % E_pot_val)

    f.write("\n\n")

    f.write("corr_func:\n")
    for corr_val in corr_func:
        f.write("%.4f " % corr_val)

    f.write("\n\n")

    f.write("r_vals:\n")
    for r_val in r_vals:
        f.write("%.4f " % r_val)

    f.close()
