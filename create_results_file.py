

def create_results_file(filename, rho, N, L, temp, h, t_max, real_sim_time, vel_rescalings, spec_heat, corr_func, r_vals):

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
    f.write("vel_rescalings = %d\n" % vel_rescalings)
    f.write("\n")
    f.write("Sim results:\n")
    f.write("------------\n")
    f.write("spec_heat = %.4f\n" % spec_heat)

    f.write("corr_func:\n")
    for corr_val in corr_func:
        f.write("%.4f " % corr_val)

    f.write("\n\n")

    f.write("r_vals:\n")
    for r_val in r_vals:
        f.write("%.4f " % r_val)

    f.close()
