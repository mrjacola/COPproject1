

def create_results_file(filename, rho, N, L, temp, h, t_max, real_sim_time):

    f = open(filename, "w+")
    f.write("Parameters\n")
    f.write("rho: %.4f\n" % rho)
    f.write("N: %.4f\n" % N)
    f.write("L: %.4f\n" % L)
    f.write("temp: %.4f\n" % temp)
    f.write("time_step: %.4f\n" % h)
    f.write("sim_time: %.4f\n" % t_max)
    f.write("real_sim_time: %.4f s\n" % real_sim_time)
    f.close()
