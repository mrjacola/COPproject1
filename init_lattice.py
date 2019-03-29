import numpy as np
from itertools import product


def init_sc(N, d, L, rand_disp=False, disp_size=0.1):
    # Function for giving initial positions in a sc lattice
    # Works for any dimension d
    if np.power(N, 1 / d) % 1 != 0:
        nr_one_d_points = int(np.power(N, 1 / d)) + 1
    else:
        nr_one_d_points = int(np.power(N, 1 / d))

    one_d_grid = np.linspace(-L / 2 * (1 - 1 / nr_one_d_points), L / 2 * (1 - 1 / nr_one_d_points), nr_one_d_points)
    x_0 = np.array(list(product(one_d_grid, repeat=d)))

    if not rand_disp:
        return x_0, N
    else:
        return x_0 + disp_size * (np.random.rand(nr_one_d_points ** 2, d) - 0.5), N


def init_3dfcc(N, L, rand_disp=False, disp_size=0.1):
    # Function for giving initial positions in a fcc lattice
    # Works only for 3D

    n = int(np.rint(np.power(N / 4, 1 / 3)))
    actual_N = int(np.power(n, 3) * 4)
    b = L / n
    lat = np.zeros([int(actual_N / 4), 4, 3])
    cell = np.asarray([[0, 0, 0], [0, 1, 1], [1, 0, 1], [1, 1, 0]])

    for c in product(range(0, n), repeat=3):
        lat[c[2] + c[1] * n + c[0] * n * n, :, :] = b / 2 * cell + [c[0] * b, c[1] * b, c[2] * b] - L / 2

    x_0 = np.reshape(lat, (actual_N, 3))

    if not rand_disp:
        return x_0
    else:
        return x_0 + disp_size * (np.random.rand(x_0.shape) - 0.5)
