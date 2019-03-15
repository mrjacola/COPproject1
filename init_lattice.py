import numpy as np
from itertools import product
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def init_sc(N, d, L, rand_disp=False, disp_size=0.1):
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
        return x_0, actual_N
    else:
        return x_0 + disp_size * (np.random.rand(x_0.shape) - 0.5), actual_N


if __name__ == "__main__":

    x = init_3dfcc(100, 10)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x[:,0],x[:,1],x[:,2])
    plt.show()

    print(x)


# L = 10
# x = init_lattice(9, 2, L)
# print(x.shape)
#
# plt.scatter(x[:,0],x[:,1])
# plt.ylim(-L/2, L/2)
# plt.xlim(-L/2, L/2)
# plt.show()






    # dx = L/np.power(N,1/d)
    # x_0 = np.zeros(N,d)
    # x_0[0, :] = (L/2)*np.ones(d)
    #
    # for i in range(N):
    #
    #     if x_0[i,0] < L/2-dx:
    #
    #         x_0[i+1,0] = x_0[i, 0] + dx
    #         x_0[i+1,1] = x_0[i,1]
    #     else:
    #         x_0[i+1,0] = -L/2
    #         x_0[i+1,1] = x_0[i, 1] + dx