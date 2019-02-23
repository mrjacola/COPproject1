import numpy as np
from itertools import product
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def init_sc(N, d, L, rand_disp=False, disp_size=0.1):
    # Works for any dimension d
    if np.sqrt(N) % 1 != 0:
        nr_one_d_points = int(np.sqrt(N)) + 1
    else:
        nr_one_d_points = int(np.sqrt(N))

    one_d_grid = np.linspace(-L/2 * (1 - 1/nr_one_d_points), L/2 * (1 - 1/nr_one_d_points), nr_one_d_points)
    x_0 = np.array(list(product(one_d_grid, repeat=d)))

    if not rand_disp:
        return x_0
    else:
        return x_0 + disp_size*(np.random.rand(nr_one_d_points**2, d) - 0.5)


def init_3dfcc(N, L, d=3, rand_disp=False, disp_size=0.1):
    # Works only for 3D

    nr_one_d_points = int(np.power(N, 1 / 3))
    x_0 = np.zeros([N, d])
    a1 = np.array([1, 1, 0])*L/(2*nr_one_d_points)
    a2 = np.array([1, 0, 1])*L/(2*nr_one_d_points)
    a3 = np.array([0, 1, 1])*L/(2*nr_one_d_points)

    for c in product(range(-nr_one_d_points//2, (nr_one_d_points + 1)//2), repeat=d):

        g = c[0]*a1 + c[1]*a2 + c[2]*a3
        x_0 = np.append(x_0, [g], 0)

    return x_0



x = init_3dfcc(10, 10)

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