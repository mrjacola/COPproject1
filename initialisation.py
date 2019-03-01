#I can't make this work in the jupyter notebook calling the function, it returns some weird errors such as 'division by zero'
 in lines where there are no divisions. However, if I copy the function directly in the notebook it works and it generates the fcc lattice
 (or at least I think so). I have uploaded it in a different file. Check if at least one of them works for you.

import numpy as np
from itertools import product
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def init_sc(N, d, L, rand_disp=False, disp_size=0.1):
    # Works for any dimension d
    if np.power(N,1/d) % 1 != 0:
        nr_one_d_points = int(np.power(N,1/d)) + 1
    else:
        nr_one_d_points = int(np.power(N,1/d))

    one_d_grid = np.linspace(-L/2 * (1 - 1/nr_one_d_points), L/2 * (1 - 1/nr_one_d_points), nr_one_d_points)
    x_0 = np.array(list(product(one_d_grid, repeat=d)))

    if not rand_disp:
        return x_0
    else:
        return x_0 + disp_size*(np.random.rand(nr_one_d_points**2, d) - 0.5)


def init_3dfcc(N, L, rand_disp=False, disp_size=0.1):
    n= int(np.rint(np.power(N/4, 1/3)))  #number of cells along one edge of the box
    actualN = int(np.power(n, 3)*4)      #actual number of particles adjusted to fit in an fcc
    b=L/n                                #side of a cell
    lat = np.zeros([int(actualN/4), 4, 3])
    cell=np.asarray([[0,0,0],[0, 1, 1],[1, 0, 1],[1, 1, 0]])

    for c in product(range(0, n), repeat=3):
        lat[c[2]+c[1]*n+c[0]*n*n,:,:] = b/2*cell+[c[0]*b, c[1]*b, c[2]*b]-(L)/2


    x = np.reshape(lat, (actualN,3))       

    return x
