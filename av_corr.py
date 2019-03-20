# This function requires the of corr and its calculation for every time step, as well as bin_edges

import numpy as np
import math


def g(corr, bin_edges, L, N):

    delta_R = bin_edges[-1] / (len(bin_edges) - 1)
    R = bin_edges[0:len(bin_edges) - 1] + delta_R / 2
    av_corr = np.divide(np.mean(corr,1), R * R)/(4 * math.pi * delta_R * N*(N - 1)) * L ** 3
    return av_corr, R
