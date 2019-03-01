#This does work for me, but I get weird errors when I try to adapt this to a function which is called when running the main program.

import numpy as np
from itertools import product
import matplotlib.pyplot as plt
L=10
N=108
n= int(np.rint(np.power(N/4, 1/3)))
actualN = int(np.power(n, 3)*4)
b=L/n
lat = np.zeros([int(actualN/4), 4, 3])
cell=np.asarray([[0,0,0],[0, 1, 1],[1, 0, 1],[1, 1, 0]])

for c in product(range(0, n), repeat=3):
    lat[c[2]+c[1]*n+c[0]*n*n,:,:] = b/2*cell+[c[0]*b, c[1]*b, c[2]*b]-(L)/2


x = np.reshape(lat, (actualN,3)) 
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x[:,0],x[:,1],x[:,2])
plt.show()
print(n)
print(actualN)
print(x)
