import numpy as np
from scipy.special import roots_legendre
import sys

nodes = sys.argv[1]
mu_min = sys.argv[2]
mu_max = sys.argv[3]

nodes = int(nodes)
a_new = float(mu_min)
b_new = float(mu_max)

print(nodes,a_new,b_new)

nodes_std, weights_std = roots_legendre(nodes)

# Transform nodes to the new interval
nodes_new = 0.5 * (b_new - a_new) * nodes_std + 0.5 * (a_new + b_new)

# Transform weights for the new interval
weights_new = 0.5 * (b_new - a_new) * weights_std

print(nodes_new, weights_new)
np.savetxt('grids/mu.dat', np.column_stack((nodes_new, weights_new)), fmt='%.6f')








