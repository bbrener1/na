#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import numpy.random as nprnd

#
# network = np.zeros((30000,30000))
#
# degrees = np.sum(network, axis = 1)
#
# x = 0
#
# while x < 80000:
#     degree_probability = (degrees + 1)/np.sum(degrees + 1)
#     i = x%30000
#     j = nprnd.choice(np.arange(30000),1,p = degree_probability)
#     if i == j:
#         continue
#     if network[i,j] > 0 or network[j,i] > 0:
#         continue
#     network[i,j] += 1
#     network[j,i] += 1
#     degrees[i] += 1
#     degrees[j] += 1
#
#     x+=1
#
#     if x%1000 == 0:
#         print x
#
# print np.sum(network)
# print np.sum(degrees)
#
# plt.figure()
# plt.hist(degrees, log=True)
# plt.show()

x = 0

for i in np.arange(1,30000):
    x += 30000.0/(i**2)

print x
