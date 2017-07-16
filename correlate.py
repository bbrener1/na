#!/usr/bin/env python

import sys
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt


def quick_correlation(observation_matrix, name = None):

    correlation_matrix = np.abs(np.nan_to_num(np.corrcoef(observation_matrix.T)))
    correlation_matrix[correlation_matrix > 10] = 10
    if name != None:
        np.save(name, np.nan_to_num(correlation_matrix))



    print "Quick Correlation Debug"
    print np.sum(correlation_matrix)
    print correlation_matrix.shape
    print "Non-zero edges: "
    print np.sum(correlation_matrix > 0)


    plt.figure()
    plt.hist(correlation_matrix.flatten(), bins=20, log=True)
    plt.savefig("quick_correlation_debug.png")

    return correlation_matrix
