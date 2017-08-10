#!/usr/bin/env python
from __future__ import division
import sys
import os
import numpy as np
import scipy

import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt

from matrix_assurance import *

from sklearn.cluster import AgglomerativeClustering
import networkx as nx



def variance_discrepancy(deviation_matrix, counts, header):

    count_variance = np.var(counts, axis=0)

    deviation_variance = np.var(deviation_matrix, axis=0)

    mean_deviation = np.mean(deviation_matrix, axis = 0)

    variance_index = np.divide(np.sqrt(count_variance),np.power(deviation_variance,3))

    mask = np.argsort(variance_index) > 4400

    np.savetxt("TF_candidates_deviation.txt", header[mask],fmt='%s')
    #
    # plt.figure()
    # plt.bar(np.arange(len(deviation_variance)),np.sort(deviation_variance))
    # plt.xticks(np.arange(len(deviation_variance)), np.arange(len(deviation_variance))[np.argsort(deviation_variance)].astype(dtype=str))
    #
    # plt.savefig("dev_variance.png")
    #
    # plt.figure()
    # plt.bar(np.arange(len(variance_index)),np.sort(variance_index))
    # plt.xticks(np.arange(len(variance_index)), np.arange(len(variance_index))[np.argsort(variance_index)].astype(dtype=str))
    # plt.savefig("variance_index.png")
    #
    # plt.figure()
    # plt.bar(np.arange(len(mean_deviation)), np.sort(mean_deviation),tick_label=np.arange(len(mean_deviation))[np.argsort(mean_deviation)].astype(dtype=str))
    # # plt.xticks(np.arange(len(mean_deviation)), np.arange(len(mean_deviation))[np.argsort(mean_deviation)].astype(dtype=str))
    # plt.savefig("mean_deviation.png")

    # variance_model = scipy.stats.f()
    #
    # variance_model.pdf()

def high_degrees(correlation_matrix, threshold, header, filename = "TF_candidates_degree.txt"):

    degrees = np.sum(correlation_matrix > threshold,axis=0)

    degree_indecies = np.arange(degrees.shape[0])[np.argsort(degrees) > (degrees.shape[0] - degrees.shape[0]/10.0)]

    np.savetxt(filename, header[degree_indecies],fmt='%s')

counts = np.load(sys.argv[1])

deviation_matrix = np.load(sys.argv[2])

correlation = np.load(sys.argv[3])

header = np.load(sys.argv[4])

threshold = float(sys.argv[5])

high_degrees(correlation, threshold, header)

variance_discrepancy(deviation_matrix, counts, header)
