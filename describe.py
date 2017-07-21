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



def describe(deviation_matrix, description = "./description/"):

    plt.figure()
    plt.hist(deviation_matrix.flatten().T,bins=21)
    plt.plot(np.linspace(-10,10,200),scipy.stats.norm.pdf(np.linspace(-10,10,200))*deviation_matrix.shape[0]*deviation_matrix.shape[1])
    plt.title("Frequency of Deviation From Neighbors (Std Deviations)")
    plt.ylabel("Number of Instances")
    plt.xlabel("Deviation (Standard Deviations)")
    plt.savefig(description + "deviation_matrix_hist.png")

    plt.figure()
    plt.hist(np.var(deviation_matrix, axis=0).T)
    plt.title("Distribution of variance by genes")
    plt.xlabel("Variance of Genes")
    plt.ylabel("Number of Genes")
    plt.savefig(description + "meta_variance.png")

    plt.figure()
    plt.imshow(deviation_matrix,cmap='seismic')
    plt.savefig(description + "unsorted_heatmap.png", dpi=300)

    plt.figure()
    sorted_matrix = deviation_matrix.T[np.argsort(np.var(deviation_matrix,axis=0))].T
    plt.imshow(sorted_matrix,cmap='seismic')
    plt.savefig(description + "singly_sorted_heatmap.png", dpi=300)

    # plt.figure()
    # sorted1 = deviation_matrix.T
    # plt.imshow(sorted1,cmap='seismic')
    # plt.savefig(description + "sorted1.png",dpi=300)
    #
    # print deviation_matrix.shape
    # print sorted1.shape
    #
    # plt.figure()
    # print "Starting Indecies"
    # indecies = np.argsort(deviation_matrix,axis = 1)
    # print "Sorted"
    # print indecies.shape

    indecies2 = np.argsort(np.sum(deviation_matrix,axis = 1))

    # print "Sorted sums"
    # indecies2.shape

    sorted2 = deviation_matrix[indecies2]
    # print sorted2.shape
    #
    # print "Picked"
    # plt.imshow(sorted2,cmap='seismic')
    # print "Painted"
    # plt.savefig(description + "sorted2.png",dpi=300)
    # print "Saved"

    plt.figure()
    sorted3 = sorted2.T
    plt.imshow(sorted3,cmap='seismic')
    plt.savefig(description + "doubly_sorted.png",dpi=300)


    cell_clustering = AgglomerativeClustering(n_clusters=15)
    plt.figure()
    sorted_singly = deviation_matrix[np.argsort(cell_clustering.fit_predict(deviation_matrix))]
    sorted_doubly = sorted_singly.T[np.argsort(np.mean(sorted_singly,axis=0))].T
    plt.imshow(sorted_doubly,cmap="seismic")
    plt.savefig(description + "clustered.png", dpi=300)

    plt.figure()
    sorted_singly = deviation_matrix[np.argsort(cell_clustering.fit_predict(deviation_matrix))]
    sorted_doubly = sorted_singly.T[np.argsort(np.var(sorted_singly,axis=0))].T
    plt.imshow(sorted_doubly,cmap="seismic")
    plt.savefig(description + "clustered_variance.png", dpi=300)

    strange = sorted_doubly[:,-50:]
    np.savetxt(description+ "strange.txt",strange)
    np.savetxt(description + "indecies_strange.txt",np.argsort(np.var(sorted_singly,axis=0)))



def main():

    prefix = sys.argv[1]

    if len(sys.argv)> 2:
        deviation_matrix = sys.argv[2]
    else:
        deviation_matrix = prefix + "deviation_matrix.npy"

    if len(sys.argv)> 3:
        description = sys.argv[3]
    else:
        if prefix[-1] == "/":
            description = "description/"
        else:
            description = "/description/"

    deviation_matrix = matrix_assurance(deviation_matrix)

    if os.path.exists(prefix+description):
        describe(deviation_matrix, description = prefix+description)
    else:
        os.makedirs(prefix+description)
        describe(deviation_matrix, description = prefix+description)

if __name__ == "__main__":
    main()
