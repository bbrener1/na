#!/usr/bin/env python
from __future__ import division
import sys
import numpy as np
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt

from sklearn.cluster import AgglomerativeClustering
import networkx as nx



def describe(deviation_matrix,):

    plt.figure()
    plt.hist(deviation_matrix.flatten().T,bins=21)
    plt.plot(np.linspace(-10,10,200),scipy.stats.norm.pdf(np.linspace(-10,10,200))*7000000)
    plt.title("Frequency of Deviation From Neighbors (Std Deviations)")
    plt.ylabel("Number of Instances")
    plt.xlabel("Deviation (Standard Deviations)")
    plt.savefig("deviation_matrix_hist.png")

    plt.figure()
    plt.hist(np.var(deviation_matrix, axis=0).T)
    plt.title("Distribution of variance by genes")
    plt.xlabel("Variance of Genes")
    plt.ylabel("Number of Genes")
    plt.savefig("meta_variance.png")

    plt.figure()
    plt.imshow(deviation_matrix,cmap='bwr')
    plt.savefig("unsorted_heatmap.png", dpi=300)

    plt.figure()
    sorted_matrix = deviation_matrix.T[np.argsort(np.var(deviation_matrix,axis=0))].T
    plt.imshow(sorted_matrix,cmap='bwr')
    plt.savefig("singly_sorted_heatmap.png", dpi=300)

    plt.figure()
    sorted1 = deviation_matrix.T
    plt.imshow(sorted1,cmap='bwr')
    plt.savefig("sorted1.png",dpi=300)

    plt.figure()
    sorted2 = sorted1[np.argsort(deviation_matrix,axis = 0)]
    plt.imshow(sorted2,cmap='bwr')
    plt.savefig("sorted2.png",dpi=300)

    plt.figure()
    sorted3 = sorted2.T
    plt.imshow(sorted3,cmap='bwr')
    plt.savefig("sorted3.png",dpi=300)


    cell_clustering = AgglomerativeClustering(n_clusters=15)
    plt.figure()
    sorted_singly = deviation_matrix[cell_clustering.fit_predict(deviation_matrix)]
    sorted_doubly = sorted_singly.T[np.argsort(np.mean(sorted_singly,axis=0))].T
    plt.imshow(sorted_doubly,cmap="seismic")
    plt.savefig("clustered.png", dpi=300)

    plt.figure()
    sorted_singly = deviation_matrix[cell_clustering.fit_predict(deviation_matrix)]
    sorted_doubly = sorted_singly.T[np.argsort(np.var(sorted_singly,axis=0))].T
    plt.imshow(sorted_doubly,cmap="seismic")
    plt.savefig("clustered_variance.png")
