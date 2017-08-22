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



def describe(deviation_matrix, correlation_matrix, cell_identity, description = "./description/"):

    plt.figure()
    plt.hist(deviation_matrix.flatten().T,bins=101)
    # plt.plot(np.linspace(-10,10,200),scipy.stats.norm.pdf(np.linspace(-10,10,200))*deviation_matrix.shape[0]*deviation_matrix.shape[1])
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


    cell_clustering = AgglomerativeClustering(n_clusters=100)
    plt.figure()
    sorted_singly = deviation_matrix[np.argsort(cell_clustering.fit_predict(deviation_matrix))]
    sorted_doubly = sorted_singly.T[np.argsort(np.mean(sorted_singly,axis=0))].T
    plt.imshow(sorted_doubly,cmap="seismic")
    plt.savefig(description + "clustered.png", dpi=300)

    plt.figure()
    sorted_singly = deviation_matrix[np.argsort(cell_clustering.fit_predict(deviation_matrix))]
    print 1656
    print cell_clustering.fit_predict(deviation_matrix).shape
    sorted_doubly = sorted_singly.T[np.argsort(np.var(sorted_singly,axis=0))].T
    plt.imshow(sorted_doubly,cmap="seismic")
    plt.savefig(description + "clustered_variance.png", dpi=300)

    gene_clustering = AgglomerativeClustering(n_clusters=100)
    plt.figure()
    ax = plt.axes()
    sorted_singly = deviation_matrix[np.argsort(cell_clustering.fit_predict(deviation_matrix))]
    print 4773
    print gene_clustering.fit_predict(sorted_singly.T).shape
    sorted_doubly = sorted_singly.T[np.argsort(gene_clustering.fit_predict(sorted_singly.T))]
    sorted_doubly = np.concatenate((sorted_doubly.T,np.ones((sorted_doubly.T.shape[0],1))*-10), axis=1)
    sorted_doubly = np.concatenate((sorted_doubly,np.ones((sorted_doubly.shape[0],1))*10), axis=1)
    im = ax.imshow(sorted_doubly, cmap='seismic')
    plt.title("Residual Expression of Genes In Cells, Clustered Hierarchically")
    plt.xlabel("Genes")
    plt.ylabel("Cells")
    plt.colorbar(im, ax=ax,shrink=.33)
    plt.savefig(description + "doubly_clustered.png",dpi=300)

    print np.max(sorted_doubly)

    plt.figure()
    double_sorting_mask = (cell_identity[:,2]>0)[np.argsort(cell_clustering.fit_predict(deviation_matrix))]
    raw = sorted_doubly[double_sorting_mask]
    minmax = np.concatenate((raw,np.ones((raw.shape[0],1))*10), axis=1)
    minmax = np.concatenate((minmax,np.ones((raw.shape[0],1))*-10), axis=1)
    plt.imshow(minmax,cmap='seismic')
    plt.title("Expression of Genes only in MPP Cells, Same Clustering")
    plt.xlabel("Genes")
    plt.ylabel("Cells")
    plt.savefig(description + "MPPcluster", dpi=300)

    gene_meta_clustering = AgglomerativeClustering(n_clusters=100, linkage='average')
    sorting_indecies = np.argsort(gene_meta_clustering.fit_predict(deviation_matrix.T))
    corr_sort_1 = correlation_matrix[sorting_indecies]
    corr_sort_2 = corr_sort_1.T[sorting_indecies]

    plt.figure()
    plt.imshow(correlation_matrix, cmap='OrRd')
    plt.savefig(description + "correlation.png",dpi=300)

    plt.figure()
    plt.imshow( 1.0/((1-corr_sort_1)**3), cmap='OrRd')
    plt.savefig(description + "coexpression_clusters1.png",dpi=300)

    plt.figure()
    plt.imshow(corr_sort_2 * 1000000, cmap='binary')
    plt.title("Patterns of co-expression among genes")
    plt.savefig(description + "coexpression_clusters2.png",dpi=300)


    strange = sorted_doubly[:,-50:]
    np.savetxt(description+ "strange.txt",strange)
    np.savetxt(description + "indecies_strange.txt",np.argsort(np.var(sorted_singly,axis=0)))



def main():

    prefix = sys.argv[1]

    if len(sys.argv)> 2:
        deviation_matrix = sys.argv[2]
    else:
        deviation_matrix = prefix + "numeric_cons_dev_matrix.npy"

    if len(sys.argv) > 3:
        correlation_matrix = sys.argv[3]
    else:
        correlation_matrix = prefix + "folded_correlation_backup.npy"

    if len(sys.argv) > 4:
        cell_identity = sys.argv[4]
    else:
        cell_identity = prefix + "cell_identity.npy"

    if len(sys.argv) > 5:
        description = sys.argv[5]
    else:
        if prefix[-1] == "/":
            description = "description/"
        else:
            description = "/description/"

    deviation_matrix = matrix_assurance(deviation_matrix)
    correlation_matrix = matrix_assurance(correlation_matrix)
    cell_identity = matrix_assurance(cell_identity)

    if os.path.exists(prefix+description):
        describe(deviation_matrix, correlation_matrix, cell_identity, description = prefix+description)
    else:
        os.makedirs(prefix+description)
        describe(deviation_matrix, correlation_matrix, cell_identity, description = prefix+description)

if __name__ == "__main__":
    main()
