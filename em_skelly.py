#!/usr/bin/env python
#
# covariance
#
# imp_obs = predict_drops(obs, covariance)
#
# covariance = deviation_matrix(imp_obs)

import numpy as np
import sys

import stripped_regression

from sklearn.decomposition import PCA
import scipy.spatial.distance as spt

def predict_gene_counts(cell_matrix, iteration, distances, linear_model):

    predicted_counts = np.zeros(cell_matrix.shape)

    for i, cell in enumerate(cell_matrix):
        predicted_counts[i] = linear_model.predict_cell(cell)

    distance_weights = np.divide(1, np.power((1+iteration/50), distances))

def main():

    prefix = sys.argv[1]

    print "Initiating"

    # counts = np.load(prefix+"/counts.npy")

    counts = prefix + "/counts.npy"

    print "Counts loaded:"

    # print counts.shape

    linear_model = stripped_regression.stripped_regression(counts, solved = "solved" in sys.argv, prefix = prefix)

    print "Model built"

    linear_model.test()


    dist_model = PCA(n_components=50)
    dist_interm = dist_model.fit_transform(counts)
    dist = spt.squareform(spt.pdist(dist_interm))

    for gene in counts












if __name__ == "__main__":
    main()
