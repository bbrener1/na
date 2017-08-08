#!/usr/bin/env python
#
# covariance
#
# imp_obs = predict_drops(obs, covariance)
#
# covariance = deviation_matrix(imp_obs)

import numpy as np

import stripped_regression

from sklearn.decomposition import PCA
import scipy.spatial.distance as spt

def predict_gene_counts(cell_matrix, iteration, distances, linear_model):

    predicted_counts = np.zeros(cell_matrix.shape)

    for i, cell in enumerate(cell_matrix):
        predicted_counts[i] = linear_model.predict_cell(cell)

    # distance_weights =

def main():

    counts = np.load(prefix+"/counts.npy")

    linear_model = stripped_regression.stripped_regression(counts, prefix = prefix)

    linear_model.test()
