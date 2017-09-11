#!/usr/bin/env python
#
# covariance
#
# imp_obs = predict_drops(obs, covariance)
#
# covariance = deviation_matrix(imp_obs)

import numpy as np
import sys

import manual_stripped_regression

from sklearn.decomposition import PCA
import scipy.spatial.distance as spt
from scipy.stats import pearsonr

import consensus_folded_deviation as cfd

def main():

    prefix = sys.argv[1]

    print "Initiating"

    # counts = np.load(prefix+"/counts.npy")

    counts = np.load(prefix + "/counts.npy")

    print "Counts loaded:"

    print counts.shape

    linear_model = manual_stripped_regression.stripped_regression(counts, solved = "solved" in sys.argv, prefix = prefix)

    print "Model built"

    linear_model.test()

    print "Testing initial self-predictive power."

    mean_matrix = np.tile(linear_model.means,(counts.shape[0],1))

    print "Masked prediction:"

    masked_imputed = np.zeros(counts.shape)

    for i, cell in enumerate(counts):

        masked_imputed[i] = linear_model.predict_cell(cell, verbose = True, masked = True)[0]

        print masked_imputed[i,:10]

        if i%100 == 0:
            print i

    print "Masked prediction self-correlation:"
    print pearsonr(counts[counts > 0], masked_imputed[counts > 0])

    print "Count correlation to means:"
    print pearsonr(counts[counts > 0], mean_matrix[counts > 0])

    print "Naive prediction:"

    naive_imputed = np.zeros(counts.shape)

    for i, cell in enumerate(counts):

        naive_imputed[i] = linear_model.predict_cell(cell, verbose = False)[0]

        if i%100 == 0:
            print i

    print "Naive linear model self-correlation:"
    print pearsonr(counts[counts > 0], naive_imputed[counts > 0])
    print "Count correlation to means:"
    print pearsonr(counts[counts > 0], mean_matrix[counts > 0])


    print "Sequential naive prediction:"
    second_naive = np.zeros(counts.shape)

    for i, cell in enumerate(naive_imputed):

        second_naive = linear_model.predict_cell(cell, verbose = False)[0]

        if i%100 == 0:
            print i

    print "Sequential naive prediction self-correlation:"
    print pearsonr(counts, second_naive)
    print "Sequential naive prediction (non-zero only):"
    print pearsonr(counts[counts > 0], second_naive[counts > 0])

    print "Computing deviation means:"

    deviation_medians = np.median(cfd.folded_deviation_matrix(counts)[3], axis=2)

    print "Prediction of non-zero values through deviation means:"
    print pearsonr(counts[counts > 0],deviation_medians[counts > 0])
    print "Predictions of all values through deviation means:"
    print pearsonr(counts, deviation_medians)

    # dist_model = PCA(n_components=50)
    # dist_interm = dist_model.fit_transform(counts)
    # dist = spt.squareform(spt.pdist(dist_interm))




def predict_by_neighbor(counts, distances, neighbors, model):

    output.write("Computing by neighbors:\n")

    for i, row in enumerate(neighbor_mean_matrix):

        integer_mask = np.random.randint(fold, size=neighbor_setting)

        neighborhood = counts[np.argsort(distance_matrix[i]) < neighbor_setting]

        for j, cycle in enumerate(folds):

            fold_mask = np.zeros(neighborhood.shape[0],dtype=bool)
            for element in cycle:
                fold_mask = np.logical_or(fold_mask,integer_mask == element)

            # print cycle
            # print integer_mask
            # print fold_mask
            # print neighbor_mean_matrix.shape
            # print std_dev_matrix.shape
            # print np.mean(neighborhood[fold_mask],axis=0)

            neighbor_mean_matrix[i,:,j] = np.mean(neighborhood[fold_mask],axis=0)

            std_dev_matrix[i,:,j] = np.std(neighborhood[fold_mask],axis=0)



        if i%100 == 0:
            output.write(str(i) + "\n")
            print i

# def predict_gene_counts(cell_matrix, iteration, distances, linear_model):
#
#     predicted_counts = np.zeros(cell_matrix.shape)
#
#     for i, cell in enumerate(cell_matrix):
#         predicted_counts[i] = linear_model.predict_cell(cell)
#
#     distance_weights = np.divide(1, np.power((1+iteration/50), distances))
#
#     np.average()










if __name__ == "__main__":
    main()
