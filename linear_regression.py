#!/usr/bin/env python

import sys

import numpy as np

import itertools

from scipy.stats import pearsonr

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

def linear_regression(observations):

    means = np.mean(observations,axis=0)

    centered = observations - np.tile(means, (observations.shape[0],1))

    slopes = np.cov(observations.T)

    print means.shape
    print means.T.shape
    print slopes.shape

    intercepts = np.tile(means,(observations.shape[1],1)).T - np.multiply(np.tile(means,(observations.shape[1],1)),slopes)


    print intercepts.shape

    correlation = np.corrcoef(observations.T)

    # correlation = np.zeros((observations.shape[1],observations.shape[1]))
    #
    # for i in range(observations.shape[1]):
    #     for j in range(i,observations.shape[1]):
    #         if j%500 == 0:
    #             print i
    #             print j
    #         corr = pearsonr(observations[:,i],observations[:,j])
    #         # correlation[i,j],correlation_p_matrix[i,j] = np.abs(corr[0]),corr[1]
    #         correlation[i,j]= np.abs(corr[0])
    #         correlation[j,i] = correlation[i,j]


    return slopes, intercepts, means, correlation

def partial_correlation(data):

    # covariance = np.cov(data.T)
    # precision = np.linalg.inv(covariance)

    precision = np.linalg.inv(np.corrcoef(data.T))

    # plt.figure()
    # plt.hist(precision.diagonal(),bins=50)
    # plt.savefig("partial.png")

    # partial_correlation_adjuster = np.multiply(np.tile(precision.diagonal(),(precision.shape[0],1)),np.tile(precision.diagonal(),(precision.shape[0],1)).T)

    # print np.min(precision)
    # print np.min(partial_correlation_adjuster)
    # print np.min(precision.diagonal())
    # print precision.diagonal()[:100]
    # print np.min(covariance.diagonal())
    # print covariance.diagonal()[:100]

    # partial_correlation_adjuster = np.sqrt(np.abs(partial_correlation_adjuster))


    partial_correlation = np.divide(precision,np.tile(precision.diagonal(),(precision.shape[0],1)).T)
    print precision[233,233]
    print np.tile(precision.diagonal(),(precision.shape[0],1)).T[233,10]

    # raise ValueError

    return partial_correlation

def predict(data, true_values, slopes, intercepts, means, index,correlation, partial):

    adj_slopes = np.multiply(slopes,correlation)
    adj_slopes[np.identity(slopes.shape[0],dtype=bool)] = 0

    weights = np.mean(adj_slopes,axis = 0)

    print "Data shape"
    print data.shape
    print np.dot(data,slopes).shape

    prediction1 = np.divide(np.dot(data,adj_slopes),weights) + intercepts[index]

    print prediction1.shape
    print true_values.shape

    print "Prediction 1"
    print pearsonr(prediction1, true_values)

    zero_mask = data != 0

    prediction2 = np.dot(data[zero_mask],adj_slopes[zero_mask].T[zero_mask].T)

    prediction2 = np.divide(prediction2,weights[zero_mask]) + intercepts[index][zero_mask]
    #
    # print prediction2.shape
    # print prediction2[:10]
    # print intercepts[index][np.logical_not(zero_mask)].shape
    # print weights[np.logical_not(zero_mask)].shape

    # print prediction2.shape
    # print prediction2[:10]
    # print true_values[np.logical_not(zero_mask)][:10]

    print "Prediction 2"
    print pearsonr(prediction2, true_values[zero_mask])

    part_slopes = np.multiply(slopes,partial)
    part_weights = np.mean(part_slopes,axis=0)

    prediction2a = np.dot(data[zero_mask],part_slopes[zero_mask].T[zero_mask].T)
    prediction2a = np.divide(prediction2,part_weights[zero_mask]) 
    # weights = np.ones(weights.shape)
    #
    # prediction2a = np.dot(data[zero_mask],adj_slopes[zero_mask].T[zero_mask].T)
    # prediction2a = np.divide(prediction2a,weights[zero_mask]) + intercepts[index][zero_mask]

    print "Prediction 2 A"
    print pearsonr(prediction2a, true_values[zero_mask])

    fold_setting = 5

    folds = list(itertools.combinations(range(fold_setting),int(fold_setting-1)))

    integer_mask = np.random.randint(fold_setting, size=data.shape[0])

    for fold in folds:

        fold_mask = np.zeros(data.shape[0],dtype=bool)

        for element in fold:
            fold_mask = np.logical_or(fold_mask,integer_mask == element)

        final_mask = np.logical_and(fold_mask, zero_mask)

        unadjusted = np.dot(data[final_mask],adj_slopes[final_mask][:,np.logical_not(final_mask)])

        final = np.divide(unadjusted,weights[np.logical_not(final_mask)]) + intercepts[index][np.logical_not(final_mask)]

        print "Prediction 3"
        print pearsonr(final,true_values[np.logical_not(final_mask)])
        print np.sum(final_mask)


def main():

    counts = np.load(sys.argv[1])

    slopes, intercepts, means, correlations = linear_regression(counts)

    partial = partial_correlation(counts)

    for pick in np.random.randint(counts.shape[0], size=50):

        print pick

        predict(counts[pick,:],counts[pick,:], slopes, intercepts, means, pick, correlations, partial)






if __name__ == "__main__":
    main()



pass
