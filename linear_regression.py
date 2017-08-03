#!/usr/bin/env python

import sys

import numpy as np

from sklearn import preprocessing as pre

import itertools

from scipy.stats import pearsonr

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

def linear_regression(observations):

    means = np.mean(observations,axis=0)

    centered = observations - np.tile(means, (observations.shape[0],1))

    slopes = np.cov(observations.T)
    #
    # print slopes[50,50]
    # print slopes[50,51]
    # print slopes.diagonal()
    # print np.tile(slopes.diagonal(),(slopes.shape[0],1)).shape
    # print np.tile(slopes.diagonal(),(slopes.shape[0],1))[:10,:10]

    # print np.tile(slopes.diagonal(),(slopes.shape[0],1))[50,10]
    # print np.tile(slopes.diagonal(),(slopes.shape[0],1))[50,11]
    #
    # print np.tile(slopes.diagonal(),(slopes.shape[0],1))[10,50]
    # print np.tile(slopes.diagonal(),(slopes.shape[0],1))[100,50]


    slopes =  np.divide(slopes,np.tile(slopes.diagonal(),(slopes.shape[0],1)).T)




    print means.shape
    print means.T.shape
    print slopes.shape

    intercepts = np.tile(means,(observations.shape[1],1)) - np.multiply(np.tile(means,(observations.shape[1],1)).T,slopes)


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

# def predict(data, true_values, slopes, intercepts, means, index,correlation, partial):
def predict(data, true_values, slopes, intercepts, means, correlation):

    # slopes[np.identity(slopes.shape[0],dtype=bool)] = 0

    zero_mask = data != 0


    correlation = np.copy(correlation)
    correlation[zero_mask] = 0
    correlation = np.multiply(correlation,np.logical_not(np.identity(correlation.shape[0],dtype=bool)))

    #
    #
    # unweighted = np.multiply(np.tile(data,(slopes.shape[0],1)).T,slopes) + intercepts
    #
    # weighted = np.multiply(unweighted,np.abs(correlation))
    #
    # # prediction1 = np.mean(unweighted,axis=0)
    #
    # prediction1 = np.divide(np.sum(weighted,axis=0),np.sum(np.abs(correlation),axis=0))
    #
    #
    #
    # print "Prediction 1"
    # print pearsonr(prediction1, true_values)
    #
    # print "Prediction 2"
    #
    # print "Quality of template"
    # sec_zero = true_values != 0
    # print pearsonr(data[np.logical_and(zero_mask,sec_zero)],true_values[np.logical_and(zero_mask,sec_zero)])

    print "Better than noise?"
    print "=================="
    print pearsonr(data, means)

    correlation = np.power(correlation,3)

    centered = data - means
    unweighted_centered = np.multiply(np.tile(centered,(slopes.shape[0],1)).T,slopes)
    weighted_centered = np.multiply(unweighted_centered,correlation)

    centered_prediction = np.divide(np.sum(weighted_centered,axis=0),np.sum(np.abs(correlation),axis=0))

    prediction2 = centered_prediction + means

    print pearsonr(prediction2, means)
    print pearsonr(prediction2, true_values)

    print "Centering"
    print "=================="
    print list(centered[100:110])
    print list(centered_prediction[100:110])

    print "Misc Correlations"
    print "=================="
    print np.sum(centered)
    print np.sum(np.abs(centered))
    print np.sum(centered_prediction)
    print np.sum(np.abs(centered_prediction))
    print pearsonr(centered,centered_prediction)
    print "\n\n"

    return prediction2

def main():

    counts = np.load(sys.argv[1])

    # counts = pre.scale(counts)

    slopes, intercepts, means, correlations = linear_regression(counts)

    # partial = partial_correlation(counts)



    for pick in np.random.randint(counts.shape[0], size=5):

        print pick


        predict(counts[pick,:],counts[pick,:], slopes, intercepts, means, correlations)






if __name__ == "__main__":
    main()



pass
