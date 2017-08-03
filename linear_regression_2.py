#!/usr/bin/env python

import sys

import numpy as np

from sklearn import preprocessing as pre

import itertools

from scipy.stats import pearsonr
from scipy.stats import linregress

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

def linear_regression(counts):

    slopes = np.zeros((counts.shape[1],counts.shape[1]))

    intercepts = np.zeros((counts.shape[1],counts.shape[1]))

    correlations = np.zeros((counts.shape[1],counts.shape[1]))

    pval = np.zeros((counts.shape[1],counts.shape[1]))

    for i in range(counts.shape[1]):
        print i
        for j in range(counts.shape[1]):

            slopes[i,j],intercepts[i,j], correlations[i,j], pval[i,j], _ = linregress(counts[:,i],counts[:,j])

    slopes[np.identity(slopes.shape[0],dtype=bool)] = 0

    means = np.mean(counts,axis=0)

    return slopes, intercepts, means, correlations, pval

    # means = np.mean(counts,axis=0)
    #
    # centered = counts - np.tile(means,(counts.shape[0],1))
    #
    # slopes = np.cov(centered.T)
    #
    # slopes =  np.divide(slopes,np.tile(slopes.diagonal(),(slopes.shape[0],1)).T)
    #
    # intercepts

def predict_cell(cell, slopes, intercepts, means, correlations, pval, truth = None ):

    raw_predicted = np.multiply(np.tile(cell,(slopes.shape[0],1)).T,slopes) + np.intercepts

    unadjusted = np.mean(raw_predicted, axis = 0)

    pvalue_adjusted = np.average(raw_predicted, axis = 0, weights = np.power(pval, -1))

    correlation_adjusted = np.average(raw_predicted, axis = 0, weights = correlations)

    if truth != None:
        print "Truth To Mean"
        print pearsonr(truth,means)
        print "Guesses"
        print "======="

        print "Unadjusted"
        print pearsonr(unadjusted,truth)
        print "Inverse P Value Adjusted"
        print pearsonr(pvalue_adjusted,truth)
        print "Correlation Adjusted"
        print pearsonr(correlation_adjusted,truth)

        print "Centered Data"
        print "======"
        print pearsonr(unadjusted-means,truth-means)
        print pearsonr(pvalue_adjusted-means,truth-means)
        print pearsonr(correlation_adjusted-means,truth-means)

        print "\n\n"

    return raw_predicted, pvalue_adjusted, correlation_adjusted

def predict_gene(gene, index,  slopes,intercepts, correlations, pval, truth = None):

    temp = np.zeros((gene.shape[0],1))
    temp[:,0]=gene
    gene = temp

    raw_predicted = np.multiply(np.tile(gene,(1,slopes.shape[0])),np.tile(slopes[index],(gene.shape[0],1)))

def main():

    counts = np.load(sys.argv[1])

    slopes, intercepts, means, correlations, pval = linear_regression(counts)

    partial = partial_correlation(counts)



    for pick in np.random.randint(counts.shape[0], size=100):

        print pick

        predict_cell(counts[pick,:], slopes, intercepts, means, correlations,partial,truth = counts[pick,:])






if __name__ == "__main__":
    main()



pass
