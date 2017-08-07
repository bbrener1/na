#!/usr/bin/env python

import multiprocessing as mlt

import sys

import numpy as np

from sklearn import preprocessing as pre

import itertools

from scipy.stats import pearsonr
from scipy.stats import linregress

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

def parallel_regression(counts):

    slopes = np.zeros((counts.shape[1],counts.shape[1]))

    intercepts = np.zeros((counts.shape[1],counts.shape[1]))

    means = np.mean(counts,axis=0)

    correlations = np.zeros((counts.shape[1],counts.shape[1]))

    pval = np.zeros((counts.shape[1],counts.shape[1]))

    # pool = mlt.Pool(processes=min(20,mlt.cpu_count))
    pool = mlt.Pool(processes=10)

    print "Parallel Regression Started"

    returns = pool.imap_unordered(compact_regression, map(lambda z: (counts[:,z[0]],counts[:,z[1]],z[0],z[1]), [(x, y) for x in range(counts.shape[1]) for y in range(counts.shape[1])] ), chunksize=100)

    for i,c in enumerate(returns):
        slopes[c[1],c[2]] = c[0][0]
        intercepts[c[1],c[2]] = c[0][1]
        correlations[c[1],c[2]] = c[0][2]
        pval[c[1],c[2]] = c[0][3]
        if i%10000==0:
            print i

    np.save("slopes_lin_reg", slopes)
    np.save("intercepts_lin_reg", intercepts)
    np.save("correlations_lin_reg", correlations)
    np.save("pval_lin_reg", pval)
    np.save("means_lin_reg", means)


    return slopes,intercepts,means,correlations,pval

def compact_regression(l):
    # print l[2:]
    # print l[:1]
    # print l[3]
    result = (linregress(l[0],l[1]),l[2],l[3])
    # print len(result[0])
    # print result[0]
    return result

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
    correlations[np.identity(correlations.shape[0],dtype=bool)] = 0
    pval[np.identity(pval.shape[0],dtype=bool)] = .99

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

def partial_correlation(data):

    precision = np.linalg.inv(np.corrcoef(data.T))

    partial_correlation = np.divide(precision,np.tile(precision.diagonal(),(precision.shape[0],1)).T)

    return partial_correlation

def predict_cell(cell, slopes, intercepts, means, correlations, pval, truth = None ):

    raw_predicted = np.multiply(np.tile(cell,(slopes.shape[0],1)).T,slopes) + intercepts

    unadjusted = np.mean(raw_predicted, axis = 0)

    pvalue_derived_weights = np.power(pval,-10)
    pvalue_derived_weights[pvalue_derived_weights > 1000000] = 1000000
    pvalue_derived_weights[pvalue_derived_weights < .000001] = .000001
    # pvalue_derived_weights = np.log10(pvalue_derived_weights)

    pvalue_adjusted = np.average(raw_predicted, axis = 0, weights = pvalue_derived_weights)

    correlation_derived_weights = np.power(1-np.abs(correlations), -1)
    correlation_derived_weights[correlation_derived_weights > 1000] = 1000

    correlation_adjusted = np.average(raw_predicted, axis = 0, weights = correlation_derived_weights)

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
        print "True sum"
        print np.sum(np.abs(truth-means))
        print "Prediction sums"
        print np.sum(np.abs(unadjusted-means))
        print np.sum(np.abs(pvalue_adjusted-means))
        print np.sum(np.abs(correlation_adjusted-means))

        print "\n\n"

    return raw_predicted, pvalue_adjusted, correlation_adjusted

def predict_gene(gene, index,  slopes,intercepts, correlations, pval, truth = None):

    temp = np.zeros((gene.shape[0],1))
    temp[:,0]=gene
    gene = temp

    raw_predicted = np.multiply(np.tile(gene,(1,slopes.shape[0])),np.tile(slopes[index],(gene.shape[0],1)))

def main():

    counts = np.load(sys.argv[1])

    print "Main successful"

    if "solved" in sys.argv:
        slopes = np.load("slopes_lin_reg.npy")
        intercepts = np.load("intercepts_lin_reg.npy")
        means = np.load("means_lin_reg.npy")
        correlations = np.load("correlations_lin_reg.npy")
        pval = np.load("pval_lin_reg.npy")
    else:
        slopes, intercepts, means, correlations, pval = parallel_regression(counts)

    partial = partial_correlation(counts)



    for pick in np.random.randint(counts.shape[0], size=20):

        print pick


        predict_cell(counts[pick,:], slopes, intercepts, means, correlations,partial,truth = counts[pick,:])






if __name__ == "__main__":
    main()



pass
