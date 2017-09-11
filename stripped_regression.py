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

import check_hash as chk
from matrix_assurance import matrix_assurance

def compact_regression(l):
    result = (linregress(l[0],l[1]),l[2],l[3])
    return result


class stripped_regression:


    def __init__(self, counts, solved=False, prefix=""):

        self.counts = matrix_assurance(counts)
        self.prefix = prefix

        print "Main successful"

        if solved or chk.check_hash(counts, "slopes_lin_reg.npy", prefix=self.prefix):
            self.slopes = np.load(prefix + "slopes_lin_reg.npy")
            self.intercepts = np.load(prefix + "intercepts_lin_reg.npy")
            self.means = np.load(prefix + "means_lin_reg.npy")
            self.correlations = np.load(prefix + "correlations_lin_reg.npy")
            self.pval = np.load(prefix + "pval_lin_reg.npy")
        else:
            self.slopes, self.intercepts, self.means, self.correlations, self.pval = self.parallel_regression(self.counts)

        self.imputed_counts = None

        self.partial = self.partial_correlation(self.counts)


    def test(self):

        for pick in np.random.randint(self.counts.shape[0], size=20):

            print pick

            self.predict_cell(self.counts[pick,:], truth = self.counts[pick,:])


    def parallel_regression(self, counts = None):

        if str(counts) == "None":
            counts = self.counts

        slopes = np.zeros((counts.shape[1],counts.shape[1]))

        intercepts = np.zeros((counts.shape[1],counts.shape[1]))

        means = np.mean(counts,axis=0)

        correlations = np.zeros((counts.shape[1],counts.shape[1]))

        pval = np.zeros((counts.shape[1],counts.shape[1]))

        pool = mlt.Pool(processes=min(mlt.cpu_count()-2,20))
        # pool = mlt.Pool(processes=10)

        print "Parallel Regression Started"

        returns = pool.imap_unordered( compact_regression, map(lambda z: (counts[:,z[0]],counts[:,z[1]],z[0],z[1]), [(x, y) for x in range(counts.shape[1]) for y in range(counts.shape[1])] ), chunksize=100)

        for i,c in enumerate(returns):
            slopes[c[1],c[2]] = c[0][0]
            intercepts[c[1],c[2]] = c[0][1]
            correlations[c[1],c[2]] = c[0][2]
            pval[c[1],c[2]] = c[0][3]
            if i%1000000==0:
                print i

        slopes[np.identity(slopes.shape[0],dtype=bool)] = 0

        correlations[np.identity(correlations.shape[0],dtype=bool)] = 0

        pval[np.identity(pval.shape[0],dtype=bool)] = .99

        np.save(self.prefix + "slopes_lin_reg", slopes)
        np.save(self.prefix + "intercepts_lin_reg", intercepts)
        np.save(self.prefix + "correlations_lin_reg", correlations)
        np.save(self.prefix + "pval_lin_reg", pval)
        np.save(self.prefix + "means_lin_reg", means)

        chk.write_hash(counts, "slopes_lin_reg.npy", prefix=self.prefix)

        self.slopes = slopes
        self.intercepts = intercepts
        self.correlations = correlations
        self.pval = pval
        self.means = means

        return slopes,intercepts,means,correlations,pval



    def partial_correlation(self, data):

        precision = np.linalg.inv(np.corrcoef(data.T))

        partial_correlation = np.divide(precision,np.tile(precision.diagonal(),(precision.shape[0],1)).T)

        return partial_correlation

    def predict_cell(self, cell, index = False, truth = None, verbose = True, masked = False, mask = None):

        if index:
            cell = self.counts[cell,:]

        ## In raw predicted, i,j is the value of gene j in the target cell, predicted based on gene i

        ## Column j of prediction matrix is all predictions of gene j in the cell. Mean of column j is
        ## mean of predictions

        ## Row i is cell state predicted by individual gene.

        raw_predicted = np.multiply(np.tile(cell,(self.slopes.shape[0],1)).T,self.slopes) + self.intercepts

        unadjusted = np.mean(raw_predicted, axis = 0)

        correlation_derived_weights = np.power(1-np.abs(self.correlations), -1)
        correlation_derived_weights[correlation_derived_weights > 1000] = 1000

        if masked:
            if str(mask) == "None":
                mask = cell == 0

            correlation_derived_weights[mask] = np.zeros(correlation_derived_weights.shape[1])

        correlation_adjusted = np.average(raw_predicted, axis = 0, weights = correlation_derived_weights)

        partial_corr_weights = np.power(1-np.abs(self.partial), -1)
        partial_corr_weights[partial_corr_weights > 1000] = 1000

        partial_adjusted = np.average(raw_predicted, axis = 0, weights = partial_corr_weights)

        # print "Computed correlation adjusted values, computing dropouts:"

        # Predictied value j given dropout of i

        ## Compute influence of each individual prediction on the weighted average:

        total_weights = np.sum(correlation_derived_weights, axis = 0)

        relative_weights = np.divide(correlation_derived_weights, np.tile(total_weights, (correlation_derived_weights.shape[0],1)))

        influence = np.multiply(raw_predicted, relative_weights)

        dropout_adjusted = np.tile(correlation_adjusted,(influence.shape[0],1)) - influence


        if str(truth) != "None" and verbose:
            print "Truth To Mean"
            print pearsonr(truth,self.means)
            print "Guess"
            print "======="

            print "Unadjusted"
            print pearsonr(unadjusted,truth)
            print "Correlation Adjusted"
            print pearsonr(correlation_adjusted,truth)
            print "Partial Correlation Adjusted"
            print pearsonr(partial_adjusted, truth)

            print "======"
            print "Centered Data"
            # print pearsonr(unadjusted-means,truth-means)
            print pearsonr(correlation_adjusted-self.means,truth-self.means)
            print "======"
            print "True sum"
            print np.sum(np.abs(truth-self.means))
            print "Prediction sum"
            # print np.sum(np.abs(unadjusted-means))
            print np.sum(np.abs(correlation_adjusted-self.means))

            print "\n\n"

        return correlation_adjusted, raw_predicted, correlation_derived_weights, dropout_adjusted

    # def masked_predict(self, cell, index = False, truth = None , mask= None, verbose = True):
    #
    #     if index:
    #         cell = self.counts[cell,:]
    #
    #     raw_predicted = np.multiply(np.tile(cell,(self.slopes.shape[0],1)).T,self.slopes) + self.intercepts
    #
    #     unadjusted = np.mean(raw_predicted, axis = 0)
    #
    #     correlation_derived_weights = np.power(1-np.abs(self.correlations), -1)
    #     correlation_derived_weights[correlation_derived_weights > 1000] = 1000
    #
    #     # Ignores any values that are 0 if the mask is not specified
    #
    #
    #
    #     correlation_adjusted = np.average(raw_predicted, axis = 0, weights = correlation_derived_weights)
    #
    #     # Predictied value i given dropout of j
    #
    #     dropout_adjusted = np.zeros((raw_predicted.shape[1],raw_predicted.shape[1]))
    #
    #     for i, column in self.counts.T:
    #         drop_weights = np.copy(correlation_derived_weights)
    #         drop_weights[:,i] = np.zeros(drop_weights.shape[1])
    #         dropout_adjusted[i] = np.average(raw_predicted, axis = 0, weights = drop_weights)
    #
    #
    #     if truth != None and verbose:
    #         print "Truth To Mean"
    #         print pearsonr(truth,means)
    #         print "Guess"
    #         print "======="
    #
    #         # print "Unadjusted"
    #         # print pearsonr(unadjusted,truth)
    #         print "Correlation Adjusted"
    #         print pearsonr(correlation_adjusted,truth)
    #
    #         print "Centered Data"
    #         print "======"
    #         # print pearsonr(unadjusted-means,truth-means)
    #         print pearsonr(correlation_adjusted-means,truth-means)
    #         print "True sum"
    #         print np.sum(np.abs(truth-means))
    #         print "Prediction sum"
    #         # print np.sum(np.abs(unadjusted-means))
    #         print np.sum(np.abs(correlation_adjusted-means))
    #
    #         print "\n\n"
    #
    #     return correlation_adjusted, raw_predicted, correlation_derived_weights, dropout_adjusted

    # def predict_gene(gene, index,  slopes,intercepts, correlations, pval, truth = None):
    #
    #     temp = np.zeros((gene.shape[0],1))
    #     temp[:,0]=gene
    #     gene = temp
    #
    #     raw_predicted = np.multiply(np.tile(gene,(1,slopes.shape[0])),np.tile(slopes[index],(gene.shape[0],1)))

    # def main():
    #
    #     counts = np.load(sys.argv[1])
    #
    #     print "Main successful"
    #
    #     if "solved" in sys.argv:
    #         slopes = np.load("slopes_lin_reg.npy")
    #         intercepts = np.load("intercepts_lin_reg.npy")
    #         means = np.load("means_lin_reg.npy")
    #         correlations = np.load("correlations_lin_reg.npy")
    #         pval = np.load("pval_lin_reg.npy")
    #     else:
    #         slopes, intercepts, means, correlations, pval = parallel_regression(counts)
    #
    #     partial = partial_correlation(counts)
    #
    #
    #
    #     for pick in np.random.randint(counts.shape[0], size=20):
    #
    #         print pick
    #
    #
    #         predict_cell(counts[pick,:], slopes, intercepts, means, correlations,partial,truth = counts[pick,:])
