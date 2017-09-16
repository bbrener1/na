#!/usr/bin/env python

import multiprocessing as mlt
from multiprocessing.managers import BaseManager, BaseProxy

from itertools import imap

import sys

import numpy as np

from sklearn import preprocessing as pre

import itertools

from scipy.stats import pearsonr
from scipy.stats import linregress
from scipy.stats.mstats import theilslopes

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

import check_hash as chk
from matrix_assurance import matrix_assurance

    # Wrapper function for prediction for the purposes of multiprocessing:

def compact_prediction(l):
    result = l[0].predict_cell(l[1], verbose=False, masked=True)[0]
    if l[2]%100==0:
        print l[2]
    return result

def compact_regression(l):
    masked = l[4]
    mask = np.ones(l[0].shape, dtype=bool)
    if masked:
        mask = np.logical_and(l[0] > 0,l[1] > 0)
    result = (linregress(l[0][mask],l[1][mask]),l[2],l[3])
    return result

def compact_ts_est(l):
    masked = l[4]
    try:
        mask = np.ones(l[0].shape, dtype=bool)
        if masked:
            mask = np.logical_and(l[0] > 0,l[1] > 0)
            if np.sum(mask) == 0:
                print "Empty mask, TS Gets to return statement"
                raw_input("Keep going?")
                return((0,0,0,0),l[2],l[3])
        result = (theilslopes(l[0][mask],l[1][mask]),l[2],l[3])
    except:
        result = ((0,0,0,0),l[2],l[3])
    return result


class stripped_regression:


    def __init__(self, counts, solved=False, prefix="", method='ols',masking=False,process_limit = False):

        self.counts = matrix_assurance(counts)
        self.prefix = prefix
        self.process_limit = process_limit

        print "Main successful"

        if solved or chk.check_hash(counts, "slopes_lin_reg.npy", prefix=self.prefix):
            self.slopes = np.load(prefix + "slopes_lin_reg.npy")
            self.intercepts = np.load(prefix + "intercepts_lin_reg.npy")
            self.means = np.load(prefix + "means_lin_reg.npy")
            self.correlations = np.load(prefix + "correlations_lin_reg.npy")
            self.pval = np.load(prefix + "pval_lin_reg.npy")
        else:
            self.slopes, self.intercepts, self.means, self.correlations, self.pval = self.parallel_regression(self.counts, method=method, process_limit=process_limit,masking = masking)

        self.imputed_counts = None

        self.partial = self.partial_correlation(self.counts)


    def test(self):

        for pick in np.random.randint(self.counts.shape[0], size=20):

            print pick

            self.predict_cell(self.counts[pick,:], truth = self.counts[pick,:], masked=True)


    def parallel_regression(self, counts = None, method = 'ols', masking = False , process_limit = False):

        if str(counts) == "None":
            counts = self.counts

        slopes = np.zeros((counts.shape[1],counts.shape[1]))

        intercepts = np.zeros((counts.shape[1],counts.shape[1]))

        means = np.mean(counts,axis=0)

        correlations = np.zeros((counts.shape[1],counts.shape[1]))

        pval = np.zeros((counts.shape[1],counts.shape[1]))

        if process_limit:
            pool = mlt.Pool(processes=min(mlt.cpu_count()-2,int(process_limit)))
        else:
            pool = mlt.Pool(processes=mlt.cpu_count()-2)
        # pool = mlt.Pool(processes=10)

        print "Processors detected:"
        print pool._processes
        print "Parallel Regression Started"

        if method == 'ols':

            returns = pool.imap_unordered( compact_regression, imap(lambda z: (counts[:,z[0]],counts[:,z[1]],z[0],z[1],masking), [(x, y) for x in range(counts.shape[1]) for y in range(counts.shape[1])] ), chunksize=100)

            for i,c in enumerate(returns):
                slopes[c[1],c[2]] = c[0][0]
                intercepts[c[1],c[2]] = c[0][1]
                correlations[c[1],c[2]] = c[0][2]
                pval[c[1],c[2]] = c[0][3]
                if i%10000==0:
                    print i

        if method == 'theil_sen':

            returns = pool.imap_unordered( compact_ts_est, map(lambda z: (counts[:,z[0]],counts[:,z[1]],z[0],z[1],masking), [(x, y) for x in range(counts.shape[1]) for y in range(counts.shape[1])] ), chunksize=100)


            for i,c in enumerate(returns):
                slopes[c[1],c[2]] = c[0][0]
                intercepts[c[1],c[2]] = c[0][1]
                if i%10000==0:
                    print i

            correlations = np.corrcoef(counts.T)

        if method != "ols" and method != "theil_sen":
            raise AttributeError("Not a legal estimator selection, use 'ols' or 'theil_sen'")


        # slopes[np.identity(slopes.shape[0],dtype=bool)] = 0

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

    def predict_cell(self, cell, index = False, truth = None, verbose = True, masked = False, mask = None, dropout = False):

        if index:
            cell = self.counts[cell,:]

        ## In raw predicted, i,j is the value of gene j in the target cell, predicted based on gene i

        ## Column j of prediction matrix is all predictions of gene j in the cell. Mean of column j is
        ## mean of predictions

        ## Row i is cell state predicted by individual gene.

        raw_predicted = np.multiply(np.tile(cell,(self.slopes.shape[0],1)).T,self.slopes) + self.intercepts

        unadjusted = np.mean(raw_predicted, axis = 0)

        correlation_derived_weights = np.power(1-np.abs(self.correlations), -1)
        correlation_derived_weights[correlation_derived_weights > 10000] = 10000


        if masked:
            if str(mask) == "None":
                mask = cell == 0

            correlation_derived_weights[mask] = np.zeros(correlation_derived_weights.shape[1])

        correlation_adjusted = np.average(raw_predicted, axis = 0, weights = correlation_derived_weights)

        slope_correlation_weights = np.power(1-np.sqrt(np.abs(np.multiply(self.correlations,self.slopes))),-1)
        slope_correlation_weights[slope_correlation_weights > 10000] = 10000

        slope_corr_adjusted = np.average(raw_predicted, axis=0, weights = slope_correlation_weights)

        # print "Computed correlation adjusted values, computing dropouts:"

        # Predictied value j given dropout of i

        ## Compute influence of each individual prediction on the weighted average:

        if dropout:
            total_weights = np.sum(correlation_derived_weights, axis = 0)

            relative_weights = np.divide(correlation_derived_weights, np.tile(total_weights, (correlation_derived_weights.shape[0],1)))

            influence = np.multiply(raw_predicted, relative_weights)

            dropout_adjusted = np.tile(correlation_adjusted,(influence.shape[0],1)) - influence

        else:
             dropout_adjusted = None

        if str(truth) != "None" and verbose:
            print "Truth To Mean"
            print pearsonr(truth,self.means)
            print "Guess"
            print "======="

            print "Unadjusted"
            print pearsonr(unadjusted,truth)
            print "Correlation Adjusted:"
            print pearsonr(correlation_adjusted,truth)
            print "Slope/Correlation Adjusted:"
            print pearsonr(slope_corr_adjusted, truth)
            print "======"
            print "Centered Data"
            print pearsonr(unadjusted-self.means,truth-self.means)
            print pearsonr(correlation_adjusted-self.means,truth-self.means)
            print pearsonr(slope_corr_adjusted-self.means,truth-self.means)
            print "======"
            print "True sum"
            print np.sum(np.abs(truth-self.means))
            print "Prediction sum"
            print np.sum(np.abs(unadjusted-self.means))
            print np.sum(np.abs(correlation_adjusted-self.means))
            print np.sum(np.abs(slope_corr_adjusted-self.means))

            print "\n\n"

        return correlation_adjusted, raw_predicted, correlation_derived_weights, dropout_adjusted



    def multi_prediction(self, counts, masked=True, override = False, process_limit = False, filename = ""):

        print "Starting Parallel Prediction"

        if len(filename) > 0:
            if chk.check_hash(counts, filename, prefix=self.prefix):
                if override:
                    return np.load(self.prefix+filename+".npy")
                else:
                    if chk.check_hash(counts, filename +"_combined", prefix=self.prefix):

                        combined = np.copy(counts)
                        combined[combined == 0] = np.load(self.prefix + filename + "_combined.npy")[combined == 0]
                        print "Returning combined array, dimensions:"
                        print combined.shape
                        return combined
                    else:
                        pass

        if process_limit:
            pool = mlt.Pool(processes=min(mlt.cpu_count()-2,int(process_limit)))
        else:
            pool = mlt.Pool(processes=mlt.cpu_count()-2)

        print "Processors detected:"
        print pool._processes

        result = pool.map(compact_prediction, zip([self]*counts.shape[0],counts,range(counts.shape[0])),chunksize=100)

        predicted_array = np.asarray(result)

        print "Computed predicted array, shape:"
        print predicted_array.shape

        if len(filename) > 0:
            np.save(self.prefix+filename,predicted_array)
            chk.write_hash(counts, filename, prefix = self.prefix)

        if override:
            return predicted_array
        else:
            combined = np.copy(counts)
            combined[combined == 0] = predicted_array[combined == 0]

            if len(filename) > 0:
                    np.save(self.prefix+filename+"_combined",combined)
                    chk.write_hash(counts, filename+"_combined", prefix = self.prefix)

            print "Returning combined array of shape:"
            print combined.shape

            print "Check combined array, match counts?"
            print np.sum(combined[counts > 0] == counts[counts > 0])
            print np.sum(counts > 0)

            return combined


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
