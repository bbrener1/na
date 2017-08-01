#!/usr/bin/env python

import sys
import os
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import check_hash as chk
from compare import *
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import scipy.stats as st


import itertools

from matrix_assurance import *

def plot_edge_certainty(correlation,connectivity,filename):
    fig = plt.figure()
    cm = plt.cm.get_cmap('magma')
    degrees = np.sum(connectivity, axis=1)
    x = np.repeat(degrees,degrees.shape[0])
    y = np.tile(degrees,degrees.shape[0])
    print "EDGE CERTAINTY DEBUG"
    print connectivity.shape
    print correlation.shape
    print filename
    ax = fig.add_subplot(111)
    # ax.set_yscale('log')
    # ax.set_xscale('log')
    ax.scatter(x.flatten()[connectivity.flatten()],y.flatten()[connectivity.flatten()], marker='x', alpha = .3, c=correlation.flatten()[connectivity.flatten()],s=1, cmap=cm)
    plt.savefig(filename,dpi=300)

def quick_correlation(observation_matrix, name = None, prefix = ""):

    name = prefix + name

    correlation_matrix = np.abs(np.nan_to_num(np.corrcoef(observation_matrix.T)))
    correlation_matrix = correlation_matrix - np.diag(np.diag(correlation_matrix))

    if name != None:
        np.save(prefix + "quick_correlation_backup", np.nan_to_num(correlation_matrix))
        chk.write_hash(observation_matrix, "quick_correlation_backup.npy", prefix)



    print "Quick Correlation Debug"
    print np.sum(correlation_matrix)
    print correlation_matrix.shape
    print "Non-zero edges: "
    print np.sum(correlation_matrix > 0)


    plt.figure()
    plt.hist(correlation_matrix.flatten(), bins=20, log=True)
    plt.savefig(prefix + "quick_correlation_debug.png")

    return correlation_matrix

def folded_correlation(observation_matrix, name = None, prefix = "", fold_number = 5):

    name = prefix + name
    folds = list(itertools.combinations(range(fold_number),int(fold_number/2)+1))

    integer_mask = np.random.randint(fold_number, size = observation_matrix.shape[0])

    correlation_matrix = np.zeros((len(folds),observation_matrix.shape[1],observation_matrix.shape[1]))

    for i, fold in enumerate(folds):
        fold_mask = np.zeros(observation_matrix.shape[0],dtype=bool)
        for element in fold:
            fold_mask = np.logical_or(fold_mask,integer_mask == element)

        print i

        correlation_matrix[i] = np.abs(np.nan_to_num(np.corrcoef(observation_matrix[fold_mask].T)))
        correlation_matrix[i] = correlation_matrix[i] - np.diag(np.diag(correlation_matrix[i]))

        print i

    print "Computing consensus"

    consensus_correlation = np.median(np.sort(correlation_matrix, axis=0)[:5,:,:],axis=0)

    print "Computing the QC array"

    # consensus_correlation = np.median(correlation_matrix, axis = 0)

    qc_array = np.nan_to_num(st.variation(correlation_matrix, axis = 0).flatten())

    print "QC array shape"
    print qc_array.shape
    print np.max(qc_array)
    print np.max(correlation_matrix)
    print np.max(consensus_correlation)
    print np.min(consensus_correlation)

    plt.figure()
    plt.hist(qc_array)
    plt.savefig(prefix + "corr_var_fold.png")

    if name != None:
        np.save(prefix + "folded_correlation_backup", np.nan_to_num(consensus_correlation))
        chk.write_hash(observation_matrix, "folded_correlation_backup.npy", prefix)


    return consensus_correlation

def quick_threshold_analysis(observations, gold, scroll = None, presolve= None, name = "", prefix = ""):

    name = prefix + name

    if scroll == None:
        scroll = map(lambda x: float(x)*.01,range(1,96,5))

    if presolve == None:
        correlation = folded_correlation(observations, name = "folded_correlation_backup", prefix= prefix)

    else:
        correlation = np.load(prefix + presolve)

    print correlation.shape


    r2_stat = np.zeros(len(scroll))
    unconnected_index = np.zeros(len(scroll))

    true_pos = np.zeros(len(scroll))
    false_pos = np.zeros(len(scroll))

    precision = np.zeros(len(scroll))
    recall = np.zeros(len(scroll))

    enrichment = np.zeros(len(scroll))

    for j, tau in enumerate(scroll):
        print tau
        # degree_ratings = np.zeros(observation_matrix.shape[1])
        # for i, gene in enumerate(correlation_strength_matrix):
        #     degree_ratings[i] = np.sum(gene > tau)
        degree_ratings = np.sum(correlation > tau, axis=1)
        gold_degree_ratings = np.sum(gold, axis=1)

        print "Unconnected nodes: " + str(np.sum(degree_ratings.flatten() < tau))
        print "Unconnected gold standard nodes:" + str(np.sum(gold_degree_ratings.flatten() < 2))


        connectivity = correlation > tau

        compare(connectivity,gold)

        true_pos[j],false_pos[j],precision[j],recall[j], enrichment[j] = restrictive_compare(connectivity,gold)



        model = LinearRegression()
        plt.figure()
        hist = plt.hist(degree_ratings, log=True, bins=50,alpha=.3,label= name)
        # print hist
        degree_hist = hist[0]
        hist_coord = hist[1][:-1]
        # print degree_hist
        # print np.log(degree_hist)
        # print hist_coord
        degree_hist = degree_hist.reshape(-1,1)
        hist_coord = hist_coord.reshape(-1,1)
        # print degree_hist
        # print np.log(degree_hist)
        # print hist_coord
        model.fit(hist_coord,np.log(degree_hist+1))
        # print model.coef_
        # print model.intercept_
        m = model.coef_[0][0]
        b = model.intercept_
        # print "Where are the residuals?"
        # print np.polyfit(degree_hist,range(degree_hist.shape[0]),1,full=True)
        # m, b = np.polyfit(hist_coord,np.log(degree_hist),1, full=True)
        # print m
        # print b
        # print m*hist_coord+b
        # print np.exp(m*hist_coord+b)
        plt.plot(hist_coord, np.exp(m*hist_coord+b),label= name)
        plt.title("Node Degree Histogram (Small World Assumption)")
        plt.xlabel("Node Degree")
        plt.legend()
        plt.savefig(name + str(tau) +".png")
        # print degree_ratings.reshape(-1,1)
        # print degree_ratings > 1
        # print np.sum((degree_ratings) > 1)
        # print tau
        r2_stat[j] = r2_score(degree_hist,np.exp(m*hist_coord+b))

        if tau == .13:
            plot_edge_certainty(correlation,connectivity, prefix+"edge_certainty.png")

    print r2_stat

    enrichment_optimum = enrichment[:22].max()
    degree_ratings = np.sum(correlation > tau, axis=1)
    degree_mask = np.argsort(degree_ratings) > (degree_ratings.shape[0]-200)
    np.save("TF_degree_mask",degree_mask)
    plt.figure()
    plt.plot(scroll,r2_stat, label = name)
    plt.ylabel("Concordance To The Small World Assumption (R^2)")
    plt.xlabel("Cutoff Value For Edge Construction")
    plt.title("Picking a cutoff value")
    plt.legend()
    plt.savefig(name + "_tau_plot.png")

    # tp_dev, fp_dev = quick_threshold_analysis(deviation_matrix, gold_std, name = "quicktest_dev")
    # tp_base, fp_base = quick_threshold_analysis(counts[dropout_mask].T, gold_std, name = "quicktest_base")
    # correlation_network_analysis(deviation_matrix,gold_std,matrix_label="outlier")
    # correlation_network_analysis(counts,gold_std,matrix_label="plain_counts")

    plt.figure()
    plt.plot(false_pos, true_pos)
    plt.savefig(prefix + "roc_compare.png")

    plt.figure()
    plt.plot(recall, precision)
    plt.xlabel("Recall (Capture rate)")
    plt.ylabel("Precision (Rate of true positives vs all positives)")
    plt.savefig(prefix + "pr_curve.png")

    plt.figure()
    plt.plot(scroll,enrichment)
    plt.savefig(prefix + "enrichment.png")

    return true_pos, false_pos

def main():

    prefix = sys.argv[1]

    if len(sys.argv)> 2:
        observations = sys.argv[2]
    else:
        observations = prefix + "folded_dev_matrix.npy"

    if len(sys.argv)> 3:
        gold = sys.argv[3]
    else:
        gold = prefix+"gold_network.npy"

    if len(sys.argv)>4:
        name = sys.argv[4]
    else:
        name = "correlation"

    observations = matrix_assurance(observations)
    gold = matrix_assurance(gold)
    custom_scroll = map(lambda x: float(x)*.01,range(1,50,1))



    if os.path.isfile(prefix + "folded_correlation_backup.npy"):
        if chk.check_hash(observations,"folded_correlation_backup.npy", prefix = prefix):
            return quick_threshold_analysis(observations, gold, scroll = custom_scroll, presolve="folded_correlation_backup.npy", name = name, prefix = prefix)
        else:
            return quick_threshold_analysis(observations, gold, scroll = custom_scroll, name = name, prefix = prefix)
    else:
        return quick_threshold_analysis(observations, gold, scroll = custom_scroll, name = name, prefix = prefix)


if __name__ == "__main__":
    main()
