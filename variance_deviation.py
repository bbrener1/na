#!/usr/bin/env python

import sys
import numpy as np
from matrix_assurance import *

import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt


def variance_discrepancy(deviation_matrix, counts, std_devs, labels, prefix = ""):
    count_variance = np.var(counts, axis=0)
    deviation_variance = np.var(deviation_matrix, axis=0)

    mean_abs_deviation = np.mean(np.abs(deviation_matrix), axis = 0)

    variance_index = np.divide(count_variance,deviation_variance)

    dev_var_mask = np.argsort(deviation_variance) < 200
    index_mask = np.argsort(variance_index) > len(variance_index) - 200
    abs_mask = np.argsort(mean_abs_deviation) > len(mean_abs_deviation) - 200
    real_mask = np.argsort(np.divide(np.std(counts,axis=0),np.mean(std_devs,axis=0))) > counts.shape[0] - 200

    plt.figure()
    plt.bar(np.arange(len(deviation_variance)),np.sort(deviation_variance))
    plt.xticks(np.arange(len(deviation_variance)), np.arange(len(deviation_variance))[np.argsort(deviation_variance)].astype(dtype=str))

    plt.savefig(prefix + "dev_variance.png")

    plt.figure()
    plt.bar(np.arange(len(variance_index)),np.sort(variance_index))
    plt.xticks(np.arange(len(variance_index)), np.arange(len(variance_index))[np.argsort(variance_index)].astype(dtype=str))
    plt.savefig(prefix + "variance_index.png")

    plt.figure()
    plt.bar(np.arange(len(mean_abs_deviation)), np.sort(mean_abs_deviation),tick_label=np.arange(len(mean_abs_deviation))[np.argsort(mean_abs_deviation)].astype(dtype=str))
    # plt.xticks(np.arange(len(mean_abs_deviation)), np.arange(len(mean_abs_deviation))[np.argsort(mean_abs_deviation)].astype(dtype=str))
    plt.savefig(prefix + "mean_abs_deviation.png")

    plt.figure()
    plt.scatter(np.std(counts,axis=0),np.mean(std_devs,axis=0),s=.1,marker='x')
    plt.savefig(prefix + "real_deviation_discrepancy.png", dpi=300)

    gene_candidates = labels[dev_var_mask]
    np.savetxt(prefix + "potential_TFs_dev_var.txt", gene_candidates, fmt="%s")

    gene_candidates = labels[index_mask]
    np.savetxt(prefix + "potential_TFs_index.txt", gene_candidates, fmt="%s")

    gene_candidates = labels[abs_mask]
    np.savetxt(prefix + "potential_TFs_abs.txt", gene_candidates, fmt="%s")

    gene_candidates = labels[real_mask]
    np.savetxt(prefix + "potential_TFs_real.txt", gene_candidates, fmt="%s")

    # variance_model = scipy.stats.f()
    #
    # variance_model.pdf()



def main():
    prefix = sys.argv[1]

    if len(sys.argv)> 2:
        deviation_matrix = sys.argv[2]
    else:
        deviation_matrix = prefix + "deviation_matrix.npy"

    if len(sys.argv)> 3:
        counts = sys.argv[3]
    else:
        counts = prefix+"counts.npy"

    if len(sys.argv)>4:
        labels = sys.argv[4]
    else:
        labels = prefix+"header_backup.npy"

    if len(sys.argv)>5:
        std_devs = sys.argv[5]
    else:
        std_devs = prefix + "std_dev_deviation_matrix.npy"

    deviation_matrix = matrix_assurance(deviation_matrix)
    counts = matrix_assurance(counts)
    labels = matrix_assurance(labels)
    std_devs = matrix_assurance(std_devs)

    variance_discrepancy(deviation_matrix, counts, std_devs, labels, prefix = prefix)


    # if os.path.isfile(prefix + "quick_correlation_backup.npy"):
    #     if chk.check_hash(observations,"quick_correlation_backup.npy", prefix = prefix):
    #
    #     else:
    # else:




if __name__ == "__main__":
    main()
