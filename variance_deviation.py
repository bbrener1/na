#!/usr/bin/env python

import sys
import numpy as np
from matrix_assurance import *

def variance_discrepancy(deviation_matrix, counts, labels):
    count_variance = np.var(counts, axis=1)
    deviation_variance = np.var(deviation_matrix, axis=1)

    mean_deviation = np.mean(deviation_matrix, axis = 1)

    variance_index = np.divide(count_variance,deviation_variance)

    mask = np.argsort(deviation_variance) < 15

    plt.figure()
    plt.bar(np.arange(len(deviation_variance)),np.sort(deviation_variance))
    plt.xticks(np.arange(len(deviation_variance)), np.arange(len(deviation_variance))[np.argsort(deviation_variance)].astype(dtype=str))

    plt.savefig("dev_variance.png")

    plt.figure()
    plt.bar(np.arange(len(variance_index)),np.sort(variance_index))
    plt.xticks(np.arange(len(variance_index)), np.arange(len(variance_index))[np.argsort(variance_index)].astype(dtype=str))
    plt.savefig("variance_index.png")

    plt.figure()
    plt.bar(np.arange(len(mean_deviation)), np.sort(mean_deviation),tick_label=np.arange(len(mean_deviation))[np.argsort(mean_deviation)].astype(dtype=str))
    # plt.xticks(np.arange(len(mean_deviation)), np.arange(len(mean_deviation))[np.argsort(mean_deviation)].astype(dtype=str))
    plt.savefig("mean_deviation.png")

    # variance_model = scipy.stats.f()
    #
    # variance_model.pdf()

    return mask

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

    deviation_matrix = matrix_assurance(deviation_matrix)
    counts = matrix_assurance(counts)
    labels = matrix_assurance(labels)


    variance_discrepancy(deviation_matrix, counts, labels)


    # if os.path.isfile(prefix + "quick_correlation_backup.npy"):
    #     if chk.check_hash(observations,"quick_correlation_backup.npy", prefix = prefix):
    #
    #     else:
    # else:




if __name__ == "__main__":
    main()
