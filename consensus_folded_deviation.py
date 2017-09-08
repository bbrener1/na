#!/usr/bin/env python

import sys
import os
import numpy as np
from sklearn.decomposition import PCA
import scipy.spatial.distance as spt
import check_hash as chk
import itertools
import scipy.stats as st
from sklearn import preprocessing as pre


import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt

from matrix_assurance import *


def folded_deviation_matrix(counts, gold_network = None, neighbor_setting = 50, pretag="", presolve=None, output = None, filename = "deviation_matrix", fold = 5, verbose = True):

    if output == None:
        output = sys.stdout
    elif isinstance(output,str):
        output = open(output,mode='w')

    if not verbose:
        output = open(os.devnull, mode='w')


    if not isinstance(output,file):
        print "OUTPUT OF DEVIATION_MATRIX.PY, FUNCTION COMPUTE_DEVIATION_MATRIX IS NOT A FILE"
        raise SyntaxError

    if presolve != None:

        if type(presolve) is str:
            deviation_matrix = np.load(pretag + presolve)
            dropout_mask = np.load(pretag+ "dropout_"+presolve)
            return deviation_matrix,dropout_mask.astype(dtype=bool)
        if isinstance(presolve,tuple):
            if isinstance(presolve[0],np.ndarray):
                return presolve[0],presolve[1]
        else:
            print "Presolve not a file address or numpy ndarray"
            raise ValueError




    #
    # print counts[:5,:5]

    # counts = np.asarray([[1,2,3],[1,2,4],[2,2,4],[3,2,1],[4,2,1],[4,2,2]])

    output.write("Deviation Matrix counts shape:\n")
    output.write(str(counts.shape) + "\n")

    int_sum = np.squeeze(np.asarray(np.sum(counts,axis=0)))

    output.write("Dropout mask shape, content, sum, and converted sum\n")
    output.write(str(int_sum.shape) + "\n")
    output.write(str(int_sum[:10]) + "\n")
    output.write(str(type(int_sum)) + "\n")

    dropout_mask = (int_sum > 0)

    output.write(str(dropout_mask.shape) + "\n")

    counts = counts[:,dropout_mask]

    #
    # distance_matrix = np.zeros((counts.shape[0],counts.shape[0]))
    #
    # for i in range(counts.shape[0]):
    #     for j in range(i,counts.shape[0]):
    #         distance_matrix[i,j] = np.linalg.norm(counts[i]-counts[j])
    #         distance_matrix[j,i] = distance_matrix[i,j]

    intermediate_model = PCA(n_components=min(50,counts.shape[1]))

    reduced_counts = intermediate_model.fit_transform(counts)

    distance_matrix = spt.squareform(spt.pdist(reduced_counts))

    output.write("Distance matrix content, shape\n")
    output.write(str(distance_matrix[:10,:10]) + "\n")
    output.write(str(distance_matrix.shape) + "\n")

    folds = list(itertools.combinations(range(fold),int(fold/2)+1))

    neighbor_mean_matrix = np.zeros((counts.shape[0],counts.shape[1],len(folds)))
    std_dev_matrix = np.zeros((counts.shape[0],counts.shape[1],len(folds)))

    output.write("Progress of deviation matrix calculation:\n")

    for i, row in enumerate(neighbor_mean_matrix):

        integer_mask = np.random.randint(fold, size=neighbor_setting)

        neighborhood = counts[np.argsort(distance_matrix[i])][:neighbor_setting]



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

            # std_dev_matrix[i,:,j] = np.std(neighborhood[fold_mask],axis=0)




        if i%100 == 0:
            output.write(str(i) + "\n")
            print i

    deviation_matrix = np.zeros(neighbor_mean_matrix.shape)

    for i in range(len(folds)):
        deviation_matrix[:,:,i] = np.divide((counts-neighbor_mean_matrix[:,:,i]),std_dev_matrix[:,:,i])
        deviation_matrix[:,:,i] = np.nan_to_num(deviation_matrix[:,:,i])
        deviation_matrix[:,:,i][deviation_matrix[:,:,i] > 10] = 10
        deviation_matrix[:,:,i][deviation_matrix[:,:,i] < -10] = -10


    numeric_deviation_matrix = np.median(deviation_matrix, axis=2)

    deviation_matrix = np.sum(np.abs(deviation_matrix) > 1.5 , axis = 2) > 5

    qc_array = np.nan_to_num(st.variation(neighbor_mean_matrix, axis = 2).flatten())

    # print "QC array shape"
    # print qc_array.shape
    # print np.max(qc_array)
    # print np.max(median_mean_matrix)
    # print np.min(median_mean_matrix)

    plt.figure()
    plt.hist(qc_array)
    plt.title("Coefficients of Inter-Fold Variation, Deviation Matrix")
    plt.savefig(pretag + "coef_var_fold.png")

    output.write("Neighbor mean matrix content:\n")
    output.write(str(neighbor_mean_matrix[:10,:10]) + "\n")


    output.write("Minimum in deviation matrix:\n")
    output.write(str( np.amin(deviation_matrix.flatten())) + "\n")

    deviation_matrix = np.nan_to_num(deviation_matrix)

    deviation_matrix[deviation_matrix < -10] = 0

    deviation_matrix[deviation_matrix > 10] = 0

    output.write("Was the deviation matrix flattened correctly? Second minimum\n")
    output.write(str( np.amin(deviation_matrix.flatten())) + "\n")

    # np.save( pretag + "reduced_" + filename, deviation_matrix[:,dropout_mask] )

    if str(gold_network) != "None":
        np.save( pretag + "reduced_gold_network", gold_network[dropout_mask].T[dropout_mask].T )

    np.save( pretag + filename,deviation_matrix)
    np.save( pretag + "numeric_" + filename, numeric_deviation_matrix)
    np.save( pretag + "dropout_" + filename, dropout_mask)
    np.save( pretag + "std_dev_" + filename, std_dev_matrix)
    chk.write_hash(counts,filename +".npy", pretag)

    output.write(str( deviation_matrix[:10,:10]) + "\n")

    return deviation_matrix, dropout_mask, numeric_deviation_matrix, neighbor_mean_matrix


def simple(folder):
    counts = np.load(folder+"counts.npy")
    folded_deviation_matrix(counts, neighbor_setting = 50, pretag=folder, filename= "deviation_matrix")

def main():

    print "Deviation matrix"

    prefix = sys.argv[1]

    if len(sys.argv)>2:
        neighbors = int(sys.argv[2])
    else:
        neighbors = 50


    if len(sys.argv)>3:
        counts = sys.argv[3]
    else:
        counts = prefix + "counts.npy"

    if len(sys.argv)>4:
        gold_network = sys.argv[4]
    else:
        gold_network = prefix + "gold_network.npy"

    if len(sys.argv)>5:
        output = sys.argv[5]
    else:
        output = prefix + "dev_matrix_log"

    if len(sys.argv)>6:
        filename = sys.argv[6]
    else:
        filename = "cons_dev_matrix"


    output = None
    counts = matrix_assurance(counts)
    gold_network = matrix_assurance(gold_network)

    print "Deviation matrix defaults initialized"

    if os.path.isfile(prefix + "cons_dev_matrix.npy"):
        if chk.check_hash(counts,"cons_dev_matrix.npy",prefix):
            return folded_deviation_matrix(counts, gold_network = gold_network, neighbor_setting = neighbors, pretag = prefix, presolve="cons_dev_matrix.npy")
        else:
            return folded_deviation_matrix(counts, gold_network = gold_network, neighbor_setting = neighbors, pretag = prefix, output = output, filename = filename)

    else:
        return folded_deviation_matrix(counts, gold_network = gold_network, neighbor_setting = neighbors, pretag=prefix,output = output, filename = filename)


        # "Invalid pretag probably? Problem in main of cons_dev_matrix.py"
        # raise ValueError


if __name__ == "__main__":
    main()

























pass
