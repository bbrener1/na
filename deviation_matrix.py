#!/usr/bin/env python

import sys
import os
import numpy as np
from sklearn.decomposition import PCA
import scipy.spatial.distance as spt
import check_hash as chk

from matrix_assurance import *


def compute_deviation_matrix(counts, neighbor_setting = 50, pretag="", presolve=None, output = None, filename = "deviation_matrix"):

    if output == None:
        output = sys.stdout
    elif isinstance(output,str):
        output = open(output,mode='w')



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

    int_sum = np.squeeze(np.asarray(np.sum(counts.T,axis=0)))

    output.write("Dropout mask shape, content, sum, and converted sum\n")
    output.write(str(int_sum.shape) + "\n")
    output.write(str(int_sum[:10]) + "\n")
    output.write(str(type(int_sum)) + "\n")

    dropout_mask = (int_sum > 0)

    output.write(str(dropout_mask.shape) + "\n")

    counts = counts[dropout_mask]

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

    neighbor_mean_matrix = np.zeros(counts.shape)
    std_dev_matrix = np.zeros(counts.shape)

    output.write("Progress of deviation matrix calculation:\n")
    for i in range(neighbor_mean_matrix.shape[0]):
        # output.write(str(distance_matrix[i]) + "\n")
        # print np.argsort(distance_matrix[i])
        # print np.argsort(distance_matrix[i]) < neighbor_setting
        # print counts[np.argsort(distance_matrix[i]) < neighbor_setting]
        # print np.mean(counts[np.argsort(distance_matrix[i]) < neighbor_setting],axis=0)
        neighbor_mean_matrix[i] = np.mean(counts[np.argsort(distance_matrix[i]) < neighbor_setting],axis=0)
        std_dev_matrix[i] = np.std(counts[np.argsort(distance_matrix[i]) < neighbor_setting],axis=0)

        if i%100 == 0:
            output.write(str(i) + "\n")


    output.write("Neighbor mean matrix content:\n")
    output.write(str(neighbor_mean_matrix[:10,:10]) + "\n")

    deviation_matrix = np.divide((counts-neighbor_mean_matrix),std_dev_matrix)

    output.write("Minimum in deviation matrix:\n")
    output.write(str( np.amin(deviation_matrix.flatten())) + "\n")

    deviation_matrix = np.nan_to_num(deviation_matrix)

    deviation_matrix[deviation_matrix < -10] = 0

    deviation_matrix[deviation_matrix > 10] = 0

    output.write("Was the deviation matrix flattened correctly? Second minimum\n")
    output.write(str( np.amin(deviation_matrix.flatten())) + "\n")

    np.save( pretag + filename,deviation_matrix)
    np.save( pretag + "dropout_" + filename,dropout_mask)
    np.save( pretag + "std_dev_" + filename, std_dev_matrix)
    chk.write_hash(counts,filename +".npy", pretag)

    output.write(str( deviation_matrix[:10,:10]) + "\n")

    return deviation_matrix, dropout_mask


def simple(folder):
    counts = np.load(folder+"counts.npy")
    compute_deviation_matrix(counts, neighbor_setting = 50, pretag=folder, filename= "deviation_matrix")

def main():

    print "Deviation matrix"

    prefix = sys.argv[1]

    if len(sys.argv)>2:
        counts = np.loadtxt(sys.argv[2]).T
    else:
        counts = prefix + "counts.npy"

    if len(sys.argv)>3:
        output = sys.argv[3]
    else:
        output = prefix + "dev_matrix_log"

    if len(sys.argv)>4:
        filename = sys.argv[4]
    else:
        filename = "deviation_matrix"


    counts = matrix_assurance(counts)

    print "Deviation matrix defaults initialized"

    if os.path.isfile(prefix + "deviation_matrix.npy"):
        if chk.check_hash(counts,"deviation_matrix.npy",prefix):
            return compute_deviation_matrix(counts, pretag = prefix, presolve="deviation_matrix.npy")
        else:
            return compute_deviation_matrix(counts, pretag = prefix, output = output, filename = filename)

    else:
        return compute_deviation_matrix(counts,pretag=prefix,output = output, filename = filename)


        # "Invalid pretag probably? Problem in main of deviation_matrix.py"
        # raise ValueError


if __name__ == "__main__":
    main()

























pass
