#!/usr/bin/env python
#
# covariance
#
# imp_obs = predict_drops(obs, covariance)
#
# covariance = deviation_matrix(imp_obs)

import numpy as np

from linear_regression import *
from sklearn.decomposition import PCA
import scipy.spatial.distance as spt

def main():

    counts = np.load(sys.argv[1])

    # counts = pre.scale(counts)

    slopes, intercepts, means, correlations = linear_regression(counts)

    # partial = partial_correlation(counts)



    for pick in np.random.randint(counts.shape[0], size=5):

        print pick

        predict(counts[pick,:],counts[pick,:], slopes, intercepts, means, correlations)
