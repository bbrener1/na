#!/usr/bin/env/ python

import numpy as np
import scipy.spatial.distance as spt
from scipy.stats import pearsonr
from sklearn.decomposition import PCA

counts = np.load("run1neigh50/counts.npy")
dist_model = PCA(n_components=50)
dist_interm = dist_model.fit_transform(counts)
dist = spt.squareform(spt.pdist(dist_interm))
import linear_regression
linear_regression
slopes, intercepts, means, correlations = linear_regression.linear_regression(counts)
partial = linear_regression.partial_correlation(counts)
deviations = np.load("run1neigh50/deviation_matrix.npy")
neighborhood_500 = np.argsort(dist[871]) 
good_predictions = np.zeros((10,4773))
bad_predictions = np.zeros((10,4773))
for i in range(10):
    print dist[871][neighborhood_500][:10][i]
    good_predictions[i] = linear_regression.predict(counts[neighborhood_500][:10][i],counts[871],slopes,intercepts,means,correlations,partial)

for i in range(10):
    print dist[871][neighborhood_500][-10:][i]
    bad_predictions[i] = linear_regression.predict(counts[neighborhood_500][-10:][i],counts[871],slopes,intercepts,means,correlations,partial)

stds = np.std(counts[neighborhood_500],axis=0)
means = np.mean(counts[neighborhood_500],axis = 0)
diff = counts[871] - np.mean(good_predictions,axis=0)
diff_std = np.divide(diff,stds)
diff_std[np.isinf(diff_std)] = 0
np.sum(np.abs(diff_std))

