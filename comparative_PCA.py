#!/usr/bin/env python

from sklearn.decomposition import PCA

import numpy as np
import sys

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

model = PCA(n_components=3)

raw_counts = np.load(sys.argv[1])

deviation_matrix = np.load(sys.argv[2])

raw_PCA = model.fit_transform(raw_counts)

deviation_PCA = model.fit_transform(deviation_matrix)

print raw_counts.shape
print deviation_matrix.shape

print raw_PCA.shape
print deviation_PCA.shape

fig1 = plt.figure()
ax = fig1.add_subplot(111,projection='3d')
ax.scatter(raw_PCA[:,0],raw_PCA[:,1],raw_PCA[:,2], alpha= 1, marker='.',s=.5)
fig1.savefig("raw_PCA.png")

fig2 = plt.figure()
ax = fig2.add_subplot(111,projection='3d')
ax.scatter(deviation_PCA[:,0],deviation_PCA[:,1],deviation_PCA[:,2],alpha = 1,marker='.', s=.5)
fig2.savefig("deviation_PCA.png")
