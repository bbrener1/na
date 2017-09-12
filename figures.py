#!/usr/bin/env python

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

import numpy as np

counts = np.load('counts.npy')

plt.figure("counts_frequency")
plt.hist(counts.ravel(),bins=101,log=True)
plt.title("Distribution of values in normalized expression matrix")
plt.xlabel("Log2 expression value relative to normalized library size (per Lun)")
plt.ylabel("Frequency of value (log10 scale)")
plt.savefig("counts_frequency.png")

plt.figure("transcript_totals_hist")
plt.hist(np.sum(np.exp2(counts),axis=1),bins=50,log=True)
plt.title("Histogram of total per-cell transcript counts")
plt.xlabel("Number of total transcripts in a cell")
plt.ylabel("Frequency (log10)")

plt.figure("trans_size_vs_exp")
plt.scatter(counts.ravel(),np.repeat(np.sum(np.exp2(counts),axis=1),counts.shape[1]),s=.1,alpha=.01)
plt.title("Expression of each gene vs size of cell transcriptome")
plt.xlabel("Log2 gene expression value")
plt.ylabel("Total per-cell transcriptome counts (log10 scale)")
plt.yscale('log')
plt.savefig('figures/trans_size_vs_exp.png', dpi=300)
