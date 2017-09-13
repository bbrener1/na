#!/usr/bin/env python

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import scipy.special

import numpy as np

counts = np.load('counts.npy')
header = np.load('header_backup.npy')

# plt.figure("counts_frequency")
# plt.hist(counts.ravel(),bins=101,log=True)
# plt.title("Distribution of values in normalized expression matrix")
# plt.xlabel("Log2 expression value relative to normalized library size (per Lun)")
# plt.ylabel("Frequency of value (log10 scale)")
# plt.savefig("counts_frequency.png")
#
# plt.figure("transcript_totals_hist")
# plt.hist(np.sum(np.exp2(counts),axis=1),bins=50,log=True)
# plt.title("Histogram of total per-cell transcript counts")
# plt.xlabel("Number of total transcripts in a cell")
# plt.ylabel("Frequency (log10)")
#
# plt.figure("trans_size_vs_exp")
# plt.scatter(counts.ravel(),np.repeat(np.sum(np.exp2(counts),axis=1),counts.shape[1]),s=.01,alpha=.3)
# plt.title("Expression of each gene vs size of cell transcriptome")
# plt.xlabel("Log2 gene expression value")
# plt.ylabel("Total per-cell transcriptome counts (log10 scale)")
# plt.yscale('log')
# plt.savefig('figures/trans_size_vs_exp.png')
#
# plt.figure("trans_size_vs_exp_alpha")
# plt.scatter(counts.ravel(),np.repeat(np.sum(np.exp2(counts),axis=1),counts.shape[1]),s=2.0*(np.random.logistic(counts.ravel())-.5),alpha=.3)
# plt.title("Expression of each gene vs size of cell transcriptome")
# plt.xlabel("Log2 gene expression value")
# plt.ylabel("Total per-cell transcriptome counts (log10 scale)")
# plt.yscale('log')
# plt.savefig('figures/trans_size_vs_exp_alpha.png')

plt.figure("gene_histogram_gigaplex")
plt.suptitle("Set of histograms of expression values for randomly chosen genes")
plt.xlabel("Log2 gene expression values")
plt.ylabel("Frequency")
for i, pick in enumerate(np.random.randint(counts.shape[1], size=20)):
    print pick
    plt.subplot(4,5,i+1)
    plt.title(header[pick],size=6)
    plt.hist(counts[:,pick], bins=21,log=True)
plt.savefig("figures/gene_histogram_gigaplex.png",dpi=300)
