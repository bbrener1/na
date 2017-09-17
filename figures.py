#!/usr/bin/env python

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import scipy.special

import numpy as np

from scipy.cluster import hierarchy as hrc
from sklearn.decomposition import PCA

from sklearn.cluster import AgglomerativeClustering


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
#
# plt.figure("gene_histogram_gigaplex")
# plt.suptitle("Set of histograms of expression values for randomly chosen genes")
# plt.xlabel("Log2 gene expression values")
# plt.ylabel("Frequency")
# for i, pick in enumerate(np.random.randint(counts.shape[1], size=20)):
#     print pick
#     plt.subplot(4,5,i+1)
#     plt.title(header[pick],size=6)
#     plt.hist(counts[:,pick], bins=21,log=True)
# plt.savefig("figures/gene_histogram_gigaplex.png",dpi=300)
#
# cfig = plt.figure("raw_expression_heatmap")
# ax1 = cfig.add_axes((.1,0,.7,1))
# ax2 = cfig.add_axes((.85,.3,.06,.4))
# im = ax1.imshow(counts, cmap='hot')
# ax1.set_xlabel("Genes")
# ax1.set_ylabel("Cells")
# ax1.set_title("Log2 expression values (Unordered)")
# cfig.colorbar(mappable=im,cax=ax2)
# plt.savefig('figures/raw_expression_heatmap', dpi=600)
#

#
#
# gene_linked = hrc.linkage(counts.T, method='average', metric='correlation')
# gene_dendrogram = hrc.dendrogram(gene_linked,no_plot=True)
#
# cell_linked = hrc.linkage(counts, method='average', metric='cosine')
# clusterization = hrc.fcluster(cell_linked, criterion='inconsistent',t=.5,)
# cell_dendrogram = hrc.dendrogram(cell_linked,no_plot=True)
#
# fig = plt.figure("doubly_clustered", figsize=(8,4))
# ax1 = fig.add_axes([.09,.1,.2,.6])
# # display_dendrogram = hrc.dendrogram(cell_linked, p=3, truncate_mode='level',orientation='left',show_contracted=True,ax=ax1)
# with plt.rc_context({'lines.linewidth':0.1}):
#     display_dendrogram = hrc.dendrogram(cell_linked,orientation='left',ax=ax1)
# # ax1.set_xlim(left=1.0,right=.75)
# # ax1.set_xscale('log')
#
# ax2 = fig.add_axes([.3,.71,.55,.2])
# # display_dendrogram = hrc.dendrogram(gene_linked, p=3, truncate_mode='level',ax=ax2)
# with plt.rc_context({'lines.linewidth':0.1}):
#     display_dendrogram = hrc.dendrogram(gene_linked, ax=ax2)
# # ax2.set_ylim(top=1,bottom=.75)
# # ax2.set_yscale('log')
#
# print counts.shape
# print len(cell_dendrogram['leaves'])
# print len(gene_dendrogram['leaves'])
#
#
# ax3 = fig.add_axes([.3,.1,.55,.6])
#
# sorted_singly = counts[np.flip(cell_dendrogram['leaves'],0)]
# sorted_doubly = sorted_singly.T[gene_dendrogram['leaves']].T
# # sorted_doubly = np.concatenate((sorted_doubly.T,np.ones((sorted_doubly.T.shape[0],1))*-10), axis=1)
# # sorted_doubly = np.concatenate((sorted_doubly,np.ones((sorted_doubly.shape[0],1))*10), axis=1)
# im = ax3.imshow(sorted_doubly, cmap='hot', aspect='auto')
# # plt.title("Residual Expression of Genes In Cells, Clustered Hierarchically")
# # plt.xlabel("Genes")
# # plt.ylabel("Cells")
# ax4 = fig.add_axes([.85,.1,.05,.6])
# # ax4.set_ylim(bottom=-10,top=10)
# fig.colorbar(mappable=im, fraction=.99, ax=ax4)
# np.save("clustered_counts",sorted_doubly)
# np.save("clustered_header",header[gene_dendrogram['leaves']])
# np.save("gene_clustering_indecies", gene_dendrogram['leaves'])
# np.save("cell_clustering_indecies", cell_dendrogram['leaves'])
# plt.savefig("figures/doubly_clustered_raw_genes.png", dpi=800)

fig = plt.figure("gene_scatter_gigaplex")
plt.suptitle("Set of scatter plots of expression values for randomly chosen gene pairs")
plt.xlabel("Log2 gene expression values")
plt.ylabel("Frequency")
for i, pick1 in enumerate(np.random.randint(counts.shape[1], size=20)):
    for j, pick2 in enumerate(np.random.randint(counts.shape[1], size=20)):
        plt.subplot(4,5,i+1)
        plt.xlabel(header[pick1],size=6)
        plt.ylabel(header[pick2],size=6)
        plt.xticks(np.arange(0,14,2),size=3)
        plt.yticks(np.arange(2,10,2),size=3)
        plt.scatter(counts[:,pick1],counts[:,pick2],marker='x',s=.1,alpha=.1,c='b')
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig("figures/gene_scatter_gigaplex.png",dpi=500)
