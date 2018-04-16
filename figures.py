#!/usr/bin/env python

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import scipy.special
from scipy.stats import linregress

import numpy as np

from scipy.cluster import hierarchy as hrc
from sklearn.decomposition import PCA
from sklearn.cluster.bicluster import SpectralBiclustering

from sklearn.cluster import AgglomerativeClustering


counts = np.load('counts.npy')
deviation_matrix = np.load('numeric_deviation_matrix.npy')
dropout_deviation_matrix = np.load('dropout_deviation_matrix.npy')
header = np.load('header_backup.npy')
cell_identities = np.load('cell_identity.npy')
correlations = np.corrcoef(counts.T)
correlations[np.identity(correlations.shape[0],dtype=bool)]=0

clustered_counts = np.load('clustered_counts.npy')
clustered_header = np.load('clustered_header.npy')
# clustered_cell_labels = np.load('clustered_cells.npy')

cell_clustering_indecies = np.load('cell_clustering_indecies.npy')
gene_clustering_indecies = np.load('gene_clustering_indecies.npy')

residual_gene_clusters = np.load('residual_gene_clusters.npy')
# np.load('residual_cell_clusters.npy')
clustered_residuals = np.load('clustered_residuals.npy')


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
# ax3.set_xticks(np.arange(0,4773,100))
# ax3.set_yticks(np.arange(0,1656,100))
# ax.tick_params(labelsize=6)
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
#
# fig = plt.figure("gene_scatter_gigaplex")
# plt.suptitle("Set of scatter plots of expression values for randomly chosen gene pairs")
# plt.xlabel("Log2 gene expression values")
# plt.ylabel("Frequency")
# for i, pick1 in enumerate(np.random.randint(counts.shape[1], size=20)):
#     for j, pick2 in enumerate(np.random.randint(counts.shape[1], size=20)):
#         plt.subplot(4,5,i+1)
#         plt.xlabel(header[pick1],size=6)
#         plt.ylabel(header[pick2],size=6)
#         plt.xticks(np.arange(0,14,2),size=3)
#         plt.yticks(np.arange(2,10,2),size=3)
#         plt.scatter(counts[:,pick1],counts[:,pick2],marker='x',s=.1,alpha=.1,c='b')
# plt.tight_layout(rect=[0, 0.03, 1, 0.95])
# plt.savefig("figures/gene_scatter_gigaplex.png",dpi=500)
#
# plt.figure("correlation_histogram",figsize=(5,4))
# plt.suptitle("Frequencies of correlations between genes")
# plt.subplot(111)
# plt.hist(correlations.ravel(),bins=21,log=True)
# plt.xlabel("Pearson correlation value, gene-gene comparison")
# plt.ylabel("Frequency (log10 scale)")
# plt.savefig("figures/correlation_histogram.png")
#
# plt.figure("random_correlates")
# plt.suptitle("Set of scatter plots for correlated gene pairs, r>.5")
# plt.xlabel("Log2 gene expression values")
# plt.ylabel("Frequency")
# x, y = np.where(np.abs(correlations) > .5)
# correlate_picks = []
# for i, pick in enumerate(np.random.randint(x.shape[0], size=20)):
#     correlate_picks.append(pick)
#     plt.subplot(4,5,i+1)
#     plt.xlabel(header[x[pick]],size=6)
#     plt.ylabel(header[y[pick]],size=6)
#     plt.xticks(np.arange(0,14,2),size=3)
#     plt.yticks(np.arange(2,10,2),size=3)
#     plt.scatter(counts[:,x[pick]],counts[:,y[pick]],marker='x',s=.1,alpha=.1,c='b')
# plt.tight_layout(rect=[0, 0.03, 1, 0.95])
# # print x.shape
# plt.savefig("figures/correlated_scatter_gigaplex.png",dpi=500)
#
# plt.suptitle("Fitted linear regressions of correlated genes")
# # plt.xlabel("Log2 gene expression values")
# # plt.ylabel("Frequency")
# # x, y = np.where(np.abs(correlations) > .5)
# for i, pick in enumerate(correlate_picks):
#     slope, intercept, rvalue, pvalue, std_err = linregress(counts[:,x[pick]], y=counts[:,y[pick]])
#     plt.subplot(4,5,i+1)
#     plt.plot(counts[:,x[pick]], intercept + counts[:,x[pick]]*slope, 'r', label = str(np.around(rvalue, decimals=3)))
#     plt.legend(fontsize=5)
# # plt.tight_layout(rect=[0, 0.03, 1, 0.95])
# # print x.shape
# plt.savefig("figures/linear_fit_gigaplex.png",dpi=500)


#
# plt.figure("random_med_correlates")
# plt.suptitle("Set of scatter plots for correlated gene pairs, .1<r<.5")
# plt.xlabel("Log2 gene expression values")
# plt.ylabel("Frequency")
# x, y = np.where(np.logical_and(np.abs(correlations) > .1,np.abs(correlations) < .5))
# for i, pick in enumerate(np.random.randint(x.shape[0], size=20)):
#     plt.subplot(4,5,i+1)
#     plt.xlabel(header[x[pick]],size=6)
#     plt.ylabel(header[y[pick]],size=6)
#     plt.xticks(np.arange(0,14,2),size=3)
#     plt.yticks(np.arange(2,10,2),size=3)
#     plt.scatter(counts[:,x[pick]],counts[:,y[pick]],marker='x',s=.1,alpha=.1,c='b')
# plt.tight_layout(rect=[0, 0.03, 1, 0.95])
# print x.shape
# plt.savefig("figures/correlated_scatter_gigaplex(intermediate).png",dpi=500)

# high_correlation_pairings_output = open("figures/high_correlations_pairings_output.txt", mode='w')
# x, y = np.where(np.abs(correlations) > .5)
# correlations_formatted_list = map(lambda z: (header[z[0]],header[z[1]],correlations[z[0],z[1]]), zip(x,y))
# for pairing in sorted(correlations_formatted_list, key=lambda x: x[2]):
#     high_correlation_pairings_output.write("\t".join(map(lambda x: str(x), pairing)) + "\n")

# gene_cluster_at_3500_output = open("figures/gene_cluster_3500.txt", mode='w')
# cell_cluster_at_700_output = open("figures/cell_cluster_700.txt", mode='w')
#
# for i, gene in enumerate(clustered_header[3600:3700]):
#     gene_cluster_at_3500_output.write(str(gene) + '\t' + str(clustered_counts[620:840,3600+i].mean(axis = 0)) + "\n")
#
# cell_clustering_indecies = np.load('cell_clustering_indecies.npy')
#
# clustered_cells = cell_identities[cell_clustering_indecies]
#
# for cell in clustered_cells[620:840]:
#     cell_cluster_at_700_output.write('\t'.join(map(lambda x: str(x),cell)) + '\n')
#
# print clustered_cells[620:840].sum(axis=0)
# print clustered_cells.sum(axis=0)
#
#
# padded_deviation = np.zeros((deviation_matrix.shape[0],deviation_matrix.shape[1]+np.sum(np.logical_not(dropout_deviation_matrix))))
# shifted_index = 0
# for i in range(len(dropout_deviation_matrix)):
#     if dropout_deviation_matrix[i]:
#         padded_deviation[:,i] = deviation_matrix[:,i-shifted_index]
#     else:
#         padded_deviation[:,i] = np.zeros(padded_deviation.shape[0])
#         shifted_index += 1

#
# print padded_deviation.shape
# rearranged_deviation = padded_deviation[cell_clustering_indecies]
# rearranged_deviation = rearranged_deviation.T[gene_clustering_indecies].T
# rearranged_deviation = np.concatenate((rearranged_deviation,np.ones((rearranged_deviation.shape[0],1))*-7.5), axis=1)
# rearranged_deviation = np.concatenate((rearranged_deviation,np.ones((rearranged_deviation.shape[0],1))*7.5), axis=1)
#
#
# plt.figure('cluster_deviation')
# plt.title("Normal Agglomerative Clustering, Residual Values")
# image = plt.imshow(rearranged_deviation,cmap='seismic')
# plt.colorbar(image, fraction=.01)
# plt.savefig('figures/neighbor_deviation.png', dpi=700)


# gene_linked = hrc.linkage(counts.T, method='average', metric='correlation')
# gene_dendrogram = hrc.dendrogram(gene_linked,no_plot=True)
#
# cell_linked = hrc.linkage(counts, method='average', metric='cosine')
# clusterization = hrc.fcluster(cell_linked, criterion='inconsistent',t=.5,)
# cell_dendrogram = hrc.dendrogram(cell_linked,no_plot=True)
#
# fig = plt.figure("neighbor_residuals", figsize=(8,4))
# fig.suptitle("Residuals of Agglomeratively Clustered Cells")
# ax1 = fig.add_axes([.09,.1,.2,.6])
# # display_dendrogram = hrc.dendrogram(cell_linked, p=3, truncate_mode='level',orientation='left',show_contracted=True,ax=ax1)
# with plt.rc_context({'lines.linewidth':0.1}):
#     display_dendrogram = hrc.dendrogram(cell_linked,orientation='left',ax=ax1)
# ax1.set_xlim(left=1.0,right=.75)
# ax1.set_xscale('log')

# ax2 = fig.add_axes([.3,.71,.55,.2])
# display_dendrogram = hrc.dendrogram(gene_linked, p=3, truncate_mode='level',ax=ax2)
# with plt.rc_context({'lines.linewidth':0.1}):
#     display_dendrogram = hrc.dendrogram(gene_linked, ax=ax2)
# ax2.set_ylim(top=1,bottom=.75)
# ax2.set_yscale('log')
#
# print counts.shape
# print len(cell_dendrogram['leaves'])
# print len(gene_dendrogram['leaves'])


# ax3 = fig.add_axes([.3,.1,.55,.6])
#
# sorted_singly = counts[np.flip(cell_dendrogram['leaves'],0)]
# sorted_doubly = sorted_singly.T[gene_dendrogram['leaves']].T
# sorted_doubly = np.concatenate((sorted_doubly.T,np.ones((sorted_doubly.T.shape[0],1))*-10), axis=1)
# sorted_doubly = np.concatenate((sorted_doubly,np.ones((sorted_doubly.shape[0],1))*10), axis=1)
# im = ax3.imshow(rearranged_deviation, cmap='seismic', aspect='auto')
# ax3.set_xticks(np.arange(0,4773,100))
# ax3.set_yticks(np.arange(0,1656,100))
# ax3.tick_params(labelsize=6)
# plt.title("Residual Expression of Genes In Cells, Clustered Hierarchically")
# plt.xlabel("Genes")
# plt.ylabel("Cells")
# ax4 = fig.add_axes([.85,.1,.05,.6])
# ax4.set_ylim(bottom=-10,top=10)
# fig.colorbar(mappable=im, fraction=.99, ax=ax4)
# np.save("clustered_counts",sorted_doubly)
# np.save("clustered_header",header[gene_dendrogram['leaves']])
# np.save("gene_clustering_indecies", gene_dendrogram['leaves'])
# np.save("cell_clustering_indecies", cell_dendrogram['leaves'])
# plt.savefig("figures/neighbor_cluster_deviation.png", dpi=800)


#

# fig = plt.figure('clustered_residuals', figsize=(8,4))
#
# cell_linked = hrc.linkage(deviation_matrix, method='average', metric='cosine')
# clusterization = hrc.fcluster(cell_linked, criterion='inconsistent',t=.5,)
# cell_dendrogram = hrc.dendrogram(cell_linked,no_plot=True)
# gene_linked = hrc.linkage(deviation_matrix.T, method='average', metric='cosine')
# gene_dendrogram = hrc.dendrogram(gene_linked,no_plot=True)
#
#
# ax1 = fig.add_axes([.05,.1,.115,.6])
# # display_dendrogram = hrc.dendrogram(cell_linked, p=3, truncate_mode='level',orientation='left',show_contracted=True,ax=ax1)
# with plt.rc_context({'lines.linewidth':0.1}):
#     display_dendrogram = hrc.dendrogram(cell_linked,orientation='left',ax=ax1)
# ax1.set_xlim(left=1.0,right=.75)
# ax1.set_xscale('log')
#
# ax2 = fig.add_axes([.18,.71,.62,.2])
# # display_dendrogram = hrc.dendrogram(gene_linked, p=3, truncate_mode='level',ax=ax2)
# with plt.rc_context({'lines.linewidth':0.1}):
#     display_dendrogram = hrc.dendrogram(gene_linked, ax=ax2)
# ax2.set_ylim(top=1,bottom=.75)
# ax2.set_yscale('log')
#
# print len(cell_dendrogram['leaves'])
# print len(gene_dendrogram['leaves'])
#
#
# ax3 = fig.add_axes([.18,.1,.62,.6])
#
# clustered_residuals = deviation_matrix.T[gene_dendrogram['leaves']]
# clustered_residuals = clustered_residuals.T[cell_dendrogram['leaves']].T
# clustered_residuals = np.concatenate((clustered_residuals.T,np.ones((clustered_residuals.T.shape[0],1))*-7), axis=1)
# clustered_residuals = np.concatenate((clustered_residuals,np.ones((clustered_residuals.shape[0],1))*7), axis=1)
# im = ax3.imshow(clustered_residuals, cmap='seismic', aspect='auto')
# # plt.title("Residual Expression of Genes In Cells, Clustered Hierarchically")
# # plt.xlabel("Genes")
# # plt.ylabel("Cells")
# ax4 = fig.add_axes([.85,.1,.05,.6])
# # ax4.set_ylim(bottom=-10,top=10)
# fig.colorbar(mappable=im, fraction=.99, ax=ax4)
# plt.savefig("figures/clustered_residuals.png",dpi=700)
#
# np.save('residual_gene_clusters', gene_dendrogram['leaves'])
# np.save('residual_cell_clusters', cell_dendrogram['leaves'])
# np.save('clustered_residuals', clustered_residuals)

# residual_gene_cluster_out = open("figures/residual_cluster_2600_1150.txt", mode='w')
# # residual_cell_cluster_out = open("figures/cell_cluster_700.txt", mode='w')
#
# residual_header_clustering = header[residual_gene_clusters]
#
# for i, gene in enumerate(clustered_header[2600:2800]):
#     residual_gene_cluster_out.write(str(gene) + '\t' + str(clustered_residuals[1150:1250,2600+i].mean(axis = 0)) + "\n")
#
# residual_gene_cluster_out = open("figures/residual_cluster_2800_100.txt", mode='w')
#
# for i, gene in enumerate(clustered_header[2800:3000]):
#     residual_gene_cluster_out.write(str(gene) + '\t' + str(clustered_residuals[100:250,2800+i].mean(axis = 0)) + "\n")
#
# residual_gene_cluster_out = open("figures/residual_cluster_400_700.txt", mode='w')
#
# for i, gene in enumerate(clustered_header[400:800]):
#     residual_gene_cluster_out.write(str(gene) + '\t' + str(clustered_residuals[700:850,400+i].mean(axis = 0)) + "\n")
#
# residual_gene_cluster_out = open("figures/residual_cluster_1500_0.txt", mode='w')
#
# for i, gene in enumerate(clustered_header[1500:2000]):
#     residual_gene_cluster_out.write(str(gene) + '\t' + str(clustered_residuals[0:100,1500+i].mean(axis = 0)) + "\n")
#
#
# residual_gene_cluster_out = open("figures/residual_cluster_4100_1100.txt", mode='w')
#
# for i, gene in enumerate(clustered_header[4100:4300]):
#     residual_gene_cluster_out.write(str(gene) + '\t' + str(clustered_residuals[1100:1400,4100+i].mean(axis = 0)) + "\n")


fig = plt.figure('biclustered', figsize=(8,4))

model = SpectralBiclustering(n_clusters=15,n_clusters=15)

ax1 = fig.add_axes([.05,.1,.115,.6])
# display_dendrogram = hrc.dendrogram(cell_linked, p=3, truncate_mode='level',orientation='left',show_contracted=True,ax=ax1)
with plt.rc_context({'lines.linewidth':0.1}):
    display_dendrogram = hrc.dendrogram(cell_linked,orientation='left',ax=ax1)
ax1.set_xlim(left=1.0,right=.75)
ax1.set_xscale('log')

ax2 = fig.add_axes([.18,.71,.62,.2])
# display_dendrogram = hrc.dendrogram(gene_linked, p=3, truncate_mode='level',ax=ax2)
with plt.rc_context({'lines.linewidth':0.1}):
    display_dendrogram = hrc.dendrogram(gene_linked, ax=ax2)
ax2.set_ylim(top=1,bottom=.75)
ax2.set_yscale('log')

print len(cell_dendrogram['leaves'])
print len(gene_dendrogram['leaves'])


ax3 = fig.add_axes([.18,.1,.62,.6])

clustered_residuals = deviation_matrix.T[gene_dendrogram['leaves']]
clustered_residuals = clustered_residuals.T[cell_dendrogram['leaves']].T
clustered_residuals = np.concatenate((clustered_residuals.T,np.ones((clustered_residuals.T.shape[0],1))*-7), axis=1)
clustered_residuals = np.concatenate((clustered_residuals,np.ones((clustered_residuals.shape[0],1))*7), axis=1)
im = ax3.imshow(clustered_residuals, cmap='seismic', aspect='auto')
# plt.title("Residual Expression of Genes In Cells, Clustered Hierarchically")
# plt.xlabel("Genes")
# plt.ylabel("Cells")
ax4 = fig.add_axes([.85,.1,.05,.6])
# ax4.set_ylim(bottom=-10,top=10)
fig.colorbar(mappable=im, fraction=.99, ax=ax4)
plt.savefig("figures/clustered_residuals.png",dpi=700)

# fig = plt.figure("Treesplits")
#
# ax1 = fig.add_subplot(211)
#
# print fig
# print ax1
#
# split_feature = 4386
#
#
# top_half = counts[top_indecies]
#
# bottom_half = counts[bottom_indecies]
#
# top_half = top_half.T[gene_clustering_indecies].T
# bottom_half = bottom_half.T[gene_clustering_indecies].T
#
# tree_split_counts = np.zeros((top_half.shape[0]+bottom_half.shape[0],top_half.shape[1]))
#
# tree_split_counts[:top_half.shape[0]] = top_half
# tree_split_counts[top_half.shape[0]:] = bottom_half
#
# ax1.imshow(tree_split_counts, cmap='hot')
#
# ax2 = fig.add_subplot(212)
#
# resorted_counts = counts[np.argsort(counts[:,split_feature])].copy()
# resorted_counts = resorted_counts[resorted_counts[:,split_feature] > 0]
#
# feature_correlations = np.cov(resorted_counts.T)[split_feature]
#
# resorted_counts = resorted_counts.T[np.argsort(feature_correlations)].T
#
# print resorted_counts[:,list(gene_clustering_indecies).index(627)]
# print resorted_counts.shape
#
# print "Split happens at " + str(len(top_indecies))
# # for i, element in resorted_counts[:,split_feature]:
# #     if i%100 == 0:
# #         print element
#
# ax2.imshow(resorted_counts,cmap='hot', aspect='auto')
#
# fig.savefig("figures/tree_split_" + str(split_feature) + ".png" ,dpi=300)
