#!/usr/bin/env python
from __future__ import division
import sys
import numpy as np
from numpy.linalg import inv
from sklearn.linear_model import LinearRegression
from sklearn.mixture import GaussianMixture
from sklearn.covariance import GraphLasso
from scipy.stats import pearsonr
import scipy.stats
import scipy.io
import pyensembl as en
from sklearn.decomposition import PCA
from sklearn.covariance import graph_lasso
from sklearn.preprocessing import scale
from sklearn.metrics import r2_score
from sklearn.cluster import AgglomerativeClustering
import scipy.spatial.distance as spt
import random
import numpy.random as nprnd
import networkx as nx



import warnings
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.animation as manimation

def synthesize_header(size, species = "mouse", out_format = "10x"):

    ens = en.EnsemblRelease(species=species)
    total_ids = ens.gene_ids()
    header_list = random.sample(total_ids,size)
    header = [ens.gene_name_of_gene_id(x) for x in header_list]
    output = open("synthetic_golden_tests/header.txt", mode = 'w')

    if out_format == "10x":
        for gene in header:
            output.write(gene + "\t")
        output.close()
        return header, output

    if out_format == "vision":
        header_str = "\t".join(header_list) + "\n"
        output.write(header_str)
        output.close()
        return header_str

def synthesize_network(genes):
    network = np.identity(genes)*.8
    for i in range(genes):
        connected = np.arange(network.shape[0])[network[i] > 0]
        j = nprnd.choice(np.delete(np.arange(genes),i),1,p = np.delete(np.sum(network, axis=1) + np.ones(genes), i)/(np.sum(network)+genes-1-np.sum(network[i])))
        while j in connected:
            j = nprnd.choice(np.delete(np.arange(genes),i),1,p = np.delete(np.sum(network, axis=1) + np.ones(genes), i)/(np.sum(network)+genes-1-np.sum(network[i])))

        coefficient = np.random.ranf()

        network[j,i] = coefficient

    for i in range(int(genes/5)):
        i = int(random.random()*genes)

        connected = np.arange(network.shape[0])[network[i] > 0]

        k = nprnd.choice(np.delete(np.arange(genes),i),1,p = np.delete(np.sum(network, axis=1) + np.ones(genes), i)/(np.sum(network)+genes-1-np.sum(network[i])))
        while k in connected:
            k = nprnd.choice(np.delete(np.arange(genes),i),1,p = np.delete(np.sum(network, axis=1) + np.ones(genes), i)/(np.sum(network)+genes-1-np.sum(network[i])))

        coefficient = np.random.ranf()

        network[k,i] = coefficient



    neg_rand_mask = nprnd.randint(0,2,network.shape,dtype=bool)
    # neg_rand_mask - np.diagonal(neg_rand_mask)
    network[neg_rand_mask] = network[neg_rand_mask] * -1
    print "Network Debug"
    print np.sum(network.flatten())
    print network
    np.savetxt("synthetic_golden_tests/synthetic_golden_network.txt", network)


    graph_a_network(network, "synthetic_golden_tests/")

    return network


def degree_histogram(network):
    plt.figure()
    plt.hist(np.sum(np.abs(network) > 0,axis=1) + np.sum(np.abs(network) > 0,axis=0),np.arange(10), log=True)
    plt.savefig("synthetic_golden_tests/synthetic_golden_degree_histogram.png")
    print "Histogram Debug"
    print np.sum(np.abs(network) > 0,axis=0)
    print np.sum(np.abs(network) > 0,axis=1)
    print np.sum(np.abs(network) > 0,axis=0) + np.sum(np.abs(network) > 0,axis=1)
    print np.sum(np.abs(network) > 0)

    plt.figure()
    plt.hist(np.sum(np.abs(network) > 0,axis=1) ,np.arange(10), log=True)
    plt.savefig("synthetic_golden_tests/synthetic_golden_degree_histogram_outbound.png")


def simulate_counts(network, cells, header = None, initial_counts = None, cycles = 200):
    if initial_counts == None:
        counts = nprnd.randn(cells,network.shape[0])*50
        # cov_counts = np.copy(counts)

    divergence = []
    network_activation = []
    summary_divergence = []

    variant_identity = np.identity(network.shape[0])*10

    FFMpegWriter = manimation.writers['ffmpeg']
    metadata = dict(title="GraphProgress",artist='Matplotlib',comment="cross your fingers")
    writer = FFMpegWriter(fps=15,metadata=metadata)

    drawing_connectivity_matrix = scipy.sparse.csc_matrix((np.abs(network * 100)).astype(dtype=int))
    graph = nx.from_scipy_sparse_matrix(drawing_connectivity_matrix, create_using=nx.DiGraph())
    layout = nx.spring_layout(graph)

    fig = plt.figure("mp4")
    plt.xlim(0,1)
    plt.ylim(0,1)
    with writer.saving(fig,"synthetic_golden_tests/graph_progress.mp4", 302):
        for i in range(cycles):
            initial = np.copy(counts)
            for j in range(cells):
                if i == 56:
                    if j%5 == 0:
                        print "Input counts"
                        print counts[j]
                        print "Input network"
                        print network
                        print "Output"
                        print np.dot(counts[j],network)

                counts[j] = nprnd.multivariate_normal(np.dot(counts[j],network).flatten(),variant_identity)
                # cov_counts[j] = nprnd.multivariate_normal(# cov_counts[j],network)

                if i == 56:
                    if j%5 == 0:

                        print "Recorded?"
                        print counts[j]


                        #     counts = np.reshape(nprnd.multivariate_normal(np.dot(counts,network*20).flatten(),np.identity(network.shape[0]*cells)*10),(cells,network.shape[0]))
            counts[counts < 0] = 0
            counts[counts > 100] = 100

            # cov_counts[# cov_counts < 0] = 0
            # cov_counts[# cov_counts > 100] = 100

            if i == 56:
                print "Do counts get flattened?"
                print counts[5]

            divergence.append(pearsonr(initial.flatten(), counts.flatten())[0])
            network_activation.append(np.sum(np.abs(counts) > 0)/(counts.shape[0]*counts.shape[1]))
            print i

            network_frame(network, np.mean(counts,axis=0)*10, writer, "mp4", layout)
            # network_frame(network, np.array([0,1,100,1000,5]),writer,'mp4',layout)

    plt.figure()
    plt.plot(divergence, alpha = .25,label="Divergence")
    # plt.plot(network_activation, alpha = .25,label="Activation")
    # plt.plot(summary_divergence, alpha = .25, label="Summation")
    plt.legend()
    plt.savefig("synthetic_golden_tests/divergence_progress.png")

    if header == None
        np.savetxt("synthetic_golden_tests/synthetic_golden_counts.txt",counts)
    else:
        counts_vision = open("synthetic_golden_counts.txt", mode='w')
        counts_vision.write(header)
        for line in counts:
            counts_vision.write(str(line))
        counts_vision.close()
    # np.savetxt("synthetic_golden_tests/synthetic_golden_cov_counts.txt",# cov_counts)


    gold_summary = open("synthetic_golden_tests/network_summary.txt",mode='w')
    gold_summary.write("Gene Count File Dimensions: ")
    gold_summary.write(str(network.shape[0]) + " genes, " + str(cells) + " cells.\n")
    gold_summary.write("Mean divergence for last 50 cycles after cycling through network (default 200 cycles): " + str(np.mean(divergence[-50:])) + "\n")
    gold_summary.write("Mean node activation for last 50 cycles: " + str(np.mean(network_activation[-50:])) + "\n")
    gold_summary.close()


    print "Counts Simulation"
    return counts

def graph_a_network(connectivity_matrix,name):
    # graph = nx.Graph()
    # graph.add_nodes_from(range(len(header))))
    # for i in range(connectivity_matrix.shape[0]):
    #     for j in range(connectivity_matrix.shape[1]):
    #         if connectivity_matrix[i,j] > 0:
    #             graph.add_edge(i,j)

    print "Trying to graph"

    drawing_connectivity_matrix = scipy.sparse.csc_matrix((np.abs(connectivity_matrix)*100).astype(dtype=int))


    print "Converted to sparse"

    graph = nx.from_scipy_sparse_matrix(drawing_connectivity_matrix, create_using=nx.DiGraph())



    edge_weight_list = []
    for edge in graph.edges_iter():
        edge_weight_list.append(connectivity_matrix[edge[0],edge[1]])

    node_degree_list = []
    degree_dict = graph.degree(weight='weight')
    for i in range(len(degree_dict)):
        node_degree_list.append(degree_dict[i])
    #
    # print "GRAPH DEBUG"
    #
    # print "\n"
    # print "\n"

    print edge_weight_list[:20]
    print node_degree_list

    print "Converted to graph"

    ed_cm = plt.get_cmap('bwr')

    plt.figure()
    nx.draw_networkx(graph, node_size = node_degree_list , edge_color=edge_weight_list, edge_cmap=ed_cm, arrows=True )
    plt.savefig(name+"network_drawing.png")

    print "Plotted graph"

    plt.figure()
    nx.draw_shell(graph)
    plt.savefig(name+"shell_network.png")

    print "Plotted second graph"

    plt.figure()
    nx.draw_spring(graph)
    plt.savefig(name+"graphviz_network.png")

    print "Plotted third graph"

def network_frame(connectivity_matrix, counts, writer, figure, layout):

        # print "Trying to graph"

        drawing_connectivity_matrix = scipy.sparse.csc_matrix((np.abs(connectivity_matrix)*100).astype(dtype=int))


        # print "Converted to sparse"

        graph = nx.from_scipy_sparse_matrix(drawing_connectivity_matrix, create_using=nx.DiGraph())

        edge_weight_list = []
        for edge in graph.edges_iter():
            edge_weight_list.append(connectivity_matrix[edge[0],edge[1]])

        # node_degree_list = []
        # degree_dict = graph.degree(weight='weight')
        # for i in range(len(degree_dict)):
        #     node_degree_list.append(degree_dict[i])
        #
        # print "GRAPH DEBUG"
        #
        # print "\n"
        # print "\n"

        # print edge_weight_list[:20]
        # print node_degree_list

        # print "Converted to graph"

        ed_cm = plt.get_cmap('bwr')



        plt.figure(figure)
        plt.clf()
        plt.cla()
        nx.draw_networkx(graph, pos = layout, node_size = (counts/10) , edge_color=edge_weight_list, font_size = 2, linewidths = .5, edge_cmap=ed_cm, arrows=True )
        writer.grab_frame()
        # print "Grabbed Frame"
        # print counts





def main():

    if len(sys.argv) <= 1:
        genes = 100
        cells = 100
    else:
        genes = int(sys.argv[1])
        cells = int(sys.argv[2])
        fmt = sys.argv[3]

    header = synthesize_header(genes, out_format = fmt)
    network = synthesize_network(genes)
    degree_histogram(network)
    simulate_counts(network,cells, header,cycles=200)



if __name__ == "__main__":
    main()
