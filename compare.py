#!/usr/bin/env python

from __future__ import division
import numpy as np
import sys


def fair_compare(inf_network, directional_matrix, output = None):
    if output = None:
        output = sys.stdout

    

def compare(network1, network2, output=None):

    if output == None:
        output = sys.stdout

    print "REGULAR======================="

    n1_edges = np.sum((np.abs(network1) > 0).flatten())
    n2_edges = np.sum((np.abs(network2) > 0).flatten())
    intersect = np.sum(np.logical_and((network1 != 0),(network2 != 0))).flatten()
    possible = network1.shape[0] * network1.shape[1]

    output.write(str(n1_edges) + " edges present in network 1, (" + str(n1_edges/possible) + "%)\n")
    output.write("Network shape: " + str(network1.shape) + '\n')
    output.write(str(n2_edges) + " edges present in network 2\n")
    output.write("Network shape: " + str(network2.shape) + '\n')
    output.write(str(intersect) + " edges shared between the networks, (" + str(np.sum(intersect)/n1_edges) + "%% of network 1 (rate of true positives),  " + str(np.sum(intersect)/n2_edges) + " %% of network 2, (capture rate)\n")

def restrictive_compare(network1,network2,output=None):
    if output == None:
        output = sys.stdout

    print "RESTRICTIVE==================="

    mask1a = np.sum(np.abs(network1), axis = 1) > 0
    mask1b = np.sum(np.abs(network1), axis = 0) > 0
    mask1 = np.logical_or(mask1a,mask1b)

    mask2a = np.sum(np.abs(network2), axis = 1) > 0
    mask2b = np.sum(np.abs(network2), axis = 0) > 0
    mask2 = np.logical_or(mask2a,mask2b)

    print "Mask debug"
    print np.sum(mask1)
    print np.sum(mask2)
    print mask1.shape
    print mask2.shape


    mask_c = mask2

    n1_edges = np.sum((np.abs(network1[mask_c].T[mask_c].T) > 0).flatten())
    n2_edges = np.sum((np.abs(network2[mask_c].T[mask_c].T) > 0).flatten())
    intersect = np.sum(np.logical_and((network1[mask_c].T[mask_c].T != 0),(network2[mask_c].T[mask_c].T != 0))).flatten()
    possible = np.sum(mask_c)**2

    output.write("When comparing only genes that the gold standard connects:\n")

    output.write(str(n1_edges) + " edges present in network 1, (" + str(n1_edges/possible) + "%)\n")
    output.write("Network shape: " + str(network1[mask_c].T[mask_c].T.shape) + '\n')
    output.write(str(n2_edges) + " edges present in network 2\n")
    output.write("Network shape: " + str(network2[mask_c].T[mask_c].T.shape) + '\n')

    output.write(str(intersect) + " edges shared between the networks, (" + str(np.sum(intersect)/n1_edges) + "%% of network 1 (rate of true positives),  " + str(np.sum(intersect)/n2_edges) + " %% of network 2, (capture rate)\n")

    true_positive = intersect/n2_edges
    false_positive = (n1_edges-intersect)/(network2.shape[0]*network2.shape[1]-n2_edges)

    precision = intersect/n1_edges
    recall = intersect/n2_edges

    enrichment = recall / (n1_edges/possible)

    return true_positive,false_positive, precision, recall, enrichment

    # output.write(str(intersect) + " edges shared between the networks, (" + str(np.sum(intersect)/n2_edges) + ")\n")
