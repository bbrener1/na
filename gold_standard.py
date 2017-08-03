#!/usr/bin/env python

import sys
import os
import numpy as np

import pyensembl as en

import check_hash as chk

def translate_trrust():

    if presolve != None:
        if type(presolve) is str:
            network_matrix = np.load(prefix + presolve)
            return network_matrix.astype(dtype=bool)
        if isinstance(presolve,np.ndarray):
            return presolve
        else:
            print "Presolve not a file address or numpy ndarray"
            raise ValueError


    print labels[:50]

    ens_obj = en.EnsemblRelease(species=species)
    error_count = 0

    # for i in range(len(labels)):
    #     try:
    #         labels[i] = ens_obj.gene_name_of_gene_id(labels[i])
    #     except ValueError:
    #         labels[i] = "error"
    #         error_count += 1

    if isinstance(labels, str):
        labels = np.load(labels)

    labels = [x.upper() for x in labels]

    if not isinstance(in_file, file):
        in_data = open(in_file)
    else:
        in_data = in_file

    set_of_interest = set(labels)
    try:
        set_of_interest.remove("error")
    except:
        pass

    network_matrix = np.zeros((len(labels),len(labels)))

    print "GOLD STANDARD DEBUG"
    print len(set_of_interest)
    print list(set_of_interest)[:50]
    print labels[:50]
    print type(in_data)
    print in_data.readline()
    in_data.seek(0)

    error_library = []

    for line in in_data:
        if line.split()[0].upper() in set_of_interest:
            try:
                network_matrix[labels.index(line.split()[0].upper()),labels.index(line.split()[1].upper())] = 1
                # print line
            except ValueError:
                # error_library.append(line)
                # print "GOLD ERROR"
                # print line
                # print line.split()[0].upper()
                # print line.split()[2].upper()
                # print labels.index(line.split()[0].upper())
                # print labels.index(line.split()[2].upper())
                # print labels[labels.index(line.split()[0].upper())]
                # print labels[labels.index(line.split()[2].upper())]
                try:
                    labels[27:].index(line.split()[0])
                except ValueError:
                    # print "Couldn't find " + line.split()[0]
                    error_library.append(str((line,line.split()[0])))
                try:
                    labels[27:].index(line.split()[1])
                except ValueError:
                    # print "Couldn't find " + line.split()[2]
                    error_library.append(str((line,line.split()[1])))
                error_count += 1
                continue

    np.save("directional_matrix", network_matrix)

    network_matrix = np.logical_or(network_matrix, network_matrix.T)
    network_matrix = network_matrix - np.diag(np.diag(network_matrix))


    np.save(prefix + "gold_network",network_matrix)
    chk.write_hash(in_file, "gold_network.npy",prefix)

    np.savetxt(prefix + "errors_gold_network.txt",np.asarray(error_library,dtype=str),fmt = '%s')

    print "Gold Standard Debug"
    print np.sum(network_matrix.flatten())
    print "GOLD STANDARD ERRORS"
    print error_count
    print error_library[:50]

    return network_matrix


def translate_gold_standard(in_file,labels,species='mouse',presolve=None, prefix= ""):



    if presolve != None:
        if type(presolve) is str:
            network_matrix = np.load(prefix + presolve)
            return network_matrix.astype(dtype=bool)
        if isinstance(presolve,np.ndarray):
            return presolve
        else:
            print "Presolve not a file address or numpy ndarray"
            raise ValueError


    print labels[:50]

    ens_obj = en.EnsemblRelease(species=species)
    error_count = 0

    # for i in range(len(labels)):
    #     try:
    #         labels[i] = ens_obj.gene_name_of_gene_id(labels[i])
    #     except ValueError:
    #         labels[i] = "error"
    #         error_count += 1

    if isinstance(labels, str):
        labels = np.load(labels)

    labels = [x.upper() for x in labels]

    if not isinstance(in_file, file):
        in_data = open(in_file)
    else:
        in_data = in_file

    set_of_interest = set(labels)
    try:
        set_of_interest.remove("error")
    except:
        pass

    network_matrix = np.zeros((len(labels),len(labels)))

    print "GOLD STANDARD DEBUG"
    print len(set_of_interest)
    print list(set_of_interest)[:50]
    print labels[:50]
    print type(in_data)
    print in_data.readline()
    in_data.seek(0)

    error_library = []

    for line in in_data:
        if line.split()[0].upper() in set_of_interest:
            try:
                network_matrix[labels.index(line.split()[0].upper()),labels.index(line.split()[2].upper())] = 1
                # print line
            except ValueError:
                # error_library.append(line)
                # print "GOLD ERROR"
                # print line
                # print line.split()[0].upper()
                # print line.split()[2].upper()
                # print labels.index(line.split()[0].upper())
                # print labels.index(line.split()[2].upper())
                # print labels[labels.index(line.split()[0].upper())]
                # print labels[labels.index(line.split()[2].upper())]
                try:
                    labels[27:].index(line.split()[0])
                except ValueError:
                    # print "Couldn't find " + line.split()[0]
                    error_library.append(str((line,line.split()[0])))
                try:
                    labels[27:].index(line.split()[2])
                except ValueError:
                    # print "Couldn't find " + line.split()[2]
                    error_library.append(str((line,line.split()[2])))
                error_count += 1
                continue

    np.save("directional_matrix", network_matrix)

    network_matrix = np.logical_or(network_matrix, network_matrix.T)
    network_matrix = network_matrix - np.diag(np.diag(network_matrix))


    np.save(prefix + "gold_network",network_matrix)
    chk.write_hash(in_file, "gold_network.npy",prefix)

    np.savetxt(prefix + "errors_gold_network.txt",np.asarray(error_library,dtype=str),fmt = '%s')

    print "Gold Standard Debug"
    print np.sum(network_matrix.flatten())
    print "GOLD STANDARD ERRORS"
    print error_count
    print error_library[:50]

    return network_matrix

def trrust_wrapper():

    prefix = sys.argv[1]

    if len(sys.argv) > 2:
        regnetworkweborg = sys.argv[2]
    else:
        regnetworkweborg = prefix+"trrust_rawdata.mouse.tsv"

    if len(sys.argv)>3:
        species = sys.argv[3]
    else:
        species = "mouse"

    if len(sys.argv)>4:
        header = np.load(sys.argv[4])
    else:
        header = np.load(prefix+"header_backup.npy")



    if os.path.isfile(prefix + "gold_network.npy"):

        if chk.check_hash(regnetworkweborg, "gold_network.npy", prefix= prefix):
            translate_gold_standard(regnetworkweborg, labels = header, species = species, presolve= "gold_network.npy", prefix = prefix)
        else:
            translate_gold_standard(regnetworkweborg, labels = header, species = species, prefix = prefix)

    else:
        translate_gold_standard(regnetworkweborg, labels = header, species = species, prefix = prefix)



def main():

    prefix = sys.argv[1]

    if len(sys.argv) > 2:
        regnetworkweborg = sys.argv[2]
    else:
        regnetworkweborg = prefix+"mouse.source"

    if len(sys.argv)>3:
        species = sys.argv[3]
    else:
        species = "mouse"

    if len(sys.argv)>4:
        header = np.load(sys.argv[4])
    else:
        header = np.load(prefix+"header_backup.npy")



    if os.path.isfile(prefix + "gold_network.npy"):

        if chk.check_hash(regnetworkweborg, "gold_network.npy", prefix= prefix):
            translate_gold_standard(regnetworkweborg, labels = header, species = species, presolve= "gold_network.npy", prefix = prefix)
        else:
            translate_gold_standard(regnetworkweborg, labels = header, species = species, prefix = prefix)

    else:
        translate_gold_standard(regnetworkweborg, labels = header, species = species, prefix = prefix)







if __name__ == "__main__":
    main()
