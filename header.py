#!/usr/bin/env python

import sys
import numpy as np
import pyensembl as en

import os

import check_hash as chk

def translate_10x_header(input_file, species="mouse", prefix = "", test = False):

    print prefix
    ens = en.EnsemblRelease(species=species)

    if type(input_file) is str:
        in_data = open(input_file)
    if type(input_file) is file:
        in_data = input_file
    if type(in_data) is not file:
        raise ValueError

    header = [x.split()[0] for x in in_data]

    print "Header debug 1"
    print header[:10]
    print len(header)


    error_count = 0
    error_list = []

    if not test:
        for i in range(len(header)):
            try:
                header[i] = ens.gene_name_of_gene_id(header[i])
            except ValueError:
                error_list.append(header[i])
                header[i] = "error"
                error_count += 1
    else:
        for i in range(len(header)):
            header[i] = ens.gene_name_of_gene_id(header[i])



    print "Header debug 2"
    print header[:10]
    print len(header)
    print error_count
    print error_list[:10]

    counts = scipy.io.mmread("/".join(in_data.name.split("/")[:-1])+"/matrix.mtx").todense()
    np.save(prefix+"counts",counts)

    try:
        float(header[0])
        raise SyntaxError
    except ValueError:

        np.save(prefix + "header_backup",np.asarray(header,dtype=str))

        print header[:28]

        return header

def translate_header(input_file, species = "mouse", prefix = "", test = False):

    ens = en.EnsemblRelease(species=species)

    if type(input_file) is str:
        in_data = open(input_file)
    if type(input_file) is file:
        in_data = input_file
    if type(in_data) is not file:
        raise ValueError

    header = in_data.readline().split()

    if len(header) < 3:
        try:
            return translate_10x_header(input_file,species=species, prefix = prefix, test = test)
        except SyntaxError:
            print "Malformed matrix file? Printing first line"
            print header
            raise IOError
    else:
        counts = np.loadtxt(input_file,skiprows=1,usecols=np.arange(27,len(header)))
        np.save(prefix+"counts.npy", counts)

    error_count = 0
    error_list = []

    if not test:

        for i in range(len(header)):
            try:
                header[i] = ens.gene_name_of_gene_id(header[i])
            except ValueError:
                error_list.append(header[i])
                header[i] = "error"
                error_count += 1

    else:
        for i in range(len(header)):
            header[i] = ens.gene_name_of_gene_id(header[i])


    trunc_header = header[27:]

    print error_count
    print header[:28]
    print header[:50]
    print trunc_header[:20]
    print len(trunc_header)


    np.save(prefix + "header_backup",np.asarray(trunc_header,dtype=str))
    chk.write_hash(input_file,"header_backup.npy", prefix)

    # header_backup = open("header_backup.txt",mode='w')
    # for element in trunc_header:
    #     header_backup.write(element + "\n")
    # header_backup.close()

    return header,trunc_header


def main():

    prefix = sys.argv[1]

    if len(sys.argv)>2:
        target = sys.argv[2]

    if len(sys.argv)>3:
        species = sys.argv[3]
    else:
        species = "mouse"

    if "-t" in sys.argv:
        test = True
    else:
        test = False


    if os.path.isfile(prefix + "header_backup.npy"):
        if not chk.check_hash(target,"header_backup.npy", prefix):

            return translate_header(target, species = species, prefix = prefix, test = test)
    else:
        return translate_header(target, species = species, prefix = prefix, test = test)








if __name__ == "__main__":
    main()
