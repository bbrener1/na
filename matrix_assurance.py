#!/usr/bin/env python

import sys
import numpy as np

def matrix_assurance(arg):

    if isinstance(arg,str):

        # print "Matrix assurance debug 1"

        # print arg.split(".")[-1]
        # print arg.split(".")[-1] == "npy"

        if arg.split(".")[-1] == "npy":

            # print "Matrix assurance debug 3"
            return np.load(arg)

        elif arg.split(".")[-1] == "txt":

            # print "Matrix assurance debug 4"

            try:
                return np.loadtxt(arg)
            except ValueError:
                first_line_reader = open(arg)
                first_line_dim = len(first_line_reader.readline().split())
                return np.loadtxt(arg,skiprows=1,usecols=np.arange(27,first_line_dim))

    if isinstance(arg,np.ndarray):
        return arg

    # print "Matrix assurance debug 2"
