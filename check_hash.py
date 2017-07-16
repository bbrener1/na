#!/usr/bin/env python

import sys
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

def check_hash(attempt_array, target, prefix = ""):

    print "Checking for previously computed matrix output"

    if isinstance(attempt_array,str):
        try:
            attempt_array = np.load(attempt_array)
        except IOError:
            attempt_array = np.loadtxt(attempt_array, dtype=str)
        except:
            print "What the hell kind of array did you pass to the presolve checker?"


    try:
        hash_file = open(prefix+"presolve_hash_dictionary.txt", mode='r')
    except IOError as e:
        if e.errno == 2:
            return False
        else:
            "Something's wrong with the hash file"

    current_presolves = {x.split()[0]:x.split()[1] for x in hash_file}
    try:
        print current_presolves[target]
        print hash(str(attempt_array))
        print int(current_presolves[target]) == hash(str(attempt_array))
    except KeyError:
        return False

    return int(current_presolves[target]) == hash(str(attempt_array))

def write_hash(write_array, target, prefix = ""):

    print "Writing hash for computation"

    if isinstance(write_array,str):
        try:
            write_array = np.load(write_array)
        except IOError:
            write_array = np.loadtxt(write_array, dtype=str)
        except:
            print "What the hell kind of array did you pass to the presolve checker?"



    try:
        hash_file = open(prefix+"presolve_hash_dictionary.txt", mode='r+')
    except IOError:
        hash_file = open(prefix+"presolve_hash_dictionary.txt", mode='w+')
    print "write hash debug"
    current_presolves = {x.split()[0]:x.split()[1] for x in hash_file}
    current_presolves[target] = hash(str(write_array))
    print current_presolves
    hash_file.seek(0)
    for key in current_presolves:
        hash_file.write(str(key) + "\t" + str(current_presolves[key]) + "\n")
