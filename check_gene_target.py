#!/usr/bin/env python

import numpy as np
import sys

print len(sys.argv)
print sys.argv

if len(sys.argv) > 5:
    threshold = float(sys.argv[5])
else:
    threshold = .2

header = np.load(sys.argv[1])



index_array = np.argwhere(np.char.find(header,sys.argv[3])+1)

print index_array
print header[index_array]

corr = np.load(sys.argv[2])
for index in index_array:
    np.savetxt(sys.argv[4] + header[index[0]] + ".txt", header[corr[index][0] > threshold] ,fmt='%s')
