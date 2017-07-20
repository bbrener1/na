#!/usr/bin/env python

import numpy as np
import sys

header = np.load(sys.argv[1])

index = np.argwhere(header == sys.argv[3])[0][0]

print index

corr = np.load(sys.argv[2])

np.savetxt(sys.argv[4], header[corr[index] > .2] ,fmt='%s')
