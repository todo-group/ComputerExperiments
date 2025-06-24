#!/usr/bin/python

# about: read matrix data (cmatrix.h format) and print it to the console
# usage: python read_matrix.py <input_file>

# requires numpy libraries

import sys
import numpy as np

f = open(sys.argv[1], "r")
line = f.readlines()
m, n = np.uint(line[0].split())
A = np.zeros((m, n))
for i in range(m):
    A[i, :] = np.float64(line[i + 1].split())
print(A)
