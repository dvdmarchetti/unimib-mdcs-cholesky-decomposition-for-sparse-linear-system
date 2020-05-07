import time
import psutil
import os

import numpy as np
import pandas as pd
from pprint import pprint
import scipy
import scipy.linalg   # SciPy Linear Algebra Library
from scipy import array, linalg, dot
from scipy.io import mmread


def cholesky(A):
    """
    Reference: https://rosettacode.org/wiki/Cholesky_decomposition#Python
    """
    L = [[0.0] * len(A) for _ in range(len(A))]
    for i, (Ai, Li) in enumerate(zip(A, L)):
        for j, Lj in enumerate(L[:i+1]):
            s = sum(Li[k] * Lj[k] for k in range(j))
            Li[j] = sqrt(Ai[i] - s) if (i == j) else \
                      (1.0 / Lj[j] * (Ai[j] - s))
    return L

# python -m pip install --user numpy scipy matplotlib ipython jupyter pandas sympy nose

matlabroot_ans = 'D:\Program Files\MATLAB\R2020a'

#import matlab.engine
#eng = matlab.engine.start_matlab()
#eng.doc(nargout=0)


A = np.array([[1,-2j],[2j,5]])
L = linalg.cholesky(A, lower=True)
#print(L)



a = mmread('matrix\ex15.mtx')

#<3x3 sparse matrix of type '<type 'numpy.float64'>'
#    with 9 stored elements in COOrdinate format>
#and one can then transform it into a dense matrix

a.todense() #it works
#print("todense: ")
#print(a)

a.toarray()
#print("toarray: ")
#print(a)

#start mesuring #####
start_time = time.time()    #time
process = psutil.Process(os.getpid()) #memory
####################
matrix = scipy.io.mmread('matrix\ex15.mtx')
#cholesky(matrix)
m = np.empty((5,5)) #The empty function creates an array. Its initial content is random and depends on the state of the memory.

eig = np.linalg.eig(m)
v = eig[0]
print(v)
ch = True   #if eigen > 0 --> you can cholesky
for element in v:
    if (element < 0):
        ch = False
if(ch):
    print(cholesky(m))
#https://stackoverflow.com/questions/16699163/what-is-the-difference-between-cholesky-in-numpy-and-scipy
#https://www.youtube.com/watch?v=4SWMzENcgSE

A = np.array([[6, 3, 4, 8], [3, 6, 5, 1], [4, 5, 10, 7], [8, 1, 7, 25]])
L = scipy.linalg.cholesky(A, lower=True)
U = scipy.linalg.cholesky(A, lower=False)

print ("A:")
pprint(A)

print ("L:")
pprint(L)

print ("U:")
pprint(U)

"""
L = np.linalg.cholesky(m)
np.dot(L, L.T.conj()) # verify that L * L.H = A
m = [[1,-2j],[2j,5]] # what happens if A is only array_like?
np.linalg.cholesky(m) # an ndarray object is returned
 # But a matrix object is returned if A is a matrix object
np.linalg.cholesky(np.matrix(m))


matrix = matrix.todense()
print(len(matrix), len(matrix[0]) )
eig = np.linalg.eigvals(matrix)

full_df = pd.DataFrame(matrix.todense()).T
print(full_df)
#L = np.linalg.cholesky(matrix)
"""


####END
print("time = --- %s seconds ---" % (time.time() - start_time));
print("mem = ", process.memory_info().rss, 'B')  # mem in bytes
#############
