import math
import os
import time
import tracemalloc


#!python -m pip install --user numpy scipy matplotlib ipython jupyter pandas sympy nose
import numpy as np
import pandas as pd
from pprint import pprint
# import psutil
import scipy
from scipy.sparse import linalg as splinalg
from scipy import io, linalg, dot, sparse
from sksparse.cholmod import cholesky


# def cholesky(A):
#     """py
#     Reference: https://rosettacode.org/wiki/Cholesky_decomposition#Python
#     """
#     L = [[0.0] * len(A) for _ in range(len(A))]
#     for i, (Ai, Li) in enumerate(zip(A, L)):
#         for j, Lj in enumerate(L[:i+1]):
#             s = sum(Li[k] * Lj[k] for k in range(j))
#             if (i == j):
#                 Li[j] = sqrt(Ai[i] - s)
#             else:
#                 Li[j] = (1.0 / Lj[j] * (Ai[j] - s))
#     return L


def sparse_cholesky(A):
    n = A.shape[0]
    LU = splinalg.splu(A, diag_pivot_thresh=0)

    # if (LU.perm_r == np.arange(n)).all() and (LU.U.diagonal() > 0).all():
    return LU.L.dot(scipy.sparse.diags(LU.U.diagonal()**0.5))
    # else:
        # raise ValueError('The matrix is not positive definite')


def solve_matrix(filename):
    print('Load matrix')
    A = scipy.io.mmread('matrix/ex15.mtx')

    # data = scipy.io.loadmat(filename)
    # A = data['Problem']['A'][0][0]
    # print(A)
    print('Calculate b')
    x_es = np.ones(A.shape[0])
    b = scipy.dot(A, x_es)

    print('Solve')
    # L = sparse_cholesky(A)
    factor = cholesky(A)
    print(factor)
    print('Calculate solution')
    x_ap = factor(b)
    # x_ap = scipy.dot(scipy.dot(L, L.T.conj()), b)
    print(x_ap)

    return {
        'size': A.shape[0],
        'proc_memory_delta': 0,
        'chol_time': 0,
        'relative_error': 0,
    }


def write_to_csv(filename, contents):
    pass

# # start mesuring
# # memory method 2
# tracemalloc.start()
# # memory and time
# start_time = time.time()    #time
# process = psutil.Process(os.getpid()) #memory
# # IN ALTERNATIVA https://pypi.org/project/memory-profiler/ ############

dirname = os.path.dirname(__file__)
solve_matrix(os.path.join(dirname, 'matrix', 'ex15.mat'))

# print("time = --- %s seconds ---" % (time.time() - start_time));
# #print(process.memory_info()) ##detailed
# print("mem = ", process.memory_info().rss, 'B')  # mem in bytes
# ###memory method 2
# current, peak = tracemalloc.get_traced_memory()
# print(f"Current memory usage is {current / 10**6}MB; Peak was {peak / 10**6}MB")
# tracemalloc.stop()
