import numpy as np
from matplotlib import pyplot as plt

from scipy.sparse.linalg import LinearOperator

from scipy.sparse import identity, diags, csr_matrix, lil_matrix
from scipy.sparse import kron as skron
from scipy.sparse import linalg as sla
from scipy import linalg as la
from numpy.polynomial.chebyshev import Chebyshev
from numpy.polynomial.chebyshev import chebinterpolate



def create_Tz(L, J, kappa, beta):
    T = lil_matrix((4**L, 4**L), dtype=float)
    for i in range(4**L):
        for j in range(2**L):
            a = ((i & (~ (2**L-1))) + j)

            a_bin = np.array([int(bin) for bin in format(a, f'0{L+2}b')])
            i_bin = np.array([int(bin) for bin in format(i, f'0{L+2}b')])

            a_bin = a_bin * 2 - 1
            i_bin = i_bin * 2 - 1

            T[i, a] = np.exp(- kappa * J * np.sum(a_bin * i_bin))

    return T


def create_Tx(L, J, beta, kappa_1):
    T = lil_matrix((4**L, 4**L), dtype=float)
    for i in range(2**L):
        i_bin = np.array([int(bin) for bin in format(i, f'0{L+2}b')]) 
        i_bin = i_bin * 2 - 1
        en_hor = 0  
        for j in range(L):
            en_hor += - kappa_1 * J * i_bin[j] * i_bin[(j + 1) % L]

        for a in range(2**L):
            a_bin = np.array([int(bin) for bin in format(a, f'0{L+2}b')])
            a_bin = a_bin * 2 - 1

            en_vert = - np.sum(a_bin * i_bin)
            
            a = (i << L) + a
            T[a, a] += np.exp(- beta * (en_hor + en_vert))

    return T


def update_matrices(L, J, K, kappa, kappa_1, h0, h1, beta):
    Tz = create_Tz(L, J, kappa, beta)
    Tx = create_Tx(L, J, beta, kappa_1)
    Tz_1 = create_Tz(L, J, kappa_1, beta)
    Tx_1 = create_Tx(L, J, beta, kappa)
    return Tz, Tx, Tx_1, Tz_1


def vector_multiplication(v, L, Tz_1, Tx_1, Tz, Tx):
    v = Tx_1 @ Tz_1 @ Tx @ Tz @ v
    return v

L = 4
J = 1
K = 1
kappa = 1
kappa_1 = 1
h0 = 0
h1 = 0
beta = 1

Tz, Tx, Tx_1, Tz_1 = update_matrices(L, J, K, kappa, kappa_1, h0, h1, beta)
T = Tx_1 @ Tz_1
print(T.todense())
print(np.allclose(Tz_1.todense(), Tz_1.T.todense()))







            

            
            


            




