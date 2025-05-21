import numpy as np
from matplotlib import pyplot as plt

from scipy.sparse.linalg import LinearOperator

from scipy.sparse import identity, diags, csr_matrix, lil_matrix
from scipy.sparse import kron as skron
from scipy.sparse import linalg as sla
from scipy import linalg as la
from numpy.polynomial.chebyshev import Chebyshev
from numpy.polynomial.chebyshev import chebinterpolate



def Dhalf_horizontal_TM(L, J, beta, h0, h1, pbc = True, fmt = 'csr'):
## site-site horizontal interaction matrix (finite space-direction L)
## only diagonal, simplified construction

    # Generate all possible spin states for L spins
    states = np.array([list(format(i, f'0{L}b')) for i in range(2 ** L)], dtype=int)
    
    # Initialize the diagonal of the transfer matrix with zeros
    size = 2 ** L
    Diagonal = np.zeros(size, dtype=float)

    
    for i in range(size):
        spins = 2 * states[i] - 1

        energy = 0.0

        for k in range(L-1):
            energy += spins[k] * spins[k+1]
            if k % 2 == 0:
                energy += h0 * spins[k]
            else:
                energy += h1 * spins[k]
    
        if ( pbc ):
            energy += spins[L-1] * spins[0]

        Diagonal[i] = energy
         
    Dhalf = diags( np.exp( beta * J/2 * Diagonal ), format = fmt  )

    return Dhalf

def create_T(L, kappa, kappa_1, K, beta):
    #T_z,1 - contribution to the interlayer interaction of first spin
    T = lil_matrix((2**L, 2**L), dtype=float)
    for i in range (2**L):
        a = (i & (~(1 << (L+1)))) + ((2 & i) << L) 

        a_bin = np.array([int(bin) for bin in format(a, f'0{L+2}b')])
        i_bin = np.array([int(bin) for bin in format(i, f'0{L+2}b')])
        
        a_bin = a_bin * 2 - 1
        i_bin = i_bin * 2 - 1
        
        T[i, a] = np.exp(beta * K * a_bin[0] * i_bin[0] - beta * kappa * K * i_bin[0] * a_bin[-2] - beta * kappa_1 * K * i_bin[0] * a_bin[1])

        a = a ^ 2
        a_bin[-2] *= -1
        
        T[i, a] = np.exp(beta * K * a_bin[0] * i_bin[0] - beta * kappa * K * i_bin[0] * a_bin[-2] - beta * kappa_1 * K * i_bin[0] * a_bin[1])
        
    return T


T = create_T(4, 1, 1, 1, 1)

print(np.allclose(T.todense(), T.T.todense()))





