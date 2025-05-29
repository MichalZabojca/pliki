import numpy as np

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
        a = i & (~(1 << (L-1)))
        a_bin = np.array([int(bin) for bin in format(a, f'0{L}b')])
        i_bin = np.array([int(bin) for bin in format(i, f'0{L}b')])
        
        a_bin = a_bin * 2 - 1
        i_bin = i_bin * 2 - 1
        
        T[i, a] =  np.exp( beta * kappa * K * i_bin[0] * a_bin[-1] - beta * kappa_1 * K * i_bin[0] * a_bin[1])
        print(i, a)
        a = a ^ (1 << (L-1))
        a_bin[0] *= -1
        print(i, a)

        T[i, a] = np.exp( beta * kappa * K * i_bin[0] * a_bin[-1] - beta * kappa_1 * K * i_bin[0] * a_bin[1])
        
    return T

def create_T_single(L, kappa, kappa_1, K, beta):
    T = lil_matrix((2**L, 2**L), dtype=float)
    for i in range (2**L):
        a = i & (~(1 << (L-1)))
        a_bin = np.array([int(bin) for bin in format(a, f'0{L}b')])
        i_bin = np.array([int(bin) for bin in format(i, f'0{L}b')])
        
        a_bin = a_bin * 2 - 1
        i_bin = i_bin * 2 - 1
        
        T[i, a] =  np.exp(beta * K * a_bin[0] * i_bin[0])

        a = a ^ (1 << (L-1))
        a_bin[0] *= -1
        
        T[i, a] = np.exp(beta * K * a_bin[0] * i_bin[0])

    return T

def cyclic_shift(L):
    P = lil_matrix((2**(L), 2**(L)))
    for i in range(2**(L)):
        a = (i << 1) + (i >> (L-1))
        a = a & (~(1 << L))
        P[a, i] = 1
    return P

def reverse_cyclic_shift(L):
    P = lil_matrix((2**L, 2**L))
    for i in range(2**L):
        a = (i >> 1) + ((i & 1) << (L - 1))
        P[a, i] = 1
    return P

def reverse(L):
    P = lil_matrix((2**L, 2**L))
    for i in range(2**L):
        b = '{:0{width}b}'.format(i, width = L) 
        a = int(b[::-1], 2)
        P[a, i] = 1
    return P
    

fmt = 'csr'
L = 4
T = create_T(L, 1, 1, 1, 1)
T_single = create_T_single(L, 1, 1, 1, 1)
P = cyclic_shift(L)
P1 = reverse_cyclic_shift(L) 
R = reverse(L)
print(R.todense())

Trans_2_site = T_single @ P1 @ T
Trans = identity(2**L, format = fmt) 

for i in range(int(L / 2)):
    Ttmp = P @ P @ T @ Trans
    Trans = Ttmp

Trans = P @ P @ R @ Trans

print(Trans.todense())
print(np.allclose(Trans.todense(), Trans.T.todense()))

