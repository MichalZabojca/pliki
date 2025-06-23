import numpy as np

from scipy.sparse.linalg import LinearOperator

from scipy.sparse import identity, diags, csr_matrix, lil_matrix
from scipy.sparse import kron as skron
from scipy.sparse import linalg as sla
from scipy import linalg as la
from numpy.polynomial.chebyshev import Chebyshev
from numpy.polynomial.chebyshev import chebinterpolate


def construct_transfer_matrix(L, K, kappa, kappa_1, beta, h):
    # Generate all possible spin states for L spins
    states = np.array([list(format(i, f'0{L}b')) for i in range(2 ** L)], dtype=int)
    
    # Initialize the transfer matrix with zeros
    size = 2 ** L
    T = np.zeros((size, size), dtype=float)
    
    for i in range(size):
        for j in range(size):
            # Convert spin states to -1 and 1 representation    
            snew = 2 * states[i] - 1
            sold = 2 * states[j] - 1 
            
            energy = 0.0
            # Nearest neighbor (horizontal) terms between old and new spins
            for k in range(L): 

                energy += - h * (snew[k] + sold[k])/2

                if (k % 2 == 1):
                    energy += (K * snew[k] * snew[(k+1) % L] + kappa * K * snew[k] * snew[(k + 2) % L]) / 2
                    energy += (K * sold[k] * sold[(k+1) % L] + kappa * K * sold[k] * sold[(k + 2) % L]) / 2
                if (k % 2 == 0):
                    energy += (snew[k] * snew[(k+1) % L] + kappa_1 * K * snew[k] * snew[(k + 2) % L]) / 2
                    energy += (sold[k] * sold[(k+1) % L] + kappa_1 * K * sold[k] * sold[(k + 2) % L]) / 2

            for k in range(L):
                if (k % 2 == 0):
                    energy += K * kappa * snew[k] * sold[k]
                if (k % 2 == 1):
                    energy += K * kappa_1 * snew[k] * sold[k]
                          
            # Left-up diagonal terms
            for k in range(1, L, 2): # sites k = 0, 1, ..., L-1
                energy += K * sold[k] * snew[(k - 1) % L]
                energy += sold[k] * snew[(k + 1) % L]

            T[i, j] = np.exp(- beta * energy)
    
    return T


def horizontal(L, kappa, kappa_1, K, beta, h):
    T = lil_matrix((2**L, 2**L), dtype=float)
    for a in range(2**L):
        a_bin = np.array([int(bin) for bin in format(a, f'0{L}b')])
        a_bin = a_bin * 2 - 1
        en = 0
        for j in range(L):
            
            en += - h * a_bin[j]

            if (j % 2 == 1):
                en += beta * K * a_bin[j] * a_bin[(j + 1) % L] + kappa * beta * K * a_bin[j] * a_bin[(j + 2) % L]
            if (j % 2 == 0):
                en += beta * a_bin[j] * a_bin[(j + 1) % L] + kappa_1 * beta * K * a_bin[j] * a_bin[(j + 2) % L]
        T[a, a] = np.exp(-en / 2)

    return T



def create_T(L, kappa, kappa_1, K, beta):
    T = lil_matrix((2**L, 2**L), dtype=float)
    for i in range (2**L):
        a = i & (~(1 << (L-1)))
        a_bin = np.array([int(bin) for bin in format(a, f'0{L}b')])
        i_bin = np.array([int(bin) for bin in format(i, f'0{L}b')])
        
        a_bin = a_bin * 2 - 1
        i_bin = i_bin * 2 - 1
        
        T[i, a] = np.exp(- kappa * beta * K * i_bin[0] * a_bin[0] - beta * i_bin[0] * a_bin[-1] - beta * K * i_bin[0] * a_bin[1])
        a = a ^ (1 << (L-1))
        a_bin[0] *= -1

        T[i, a] = np.exp(- kappa * beta * K * i_bin[0] * a_bin[0] - beta * i_bin[0] * a_bin[-1] - beta * K * i_bin[0] * a_bin[1])
        
    return T

def create_T_single(L, kappa, kappa_1, K, beta):
    T = lil_matrix((2**L, 2**L), dtype=float)
    for i in range (2**L):
        a = i & (~(1 << (L-1)))
        a_bin = np.array([int(bin) for bin in format(a, f'0{L}b')])
        i_bin = np.array([int(bin) for bin in format(i, f'0{L}b')])
        
        a_bin = a_bin * 2 - 1
        i_bin = i_bin * 2 - 1
        
        T[i, a] =  np.exp(- kappa_1 * beta * K * a_bin[0] * i_bin[0])

        a = a ^ (1 << (L-1))
        a_bin[0] *= -1
        
        T[i, a] = np.exp(- kappa_1 * beta * K * a_bin[0] * i_bin[0])

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
kappa = 1
kappa_1 = 1
K = 1
beta = 1
h = 1 
T = create_T(L, kappa, kappa_1, K, beta)
T_single = create_T_single(L, kappa, kappa_1, K, beta)
Dh_2 = horizontal(L, kappa, kappa_1, K, beta, h)
P = cyclic_shift(L)
P1 = reverse_cyclic_shift(L) 
R = reverse(L)

Trans_2_site = T_single @ P1 @ T
Trans = identity(2**L, format = fmt) 
Trans = P @ P @ T @ Trans

for i in range(int(L / 2 - 1)):
    Ttmp = P @ P @ P @ Trans_2_site @ Trans
    Trans = Ttmp

Trans = Dh_2 @ P @ T_single @ P1 @ Trans @ Dh_2

T_reference = construct_transfer_matrix(L, K, kappa, kappa_1, beta, h)

print(np.allclose(Trans.todense(), T_reference))




def update_matrices(L, kappa, kappa_1, K, beta, h):
    T = create_T(L, kappa, kappa_1, K, beta)
    T_single = create_T_single(L, kappa, kappa_1, K, beta)
    Dh_2 = horizontal(L, kappa, kappa_1, K, beta, h)
    P = cyclic_shift(L)
    P1 = reverse_cyclic_shift(L) 
    R = reverse(L)
        
    Trans_2_site = T_single @ P1 @ T

    return T, T_single, Dh_2, P, P1, R, Trans_2_site



def vector_multiplication(v, L, T, T_single, Dh_2, P, P1, R, Trans_2_site):
    v = Dh_2 @ v

    for i in range(int(L/2-1)):
        v = P @ P @ P @ Trans_2_site @ v
    
    v = Dh_2 @ P @ T_single @ P1 @ v

    return v
    
lallalalalalala
