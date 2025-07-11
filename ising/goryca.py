import numpy as np

from scipy.sparse.linalg import LinearOperator

from scipy.sparse import identity, diags, csr_matrix, lil_matrix
from scipy.sparse import kron as skron
from scipy.sparse import linalg as sla
from scipy import linalg as la
from numpy.polynomial.chebyshev import Chebyshev
from numpy.polynomial.chebyshev import chebinterpolate


def construct_transfer_matrix(L, K, kappa, kappa_1, beta, h, h1):
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


                if (k % 2 == 1):
                    energy += - h * (snew[k] + sold[k])/2
                    energy += (- K * snew[k] * snew[(k+1) % L] + kappa * K * snew[k] * snew[(k + 2) % L]) / 2
                    energy += (- K * sold[k] * sold[(k+1) % L] + kappa * K * sold[k] * sold[(k + 2) % L]) / 2
                if (k % 2 == 0):
                    energy += - h1 * (snew[k] + sold[k])/2
                    energy += (snew[k] * snew[(k+1) % L] + kappa_1 * K * snew[k] * snew[(k + 2) % L]) / 2
                    energy += (sold[k] * sold[(k+1) % L] + kappa_1 * K * sold[k] * sold[(k + 2) % L]) / 2

            for k in range(L):
                if (k % 2 == 0):
                    energy += K * kappa * snew[k] * sold[k]
                if (k % 2 == 1):
                    energy += K * kappa_1 * snew[k] * sold[k]
                          
            # Left-up diagonal terms
            for k in range(1, L, 2): # sites k = 0, 1, ..., L-1
                energy += - K * sold[k] * snew[(k - 1) % L]
                energy += sold[k] * snew[(k + 1) % L]

            T[i, j] = np.exp(- beta * energy)
    
    return T


def horizontal(L, kappa, kappa_1, K, beta, h, h1):
    T = lil_matrix((2**L, 2**L), dtype=float)
    for a in range(2**L):
        a_bin = np.array([int(bin) for bin in format(a, f'0{L}b')])
        a_bin = a_bin * 2 - 1
        en = 0
        for j in range(L):
            if (j % 2 == 1):
                en += -beta * h * a_bin[j]
                en += -beta * K * a_bin[j] * a_bin[(j + 1) % L] + kappa * beta * K * a_bin[j] * a_bin[(j + 2) % L]
            if (j % 2 == 0):
                en += - beta * h1 * a_bin[j]
                en += beta * a_bin[j] * a_bin[(j + 1) % L] + kappa_1 * beta * K * a_bin[j] * a_bin[(j + 2) % L]
        T[a, a] = np.exp(-en / 2)

    return T




def horizontal_v2(L, kappa, kappa_1, K, beta, h, h1):
    T = lil_matrix((2**L, 2**L), dtype=float)
    for a in range(2**L):
        a_bin = np.array([int(bin) for bin in format(a, f'0{L}b')])
        a_bin = a_bin * 2 - 1
        en = 0
        for j in range(L):
            if (j % 2 == 1):
                en += -K * a_bin[j] * a_bin[(j + 1) % L] + kappa * K * a_bin[j] * a_bin[(j + 2) % L]
            if (j % 2 == 0):
                en += a_bin[j] * a_bin[(j + 1) % L] + kappa_1 * K * a_bin[j] * a_bin[(j + 2) % L]
        T[a, a] = np.exp(- beta * en / 2)

    return T

def horizontal_h(L, beta, h, h1):
    T = lil_matrix((2**L, 2**L), dtype=float)
    for a in range(2**L):
        a_bin = np.array([int(bin) for bin in format(a, f'0{L}b')])
        a_bin = a_bin * 2 - 1
        en = 0
        for j in range(L):
            if (j % 2 == 1):
                en += - h * a_bin[j]

        T[a, a] = np.exp(-beta * en / 2)
    return T

def horizontal_h1(L, beta, h, h1):
    T = lil_matrix((2**L, 2**L), dtype=float)
    for a in range(2**L):
        a_bin = np.array([int(bin) for bin in format(a, f'0{L}b')])
        a_bin = a_bin * 2 - 1
        en = 0
        for j in range(L):
            if (j % 2 == 0):
                en += - h1 * a_bin[j]

        T[a, a] = np.exp(-beta * en / 2)
    return T

def create_T(L, kappa, kappa_1, K, beta):
    T = lil_matrix((2**L, 2**L), dtype=float)
    for i in range (2**L):
        a = i & (~(1 << (L-1)))
        a_bin = np.array([int(bin) for bin in format(a, f'0{L}b')])
        i_bin = np.array([int(bin) for bin in format(i, f'0{L}b')])
        
        a_bin = a_bin * 2 - 1
        i_bin = i_bin * 2 - 1
        
        T[i, a] = np.exp(-kappa * beta * K * i_bin[0] * a_bin[0] - beta * i_bin[0] * a_bin[-1] + beta * K * i_bin[0] * a_bin[1])
        a = a ^ (1 << (L-1))
        a_bin[0] *= -1

        T[i, a] = np.exp(-kappa * beta * K * i_bin[0] * a_bin[0] - beta * i_bin[0] * a_bin[-1] + beta * K * i_bin[0] * a_bin[1])
        
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


def negate(L):
    N = lil_matrix((2**L, 2**L))
    for i in range(2**L):
        for j in range(L):
            if (j % 2 == 1):
                a = i ^ (1 << j)
        N[a, i] = 1
    return N

    
'''
fmt = 'csr'
L = 2
kappa = 0
kappa_1 = 0
K = 1
beta = 15


B = 1
phi = np.pi/4
h = -5
h1 = -6



T = create_T(L, kappa, kappa_1, K, beta)
print("created T")
T_single = create_T_single(L, kappa, kappa_1, K, beta)
print("created T_single")
Dh_2 = horizontal_v2(L, kappa, kappa_1, K, beta, h, h1)
print("created Dh")
H = horizontal_h(L, beta, h, h1)
H1 = horizontal_h1(L, beta, h, h1)
print("created_H")
P = cyclic_shift(L)
print("created P")
P1 = reverse_cyclic_shift(L) 
print("created P1")
R = reverse(L)
print("created R")

Trans_2_site = T_single @ P1 @ T
Trans = identity(2**L, format = fmt) 
Trans = P @ P @ T @ Trans
print("before for loop")

for i in range(int(L / 2 - 1)):
    print(i)
    Ttmp = P @ P @ P @ Trans_2_site @ Trans
    Trans = Ttmp

print("after for loop")

Trans_final = H1 @ H @ Dh_2 @ P @ T_single @ P1 @ Trans @ H1 @ H @ Dh_2

T_reference = construct_transfer_matrix(L, K, kappa, kappa_1, beta, h, h1)


print("created Trans")

print(np.allclose(T_reference, Trans_final.todense()))


H_1 = horizontal_h(L, beta, h, h1)
H1_1 = horizontal_h1(L, beta, h, h1)

h1 =  6

H = horizontal_h(L, beta, h, h1)
H1 = horizontal_h1(L, beta, h, h1)

N = negate(L)

print(H1_1.todense())
print(H1.todense())

print(np.allclose((N @ H1_1 @ N).todense(), H1.todense()))
print(np.allclose((N @ H_1 @ N).todense(), H.todense()))

Trans_1 = H1_1 @ H @ Dh_2 @ P @ T_single @ P1 @ Trans @ H1_1 @ H @ Dh_2
Trans = N @ H1 @ H @ Dh_2 @ P @ T_single @ P1 @ Trans @ H1 @ H @ Dh_2 @ N

print(N.todense())
print(np.allclose(Trans_1.todense(), Trans.todense()))
print(Trans.todense())
print(Trans_1.todense())
'''





def update_matrices(L, kappa, kappa_1, K, beta, h, h1):
    print("update_start")
    T = create_T(L, kappa, kappa_1, K, beta)
    print("created T")
    T_single = create_T_single(L, kappa, kappa_1, K, beta)
    print("created T_single")
    Dh_2 = horizontal_v2(L, kappa, kappa_1, K, beta, h, h1)
    print("created Dh")
    H_h = horizontal_h(L, beta, h, h1)
    print("created H_h")
    H_h1 = horizontal_h1(L, beta, h, h1)
    print("created H_h1")
    P = cyclic_shift(L)
    print("created P")
    P1 = reverse_cyclic_shift(L) 
    print("created P1")
    R = reverse(L)
    print("created R")
        
    Trans_2_site = T_single @ P1 @ T
    print("update_end")
    return T, T_single, Dh_2, P, P1, R, Trans_2_site, H_h, H_h1



def vector_multiplication(v, L, T, T_single, Dh_2, P, P1, R, Trans_2_site, H_h, H_h1):
    v = H_h1 @ H_h @ Dh_2 @ v

    v = P @ P @ T @ v

    for i in range(int(L/2-1)):
        v = P @ P @ P @ Trans_2_site @ v
    
    v = H_h1 @ H_h @ Dh_2 @ P @ T_single @ P1 @ v

    return v


def update_H(H_h, H_h1, new_h, new_h1, old_h, old_h1, L):
    '''
    H_h = np.power(H_h, new_h/old_h)
    H_h1 = np.power(H_h1, new_h1/old_h1)
    '''
    for i in range(2**L):
        H_h[i, i] = np.power(H_h[i, i], new_h/old_h)
        H_h1[i, i] = np.power(H_h1[i, i], new_h1/old_h1)
    return H_h, H_h1


