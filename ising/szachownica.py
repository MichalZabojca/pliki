import numpy as np
from matplotlib import pyplot as plt

from scipy.sparse.linalg import LinearOperator

from scipy.sparse import identity, diags, csr_matrix, lil_matrix
from scipy.sparse import kron as skron
from scipy.sparse import linalg as sla
from scipy import linalg as la
from numpy.polynomial.chebyshev import Chebyshev
from numpy.polynomial.chebyshev import chebinterpolate


from scipy import special

def construct_transfer_matrix_dnni(L, J, K, kappa, kappa_1, rightup, pbc):
    """
    Constructs the transfer matrix for the DNNI model with nearest neighbor couplings J and K,
    and diagonal terms with prefactor (- kappa), on a strip of width L.
    
    The default parameteres: rightup=False, pbc=False give Triangular NN Ising (TNNI) model; open bc

    Parameters:
    L (int): Width of the strip (number of spins).
    J (float): Interaction strength for horizontal bonds.
    K (float): Interaction strength for vertical bonds.
    kappa (float): Prefactor for the diagonal interaction terms.
    rightup (bool): If True, includes the right-up diagonal terms.
    pbc (bool): If True, applies periodic boundary conditions.
    Returns:
    T (numpy array): The transfer matrix of size 2^L x 2^L.
    """
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
            for k in range(L-1): 
                energy += J/2 * snew[k] * snew[k+1]
                energy += J/2 * sold[k] * sold[k+1]

            # Vertical bond terms between old and new spins
            for k in range(L):
                energy += K * snew[k] * sold[k]
                          
            # Left-up diagonal terms
            for k in range(0, L-2, 2): # sites k = 0, 1, ..., L-1
                energy += - kappa_1 * K * snew[k] * sold[k+1]
                energy += - kappa * K * snew[k+1] * sold[k+2]
            energy += -kappa_1 * K *snew[L-2] * sold[L-1]

            # Right-up diagonal terms if enabled
            if ( rightup ):
                for k in range(0, L-2, 2): # right up diagonal
                    energy += - kappa_1 * K * snew[k+1] * sold[k]
                    energy += - kappa * K * snew[k+2] * sold[k+1]
                energy += -kappa_1 * K * snew[L-1] * sold[L-2]

            # Apply periodic boundary conditions if enabled
            if ( pbc ): # not tested yet ???
                # Horizontal terms for periodic boundaries
                energy += J/2 * snew[L-1] * snew[0]
                energy += J/2 * sold[L-1] * sold[0]
                
                # Diagonal terms for periodic boundaries
                energy += - kappa * K * snew[L-1] * sold[0]
                if ( rightup ):
                    energy += - kappa * K * sold[L-1] * snew[0]
                        
            # Set the transfer matrix element            
            T[i, j] = np.exp(energy)
    
    return T


def Dhalf_horizontal_TM(L, J, beta, h0, h1, pbc = True, fmt = 'csr'):
## site-site horizontal interaction matrix (finite space-direction L)
## only diagonal, simplified construction

    # Generate all possible spin states for L spins
    states = np.array([list(format(i, f'0{L}b')) for i in range(2 ** L)], dtype=int)
    
    # Initialize the diagonal of the transfer matrix with zeros
    size = 2 ** L
    Diagonal = np.zeros(size, dtype=float)

    
    for i in range(size):
        # Convert states to spin +/-1 representation    
        spins = 2 * states[i] - 1

        energy = 0.0

        # Nearest neighbor (horizontal) terms between all spins
        for k in range(L-1):
            energy += spins[k] * spins[k+1]
            if k % 2 == 0:
                energy += h0 * spins[k]
            else:
                energy += h1 * spins[k]
    
        # Apply periodic boundary conditions if enabled
        if ( pbc ): # not tested yet ???
            # Horizontal terms for periodic boundaries
            energy += spins[L-1] * spins[0]

        Diagonal[i] = energy
         
    # exponentiate the diagonal, bring back to sparse, take sqrt()
    Dhalf = diags( np.exp( beta * J/2 * Diagonal ), format = fmt  )

    return Dhalf


def create_T(L, kappa, kappa_1, K, beta):
    #T_z,1 - contribution to the interlayer interaction of first spin
    T = lil_matrix((2**(L+2), 2**(L+2)), dtype=float)
    for i in range (2**(L+2)):
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


def create_Tzr1(L, kappa, kappa_1,  K, beta):
    #T_z,1|r1 - same as T_z,1, but takes r1' instead of s2' for computing the right-down interaction
    T = lil_matrix((2**(L+2), 2**(L+2)), dtype=float)
    for i in range (2**(L+2)):
        a = (i & (~(1 << (L+1)))) + ((2 & i) << L) 

        a_bin = np.array([int(bin) for bin in format(a, f'0{L+2}b')])
        i_bin = np.array([int(bin) for bin in format(i, f'0{L+2}b')])
        
        a_bin = a_bin * 2 - 1
        i_bin = i_bin * 2 - 1
        
        T[i, a] = np.exp(beta * K * a_bin[0] * i_bin[0] - beta * kappa * K * i_bin[0] * a_bin[-2] - beta * kappa_1 * K * i_bin[0] * a_bin[-1])

        a = a ^ 2
        a_bin[-2] *= -1
        
        T[i, a] = np.exp(beta * K * a_bin[0] * i_bin[0] - beta * kappa * K * i_bin[0] * a_bin[-2] - beta * kappa_1 * K * i_bin[0] * a_bin[-1])
        
    return T


def create_shift_matrix(L):
    #P - cyclic spin shift, doesn't change auxiliary spins
    P = lil_matrix((2**(L+2), 2**(L+2)))
    for i in range(2**(L+2)):
        a = (i - i % 4) << 1
        a = ((a >> (L)) & 4) + (a & (~(1 << (L+2)))) + i % 4
        P[a, i] = 1
    return P


def cyclic_1bit_shift(L):
    P = lil_matrix((2**(L), 2**(L)))
    for i in range(2**(L)):
        a = (i << 1) + (i >> (L-1))
        a = a & (~(1 << L))
        P[a, i] = 1
    return P

def cyclic_1bit_shift_reverse(L):
    P = lil_matrix((2**(L), 2**(L)))
    for i in range(2**(L)):
        a = (i >> 1) + ((i & 1) << (L-1))
        P[a, i] = 1
    return P




def create_S(L):
    #S - maps vector of size 2^L to vector of size 2^(L+2)
    S = lil_matrix((2**(L+2), 2**L))
    for i in range (2**L):
        a = (i << 2) + (i >> (L - 1)) + ((i << 1) & 2)
        S[a, i] = 1
    
    return S


def create_S1(L):
    #S^(-1) - recovers vector of size 2^L from vector of size 2^(L+2)
    S1 = lil_matrix((2**L, 2**(L+2)))
    for i in range(2**(L+2)):
        a = i >> 2
        S1[a, i] = 1
        
    return S1

fmt = "csr"
L = 4
J = 1
K = 1
kappa = 0.8
beta = 1
kappa_1 = 0.6

B = 0
phi = np.pi/8
h0 = np.cos(phi) * B
h1 = np.sin(phi) * B

Dh_1 = Dhalf_horizontal_TM(L, J, beta, h0, h1, pbc=True, fmt='csr')
Dh_2 = Dhalf_horizontal_TM(L, J, beta, h1, h0, pbc=True, fmt='csr')
T_left = create_T(L, kappa, kappa_1, K, beta)
T_right = create_T(L, kappa_1, kappa, K, beta)
Tzr1_left = create_Tzr1(L, kappa, kappa_1, K, beta)
Tzr1_right = create_Tzr1(L, kappa_1, kappa, K, beta)
P = create_shift_matrix(L)
S1 = create_S1(L)
S = create_S(L)

def update_matrices(L, J, K, kappa, kappa_1, h0, h1, beta):
    Dh_1 = Dhalf_horizontal_TM(L, J, beta, h0, h1, pbc=True, fmt='csr')
    Dh_2 = Dhalf_horizontal_TM(L, J, beta, h1, h0, pbc=True, fmt='csr')
    T_left = create_T(L, kappa, kappa_1, K, beta)
    T_right = create_T(L, kappa_1, kappa, K, beta)
    Tzr1_left = create_Tzr1(L, kappa, kappa_1, K, beta)
    Tzr1_right = create_Tzr1(L, kappa_1, kappa, K, beta)
    P = create_shift_matrix(L)
    S1 = create_S1(L)
    S = create_S(L)
    C = cyclic_1bit_shift(L)
    C1 = cyclic_1bit_shift_reverse(L)
    return Dh_1, Dh_2, T_left, T_right, Tzr1_left, Tzr1_right, P, S1, S, C, C1

def vector_multiplication(v, L, Dh_1, Dh_2, T_left, T_right, Tzr1_left, Tzr1_right, P, S1, S, C, C1):
    v = S @ Dh_1 @ v
    for i in range(int((L-2)/2)):
        v = P @ T_right @ P @ T_left @ v
    v = Dh_1 @ S1 @ P @ Tzr1_right @ P @ T_left @ v

    v = S @ Dh_2 @ v
    for i in range(int((L-2)/2)):
        v = P @ T_left @ P @ T_right @ v
    v = Dh_2 @ S1 @ P @ Tzr1_left @ P @ T_right @ v

    return v


def reverse(L):
    P = lil_matrix((2**L, 2**L))
    for i in range(2**L):
        b = '{:0{width}b}'.format(i, width = L) 
        a = int(b[::-1], 2)
        P[a, i] = 1
    return P

R = reverse(L)

Trans = identity(2**(L+2), format = fmt) 
for i in range(int((L-2)/2)):
    Ttmp = P @ T_right @ P @ T_left @ Trans
    Trans = Ttmp

T_layer_1 = Dh_2 @ S1 @ P @ Tzr1_right @ P @ T_left @ Trans @ S @ Dh_1


Trans = identity(2**(L+2), format = fmt)
for i in range(int((L-2)/2)):
    Ttmp = P @ T_left @ P @ T_right @ Trans
    Trans = Ttmp

T_layer_2 = Dh_1 @ S1 @ P @ Tzr1_left @ P @ T_right @ Trans @ S @ Dh_2

C = cyclic_1bit_shift(L)
C1 = cyclic_1bit_shift_reverse(L)
Trans_final = C1 @ R @ T_layer_1 @ R @ C @ T_layer_1 @ C1 @ R


print(np.allclose((Trans_final).todense(), (Trans_final).todense().T))



'''

/// - kappa
... - kappa_1
 o------o------o------o
 |\    /|.    .|\    /|
 | \  / | .  . | \  / |
 |  \/  |  ..  |  \/  |
 |  /\  |  ..  |  /\  |
 | /  \ | .  . | /  \ |
 |/    \|.    .|/    \|
 o------o------o------o  
 
T_left:
        o
       /|.
      / | .
     /  |  .
    /   |   . 
   /    |    .
  /     |     .
 o      o      o

T_right:
        o
       .|\
      . | \
     .  |  \
    .   |   \
   .    |    \
  .     |     \
 o      o      o

'''

'''
Trans_final = T_layer_1 @ T_layer_2

Ttransfer = construct_transfer_matrix_dnni(L, J, K, kappa, kappa_1, rightup = True, pbc = True)

print("Is Transfer Matrix product for DNNI all right? ",np.allclose(T_layer_1.todense(),Ttransfer))'
'''


