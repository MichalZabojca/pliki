import numpy as np
import matplotlib.pylab as plt

#from scipy.sparse.linalg import LinearOperator

from scipy.sparse import identity, diags, csr_matrix, lil_matrix
from scipy.sparse import kron as skron
from scipy.sparse import linalg as sla
from scipy import linalg as la


from scipy import special

# Pauli matrices
s0 =  np.array([[1, 0], [0, 1]], float) 
sx =  np.array([[0, 1], [1, 0]], float)
sz =  np.array([[0, 1], [0, -1]], float)


def construct_transfer_matrix_dnni(L, J, K, kappa, rightup=False, pbc=False):
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
            for k in range(L-1): # sites k = 0, 1, ..., L-1
                energy += - kappa * K * snew[k] * sold[k+1]
            
            # Right-up diagonal terms if enabled
            if ( rightup ):
                for k in range(L-1): # right up diagonal
                    energy += - kappa * K * snew[k+1] * sold[k]
            
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


def bit_shift(n, d, L_BITS=8):
    """
    Bitwise shift of bit representation of n by d bits to the LEFT, 
    cyclic on L lowest bits (rewritten with DeepSeek)

    Returns:
      integer representation of shifted n
    """
    d = d % L_BITS
    mask = (1 << L_BITS) - 1
    return ((n << d) | (n >> (L_BITS - d))) & mask

def create_shift_matrix(L, d=1):
    """
    Creates a sparse shift matrix.
  
    Args:
      d: The amount of shift.
      L: The number of bits.

    Returns:
      A sparse shift matrix in CSR format.
    """
    N = 2**L
    sparse_matrix = lil_matrix((N, N))

    for col in range(N):
        shifted_col = bit_shift(col, d, L)
        sparse_matrix[shifted_col, col] = 1

    return sparse_matrix.tocsr()


def two_site_TM(K, kappa):
    """
    Constructs the transfer matrix for two sites, single diagonal
    
    Parameters:
    J (float) : Interaction term for horizontal bonds (finite direction)
    K (float) : Interaction term for vertical bonds 
    pbc (bool): use of periodic boundary condition
    
    Returns:
    T (scipy.sparse.csr_matrix): The transfer matrix.
    """
    
    # (all) basis states for a chain of L=2 spins
    states = np.array([[0, 0],[0, 1],[1, 0],[1, 1]])
    
    L=2
    # Initialize the transfer matrix, sparse format
    T = lil_matrix((2**L, 2**L), dtype=float)

    for i in range(2**L):
        for j in range(2**L):
            snew = 2 * states[i] - 1
            sold = 2 * states[j] - 1
            if ( snew[1] == sold[1] ) :
                # Interaction energy between spins in states[i] and states[j]
                interaction_energy = 0
                # vertical T matrix
                interaction_energy += K * snew[0] * sold[0]
                # left up diagonal
                interaction_energy += - kappa * K * snew[0] * sold[1]           
                # Transfer matrix element
                T[i, j] = np.exp(interaction_energy)
    
    return T.tocsr() ## warning: this matrix is not sparse in practice


def Dhalf_horizontal_TM(L, J, pbc = True, fmt = 'csr'):
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
    
        # Apply periodic boundary conditions if enabled
        if ( pbc ): # not tested yet ???
            # Horizontal terms for periodic boundaries
            energy += spins[L-1] * spins[0]

        Diagonal[i] = energy
         
    # exponentiate the diagonal, bring back to sparse, take sqrt()
    Dhalf = diags( np.exp( J/2 * Diagonal ), format = fmt  )

    return Dhalf

from scipy.sparse import csr_array

fmt = 'csr'

# check element-wise agreement for the parameters:
L = 3
kappa = 0
K= np.log(2)
J = 1

# calculate benchmark
Ttransfer = construct_transfer_matrix_dnni(L, J, K, kappa, rightup = False, pbc = False)


################################################################
# general L, rightup = False, DNNI model, open bc              #
#                                                              #
################################################################

# one-site T-matrix extended to a chain
T_site = skron( np.exp(K) * s0 + np.exp(-K) * sx, identity( 2**(L-1),  format = fmt ), format = fmt)
# two-site T-matrix extended to a chain
T_2site = skron( two_site_TM(K,kappa), identity( 2**(L-2),  format = fmt ), format = fmt)


# shift matrix
P  = create_shift_matrix(L, d=1)

Trans = identity(  2**(L), format = fmt )
for ind in range(L-1): # calculate ( P T_2site)^(L-1)
    Ttmp = P @ T_2site @ Trans
    Trans = Ttmp

# diagonal matrix
Dh = Dhalf_horizontal_TM(L, J, pbc=False)

# final calculation
Trans_final = (Dh @ P @ T_site @ Trans @ Dh)

def create_T(K, J, L):
    T = lil_matrix((2**L, 2**L))
    for i in range(2**L):
        T[i, i] = np.exp(K*J)
        a = i ^ (1 << (L - 1))
        T[a, i] = np.exp(-K*J)
    return T

def create_P(L):
    P = lil_matrix((2**L, 2**L))
    for i in range(2**L):
        a = ((i << 1) & (~(1 << L))) + (i >> (L - 1))
        P[a, i] = 1
    return P

P_bitwise = create_P(L)
T_bitwise = create_T(K, J, L)
Trans_bitwise = identity(  2**(L), format = fmt )
for ind in range(L): # calculate ( P T_2site)^(L-1)
    Ttmp = P_bitwise @ T_bitwise @ Trans_bitwise
    Trans_bitwise = Ttmp

Trans_bitwise_final = Dh @ Trans_bitwise @ Dh

print(Trans_bitwise.todense())
print(P_bitwise.todense())
print(Trans_bitwise_final.todense())
print(Ttransfer)

'''
def mult(x):
    Trans_final
    return Trans_final @ x


A = sla.LinearOperator((2**L, 2**L), matvec = mult)
print(sla.eigsh(A, k=1, return_eigenvectors = False))
print(sla.eigsh(Trans_final, k=1, return_eigenvectors = False))'''

print("Is Transfer Matrix product for TNNI all right? ",np.allclose(Trans_final.todense(),Ttransfer))

