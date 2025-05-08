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


def Dhalf_horizontal_TM(L, J, beta, pbc = True, fmt = 'csr'):
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
    Dhalf = diags( np.exp( beta * J/2 * Diagonal ), format = fmt  )

    return Dhalf


def create_T(L, kappa, K, beta):
    #T_z,1 - contribution to the interlayer interaction of first spin
    T = lil_matrix((2**(L+2), 2**(L+2)), dtype=float)
    for i in range (2**(L+2)):
        a = (i & (~(1 << (L+1)))) + ((2 & i) << L) 

        a_bin = np.array([int(bin) for bin in format(a, f'0{L+2}b')])
        i_bin = np.array([int(bin) for bin in format(i, f'0{L+2}b')])
        
        a_bin = a_bin * 2 - 1
        i_bin = i_bin * 2 - 1
        
        T[i, a] = np.exp(beta * K * a_bin[0] * i_bin[0] - beta * kappa * K * i_bin[0] * (a_bin[1] + a_bin[-2]))

        a = a ^ 2
        a_bin[-2] *= -1
        
        T[i, a] = np.exp(beta * K * a_bin[0] * i_bin[0] - beta * kappa * K * i_bin[0] * (a_bin[1] + a_bin[-2]))
        
    return T


def create_Tzr1(L, kappa, K, beta):
    #T_z,1|r1 - same as T_z,1, but takes r1' instead of s2' for computing the right-down interaction
    T = lil_matrix((2**(L+2), 2**(L+2)), dtype=float)
    for i in range (2**(L+2)):
        a = (i & (~(1 << (L+1)))) + ((2 & i) << L) 

        a_bin = np.array([int(bin) for bin in format(a, f'0{L+2}b')])
        i_bin = np.array([int(bin) for bin in format(i, f'0{L+2}b')])
        
        a_bin = a_bin * 2 - 1
        i_bin = i_bin * 2 - 1
        
        T[i, a] = np.exp(beta * K * a_bin[0] * i_bin[0] - beta * kappa * K * i_bin[0] * (a_bin[-1] + a_bin[-2]))

        a = a ^ 2
        a_bin[-2] *= -1
        
        T[i, a] = np.exp(beta * K * a_bin[0] * i_bin[0] - beta * kappa * K * i_bin[0] * (a_bin[-1] + a_bin[-2]))
        
    return T


def create_shift_matrix(L):
    #P - cyclic spin shift, doesn't change auxiliary spins
    P = lil_matrix((2**(L+2), 2**(L+2)))
    for i in range(2**(L+2)):
        a = (i - i % 4) << 1
        a = ((a >> (L)) & 4) + (a & (~(1 << (L+2)))) + i % 4
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


def create_transfer_matrix(L, beta):

    K=np.log(2)
    kappa = 1
    J = 1
    fmt = 'csr'


    #create auxiliary matrices S and S^(-1)
    S = create_S(L)
    S1 = create_S1(L)

    #create bit shift matrix P
    P = create_shift_matrix(L)

    #create T_z,1
    T = create_T(L, kappa, K, beta)

    #create horizontal interaction matrix
    Dh = Dhalf_horizontal_TM(L, J, beta, pbc=True)

    #create T_z,1|r1 (for pbc)
    Tzr1 = create_Tzr1(L, kappa, K, beta)


    #calculate Trans = (P T_z,1)^(L-1)
    Trans = identity( 2**(L+2), format = fmt )
    for ind in range(L-1):
        Ttmp = P @ T @ Trans
        Trans = Ttmp

    #calculate transfer matrix
    Trans_final = Dh @ S1 @ P @ Tzr1 @ Trans @ S @ Dh

    return Trans_final



#Trans_final = create_transfer_matrix(2, 1)
#Ttransfer = construct_transfer_matrix_dnni(L, J, K, kappa, rightup = True, pbc = True)

#print(Trans_final.todense())
#print(Ttransfer)
#print("Is Transfer Matrix product for DNNI all right? ",np.allclose(Trans_final.todense(),Ttransfer))

#print( np.sum(np.power((la.eigvalsh(Trans_final.todense())), N) ))


def vector_multiplication(v, beta, L):
    K = 1
    kappa = 0
    J = 1
    S = create_S(L)
    S1 = create_S1(L)
    P = create_shift_matrix(L)
    T = create_T(L, kappa, K, beta)
    Dh = Dhalf_horizontal_TM(L, J, beta, pbc=True)
    Tzr1 = create_Tzr1(L, kappa, K, beta)

    v = Dh @ v
    v = S @ v
    for i in range(L-1):
        v = T @ v
        v = P @ v
    v = Tzr1 @ v
    v = P @ v
    v = S1 @ v
    v = Dh @ v
    return v



L=2
N=100
k=1


n = 20
start = 0.1
end = 8
temps = np.cos((2 * np.arange(n) + 1) / (2 * n))
temps *= (end - start) / (np.max(temps) - np.min(temps))
temps = temps - np.min(temps) + start
print(np.max(temps))

betas = 1/(k*temps)
f = np.zeros(betas.size)

for i in range(betas.size):
    print(i)
    Transfer_matrix = sla.LinearOperator((2**L, 2**L), matvec = lambda v: vector_multiplication(v, betas[i], L))
    f[i] = -k*temps[i]*np.log(sla.eigsh(Transfer_matrix, k=1, return_eigenvectors = False))

fig = plt.figure(dpi=300)
wyk = fig.add_subplot(121)
wyk_cv = fig.add_subplot(122)
wyk.plot(temps, f)

deg=10
cheb_series = Chebyshev.interpolate(, deg, domain = [start, end])
cheb_deriv = cheb_series.deriv(m = 2)
cv_cheb = cheb_deriv(temps)
wyk_cv.plot(temps, cv_cheb)


fig.savefig("f_plot.png")
plt.show()




             