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

def horizontal_TM(L, J, beta, h0, h1, pbc = True, fmt = 'csr'):









