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

from goryca import vector_multiplication, update_matrices

L = 4
J = 1
K = 1
kappa = 0.6
kappa_1 = 0
h = 10
h1 = 1

def calculate_f(temp, L, k):
    beta = 1/(k*temp)
    matrices = update_matrices(L, kappa, kappa_1, K, beta, h, h1)
    Transfer_matrix = sla.LinearOperator((2**L, 2**L), matvec = lambda v: vector_multiplication(v, L, *matrices))
    return k*temp*np.log(sla.eigs(Transfer_matrix, k=1, return_eigenvectors = False)[0])

def calculate_f_array(temp, L, k, h):
    f = np.zeros(temp.size)
    for i in range(temp.size):
        beta = 1/(k*temp[i])
        matrices = update_matrices(L, kappa, kappa_1, K, beta, h, h1)
        Transfer_matrix = sla.LinearOperator((2**L, 2**L), matvec = lambda v: vector_multiplication(v, L, *matrices))
        f[i] = k*temp[i]*np.log(sla.eigs(Transfer_matrix, k=1, return_eigenvectors = False)[0])
    return f

k=1
n = 30
start = 0.3
end = 5
temps = np.linspace(start, end, n)

f = np.zeros(temps.size)

for i in range(temps.size):
    print(i)
    f[i] = calculate_f(temps[i], L, k)

fig = plt.figure(dpi=200)
wyk = fig.add_subplot(121)
wyk_cv = fig.add_subplot(122)
wyk.plot(temps, f)

h_array = [10]
for i in h_array:
    cheb = Chebyshev.interpolate(calculate_f_array, 30, domain = [start, end],args = (L, k, i))
    cheb_der = cheb.deriv(2)
    cv_cheb = cheb_der.linspace(n, [start, end])
    wyk_cv.plot(cv_cheb[0], cv_cheb[1])

plt.show()
