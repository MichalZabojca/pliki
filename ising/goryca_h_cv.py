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


L = 8
J = 1
K = 1
kappa = 0.5
kappa_1 = 0.2
h = 1


fig_mag = plt.figure(dpi = 200)
wyk = fig_mag.add_subplot(111)

temp = 2.1
k = 1
n = 100
start = 0
end = 5
mag = np.linspace(start, end, n)
f = np.zeros(mag.size)
beta = 1/(k * temp)


def calculate_f(mag, L, k, temp):
    beta = 1/(k * temp)
    matrices = update_matrices(L, kappa, kappa_1, K, beta, mag)
    Transfer_matrix = sla.LinearOperator((2**L, 2**L), matvec = lambda v: vector_multiplication(v, L, *matrices))
    return k*temp*np.log(sla.eigs(Transfer_matrix, k=1, return_eigenvectors = False)[0])

for i in range(mag.size):
    print(i)
    f[i] = calculate_f(mag[i], L, k, temp)


wyk.plot(mag, f)

kappa = 1
kappa_1 = 1
temp = 10

for i in range(mag.size):
    print(i)
    f[i] = calculate_f(mag[i], L, k, temp)

wyk.plot(mag, f)

plt.show()
