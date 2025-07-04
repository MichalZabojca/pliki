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

from goryca import vector_multiplication, update_matrices, update_H, horizontal_h, horizontal_h1


L = 4
J = 1
K = 1
kappa = 0
kappa_1 = 0


fig_mag = plt.figure(dpi = 200)
wyk = fig_mag.add_subplot(111)

temp = 20
k = 1
n = 10
start = -2
end = 2
hy = np.linspace(start, end, n)
hx = np.linspace(start, end, n)

f = np.zeros(shape = (n, n))
beta = 1/(k * temp)


def calculate_f(hx, hy, L, k, temp):
    beta = 1/(k * temp)
    matrices = update_matrices(L, kappa, kappa_1, K, beta, hx, hy)
    Transfer_matrix = sla.LinearOperator((2**L, 2**L), matvec = lambda v: vector_multiplication(v, L, *matrices))
    return k*temp*np.log(sla.eigs(Transfer_matrix, k=1, return_eigenvectors = False)[0])

def calculate_f_update_h(matrices, h_matrices, L, temp):
    Transfer_matrix = sla.LinearOperator((2**L, 2**L), matvec = lambda v: vector_multiplication(v, L, *matrices, *h_matrices))
    return k*temp*np.log(sla.eigs(Transfer_matrix, k=1, return_eigenvectors = False)[0]) 

'''
for i in range(n):
    for j in range(n):
        print(i, j)
        f[i, j] = calculate_f(hy[i], hx[j], L, k, temp)
'''

def plot_cv(temp):
    n = 30
    cv = np.zeros(n)
    f = np.zeros(n)
    u = np.zeros(n)
    start = -20
    stop = 20
    h = np.linspace(start, stop, n)
    dt = 0.05

    hx= h[0]
    hy= h[0]

    beta = 1/(k*temp)
    matrices = update_matrices(L, kappa, kappa_1, K, beta, hx, hy)[:-2]
    h_matrices = horizontal_h(L, beta, hx, hy), horizontal_h1(L, beta, hx, hy)

    beta = 1/(k*(temp+dt))
    matrices_dt = update_matrices(L, kappa, kappa_1, K, beta, hx, hy)[:-2]
    h_matrices_dt = horizontal_h(L, beta, hx, hy), horizontal_h1(L, beta, hx, hy)

    beta = 1/(k*(temp+2*dt))
    matrices_2dt = update_matrices(L, kappa, kappa_1, K, beta, hx, hy)[:-2]
    h_matrices_2dt = horizontal_h(L, beta, hx, hy), horizontal_h1(L, beta, hx, hy)


    hx_old = h[0]
    hy_old = h[0]

    for j in range(1):
        for i in range(n):
            print(j)
            hx = h[i]
            hy = h[i]
            
            h_matrices = update_H(*h_matrices, hx, hy, hx_old, hy_old, L)
            h_matrices_dt = update_H(*h_matrices_dt, hx, hy, hx_old, hy_old, L)
            h_matrices_2dt = update_H(*h_matrices_2dt, hx, hy, hx_old, hy_old, L)
            f[i] = calculate_f_update_h(matrices, h_matrices, L, temp)
            print("h_matrices_updated")

            der_low = (calculate_f_update_h(matrices_dt, h_matrices_dt, L, temp + dt) - f[i]) / dt
            der_high = (calculate_f_update_h(matrices_2dt, h_matrices_2dt, L, temp + 2 * dt) - calculate_f_update_h(matrices_dt, h_matrices_dt,  L, temp + dt))/dt 
            u[i] = der_low
            cv[i] = (der_high - der_low)/dt
            print("cv_calculated")

            hy_old = hy
            hx_old = hx
        
    #wyk.imshow(cv)
        
    wyk.plot(h * 2, cv, label = f"{temp}")


'''
temps = [0.8, 1.3, 1.8, 2.1, 5]

for i in temps:
    plot_cv(i)
'''


def imshow_h(temp):
    n = 20
    cv = np.zeros(shape = (n, n))
    f = np.zeros(shape = (n, n))
    u = np.zeros(shape = (n, n))
    start = -50
    stop = 50
    h = np.linspace(start, stop, n)
    dt = 0.01

    hx= h[0]
    hy= h[0]

    beta = 1/(k*temp)
    matrices = update_matrices(L, kappa, kappa_1, K, beta, hx, hy)[:-2]
    h_matrices = horizontal_h(L, beta, hx, hy), horizontal_h1(L, beta, hx, hy)

    beta = 1/(k*(temp+dt))
    matrices_dt = update_matrices(L, kappa, kappa_1, K, beta, hx, hy)[:-2]
    h_matrices_dt = horizontal_h(L, beta, hx, hy), horizontal_h1(L, beta, hx, hy)

    beta = 1/(k*(temp+2*dt))
    matrices_2dt = update_matrices(L, kappa, kappa_1, K, beta, hx, hy)[:-2]
    h_matrices_2dt = horizontal_h(L, beta, hx, hy), horizontal_h1(L, beta, hx, hy)


    hx_old = h[0]
    hy_old = h[0]

    for j in range(n):
        print(j)
        for i in range(n):
            hx = h[j]
            hy = h[i]
            
            h_matrices = update_H(*h_matrices, hx, hy, hx_old, hy_old, L)
            h_matrices_dt = update_H(*h_matrices_dt, hx, hy, hx_old, hy_old, L)
            h_matrices_2dt = update_H(*h_matrices_2dt, hx, hy, hx_old, hy_old, L)
            f[j, i] = calculate_f_update_h(matrices, h_matrices, L, temp)
            print("h_matrices_updated")

            der_low = (calculate_f_update_h(matrices_dt, h_matrices_dt, L, temp + dt) - f[j, i]) / dt
            der_high = (calculate_f_update_h(matrices_2dt, h_matrices_2dt, L, temp + 2 * dt) - calculate_f_update_h(matrices_dt, h_matrices_dt,  L, temp + dt))/dt 
            u[j, i] = der_low
            cv[j, i] = (der_high - der_low)/dt
            print("cv_calculated")

            hy_old = hy
            hx_old = hx
        
    #wyk.imshow(cv)
        
    wyk.imshow(cv)

imshow_h(10)

plt.show()
