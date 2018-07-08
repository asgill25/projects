# -*- coding: utf-8 -*-
import numpy as np
from scipy import ndimage
from scipy import misc
import scipy.sparse as sp # spdiags , identity
import scipy.sparse.linalg
import matplotlib.pyplot as plt
from functools import reduce


def Neighbours(msk):

    non_zero_elts = np.nonzero(msk)
    row_length = msk.shape[1]
    num = len(non_zero_elts[0])
    i_c = np.zeros(num)
    i_n = np.zeros(num)
    i_s = np.zeros(num)
    i_w = np.zeros(num)
    i_e = np.zeros(num)

    for i in range(0, num):
        [j, k] = [non_zero_elts[0][i], non_zero_elts[1][i]]
        i_c[i] = row_length*j + k
        i_n[i] = row_length*(j - 1) + k
        i_s[i] = row_length*(j + 1) + k
        i_w[i] = row_length*j + k - 1
        i_e[i] = row_length*j + k + 1

    neigh = np.array([i_c, i_n, i_s, i_w, i_e])
    return neigh
 

def DiscFun(Neigh, f):
    no_points = len(Neigh[0])
    D = (4 * sp.eye(no_points)).tolil()
    bc = np.zeros(no_points)

    for i in range(1, 5):
        col_vals = np.isin(Neigh[i], Neigh[0])
        row_vals = np.isin(Neigh[0], Neigh[i])
        boundary_vals = np.isin(Neigh[i], Neigh[0], invert=True)
        D[row_vals, col_vals] = -1
        # print(Neigh[i][boundary_vals])
        f_vals = f[Neigh[i][boundary_vals].astype(int)]
        bc[boundary_vals] = bc[boundary_vals] + f_vals

    return D.tocsr(), bc


def Deposition(IB, p, alpha):
    # IB := indices of the neighbours of the fractals
    # p := solution field everywhere
    #alpha  := exponent
    # −−−-----------------------
    # Idx   := index to place next particle

    u = p[IB.astype(int)]**alpha
    some_sum = sum(u)

    u = u/some_sum
    v = np.cumsum(u)
    r = np.random.random()
    Idx = IB[v > r][0]
    return Idx


if __name__ == '__main__':
    grid_len = 256
    alpha = 3.0/2.0
    grid_shape = (grid_len, grid_len)
    domain = np.zeros(grid_shape)
    omega = np.zeros(grid_shape)
    omega[1:grid_len-1, 1:grid_len-1] = np.ones(np.subtract(grid_shape, (2, 2)))
    centre_point = (int(grid_len/2), int(grid_len/2))
    domain[centre_point] = 1
    X = np.arange(0., grid_len, 1.)
    Y = np.arange(0., grid_len, 1.)
    [X, Y] = np.meshgrid(X, Y)
    boundary_condition = np.log(np.sqrt((np.absolute(X - centre_point[0])**2 + np.abs(Y - centre_point[1])**2))) / (2 ** np.pi)
    boundary_condition[1:grid_len-1, 1:grid_len-1] = np.zeros(np.subtract(grid_shape, (2, 2)))
    # domain[int(grid_len / 2 + 1), int(grid_len / 2)] = 1
    # domain[int(grid_len / 2 + 2), int(grid_len / 2)] = 1
    # domain[int(grid_len / 2 + 3), int(grid_len / 2)] = 1
    # domain[int(grid_len / 2 + 4), int(grid_len / 2)] = 1
    # domain[int(grid_len / 2 + 4), int(grid_len / 2 + 1)] = 1
    # domain[int(grid_len / 2 + 4), int(grid_len / 2 + 2)] = 1
    # domain[int(grid_len / 2 + 4), int(grid_len / 2 + 3)] = 1
    # domain[int(grid_len / 2 + 4), int(grid_len / 2 + 4)] = 1
    # domain[int(grid_len / 2 + 3), int(grid_len / 2 + 4)] = 1
    # domain[int(grid_len / 2 + 2), int(grid_len / 2 + 4)] = 1
    # domain[int(grid_len / 2 + 1), int(grid_len / 2 + 4)] = 1
    # domain[int(grid_len / 2), int(grid_len / 2 + 4)] = 1
    # domain[int(grid_len / 2), int(grid_len / 2 + 3)] = 1
    # domain[int(grid_len / 2), int(grid_len / 2 + 2)] = 1
    # domain[int(grid_len / 2), int(grid_len / 2 + 1)] = 1
    for i in range(0, 4000):
        print(i)
        Neigh = Neighbours(domain)
        flat_domain = domain.flatten()
        flat_domain[Neigh[1].astype(int)] = 1
        flat_domain[Neigh[2].astype(int)] = 1
        flat_domain[Neigh[3].astype(int)] = 1
        flat_domain[Neigh[4].astype(int)] = 1
        flat_domain[Neigh[0].astype(int)] = 0
        # interior_vals = reduce(np.intersect1d, (Neigh[1], Neigh[2], Neigh[3], Neigh[4]))
        # flat_domain[interior_vals.astype(int)] = 0
        boundary_domain = flat_domain.reshape(domain.shape)
        boundary_indices = np.arange(0, grid_len**2, 1)[np.nonzero(flat_domain)]
        disc = DiscFun(Neighbours(omega - domain), (domain + boundary_condition).flatten())
        u_sol = np.zeros(domain.shape).flatten()
        u_sol[np.nonzero((omega - domain).flatten())] = sp.linalg.spsolve(disc[0], disc[1])
        # initial_domain[np.nonzero(boundary_domain)] = u_sol[np.nonzero(boundary_domain)]
        field_sol = u_sol.reshape(domain.shape) + domain

        new_index = Deposition(boundary_indices, field_sol.flatten(), alpha)
        domain = domain.flatten()
        domain[new_index] = 1
        domain = domain.reshape(grid_shape)
    # An example to plot the fractals ( in the matrix) as well as the function f
    # THIS IS  JUST AN EXAMPLE WITH RANDOM DATA
    # to be replaced with the solution field and matrix

    # # Make data.
    # X = np.arange(-1, 1, 0.25)
    # Y = np.arange(-1, 1, 0.25)
    # X, Y = np.meshgrid(X, Y)
    # A = np.sqrt(X**2 + Y**2)
    # # matrix containing the Fractals positions
    # B = np.random.randint(2, size=np.shape(A))

    # example figure & save
    fig, axes = plt.subplots(1, 1)
    axes.matshow(domain, cmap=plt.cm.Blues)
    axes.contour(field_sol, cmap=plt.cm.hot)
    # DON'T COMMENT THIS OUT!!! JUST IN CASE IT WORKS!
    fig.savefig('Fractal.eps', format='eps')   # save the figure to file
    plt.show()
