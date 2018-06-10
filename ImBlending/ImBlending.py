import numpy as np
from scipy import ndimage
from scipy import misc
import scipy.sparse as sp  # spdiags , identity
import scipy.sparse.linalg
import matplotlib.pyplot as plt
import xlrd
from functools import reduce


def divergence2d(F):
    return np.gradient(F[0])[0] + np.gradient(F[1])[1]


def mixed_gradient(grad1, grad2):
    norm1 = np.square(grad1[0]) + np.square(grad1[1])
    norm2 = np.square(grad2[0]) + np.square(grad2[1])
    higher_vals = np.greater(norm2, norm1)
    lower_vals = np.invert(higher_vals)
    mixed = [np.zeros(grad1[0].shape), np.zeros(grad1[0].shape)]
    mixed[0][higher_vals] = grad2[0][higher_vals]
    mixed[1][higher_vals] = grad2[1][higher_vals]
    mixed[0][lower_vals] = grad1[0][lower_vals]
    mixed[1][lower_vals] = grad1[1][lower_vals]
    return mixed

def im2double(im):
    min_val = np.min(np.ravel(im))
    max_val = np.max(np.ravel(im))
    # print(min_val)
    out = (im.astype('float') - min_val) / (max_val - min_val)

    return out


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


def DiscFun_RGB(Neigh, f):
    # NEED TO FIX THIS THING!!!!
    no_points = len(Neigh[0])
    D = (4 * sp.eye(no_points)).tolil()
    bc = np.zeros((no_points, 3))

    for i in range(1, 5):
        col_vals = np.isin(Neigh[i], Neigh[0])
        row_vals = np.isin(Neigh[0], Neigh[i])
        boundary_vals = np.isin(Neigh[i], Neigh[0], invert=True)
        D[row_vals, col_vals] = -1
        for j in range(0, 3):
            bc[boundary_vals, j] = bc[boundary_vals, j] + f[Neigh[i][boundary_vals].astype(int), j]

    return D.tocsr(), bc


if __name__ == '__main__':

    # plt.close('all')

    # for colour images set RGB=1
    # the first part of the CW is enough to leave the image in BW
    RGB = 1

    # mask always loaded as BW
    msk = misc.imread('mask.jpg', mode='1')
    msk = im2double(msk)
    msk_pts = np.equal(msk.flatten(), 1)

    if RGB == 1:
        fstar = misc.imread('destination.jpg')  # fstar is the destination image
        gstar = misc.imread('source.jpg')  # gstar is the source image made the same size as fstar
    else:
        fstar = misc.imread('destination.jpg', mode='L')
        gstar = misc.imread('source.jpg', mode='L')

    # Convert image fields to float points
    fstar = im2double(fstar)
    gstar = im2double(gstar)

    if RGB == 1:
        disc = DiscFun_RGB(Neighbours(msk), fstar.reshape(-1, 3))
        grad_f = [np.gradient(fstar[:, :, i]) for i in range(0, 3)]
        grad_g = [np.gradient(gstar[:, :, i]) for i in range(0, 3)]

        rhs = [- divergence2d(mixed_gradient(grad_g[i], grad_f[i])).flatten()[msk_pts] for i in range(0, 3)]

        u_sol = [sp.linalg.spsolve(disc[0], rhs[i] + disc[1][:, i]) for i in range(0, 3)]

        u = fstar.reshape(-1, 3)
        for i in range(0, 3):
            u[msk_pts, i] = u_sol[i]

    else:
        disc = DiscFun(Neighbours(msk), fstar.flatten())
        grad_f = np.gradient(fstar)
        grad_g = np.gradient(gstar)
        rhs = - divergence2d(mixed_gradient(grad_g, grad_f)).flatten()[msk_pts]
        rhs = - divergence2d(grad_f).flatten()
        u_sol = sp.linalg.spsolve(disc[0], disc[1] + rhs)
        u = fstar.flatten()
        u[msk_pts] = u_sol

    u = im2double(u)
    # plotting
    plt.figure(1)
    # plt.subplot(121)
    plt.imshow(u.reshape(fstar.shape))

    # plt.subplot(122)
    # plt.imshow(u2.reshape(fstar.shape), cmap='gray')
    plt.show()

