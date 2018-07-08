import numpy as np
import scipy.sparse as sp  # spdiags , identity
import scipy.sparse.linalg
import matplotlib.pyplot as plt


def DiscrtMat(nx, ny, h):
    # nx := domain dimension in x
    # ny := domain dimension in y
    # h  := discretization (uniform)
    # ----------
    # OUTPUTS
    # L := discretization matrix of the Laplacian

    block = (-4/(h**2))*sp.eye(nx) + sp.eye(nx, k=1)/(h**2) + sp.eye(nx, k=-1)/(h**2)
    b = [block]*ny
    L = sp.block_diag(b, format="csr")

    ones_matrix = sp.eye(nx*ny, k=nx)/(h**2) + sp.eye(nx*ny, k=-nx)/(h**2)
    L = L + ones_matrix
    return L


def gs(A, b, x0, its):
    # INPUTS
    # A := Discretization Matrix
    # b := RHS vector
    # x0 := initial guess
    # its := max number of iteration to perform
    # ----------
    # OUTPUTS
    # x  := solution after its iterations
    # er := necessary only for problem 1.1

    x = x0
    err = np.zeros(its)

    C = A[1,1]
    nm = len(x0)

    for i in range(1, its + 1):
        for j in range(0, nm):
            laplacian = (A*x).flatten()
            x[j] = (b[j]-laplacian[j])/C + x[j]
        res = (A * x).flatten() - b
        err[i-1] = np.linalg.norm(res)
        print(i/its)

    return x, err


def RestrictionMat(m, n):
    # m, n := rows & cols @ level h
    # OUTPUT----
    # R  := Restriction matrix

    new_m = int((m-1)/2)
    new_n = int((n-1)/2)
    off_vec = np.array([1/16, 1/8, 1/16])
    on_vec = np.array([1/8, 1/4, 1/8])
    off_block = np.zeros((new_m, m))
    on_block = np.zeros((new_m, m))

    for i in range(0, new_m):
        off_block[i,2*i:2*i + 3] = off_vec

    for i in range(0, new_m):
        on_block[i,2*i:2*i + 3] = on_vec

    off_block = sp.csr_matrix(off_block)
    on_block = sp.csr_matrix(on_block)

    R = sp.csr_matrix((new_m*new_n, m * n))

    long_block = sp.hstack([off_block, on_block, off_block])

    for i in range(0, new_n):
        R[new_m*i : new_m*(i+1), 2*m*i : 2*m*i + 3*m] = long_block
    return sp.csr_matrix(R)


########################################################################################################
### YOU MAY FIND THE FOLLOWING USEFUL, x, y is the grid, k is a counter

def C_initial(x, y, k):
    f = np.sin(k * np.pi * x[None, :]) * np.sin(k * np.pi * y[:, None])
    return f


def FCircle(x, y):
    f = y[:, None] ** 2 + x[None, :] ** 2
    return f


def F_RHS(x, y, k):
    if k == 0:
        f = -2 * ((1 - 6 * x[None, :] ** 2) * (y[:, None] ** 2) * (1 - y[:, None] ** 2) + (1 - 6 * y[:, None] ** 2) * (
                x[None, :] ** 2) * (1 - x[None, :] ** 2));
    else:
        f = np.zeros((len(y), len(x)))
    return f


def Sol_Analitical(x, y, k):
    if k == 0:
        S = (x[None, :] ** 2 - x[None, :] ** 4) * (y[:, None] ** 4 - y[:, None] ** 2)
    else:
        S = 10 * np.sin(np.pi * x[None, :]) * np.sinh(np.pi * y[:, None]) / np.sinh(np.pi)
    return S


#########################################################################################################


if __name__ == '__main__':
    # Problem set up
    len_y = 2.
    len_x = 1.
    kk = 9  # exponent of 2 to be used
    h = 1. / (2. ** kk)  # grid spacing, in this case 128 intervals
    yh = np.arange(0., len_y + h, h)
    xh = np.arange(0., len_x + h, h)
    [X, Y] = np.meshgrid(xh, yh)
    n = len(xh) - 2  # x dir = cols
    m = len(yh) - 2  # y dir = rows

    D = DiscrtMat(n, m, h)  # generates the discretization matrix of the Laplacian
    # D_inv = sp.linalg.inv(D)

    # ----------------------------------------------------------------------------
    # ------------------------ first function inversion --------------------------
    # ----------------------------------------------------------------------------

    # f_vals = F_RHS(xh[1:n+1], yh[1:m+1], 0).flatten()
    # my_sol1 = np.zeros((n + 2, m + 2))
    # my_sol1[1:n + 1, 1:m + 1] = sp.linalg.spsolve(D, f_vals).reshape(n, m)

    # actual_sol1 = Sol_Analitical(xh, yh, 0)

    # plt.figure()
    # actual_contour1 = plt.contourf(X, Y, actual_sol1, cmap=plt.cm.jet)
    # # plt.clabel(actual_contour, inline=1, fontsize=7)
    # plt.title('Analytical Solution (1)')
    # plt.colorbar()
    # # plt.savefig('1_analytical.png')

    # plt.figure()
    # sol_contour1 = plt.contourf(X, Y, my_sol1, actual_contour1.levels, cmap=plt.cm.jet)
    # # plt.clabel(sol_contour, inline=1, fontsize=7)
    # plt.title('Inverted Solution (1)')
    # plt.colorbar()
    # # plt.savefig('1_inverted.png')

    # norm_l2_1 = np.linalg.norm(my_sol1 - actual_sol1, 2)
    # norm_l_inf_1 = np.linalg.norm(my_sol1 - actual_sol1, np.inf)

    # print('L_2 norm (1): ' + str(norm_l2_1))
    # print('L_\u221E norm (1): ' + str(norm_l_inf_1))

    # ----------------------------------------------------------------------------
    # ------------------------- second function inversion ------------------------
    # ----------------------------------------------------------------------------

    # boundary_vals = 10 * np.sin(np.pi * xh)
    # C_vec = np.append(np.zeros(n*(m-1)), boundary_vals[1:n+1])
    # my_sol2 = np.zeros((n+2, m+2))
    # my_sol2[n+1, :] = boundary_vals
    # u_sol = sp.linalg.spsolve(D, C_vec/(-(h**2)))
    # my_sol2[1:n+1, 1:m+1] = u_sol.reshape(n, m)

    # actual_sol2 = Sol_Analitical(xh, yh, 1)

    # plt.figure()
    # actual_contour2 = plt.contourf(X, Y, actual_sol2, cmap=plt.cm.jet)
    # # plt.clabel(actual_contour, inline=1, fontsize=7)
    # plt.title('Analytical Solution (2)')
    # plt.colorbar()
    # # plt.savefig('2_analytical.png')

    # plt.figure()
    # sol_contour2 = plt.contourf(X, Y, my_sol2, actual_contour2.levels, cmap=plt.cm.jet)
    # # plt.clabel(sol_contour, inline=1, fontsize=7)
    # plt.title('Inverted Solution (2)')
    # plt.colorbar()
    # # plt.savefig('2_inverted.png')

    # plt.show()

    # norm_l2_2 = np.linalg.norm(my_sol2 - actual_sol2, 2)
    # norm_l_inf_2 = np.linalg.norm(my_sol2 - actual_sol2, np.inf)

    # print('L_2 norm (2): ' + str(norm_l2_2))
    # print('L_\u221E norm (2): ' + str(norm_l_inf_2))

    # ----------------------------------------------------------------------------
    # ------------------------- gauss-sidel 1st function -------------------------
    # ----------------------------------------------------------------------------

    # u_gs_1 = gs(D, f_vals, np.zeros(n*m), 50)[0]
    # sol_gs_1 = np.zeros((n+2, m+2))
    # sol_gs_1[1:n+1, 1:m+1] = u_gs_1.reshape(n,m)

    # plt.figure()
    # plt.contourf(X, Y, sol_gs_1, cmap=plt.cm.jet)
    # plt.title('Gauss-Seidel Solution, 50 iterates (1)')
    # plt.colorbar()
    # # plt.savefig('sol_gs_1.png')
    # plt.show()

    # ----------------------------------------------------------------------------
    # ------------------------- gauss-sidel 2nd function -------------------------
    # ----------------------------------------------------------------------------

    # u_gs_2 = gs(D, C_vec/(-(h**2)), np.zeros(n*m), 50)[0]
    # sol_gs_2 = np.zeros((n+2, m+2))
    # sol_gs_2[n + 1, :] = boundary_vals
    # sol_gs_2[1:n+1, 1:m+1] = u_gs_2.reshape(n,m)

    # plt.figure()
    # plt.contourf(X, Y, sol_gs_2, cmap=plt.cm.jet)
    # plt.title('Gauss-Seidel Solution, 50 iterates (2)')
    # plt.colorbar()
    # plt.savefig('sol_gs_2.png')
    # plt.show()

    # ----------------------------------------------------------------------------
    # ----------------------- comparing start-vals with gs -----------------------
    # ----------------------------------------------------------------------------

    # gs_sol = np.zeros((n+2, m+2))
    #
    # start_val_2 = C_initial(xh, yh, 2)[1:n+1, 1:m+1].flatten()
    # start_val_4 = C_initial(xh, yh, 4)[1:n+1, 1:m+1].flatten()
    # start_val_8 = C_initial(xh, yh, 8)[1:n+1, 1:m+1].flatten()
    #
    # err_2 = gs(D, np.zeros(n * m), start_val_2, 100)[1]
    # err_4 = gs(D, np.zeros(n * m), start_val_4, 100)[1]
    # err_8 = gs(D, np.zeros(n * m), start_val_8, 100)[1]
    # its_list = range(0, 100)
    #
    # plt.figure()
    # plt.plot(its_list, err_2, 'r', its_list, err_4, 'b', its_list, err_8, 'g')
    # plt.legend(['k = 2', 'k = 4', 'k = 8'])
    # plt.ylim([0., 30000])
    # plt.show()

    #################################################
    # Testing Restriction Matrix and Prolongation
    #################################################

    # res_mat = RestrictionMat(n, m)
    # pro_mat = 4 * res_mat.transpose()
    # original = FCircle(xh[1:n+1],yh[1:m+1])
    # restricted = res_mat * original.flatten()
    # prolonged = pro_mat * restricted
    # res_vals = restricted.reshape(int((m-1)/2),int((n-1)/2))
    # pro_vals = prolonged.reshape(m, n)
    # [X, Y] = np.meshgrid(np.arange(1/(2**(kk-1)),1,1/(2**(kk-1))), np.arange(1/(2**(kk-1)),2,1/(2**(kk-1))))
    # plt.contourf(X,Y,res_vals, cmap=plt.cm.jet)
    # plt.title('Restricted')
    # plt.colorbar()
    # plt.savefig('restricted_9.png')
    # plt.show()
    #
    # print(pro_mat*res_mat)
