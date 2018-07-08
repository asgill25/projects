import numpy as np


def PentRot(P, TransList):
    for pent in TransList:
        transforms = [P[pent]]
        no_rots = 0

        my_pent = P[pent]
        while no_rots < 4:
            flipped_pent = np.fliplr(my_pent)
            if not any([np.array_equal(i, flipped_pent) for i in transforms]):
                P[pent + str(no_rots) + 'fliprot'] = flipped_pent
                transforms.append(flipped_pent)
            rotated_pent = np.rot90(my_pent)
            no_rots += 1
            if not any([np.array_equal(i, rotated_pent) for i in transforms]):
                P[pent + str(no_rots) + 'rot'] = rotated_pent
                transforms.append(rotated_pent)
            my_pent = rotated_pent

    return P


def buildY(N, M, P, TransList=[]):

    if TransList:
        P = PentRot(P, TransList)
        TransList = [key for key in P]

    Y = {}
    for pent in TransList:
        pent_shape = P[pent].shape
        for row_start in range(0, (N + 1 - pent_shape[0])*M, M):
            for pos in range(row_start + 1, row_start + M + 2 - pent_shape[1]):
                marched_pent = np.array([range(pos + i, pos + i + pent_shape[1]) for i in range(0, M * pent_shape[0], M)])
                marched_pent[P[pent] == 0] = 0
                Y[pent + str(pos)] = marched_pent[np.nonzero(marched_pent)]

    return Y


def solveExactCover(X, Y, xcov, SolNum):
    if not X:
        print('FINAL SOLUTION = ', xcov)
        SolNum = SolNum + 1
        print('#', SolNum)
        return xcov, SolNum
    else:
        c = min(X, key=lambda c: len(X[c]))
        rr = X[c].copy()
        for r in rr:
            if not any(r[0] in t for t in xcov):
                # Above line makes sure that pentominos in different
                # locations are not counted as different shapes
                xcov.append(r)
                rows = cover(X, Y, r)
                xcov, SolNum = solveExactCover(X, Y, xcov, SolNum)
                uncover(X, Y, rows)
                xcov.pop()
        return xcov, SolNum


def cover(X, Y, r):
    # routine covers (eliminates) a row 'r'
    rows = set()
    for i in Y[r]:
        rows = rows | X[i]
    for j in rows:
        for k in Y[j]:
            X[k].remove(j)
            if (len(X[k]) == 0) and (k in Y[r]):
                del X[k]
    return rows


def uncover(X, Y, rows):
    # routine uncovers (restores) a set of rows
    for j in rows:
        for k in Y[j]:
            if k in X:
                X[k] |= {j}
            else:
                X[k] = {j}


if __name__ == '__main__':

    # These are the dimensions of our board
    nB = 3
    mB = 20

    # Initializing the pentomino dictionary here
    P = {}
    P['A'] = np.array([[1, 1, 0], [0, 1, 1], [0, 1, 0]])
    P['B'] = np.array([[1], [1], [1], [1], [1]])
    P['C'] = np.array([[0, 1, 1, 1], [1, 1, 0, 0]])
    P['D'] = np.array([[0, 1, 1], [0, 1, 0], [1, 1, 0]])
    P['E'] = np.array([[1, 0], [1, 0], [1, 0], [1, 1]])
    P['F'] = np.array([[1, 1, 1], [0, 0, 1], [0, 0, 1]])
    P['G'] = np.array([[0, 1, 1], [1, 1, 0], [1, 0, 0]])
    P['H'] = np.array([[0, 1], [0, 1], [1, 1], [0, 1]])
    P['I'] = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]])
    P['J'] = np.array([[0, 1, 0], [0, 1, 0], [1, 1, 1]])
    P['K'] = np.array([[1, 1, 0], [1, 1, 1]])
    P['L'] = np.array([[1, 1], [0, 1], [1, 1]])

    Y = buildY(nB, mB, P, [key for key in P])
    N = nB * mB
    K = range(1, N + 1)

    # build the X-dictionary from the Y-dictionary
    # You may use my approach or your own here
    X = {j: set() for j in K}
    for i in Y:
        for j in Y[i]:
            X[j].add(i)

    xcov = solveExactCover(X, Y, [], 0)