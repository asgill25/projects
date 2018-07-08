import numpy as np
import math
import time
import random
from itertools import zip_longest
import csv
import wave, struct


def gf_multTable(p):
    # most common irreducible generating polynomials gpoly
    if p == 2:
        gpoly = int('0b111', 2)
    elif p == 3:
        gpoly = int('0b1011', 2)
    elif p == 4:
        gpoly = int('0b10011', 2)
    elif p == 5:
        gpoly = int('0b100101', 2)
    elif p == 6:
        gpoly = int('0b1000011', 2)
    elif p == 7:
        gpoly = int('0b10001001', 2)
    elif p == 8:
        gpoly = int('0b100011101', 2)
    elif p == 9:
        gpoly = int('0b1000010001', 2)
    elif p == 10:
        gpoly = int('0b10000001001', 2)
    elif p == 11:
        gpoly = int('0b100000000101', 2)
    elif p == 12:
        gpoly = int('0b1000001010011', 2)
    elif p == 13:
        gpoly = int('0b10000000011011', 2)
    elif p == 14:
        gpoly = int('0b100010001000011', 2)
    elif p == 15:
        gpoly = int('0b1000000000000011', 2)
    elif p == 16:
        gpoly = int('0b10001000000001011', 2)
    elif p == 17:
        gpoly = int('0b100000000000001001', 2)
    elif p == 18:
        gpoly = int('0b1000000000010000001', 2)
    elif p == 19:
        gpoly = int('0b10000000000000100111', 2)
    elif p == 20:
        gpoly = int('0b100000000000000001001', 2)
    elif p == 21:
        gpoly = int('0b1000000000000000000101', 2)
    elif p == 22:
        gpoly = int('0b10000000000000000000011', 2)
    elif p == 23:
        gpoly = int('0b100000000000000000100001', 2)
    elif p == 24:
        gpoly = int('0b1000000000000000010000111', 2)
    else:
        return -1

    # initialization
    GF = [0]*(2**p-1)
    iGF = [0]*(2**p)
    GF[0], iGF[1], b = 1, 0, 1

    for i in range(1, 2**p-1):
        b = b << 1        # multiply by alpha
        if b.bit_length() > p:
            b = b ^ gpoly  # remap (xor) using gpoly
        GF[i] = b
        iGF[b] = i
        
    return GF, iGF


def gf_mult(a, b, GF, iGF):
    # scalar multiplication over GF(2^p)
    # returns a*b
    if a == 0 or b == 0:
        return 0
    k = len(GF)
    ab = (iGF[a] + iGF[b]) % k
    return GF[ab]


def gf_div(a, b, GF, iGF):
    # scalar division over GF(2^p)
    # Returns c = a*b^(-1)
    if b != 0:
        if a == 0:
            return 0
        k = len(GF)
        ab = (iGF[a] - iGF[b]) % k
        return GF[ab]
    else:
        print('divide by zero')
    return c


def gf_pmult(p1, p2, GF, iGF):
    # polynomial multiplication over GF(2^p)
    # Takes p1 and p2 at returns pp=p1*p2
    l1, l2 = len(p1), len(p2)
    pp = np.zeros(l2+l1-1, dtype=int)
    for i in range(l2):
        for j in range(i, i+l1):
            m = gf_mult(p1[j-i], p2[i], GF, iGF)
            pp[j] ^= m
    return pp


def gf_pdiv(pnum, pden, GF, iGF):
    # polynomial division over GF(2^p)
    # takes pnum and pden and returns q and r such that
    # pnum=pden*q+r
    lnum, lden = len(pnum), len(pden)
    k = lnum - lden
    if k < 0:
        return np.array([0]), pnum
    d = np.zeros(lden, dtype=int)
    q = np.zeros(k+1, dtype=int)
    r = pnum[0:lden-1]
    for i in range(k+1):
        r = np.append(r, pnum[lden-1+i])
        fac = gf_div(r[0], pden[0], GF, iGF)
        for j in range(lden):
            dd = gf_mult(fac, pden[j], GF, iGF)
            d[j] = r[j] ^ gf_mult(fac, pden[j], GF, iGF)
        r = d[1:]
        q[i] = fac
    return q, r


# def evalPoly(P, a, GF, iGF):
#     # evaluate polynomial (Horner)
#     # takes the polynomial P and returns v=P(a)
#     result = 0
#     order = len(P) - 1
#     if order == 0:
#         return P[0]
#     else:
#         return P[order] ^ gf_mult(a, evalPoly(P[:order], a, GF, iGF), GF, iGF)

def evalPoly(P, a, GF, iGF):
    # evaluate polynomial (Horner)
    v = 0
    for r in P:
        v = (gf_mult(v, a, GF, iGF)) ^ r
    return v


def encodingPoly(n, GF, iGF):
    # Compute encoding polynomial G for n parity bits.
    G = np.array([1])
    for i in range(0, n):
        G = gf_pmult(np.array([1, GF[i]]), G, GF, iGF)
    return G


def RS_encode(m, npar, GF, iGF):
    # Takes the message polynomial m and returns the
    # RS encoded message C with npar parity bits added.

    num = np.zeros(len(m) + npar).astype(int)
    num[0:len(m)] = m
    den = encodingPoly(npar, GF, iGF)
    Q, R = gf_pdiv(num, den, GF, iGF)

    C = num
    C[len(m):] = R
    return C


def compSyndromePoly(R, n, GF, iGF):
    # computes syndrome polynomial for the RS encoded polynomial
    # R with n parity bits.

    # R = np.array([4, 0, 0, 0, 0, 0, 1, 0])
    alpha = GF[0:n]
    syndromePoly = np.array([evalPoly(R, i, GF, iGF) for i in reversed(alpha)])
    return syndromePoly


def trim(a):
    # trim off zeros
    ia = 0
    while (a[0] == 0):
        a, ia = a[1:], 1
    return a, ia


def writeAudio(audio, fileName):
    audioLength = len(audio)
    params = wave._wave_params(nchannels=1, sampwidth=1, framerate=11025/3, nframes=audioLength, comptype='NONE', compname='not compressed')
    audioFile = wave.open(fileName, 'wb')
    audioFile.setparams(params)
    # wavef.writeframes(audio)
    for value in audio:
        data = struct.pack('<B', value)
        audioFile.writeframesraw(data)
    audioFile.close()
    return 0


def xEuclid(x2t, S, GF, iGF):
    # extended Euclid algorithm for polynomials
    S, iS = trim(S)  # trim S
    tdeg = (len(x2t)-1)/2
    a1, a2 = x2t, S    # initialization
    v1, v2 = [0], [1]
    while len(a2) >= tdeg+1:  # xEuclid
        q, a = gf_pdiv(a1, a2, GF, iGF)
        v3 = gf_pmult(v2, q, GF, iGF)
        z0 = np.zeros(len(v3)-len(v1), dtype=int)
        v = v3 ^ np.append(z0, v1)
        a1, a2 = a2, a   # update
        v1, v2 = v2, v
        a2, ia = trim(a2)  # trim a
    if (ia == 1) or (iS == 1):  # 'degenerate' case
        return v, a
    else:                  # standard case
        vf = gf_div(1, v[-1], GF, iGF)
        v = gf_pmult(v, [vf], GF, iGF)
        a = gf_pmult(a, [vf], GF, iGF)
        return v, a


def locateError(Lam, GF, iGF):
    # error locator from Lambda (Chien)

    # we go through every power of alpha to determine the roots
    eval_vals = np.array([evalPoly(Lam, a, GF, iGF) for a in GF])
    error_locs = np.where(eval_vals == 0)[0]

    # now take inverses of these
    error_locs = -error_locs % len(GF)
    return error_locs


def determineError(Lam, Om, loc, GF, iGF):
    # takes the polynomials Lam, Om, a list of the error locations loc
    # GF and iGF and returns a list of the corresponding error locations ev.
    ev = np.zeros(loc.shape).astype(int)
    p = len(GF)

    # compute derivative of Lambda
    dLam = Lam[:]
    for i in range(0, len(Lam), 2):
        dLam[len(Lam) - 1 - i] = 0
    dLam = trim(dLam[:len(dLam)-1])[0]

    for i, j in zip(range(len(loc)), loc):
        ev[i] = gf_div(evalPoly(Om, GF[-j % p], GF, iGF), evalPoly(dLam, GF[-j % p], GF, iGF), GF, iGF)
        ev[i] = gf_mult(ev[i], GF[j], GF, iGF)
    return ev


def RS_decode(R, npar, GF, iGF):
    # Takes the encoded polynomial R, number of parity bits npar
    # GF and iGF and returns the decoded polynomial with no errors CC
    x2t = np.zeros(npar + 1).astype(int)
    x2t[0] = 1
    syndromePoly = compSyndromePoly(R, npar, GF, iGF)
    if syndromePoly.max() > 0:
        lam, om = xEuclid(x2t, syndromePoly, GF, iGF)
        lam = trim(lam)[0]
        om = trim(om)[0]
        error_locs = locateError(lam, GF, iGF)
        error = determineError(lam, om, error_locs, GF, iGF)
        CC = R.copy()
        CC[len(CC) - error_locs - 1] = CC[len(CC) - error_locs - 1] ^ error
    else:
        CC = R.copy()
    return CC


if __name__ == '__main__':
    # code used outside of the above functions can be found in extraCode.py
    print('code used outside of the above functions can be found in extraCode.py')
