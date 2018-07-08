import numpy as np
import scipy as sp
import random
import time
from math import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm


def SSF(x):
    # fitness function
    o = np.ones(np.shape(x)[0])
    f_bias = -450

    return sum((x-o)**2) + f_bias

# =================================================================
# Exercise 1
# Test functions


def SRF(x):
    o = np.ones(np.shape(x)[0])
    z = x - o
    f_bias = -330
    f = f_bias
    for i in range(0, len(x)):
        f += z[i]**2 - 10*np.cos(2 * np.pi * z[i]) + 10

    return f


def SSP(x):
    o = np.ones(np.shape(x)[0])
    f_bias = -450
    f = f_bias

    for i in range(len(x)):
        f += sum((x - o)[0:i + 1])**2

    return f


def SRosF(x):
    o = np.ones(np.shape(x)[0])
    z = x - o + 1
    f_bias = 390
    f = f_bias
    for i in range(0, len(x) - 1):
        f += 100*(z[i]**2 - z[i+1])**2 + (z[i] - 1)**2

    return f
# ==================================================================


# ==================================================================
# Exercise 2


def Strategy(Pop, F, Rlist, strat=1, iBest=-1):
    """ Function that takes in Population matrix Pop,
    parameters F, list of random numbers Rlist and variable strat
    to indicate which strategy used """
    
    if strat == 1:
        return Pop[Rlist[0]] + F*(Pop[Rlist[1]] - Pop[Rlist[2]])
    elif strat == 2:
        return Pop[Rlist[0]] + F * (Pop[Rlist[1]] - Pop[Rlist[2]]) + F * (Pop[Rlist[3]] - Pop[Rlist[4]])
    elif strat == 3:
        return Pop[iBest] + F * (Pop[Rlist[0]] - Pop[Rlist[1]]) + F * (Pop[Rlist[2]] - Pop[Rlist[3]])
# ==================================================================


def DiffEvo(x, vMin, vMax, Ftup, nGen, strat=1):

    FnC = Ftup[0]
    o = Ftup[1]
    Fo = FnC(o)
    nPop, nV = np.shape(x)

    # differential evolution parameters
    F = 0.8     # differentiation constant
    CR = 0.25   # crossover constant

    # initialization
    X = np.zeros(nV)            # trial vector
    Pop = np.zeros((nPop, nV))     # population
    Fit = np.zeros(nPop)          # fitness of the population
    iBest = 1                       # index of the best solution

    # create population and evaluate fitness
    Pop = np.random.uniform(vMin, vMax, (nPop, nV))
    Fit = np.apply_along_axis(FnC, 1, Pop)

    # optimization

    for g in range(nGen):   # for each generation
        if g % 100 == 0:
            print('We are at iteration: ', g)
            print('Best solution so far: ', Pop[iBest, :])
        
        for j in range(nPop):
            # randomly sample 3 times from population
            jj = [k for k in range(nPop)]
            jj.pop(j)
            Rparam = random.sample(jj, 5)
            Rnd = int(np.floor(random.random()*nV))
            for i in range(nV):
                if (random.random() < CR) or (Rnd == i):
                    X[i] = Strategy(Pop[:, i], F, Rparam, strat, iBest)
                    if (X[i] < vMin) or (X[i] > vMax):
                        X[i] = vMin + (vMax-vMin)*random.random()
                else:
                    X[i] = Pop[j,i]

            f = np.apply_along_axis(FnC, 0, X)

            # if trial is better or equal than current
            if f <= Fit[j]:
                Pop[j,:] = X    # replace current by trial
                Fit[j] = f
                # if trial is better than the best
                if f <= Fit[iBest]:
                    iBest = j   # update the best's index
        
        # Stopping condition 
        if FnC(Pop[iBest, :]) == Fo or np.max(np.abs(Pop[iBest,:]-o)) < 10**(-7):
            print('DiffEvo stopped at: ', g)
            break
    
    return Pop[iBest,:]


# =================================================================
# Exercise 3
def DiffEvoSA(x, vMin, vMax, Ftup, nGen, strat=1):
    # Self-adaptive version of DiffEvo
    FnC = Ftup[0]
    o = Ftup[1]
    Fo = FnC(o)
    nPop, nV = np.shape(x)

    # differential evolution parameters
    F = 0.8 * np.ones(nPop)  # differentiation constant
    CR = 0.25 * np.ones(nPop)  # crossover constant

    # initialization
    X = np.zeros(nV)  # trial vector
    Pop = np.zeros((nPop, nV))  # population
    Fit = np.zeros(nPop)  # fitness of the population
    iBest = 1  # index of the best solution

    # create population and evaluate fitness
    Pop = np.random.uniform(vMin, vMax, (nPop, nV))
    Fit = np.apply_along_axis(FnC, 1, Pop)

    # optimization

    for g in range(nGen):  # for each generation
        if g % 100 == 0:
            print('We are at iteration: ', g)
            print('Best solution so far: ', Pop[iBest, :])

        for j in range(nPop):
            # randomly sample 3 times from population
            jj = [k for k in range(nPop)]
            jj.pop(j)
            Rparam = random.sample(jj, 5)
            Rnd = int(np.floor(random.random() * nV))
            for i in range(nV):
                if (random.random() < CR[j]) or (Rnd == i):
                    X[i] = Strategy(Pop[:, i], F[j], Rparam, strat, iBest)
                    if (X[i] < vMin) or (X[i] > vMax):
                        X[i] = vMin + (vMax - vMin) * random.random()
                else:
                    X[i] = Pop[j, i]

            f = np.apply_along_axis(FnC, 0, X)

            # if trial is better or equal than current
            if f <= Fit[j]:
                Pop[j, :] = X  # replace current by trial
                Fit[j] = f
                # if trial is better than the best
                if f <= Fit[iBest]:
                    iBest = j  # update the best's index

            # self-adaptive differential evolution
            if random.random() < 0.1:
                F[j] = 0.1 + 0.9 * random.random()
            if random.random() < 0.1:
                CR[j] = random.random()

        # Stopping condition
        if FnC(Pop[iBest, :]) == Fo or np.max(np.abs(Pop[iBest, :] - o)) < 10 ** (-7):
            print('DiffEvoSA stopped at: ', g)
            break

    return Pop[iBest, :]
# =================================================================


# =================================================================
# Exercise 4
def DiffEvoSAPD(x, vMin, vMax, Ftup, nGen, GR, strat=1):
    """ 
    Self-adaptive differential evolution algorithm with population size reduction
    takes in variable array x, optimisation domain vmin, vmax, test function FnC
    and a list GR, containing the generations where we reduce the population """
    FnC = Ftup[0]
    o = Ftup[1]
    Fo = FnC(o)
    nPop, nV = np.shape(x)

    # differential evolution parameters
    F = 0.8 * np.ones(nPop)  # differentiation constant
    CR = 0.25 * np.ones(nPop)  # crossover constant

    # initialization
    X = np.zeros(nV)  # trial vector
    Pop = np.zeros((nPop, nV))  # population
    Fit = np.zeros(nPop)  # fitness of the population
    iBest = 1  # index of the best solution

    # create population and evaluate fitness
    Pop = np.random.uniform(vMin, vMax, (nPop, nV))
    Fit = np.apply_along_axis(FnC, 1, Pop)

    # optimization

    for g in range(nGen):  # for each generation
        if g % 100 == 0:
            print('We are at iteration: ', g)
            print('Best solution so far: ', Pop[iBest, :])

        for j in range(nPop):
            # randomly sample 3 times from population
            jj = [k for k in range(nPop)]
            jj.pop(j)
            Rparam = random.sample(jj, 5)
            Rnd = int(np.floor(random.random() * nV))
            for i in range(nV):
                if (random.random() < CR[j]) or (Rnd == i):
                    X[i] = Strategy(Pop[:, i], F[j], Rparam, strat, iBest)
                    if (X[i] < vMin) or (X[i] > vMax):
                        X[i] = vMin + (vMax - vMin) * random.random()
                else:
                    X[i] = Pop[j, i]

            f = np.apply_along_axis(FnC, 0, X)

            # if trial is better or equal than current
            if f <= Fit[j]:
                Pop[j, :] = X  # replace current by trial
                Fit[j] = f
                # if trial is better than the best
                if f <= Fit[iBest]:
                    iBest = j  # update the best's index

            # self-adaptive differential evolution
            if random.random() < 0.1:
                F[j] = 0.1 + 0.9 * random.random()
            if random.random() < 0.1:
                CR[j] = random.random()

        # Stopping condition
        if FnC(Pop[iBest, :]) == Fo or np.max(np.abs(Pop[iBest, :] - o)) < 10 ** (-7):
            print('DiffEvoSAPD stopped at: ', g)
            break

        # population reduction
        if g in GR:
            nPop = int((nPop + 1) / 2)
            pop1 = Pop[0:nPop, :]
            pop2 = Pop[nPop:, :]
            if iBest < nPop and iBest >= int((len(pop1))/2):
                iBest = iBest - int((len(pop1))/2)
            if nPop <= iBest < nPop + int((len(pop2))/2):
                iBest = iBest - int((len(pop1))/2)
            if iBest >= nPop + int((len(pop2)) / 2):
                iBest = iBest - int((len(pop2))/2)
                iBest = iBest - int((len(pop1))/2)
            for i in range(0, int((len(pop1) + 1)/2)):
                if FnC(pop1[i + int((len(pop1))/2), :]) < FnC(pop1[i, :]):
                    pop1[i, :] = pop1[i + int((len(pop1))/2), :]
            for i in range(0, int((len(pop2) + 1)/2)):
                if FnC(pop2[i + int((len(pop2))/2), :]) < FnC(pop2[i, :]):
                    pop2[i, :] = pop2[i + int((len(pop2))/2), :]
            pop1 = pop1[0:int((len(pop1) + 1)/2), :]
            pop2 = pop2[0:int((len(pop2) + 1) / 2), :]
            Pop = np.concatenate((pop1, pop2))
    return Pop[iBest, :]
# =================================================================


if __name__ == '__main__':

    # Marker main code
    # Code used by marker to test code

    ##########################################
    # PLEASE NOTE!!!!! I made small correction to the marker main code:
    #                                   Line 327. Replaced y1 with y2 to print correct value
    # Also added strat as parameter to DiffEvo functions and Python 3 style print
    ##########################################

    nGen = 1000
    vMin, vMax = -3, 3
    x     = np.random.uniform(vMin,vMax,(1000,10))

    # Pick test function and strategy
    FnC   = SRosF
    Ftup  = (FnC,np.ones(np.shape(x)[1]))
    strat = 3

    print('Differential evolution: ')
    t0 = time.time()
    y1 = DiffEvo(x,vMin,vMax,Ftup,nGen,strat)
    t1 = time.time()
    T1 = t1-t0
    print('time it took = ', T1)
    print('Best solution is given by: ', y1, ' with function value: ', FnC(y1))

    print('Self-adaptive Differential evolution: ')
    t0 = time.time()
    y2 = DiffEvoSA(x,vMin,vMax,Ftup,nGen,strat)
    t1 = time.time()
    T2 = t1-t0
    print('time it took = ', T2)
    print('Best solution is given by: ', y2, ' with function value: ', FnC(y2))

    print('Self-adaptive differential evolution with Population reduction: ')
    t0 = time.time()
    GR = [125,250,500,961]
    y3 = DiffEvoSAPD(x,vMin,vMax,Ftup,nGen,GR,strat)
    t1 = time.time()
    T3 = t1-t0
    print('time it took = ', T3)
    print('Best solution is given by: ', y3, ' with function value: ', FnC(y3))
    print()
    print()
    print('Differential evolution: ', T1, FnC(y1))
    print('Self-adaptive Differential evolution: ', T2, FnC(y2))
    print('Self-adaptive Differential evolution with population reduction: ', T3, FnC(y3))