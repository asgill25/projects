    """
    
    nGen = 1000
    vMin,vMax = -100,100
    x     = np.random.uniform(vMin,vMax,(1000,10))

    # Pick test function and strategy
    FnC   = SSF
    Ftup  = (FnC,np.ones(np.shape(x)[1]))
    strat = 3
    
    print'Differential evoultuion: '
    t0 = time.time()
    y1 = DiffEvo(x,vMin,vMax,Ftup,nGen)
    t1 = time.time()
    T1 = t1-t0
    print 'time it took = ', T1
    print 'Best solution is given by: ', y1, ' with function value: ', FnC(y1)

    print 'Self-adaptive Differential evoultuion: '
    t0 = time.time()
    y2 = DiffEvoSA(x,vMin,vMax,Ftup,nGen)
    t1 = time.time()
    T2 = t1-t0
    print 'time it took = ', T2
    print 'Best solution is given by: ', y1, ' with function value: ', FnC(y2)


    print 'Self-adaptive differnatial evolution with Population reduction: '
    t0 = time.time()
    GR = [125,250,500,961]
    y3 = DiffEvoSAPD(x,vMin,vMax,Ftup,nGen,GR)
    t1 = time.time()
    T3 = t1-t0
    print 'time it took = ', T3
    print 'Best solution is given by: ', y3, ' with function value: ', FnC(y3)
    print
    print
    print 'Differential evolution: ', T1, FnC(y1)
    print 'Self-adaptive Differential evolution: ', T2, FnC(y2)
    print 'Self-adaptive Differential evolution with population reduction: ', T3, FnC(y3) """