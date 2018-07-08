# simple Voronoi demo with matplotlib

import numpy as np
from Voronoi2Dsubmission import *

# visualization
import matplotlib.pyplot as plt
import matplotlib.tri
import matplotlib.collections

if __name__ == '__main__':

    numPoints = 30
    radius   = 1

    # angle
    th  = np.linspace(0,1,numPoints).reshape(numPoints,1) #np.random.random((numSeeds,1))

    # radius
    rr  = np.linspace(0,1,numPoints).reshape(numPoints,1) #np.random.random((numSeeds,1))

    # set up spiral
    xs = radius/2 + 0.85*radius*rr*np.cos(4.*np.pi*th)
    ys = radius/2 + 0.85*radius*rr*np.sin(4.*np.pi*th)

    seeds = np.hstack((xs,ys))

    # tesselate
    center = np.mean(seeds,axis=0)
    dt     = Tesselate2D(center, 50*radius)

    # insert all seeds one by one
    for s in seeds:
        dt.addPoint(s)

    # build Voronoi diagram as a list of coordinates and regions
    vc,vr = dt.exportVoronoiRegions()

    print(vc[12])
    print(vr[12])

    # visualize with matplotlib.pyplot
    fig,ax = plt.subplots()
    ax.margins(0.1)
    ax.set_aspect('equal')
    plt.axis([-1.3,1.3,-1.3,1.3])

    # plot Voronoi diagram edges
    for r in vr:
        polygon = [vc[i] for i in vr[r]]       # build polygon for each region
        plt.plot(*zip(*polygon),color="green") # plot polygon edges in green
    plt.plot(seeds[:,0],seeds[:,1],'b.')       # plot the center points in blue


    # Dump plot to file
    plt.savefig('voronoi2D.png', dpi=96)
    plt.show()
