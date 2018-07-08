# simple Voronoi demo with matplotlib

import numpy as np
from voronoi2D import *

# visualization
import matplotlib.pyplot as plt
import matplotlib.tri
import matplotlib.collections

if __name__ == '__main__':

    ###########################################################
    # Generate 'numSeeds' seeds
    numSeeds  = 100
    numEdge   = 60

    centSteps = 51
    radius    = 1.0

    # angle
    th  = np.linspace(0,1,numSeeds).reshape(numSeeds,1)
    # radius
    rr  = np.linspace(0,1,numSeeds).reshape(numSeeds,1)

    the = np.linspace(0,2*np.pi,numEdge+1)[0:-1]

    xs = 0.85*radius*rr*np.cos(4.*np.pi*th)
    ys = 0.85*radius*rr*np.sin(4.*np.pi*th)

    xse = radius*np.cos(the)
    yse = radius*np.sin(the)

    xs = np.vstack((xs,np.transpose([xse])))
    ys = np.vstack((ys,np.transpose([yse])))
    seeds = np.hstack((xs,ys))

    # iteration
    for iii in range(centSteps):

        # tesselate
        center = np.mean(seeds,axis=0)
        dt     = Tesselate2D(center, 50*radius)

        # insert all seeds one by one
        for s in seeds:
            dt.addPoint(s)

        # build Voronoi diagram as a list of coordinates and regions
        vc, vr = dt.exportVoronoiRegions()

        if (iii%5 == 0):
            print 'iteration ',iii,' out of ',centSteps
            # visualize with matplotlib.pyplot
            fig,ax = plt.subplots()
            ax.margins(0.1)
            ax.set_aspect('equal')
            plt.axis([-1.3,1.3,-1.3,1.3])

            # plot Voronoi diagram edges (in red)
            for r in vr:
                polygon = [vc[i] for i in vr[r]]       # build polygon for each region
                plt.plot(*zip(*polygon),color="green") # plot polygon edges in green
            plt.plot(seeds[:,0],seeds[:,1],'b.')       # plot the center points in blue
            plt.show()

        nn = len(vr)  # number of polygons
        for kk in range(nn):
            xx = np.asarray([vc[i] for i in vr[kk][:]])[:,0]
            yy = np.asarray([vc[i] for i in vr[kk][:]])[:,1]
            if (abs(xx).max() < 1) and (abs(yy).max() < 1):  # check not on boundary
                seeds[kk, :] = CompCentroid(xx, yy)


    print 'output check = ', seeds[0], seeds[-1]
    # Dump plot to file
    plt.savefig('centroidalvoronoi2D.png', dpi=96)
