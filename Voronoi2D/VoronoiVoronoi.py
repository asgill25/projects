from PIL import Image

import numpy as np
from voronoi2D import *
import matplotlib.pyplot as plt
import matplotlib.tri
import matplotlib.collections
from matplotlib.patches import Polygon
from scipy import signal


def convert2BW(img_in):
    # open the image file
    im_file = Image.open(img_in)

    # convert image to monochrome
    im_BW = im_file.convert('L')

    # convert to numpy array
    im_array = np.array(im_BW)

    return im_array/255. #convert to interval [0,1]


def GradientSquared(A):
    # use central differencing to compute the square of the gradient
    n,m = np.shape(A)
    gradx = np.zeros((n-2,m-2))
    grady = np.zeros((n-2,m-2))

    gradx =  A[0:-2,1:-1] - 2*A[1:-1:,1:-1] + A[2:,1:-1]
    grady =  A[1:-1,0:-2] - 2*A[1:-1:,1:-1] + A[1:-1,2:]

    return gradx**2 + grady**2


def smooth(A,ntimes=1):
    # smooth A by applying moving average filter ntimes
    kernel = np.ones([5,5])/25#np.array([[0., 1., 0.],[0., 0., 0.],[0., 0., 0.]])
    As = A
    for i in range(ntimes):
        As = signal.convolve2d(As, kernel, boundary='symm', mode='same')

    return As


def PlotDiagram(vr, colour, saveName="georgytest.png"):
    # create a Voronoi diagram plot with matplotlib.pyplot
    fig,ax = plt.subplots()
    ax.margins(0.1)
    ax.set_aspect('equal')
    plt.axis([-0.1*Nx,1.1*Nx,-0.1*Ny,1.1*Ny])

    # Plot Voronoi diagram edges (in red)
    for ir,r in enumerate(vr):
        if not ifix[ir] and (ifix[ir] == 0):

            c = colors[ir]

            polygon = [vc[i] for i in vr[r]]           # Build polygon for each region
            # plt.plot(*zip(*polygon),color="green")   # Plot polygon edges in red

            # without polygon boundaries
            ax.add_patch(Polygon(np.transpose(zip(*polygon)), closed=True,fill=True,color=c))

            # with polygon boundaries
            # ax.add_patch(Polygon(np.transpose(zip(*polygon)), closed=True,fill=True,facecolor=c,edgecolor='k'))

    plt.axis('off')
    plt.savefig(saveName, bbox_inches='tight')
    plt.show()
    return



if __name__ == '__main__':

    # import array, flip and transpose for polygon plot
    ImageArray = np.fliplr(np.transpose(convert2BW('GeorgyVoronoy.jpg')))
    Nx,Ny = np.shape(ImageArray)

    im_f = Image.fromarray(ImageArray.astype(np.uint8))

    radius = 4.*max([Nx,Ny])

    # DO NOT MODIFY INPUT PARAMETERS
    numPts    = 400      # seed points
    numPtsEdge = 30    # number of boundary points per edge at
    centSteps =  20    # number of centroid iterations


    # random seed points (with deterministic pseudo random numbers)
    np.random.seed(1) # do not modify the Mersenne Twister seed!
    xs  = Nx*.95*np.random.random((numPts,1)) + .025*Nx;
    np.random.seed(2) # do not modify the Mersenne Twister seed!
    ys  = Ny*.95*np.random.random((numPts,1)) + .025*Ny;


    # edge values
    xse = np.zeros(4*numPtsEdge)
    yse = np.zeros(4*numPtsEdge)
    xse[0:numPtsEdge]  = np.linspace(0.,Nx,numPtsEdge)
    yse[0:numPtsEdge]  = -.05*Ny
    xse[numPtsEdge:2*numPtsEdge] = 1.05*Nx
    yse[numPtsEdge:2*numPtsEdge] = np.linspace(0.,Ny,numPtsEdge)
    xse[2*numPtsEdge:3*numPtsEdge] = np.linspace(0.,Nx,numPtsEdge)
    yse[2*numPtsEdge:3*numPtsEdge] = 1.05*Ny
    xse[3*numPtsEdge:4*numPtsEdge] = -.05*Nx
    yse[3*numPtsEdge:4*numPtsEdge] = np.linspace(0.,Ny,numPtsEdge)


    ne  = len(xse)
    xs  = np.vstack((xs,np.transpose([xse])))
    ys  = np.vstack((ys,np.transpose([yse])))

    # plt.plot(xs,ys,'k.')
    # plt.xlim([-100, Nx+100])
    # plt.ylim([-100, Ny+100])
    # plt.show()

    seeds = np.hstack((xs,ys))
    ifix  = np.zeros(np.shape(seeds)[0],dtype=int)
    ifix[-ne:] = 1

    # density function for weighted version ===================================
    # compute euclidean norm of the gradient of the image
    gradient = np.sqrt(GradientSquared(ImageArray))

    # smooth eulcidean norm of the gradient (a lot)
    gradient = smooth(gradient,15)
    gradient /= np.amax(gradient)

    # high values at boundary
    gradient[:,0]  = 1.
    gradient[:,-1] = 1.
    gradient[0,:]  = 1.
    gradient[-1,:] = 1.

    # plt.imshow(np.transpose(np.fliplr(gradient)), cmap='hot', interpolation='nearest')
    # plt.colorbar()
    # plt.show()

    # function to evaluate gradient
    GradientFunction = \
    lambda x, y: gradient[max([0,min([np.int(x)-1,Nx-3])]),max([0,min([np.int(y)-1,Ny-3])])]

    # function approximating integral over polygon of density gradient by density at centroid times area
    CompCellMass = \
    lambda x_polygon, y_polygon, xc,yc: CompArea(x_polygon,y_polygon)*GradientFunction(xc,yc)


    ifix  = np.zeros(np.shape(seeds)[0],dtype=int)
    ifix[-ne:] = 1

    # centroidalization
    for iii in range(centSteps):

        if (iii%2 == 0):
            print 'iteration ',iii,' out of ',centSteps

        center = np.mean(seeds, axis=0)
        dt     = Tesselate2D(center, 50*radius)

        # insert all seeds one by one
        for s in seeds:
            dt.addPoint(s)

        # Build Voronoi diagram as a list of coordinates and regions
        vc,vr = dt.exportVoronoiRegions()

        nn = len(vr)
        colors = {}
        for kk in range(nn):
            xx = np.asarray([vc[i] for i in vr[kk][:]])[:,0]
            yy = np.asarray([vc[i] for i in vr[kk][:]])[:,1]
            if (ifix[kk] == 0):
                cx, cy = seeds[kk, :] = CompCentroid(xx, yy)
                # Exercise 1 & 2 =============
                # your code
                # ============================
                if iii == centSteps-1:
                    ix = max(min(np.int(cx),Nx-1),0) # continuation at boundary
                    iy = max(min(np.int(cy),Ny-1),0)
                    c = ImageArray[ix,iy]
                    colors[kk] = [c, c, c]


    PlotDiagram(vr, colors)
