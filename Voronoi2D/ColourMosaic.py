from PIL import Image
from scipy import misc
from scipy import signal

import numpy as np
from voronoi2D import *
import matplotlib.pyplot as plt
import matplotlib.tri
import matplotlib.collections
from matplotlib.patches import Polygon


def convert2BW(img_in):
    # open the image file
    im_file = Image.open(img_in)

    # convert image to monochrome
    im_BW = im_file.convert('L')

    # convert to num-array
    im_array = np.array(im_BW)

    return im_array/255.

def getRGB(img_in):
    # open the image file
    pixels = misc.imread(img_in)

    r = pixels[:,:,0]/255.
    g = pixels[:,:,1]/255.
    b = pixels[:,:,2]/255.
    return r, g, b

def GradientSquared(A):
    # use central differencing to compute the square of the gradient
    n,m = np.shape(A)
    gradx = np.zeros((n-2,m-2))
    grady = np.zeros((n-2,m-2))

    gradx =  A[0:-2,1:-1] - 2*A[1:-1:,1:-1] + A[2:,1:-1]
    grady =  A[1:-1,0:-2] - 2*A[1:-1:,1:-1] + A[1:-1,2:]

    return gradx**2 + grady**2


def smooth(A, ntimes=1):
    # smooth A using moving average filter
    kernel = np.ones([5,5])/25
    As = A
    for i in range(ntimes):
        As = signal.convolve2d(As, kernel, boundary='symm', mode='same')

    return As


def PlotDiagram(vr, colour,Nx, Ny, saveName="test.png"):
    # create a Voronoi diagram plot with matplotlib.pyplot
    fig,ax = plt.subplots()
    ax.margins(0.1)
    ax.set_aspect('equal')
    plt.axis([-.05*Nx,1.05*Nx,-.05*Ny,1.05*Ny])

    # Plot Voronoi diagram edges (in red)
    for ir,r in enumerate(vr):
        if not ifix[ir] and (ifix[ir] == 0):

            c = colors[ir]

            polygon = [vc[i] for i in vr[r]]           # Build polygon for each region
            # plt.plot(*zip(*polygon),color="green")   # Plot polygon edges in red

            # without polygon boundaries
            #ax.add_patch(Polygon(np.transpose(zip(*polygon)), closed=True,fill=True,color=c))

            # with polygon boundaries
            ax.add_patch(Polygon(np.transpose(zip(*polygon)), closed=True,fill=True,facecolor=c, edgecolor='k'))
    plt.axis('off')
    plt.savefig(saveName, bbox_inches='tight')
    plt.show()
    plt.show()
    return




if __name__ == '__main__':

    # these parameters should work well for many images.
    numPtsX = 6       # number of grid points in x dir to start with
    numPtsY = 6       # number of grid points in y dir to start with
    numPtsEdge = 20   # number of constraining points at boundary
    centSteps = 12    # number of centroid iterations (16 should be ok), including seeding steps
    numSeedSteps = 10  # number of initial iteratsion after which additional points are added
    maxCells = 3000   # maximum number of Voronoi cells

    # filename = 'GeorgyVoronoy.jpg'
    # filename = 'QueensTower.jpg'
    # filename = 'duck.jpg'
    # filename = 'giraffe.jpg'
    # filename = 'infinitywar.png'
    filename = 'qtatnight.jpg'


    # import pixel colour arrays and flip them to make them suitable for matplotlib
    RR, GG, BB = getRGB(filename)
    RR = np.fliplr(np.transpose(RR))
    GG = np.fliplr(np.transpose(GG))
    BB = np.fliplr(np.transpose(BB))

    ImageArray = np.fliplr(np.transpose(convert2BW(filename)))


    Nx,Ny = np.shape(ImageArray)
    radius = 4.*max([Nx,Ny])

    # compute euclidean norm of the gradient of the image
    gradient = np.sqrt(GradientSquared(ImageArray))

    # smooth eulcidean norm of the gradient (a lot)
    gradient = smooth(gradient, 15)  # maybe change this
    gradient /= np.amax(gradient)

    # high values at boundary --> seed more cells there.
    gradient[:,0]  = 1.
    gradient[:,-1] = 1.
    gradient[0,:]  = 1.
    gradient[-1,:] = 1.


    # function to evaluate gradient
    GradientFunction = \
    lambda x, y: gradient[max([0,min([np.int(x)-1,Nx-3])]),max([0,min([np.int(y)-1,Ny-3])])]

    # function approximating integral over polygon of density gradient by density at centroid times area
    CompCellMass = \
    lambda x_polygon, y_polygon, xc, yc: CompArea(x_polygon,y_polygon)*GradientFunction(xc, yc)


    # initialise seeds as grid
    x_tmp = Nx*np.linspace(.1 , 0.9, numPtsX)
    y_tmp = Ny*np.linspace(.1 , 0.9, numPtsY)
    xs,ys = np.meshgrid(x_tmp,y_tmp)
    xs = np.reshape(xs,[numPtsX*numPtsX,1])
    ys = np.reshape(ys,[numPtsX*numPtsX,1])


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

    seeds = np.hstack((xs,ys))

    ifix  = np.zeros(np.shape(seeds)[0],dtype=int)
    ifix[-ne:] = 1


    ifix  = np.zeros(np.shape(seeds)[0],dtype=int)
    ifix[-ne:] = 1

    # centroidalization and addition of points at cell centres
    for iii in range(centSteps):

        print 'iteration ',iii+1,' of ',centSteps

        center = np.mean(seeds,axis=0)
        dt     = Tesselate2D(center, 50*radius)

        # insert all seeds one by one
        for s in seeds:
            dt.addPoint(s)

        # Build Voronoi diagram as a list of coordinates and regions
        vc,vr = dt.exportVoronoiRegions()

        nn = len(vr)

        # cell colouring
        colors = {}

        # cell masses
        mass = np.zeros(nn)

        new_seeds = np.copy(seeds)

        for kk in range(nn):

            xx = np.asarray([vc[i] for i in vr[kk][:]])[:, 0]
            yy = np.asarray([vc[i] for i in vr[kk][:]])[:, 1]

            # begin Ex 3 ======================================================

            if not ifix[kk] or (ifix[kk] == 0):

                cx, cy = CompCentroid(xx, yy)
                mx, my = CompCentreOfMass(xx, yy, GradientFunction)
                # which one when?
                if iii > numSeedSteps-1:
                    seeds[kk, :] = cx, cy
                else:
                    seeds[kk, :] = mx, my  # density weighted CVT
                    new_seeds[kk, :] = cx, cy

                if iii == centSteps-1:
                    # at last step of iteration, save cell colour in colors.
                    ix = max(min(np.int(cx), Nx-1), 0) # continuation at boundary
                    iy = max(min(np.int(cy), Ny-1), 0)
                    r = RR[ix, iy]
                    g = GG[ix, iy]
                    b = BB[ix, iy]
                    colors[kk] = [r, g, b]

                if -1 < np.int(cx) < Nx and -1 < np.int(cy) < Ny:
                    # save cell mass in mass array
                    mass[kk] = CompCellMass(xx, yy, mx, my)

                # =========================================================


        # Ex 3 ===========================================================
        # add additional points (the ones in new_seeds) to seeds
        if np.shape(seeds)[0] < maxCells and iii < numSeedSteps:
            median = np.median(mass[np.nonzero(mass)])
            for m, new_coord in zip(mass, new_seeds):  # if mass is above median
                if m > median:
                    seeds = np.vstack((seeds, new_coord))
                    ifix = np.append(ifix, 0)

        #if iii >= numSeedSteps: #do something different
        #
        # end Ex 3 =========================================================

        print 'current number of cells = ', np.shape(seeds)[0]


    PlotDiagram(vr, colors, Nx, Ny)
