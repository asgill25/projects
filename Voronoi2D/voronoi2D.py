"""
Delaunay and Voronoi tesselation in two dimensions using the Bowyer-Watson algorithm
"""

import numpy as np
from math import sqrt


def CompCentroid(x,y):
    M = len(x)
    if M != len(y):
        return None, None

    c_x = 0.
    c_y = 0.

    for j in range(0, M):
        c_x += (x[j] + x[(j+1) % M])*(x[j]*y[(j+1) % M] - x[(j+1) % M]*y[j])
        c_y += (y[j] + y[(j+1) % M]) * (x[j] * y[(j+1) % M] - x[(j+1) % M] * y[j])

    A = CompArea(x, y)

    c_x = c_x/(6. * A)
    c_y = c_y / (6. * A)

    return c_x, c_y

def CompArea(x, y):
    M = len(x)
    A = 0.

    for j in range(0, M):
        A += (x[j] * y[(j+1) % M] - x[(j+1) % M] * y[j])

    A = A/2.

    return A


def CompCentreOfMass(x, y, massfunction):
    M = len(x)
    if M != len(y):
        return None, None

    mx = 0.
    my = 0.
    denom = 0.

    for i in range(0, M):
        denom += massfunction(x[i], y[i])
        mx += massfunction(x[i], y[i]) * x[i]
        my += massfunction(x[i], y[i]) * y[i]

    mx = mx/denom
    my = my/denom

    return mx, my


class Tesselate2D:
    # class to tesselate the 2D plane given a set of points

    def __init__(self,ctr=(0,0),r=9999):
        # initialize the domain

        ctr = np.asarray(ctr)
        # create coordinates for the bounding square
        self.coords = [ctr + r*np.array((-1,-1)),
                       ctr + r*np.array((+1,-1)),
                       ctr + r*np.array((+1,+1)),
                       ctr + r*np.array((-1,+1))]

        # create two dictionaries to store triangle neighbors and circumcircles
        self.triangles = {}      # corner indices of (up to) three neighboring triangles
        self.circles   = {}      # coordinates of center and radius

        # create two counter-clockwise triangles for the frame
        T1 = (0,1,3)
        T2 = (2,3,1)
        self.triangles[T1] = [T2,None,None]
        self.triangles[T2] = [T1,None,None]

        # compute center coordinates and radius of circumcircle for each triangle
        for t in self.triangles:
            self.circles[t] = self.circumcenter(t)

    def circumcenter(self,tri):
        # compute the circumcenter and circumradius of a triangle

        pts  = np.asarray([self.coords[v] for v in tri])
        pts2 = np.dot(pts,pts.T)
        A    = np.bmat([[2*pts2, [[1],
                                  [1],
                                  [1]]],
                       [[[1, 1, 1, 0]]]])

        #print pts
        #print pts2
        #print A
        #aaa = raw_input('stopped here')

        b           = np.hstack((np.sum(pts*pts,axis=1), [1]))
        x           = np.linalg.solve(A,b)
        bary_coords = x[:-1]
        center      = np.dot(bary_coords,pts)

        # radius = np.linalg.norm(pts[0] - center) # euclidean distance
        radius = np.sum(np.square(pts[0] - center))  # squared distance
        return (center,radius)

    def inCircleFast(self,tri,p):
        # check if point p is inside of precomputed circumcircle of tri

        center,radius = self.circles[tri]
        return np.sum(np.square(center - p)) <= radius

    def inCircleRobust(self,tri,p):
        # check if point p is inside of precomputed circumcircle of tri
        # robust (slower) version
        # ref: http://www.cs.cmu.edu/~quake/robust.html

        m1 = np.asarray([self.coords[v]-p for v in tri])
        m2 = np.sum(np.square(m1), axis=1).reshape((3,1))
        m  = np.hstack((m1,m2))    # The 3x3 matrix to check
        return np.linalg.det(m) <= 0

    def addPoint(self,p):
        # add a point to the current tesselation, refine by Bowyer-Watson

        p   = np.asarray(p)
        idx = len(self.coords)
        # print("coords[", idx,"] ->",p)
        self.coords.append(p)

        # search the triangle(s) whose circumcircle contains p
        bad_triangles = []

        # THIS NEXT LOOP SCANS OVER >>> ALL <<< TRIANGLES !! (SLOW)
        for T in self.triangles:
            # choose one method: inCircleRobust(T, p) or inCircleFast(T, p)
            if self.inCircleFast(T,p):
            #if self.inCircleRobust(T,p):
                bad_triangles.append(T)

        # find the counter-clockwise boundary (star shape) of the
        # bad triangles, expressed as a list of edges (point pairs)
        # and the opposite triangle to each edge
        boundary = []
        # choose a "random" triangle and edge
        if bad_triangles:
            T        = bad_triangles[0]
            edge     = 0
            # get the opposite triangle of this edge
            while True:
                # check if edge of triangle T is on the boundary...
                # if opposite triangle of this edge is external to the list
                tri_op = self.triangles[T][edge]
                if tri_op not in bad_triangles:
                    # insert edge and external triangle into boundary list
                    boundary.append((T[(edge+1)%3],T[(edge-1)%3],tri_op))

                    # move to next CCW edge in this triangle
                    edge = (edge+1)%3

                    # check if boundary is a closed loop
                    if boundary[0][0] == boundary[-1][1]:
                        break
                else:
                    # move to next CCW edge in opposite triangle
                    edge = (self.triangles[tri_op].index(T)+1)%3
                    T    = tri_op

            # remove triangles too near of point p of our solution
            for T in bad_triangles:
                del self.triangles[T]
                del self.circles[T]

        # retriangulate the hole left by bad_triangles
        new_triangles = []
        for (e0,e1,tri_op) in boundary:
            # create a new triangle using point p and edge extremes
            T = (idx,e0,e1)

            # store circumcenter and circumradius of the triangle
            self.circles[T] = self.circumcenter(T)

            # set opposite triangle of the edge as neighbour of T
            self.triangles[T] = [tri_op,None,None]

            # try to set T as neighbor of the opposite triangle
            if tri_op:
                # search the neighbor of tri_op that use edge (e1,e0)
                for i,neigh in enumerate(self.triangles[tri_op]):
                    if neigh:
                        if e1 in neigh and e0 in neigh:
                            # change link to use our new triangle
                            self.triangles[tri_op][i] = T

            # add triangle to a temporal list
            new_triangles.append(T)

        # link the new triangles each another
        N = len(new_triangles)
        for i, T in enumerate(new_triangles):
            self.triangles[T][1] = new_triangles[(i+1)%N]   # next
            self.triangles[T][2] = new_triangles[(i-1)%N]   # previous

    def exportVoronoiRegions(self):
        # export coordinates and regions of Voronoi diagram as indexed data

        # remember to compute circumcircles if not done before
        # for t in self.triangles:
        #     self.circles[t] = self.circumcenter(t)

        useVertex = {i: [] for i in range(len(self.coords))}
        vor_coors = []
        index     = {}
        # build a list of coordinates and an index per triangle/region
        for tidx, (a,b,c) in enumerate(self.triangles):
            vor_coors.append(self.circles[(a,b,c)][0])
            # insert triangle, rotating it so the key is the "last" vertex
            useVertex[a] += [(b,c,a)]
            useVertex[b] += [(c,a,b)]
            useVertex[c] += [(a,b,c)]
            # set tidx as the index to use with this triangles
            index[(a,b,c)] = tidx
            index[(c,a,b)] = tidx
            index[(b,c,a)] = tidx

        # init regions per coordinate dictionary
        regions = {}
        # sort each region in a coherent order, and substitute each triangle
        # by its index
        for i in range(4,len(self.coords)):
            v = useVertex[i][0][0]  # get a vertex of a triangle
            r = []
            for _ in range(len(useVertex[i])):
                # search the triangle beginning with vertex v
                t = [t for t in useVertex[i] if t[0] == v][0]
                r.append(index[t])  # add the index of this triangle to region
                v = t[1]            # choose the next vertex to search
            regions[i-4] = r        # store region

        return vor_coors, regions

    def exportTriangles(self):
        """Export the current list of Delaunay triangles
        """
        # Filter out triangles with any vertex in the extended BBox
        return [(a-4,b-4,c-4)
                for (a,b,c) in self.triangles if a > 3 and b > 3 and c > 3]
