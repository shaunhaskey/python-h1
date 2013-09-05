
"""
Some routines for math, geometry,...
Probably useful...

14.Apr.2013 bhs
* some routines added: lineToFluxSurf(), get_kd(), intersectIn()

10.Aug.2012 bhs
* make use of H1ErrorClasses

29.Jun.2012 bhs
* start 

"""


import numpy

from H1ErrorClasses import *



def pointsOfLine(x0, x1, np):
    """
    Compute n points of a straight line defines by \vec{x0} and \vec{x1}
    in Cartesian coordinates. The points are equally spaced.

    ATTENTION: The points x0 and x1 are part of the solution. 
    Therefore, it must hold: np >= 2.

    INPUT:
      x0, x1 ... Cartesian coordinates of the two end points 
                 defining the line
      np ....... number of points to be computed.

    RETURN:
      parr ..... np x 3 array with coordinates of points
    """

    if(type(np) != int):
        raise GeometryError("Number of points, np, must be of type int!")

    if(np < 2):
        raise GeometryError("Number of points, np, must be greater than two!")

    x0 = numpy.array(x0)
    x1 = numpy.array(x1)

    if((x0==x1).tolist().count(True) == 3):
        raise GeometryError("".join(["x0 and x1 are the same point. ",
                                     "Must differ to define a line!"]))

    parr = numpy.zeros([np,3], numpy.float)
    llen = numpy.sqrt(numpy.sum((x1-x0)**2))
    larr = numpy.linspace(0., llen, np)
    xx = (x1-x0)/llen

    for i,l in enumerate(larr):
        parr[i] = x0 + l*xx

    return(parr)


def cart2Cyl(xyz):
    """
    Coordinate transformation (points) from Cartesian to cylindrical system.

    INPUT:
      xyz ... Cartesian coordinates, array like (n,3)

    RETURN:
      RpZ ... cylindrical coordinates, array like (n,3)
    """
    x = numpy.array(xyz[:,0])
    y = numpy.array(xyz[:,1])

    RpZ = numpy.empty(xyz.shape)
    RpZ[:,0] = numpy.sqrt(x**2 + y**2)
    RpZ[:,1] = numpy.arctan2(y,x) 
    RpZ[:,2] = xyz[:,2].copy()

    return(RpZ)

def cyl2Cart(RpZ):
    """
    Coordinate transformation (points) from cylindrical to Cartesian system.

    INPUT:
      RpZ ... cylindrical coordinates, array like (n,3)

    RETURN:
      xyz ... Cartesian coordinates, array like (n,3)
    """
    R = numpy.array(RpZ[:,0])
    phi = numpy.array(RpZ[:,1])

    xyz = numpy.empty(RpZ.shape)
    xyz[:,0] = R*numpy.cos(phi)
    xyz[:,1] = R*numpy.sin(phi)
    xyz[:,2] = RpZ[:,2].copy()
    return(xyz)


def contained_in(arr, val):

    """ True if val>= min(arr) and val <= max(arr)
    """
    return((val >= min(arr.flatten())) and (val <= max(arr.flatten())))



def nint(xin):
    """
    Compute nearest integer.

    INPUT:
      xin ... input values, array type

    RETURN:
      ret ... nearest integer; numpy array
    """
    xs = numpy.sign(numpy.array(xin)).astype(int)
    xv = numpy.fabs(numpy.array(xin))
    ret = numpy.empty(xv.shape).astype(int)
    l = xv-numpy.floor(xv) < 0.5
    ret[l] = numpy.floor(xv[l]).astype(int)
    ret[~l] = numpy.ceil(xv[~l]).astype(int)
    
    return(ret*xs)



def isLeft(P0, P1, P2):
    """
    Tests if a point is Left|On|Right of an infinite line.

    Taken from:
    http://softsurfer.com/Archive/algorithm_0103/algorithm_0103.htm

    INPUT:  
      P0, P1, P2 . three points 
                   dimension:  Px[2]    (just x,y coordinate)

    RETURN: 
      >0 for P2 left of the line through P0 and P1
      =0 for P2 on the line
      <0 for P2 right of the line


    See: the January 2001 Algorithm "Area of 2D and 3D Triangles and Polygons"

    """
    return ( (P1[0] - P0[0]) * (P2[1] - P0[1])
            - (P2[0] - P0[0]) * (P1[1] - P0[1]) )


def wn_PnPoly(P, V):
    """
    Winding number test for a point in a polygon.
    The polygone is in a plane.

    Taken from:
    http://softsurfer.com/Archive/algorithm_0103/algorithm_0103.htm

    INPUT:   
      P ... point to be checked
            dimension:  P[2]    (just x,y coordinate)
      V ... vertex points of a polygon V[n+1] with V[n]=V[0]
            dimension:  V[n,2]

    RETURN:  
      wn ... the winding number (=0 only if P is outside V[])
    """

    # the winding number counter
    wn = 0
    n = V.shape[0]-1

    # loop through all edges of the polygon
    for i in range(n):                  # edge from V[i] to V[i+1]
        if (V[i,1] <= P[1]):            # start y <= P.y
            if (V[i+1,1] > P[1]):       # an upward crossing
                if (isLeft(V[i], V[i+1], P) > 0):  # P left of edge
                    wn = wn + 1         # have a valid up intersect
                # end   if (isLeft(V[i], V[i+1], P) > 0):
            # end   if (V[i+1,1] > P[1]): 
        else:                           # start y > P.y (no test needed)
            if (V[i+1,1] <= P[1]):      # a downward crossing
                if (isLeft( V[i], V[i+1], P) < 0):  # P right of edge
                    wn = wn - 1         # have a valid down intersect
                # end   if (isLeft( V[i], V[i+1], P) < 0): 
            # end   if (V[i+1,1] <= P[1]):
        # end   if (V[i,1] <= P[1]):
    # end   for i in range(n):
    return(wn)




def lineToFluxSurf(MA, P, CRS):
    """
    Compute intersection of  a straight line defined by the 
    magnetic axis, MA, and a point, P, inside the magnetic 
    surface, CRS, with the last closed magnetic surface, lcms.
    Intersection must be next to P rather than next to MA.

    The cross section data, CRS, are stored in a way, that the 
    intersection should be close to the centre of CRS.

    INPUT:
      MA .. (R,Z) coordinates of magnetic axis
      P ... (R,Z) coordinates of point inside lcms
      CRS . (R,Z) coordinates of points on lcms;
                  type: numpy array n,2

    RETURN:
      RZ ... (R,Z) coordinates of point on flux surface
    """

    MA = numpy.array(MA)
    P = numpy.array(P)
    CRS = numpy.array(CRS)

    (ma_k, ma_d) = get_kd(MA, P)
    if(not numpy.isnan(ma_k)):
        # use y = k*x + d
        ma_f = True
    else:
        # use x = x1
        ma_f = False
    # end   if(ma_k != numpy.nan):

    ind = numpy.arange(len(CRS),0,-1) - 1
    for i in ind:
        #print "----- A ", i
        crs1 = CRS[i]
        crs2 = CRS[i-1]
        (crs_k, crs_d) = get_kd(crs1, crs2)
        #print "crs_k, ma_f ", crs_k, ma_f
        if(not numpy.isnan(crs_k)):
            #print "crs_k != numpy.nan"
            if(ma_f):
                x = (crs_d - ma_d) / (ma_k - crs_k)
            else:
                x = MA[0]
            # end   if(ma_f):
            y = ma_k * x + ma_d
        else:
            #print "crs_k == numpy.nan"
            if(ma_f):
                x = crs1[0]
                if(ma_k != 0.0):
                    y = crs_k * x + crs_d
                else:
                    y = ma_k * x + ma_d
                # end   if(ma_k != 0.0):
                print "lets see..."
#            else:
#                if(MA[0] != crs1[0]):
#                    return(numpy.array([numpy.nan,numpy.nan]))
#                else:
#                    raise Exception, "".join(["Problem: magnetic axis and ",
#                                              "point P are both on a line ",
#                                              "which is defined by the ",
#                                              "considered part of the ",
#                                              "cross section"])
                # end   if(MA[0] != crs1[0]):
            # end   if(ma_f):
        # end   if(crs_k != numpy.nan):
        # now check whether the point of intersection is the right one
        RZ = numpy.array([x,y])
        #print "ma_k, ma_d, crs_k, crs_d, ma_f ", ma_k, ma_d, crs_k, crs_d, ma_f
        #print "crs1, crs2, RZ ", crs1, crs2, RZ
        isOnLine1 = intersectIn(crs1, crs2, RZ)
        isOnLine2 = intersectIn(MA, RZ, P)
        #print "isOnLine1, isOnLine2 ", isOnLine1, isOnLine2
        if(isOnLine1 and isOnLine2):
            return(RZ)
        # end   if(isOnLine1 and isOnLine2):
        #print "----- B"
    # end   for p in P:

    print "bad luck, no start value for inversion found - re-consider coding"

    return(numpy.array([numpy.nan,numpy.nan]))


def get_kd(X1, X2):
    """
    Starting with two points, X1 and X2, compute slope, k, and
    offset, d, for the equation y = k*x + d.
    In case X1[0] == X2[0] use the form x = x1 and return 
    k = numpy.nan and d = X1[0].

    INPUT:
      X1, X2 ... two points of the straight line (R,Z)

    RETURN:
      k, d ... for y = k*x + d  or 
               k = numpy.nan and d = X1[0]   for x = x1
    """
    X1 = numpy.array(X1)
    X2 = numpy.array(X2)
    dd = X1 - X2
    if(dd[0] == 0.0):
        k = numpy.nan
        d = X1[0]
    else:
        k = dd[1]/dd[0]
        d = X1[1] - k*X1[0]
    # end   if(d[0] == 0.0):

    return(k, d)


def intersectIn(X1, X2, P):
    """
    Check whether point, P=(x,y), is on line within X1 and X2.

    INPUT:
      X1, X2 ... points X1=(x1,y1) and X2=(x2,y2) defining a line
      P ........ point P=(x,y) 

    RETURN:
      isOnLine ... True if P is between X1 and X2
    """
    acc = 0.002

    Xl2 = X2 - X1
    Pl = P - X1
    #print "Xl2, Pl ", Xl2, Pl

    if(numpy.signbit(Pl[0]) != numpy.signbit(Xl2[0])):
        #print '1'
        return(False)
    # end   if(...

    if(numpy.signbit(Pl[1]) != numpy.signbit(Xl2[1])):
        #print '2'
        return(False)
    # end   if(...

    if(not numpy.signbit(Pl[0]) and (Pl[0] > Xl2[0])
       and (abs(Pl[0] - Xl2[0]) > acc)):
        #print '3'
        return(False)
    # end   if(...

    if(not numpy.signbit(Pl[1]) and (Pl[1] > Xl2[1])
       and (abs(Pl[1] - Xl2[1]) > acc)):
        #print '4'
        return(False)
    # end   if(...

    if(numpy.signbit(Pl[0]) and (Pl[0] < Xl2[0])
       and (abs(Pl[0] - Xl2[0]) > acc)):
        #print '5'
        return(False)
    # end   if(...

    if(numpy.signbit(Pl[1]) and (Pl[1] < Xl2[1])
       and (abs(Pl[1] - Xl2[1]) > acc)):
        #print '6'
        return(False)
    # end   if(...

    return(True)


def rotCoords(c_in, angle):
    """
    Rotate point coordinates clockwise by specified angle.

    INPUT:
      c_in ... point coordinates; type: numpy array; shape: n,2
      angle .. rotation angle in rad

    RETURN:
      coords . points coordinates after rotation; 
               type: numpy array; shape: n,2
    """
    angle = -angle
    c_in = numpy.matrix(c_in)
    coords = numpy.empty(c_in.shape)
    m = numpy.matrix([[numpy.cos(angle), -numpy.sin(angle)],
                      [numpy.sin(angle),  numpy.cos(angle)]])
    for i,c in enumerate(c_in):
        coords[i,:] = numpy.array(m*c.transpose()).transpose()[0]
    # end   for i,c in enumerate(c_in):
    return(coords)
