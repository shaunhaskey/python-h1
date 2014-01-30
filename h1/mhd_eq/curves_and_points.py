from numpy import sqrt, arcsin, arccos, pi, degrees, sin, cos
from numpy import sort, shape, array, argsort, transpose, arange, linspace
# this is an unsigned distance
from matplotlib.mlab import dist_point_to_segment
import pylab as pl

verbose = 0  # can set verbose with a function
def set_verbose(level=1):
    """ set the verbosity level - 2 and higher generate plots, may need
    to show() to seem them!
    """
    global verbose
    verbose=level

def get_line_abc(p1, p2, normalize=False):
    """ Calculate coefficients for a line through p1 and p2 of the form ax+by+c=0
    Optionally normalize so that a^2+b^2=1
    >>> print(get_line_abc([1,0],[0,1]))
    (1, 1, -1)
    >>> print("%4.4f %4.4f %4.4f" % (get_line_abc([3,0],[0,4], normalize=True)))
    0.8000 0.6000 -2.4000
    """
    a = p2[1] - p1[1]
    b = p1[0] - p2[0]
    c = -p1[0]*(p2[1]-p1[1]) + p1[1]*(p2[0]-p1[0])
    if (normalize):
        n=sqrt(a**2+b**2)
        return (a/n, b/n, c/n)
    else: return (a, b, c)

def signed_dist_point_line(p1, p2, p3):
    """ Return the signed distance from p1 to the line (not segment)
    defined nu p2, p3.  If rotating from p2 to p3 with p1 as the axis
    is +ve (CCW), distance is positive
    >>> print("%6.6f" % signed_dist_point_line([0,0],[1,0],[0,1]))
    0.707107
    >>> print("%6.6f" % signed_dist_point_line([1,1],[2,1],[1,2]))
    0.707107
    """
    abc = get_line_abc(p2, p3, normalize=True)
    return -(abc[0]*p1[0]+abc[1]*p1[1]+abc[2])

def check_curve_segments(curve):
    """ Return the maximum change in angle between successive segments
    in a piecewise linear curve, and a Boolean showing whether the
    rotation sense changes sign.
    This can be the most time consuming step of all - use at higher level
    The curve is 2D - shape(N, 2) so curve[0] -> [x,y] - see examples

    >>> (sinang,flipped) = check_curve_segments([[3,0],[2,1],[1,1],[0,0]])
    >>> print("%4.4f, %s" % (degrees(sinang), flipped))
    45.0000, False
    >>> (sinang,flipped) = check_curve_segments([[3,0],[2,1],[1,1],[1,2]])
    >>> print("%4.4f, %s" % (degrees(sinang), flipped))
    90.0000, True
    """
    sc = shape(curve) 
    err = False
    if len(sc) != 2: err=True
    elif sc[1]!= 2: err=True
    
    if err == 1: raise ValueError, str(
        'curve is shape %s, must be shape (N,2)' % (repr(sc)))

    if verbose>4: 
        pl.figure()
        pl.plot(transpose(curve)[0],transpose(curve)[1],'-+r')
        pl.gca().set_aspect('equal')
    for (i, pt) in enumerate(curve):
        if i == 0: 
            angle = 0
            flipped = False
        else:
            (a,b,c) = get_line_abc(last_point, pt, normalize=True)
            #print(i, a,b,c)
            if i == 1: 
                max_sin = 0.0
            else:
                # sin folds 0-pi/2 with pi/2 0 so cannot detect > 90!
                # should consider cosine as range is 0-pi (or in fact atan2)
                # but cos wont deptect flips!
                sin_angle = last_a * b - last_b * a # do we need to divide?
                #cos_angle = last_a * a + last_b * b # do we need to divide?
                #print a,b,last_a, last_b, sin_angle
                if i == 2: first_sin = sin_angle
                if sin_angle*first_sin < 0: 
                    flipped = True
                    if verbose>4: print('flipped at %d' %i )
                if abs(sin_angle)>max_sin: max_sin = abs(sin_angle)
            (last_a, last_b, last_c) = (a, b, c)
        last_point = pt

    return(arcsin(max_sin), flipped)         

def dist_to_curve(point=None, curve=None, debug=True, eps=None, check_curve=True):
    """ return the signed distance from a point to a piecewise linear
    curve defined as a an array.  The sign of the distance is positive
    if the curve traces a CCW path relative to the point.  

    The check curve option will reject curves with angles between segments
    in excess of 60 degrees (limitation of the sort step in this
    algorithm.
    
    check_curve can be the most time consuming step of all - use at higher level
    ;doesn't work  print("%6.6f" % dist_to_curve([0,0], [[2,0],[1,1],[0,1]]))
    1.000000
    >>> print("%6.6f" % dist_to_curve([1.5,0], [[3,0],[2,1],[1,1],[0,0]]))
    1.000000
    """
    if eps == None: eps=1e-8
    if shape(point) != (2,): raise ValueError, str(
        'point is shape %s, must be shape (2,)' % (repr(shape(point))))    
    sc = shape(curve) 
    err = False
    if len(sc) != 2: err=True
    elif sc[1]!= 2: err=True
    
    if err == 1: raise ValueError, str(
        'curve is shape %s, must be shape (N,2)' % (repr(sc)))

    if check_curve:
        (max_angle, sign_flip) = check_curve_segments(curve)
        if (max_angle > pi/3):# or sign_flip:
            raise ValueError,str(
                'curve rotates too quickly (%2.2f deg)'
                ' or changes rotation sense' % (degrees(max_angle)))

    # order by square of distance 
    distsq = ((array(curve)[:,0] - point[0])**2 + 
              (array(curve)[:,1] - point[1])**2)
    idx = argsort(distsq)
    idxnext = 1
    while (abs(curve[0][0] - curve[idx[idxnext]][0]) +  
           abs(curve[0][1] - curve[idx[idxnext]][1])) < eps:
        idxnext += 1

    #print(point, curve[idx[0]],curve[idx[idxnext]])
    #return(-dist_point_to_segment(point, curve[idx[0]], curve[idx[idxnext]]))
    s_dist = signed_dist_point_line(point, curve[idx[0]],curve[idx[idxnext]])
    if verbose>3: 
        pl.plot(point[0]+abs(s_dist)*cos(linspace(0,2*pi,30)),
                point[1]+abs(s_dist)*sin(linspace(0,2*pi,30)))
        pl.gca().set_aspect('equal')
    if idx[idxnext] > idx[0]: 
        return(s_dist)
    else:
        return(-s_dist)

def rms_distance_points_curve(points, curve, check_curve=True):
    """ return the RMS distance from of the points from the piecewise linear
    curve
    >>> d=rms_distance_points_curve([[-.5,1.5],[1,1]], [[3,0],[2,1],[1,1],[0,0]])
    >>> print("%8.8f" % d)
    1.00000000
    """
    sc = shape(curve) 
    err = False
    if len(sc) != 2: err=True
    elif sc[1]!= 2: err=True
    
    if err == 1: raise ValueError, str(
        'curve is shape %s, must be shape (N,2)' % (repr(sc)))

    sc = shape(points) 
    err = False
    if len(sc) != 2: err=True
    elif sc[1]!= 2: err=True
    
    if err == 1: raise ValueError, str(
        'points is shape %s, must be shape (N,2)' % (repr(sc)))
    if verbose>0: 
        print("Finding RMS distance of %d points to a curve of %d segments" %
              (len(points), len(curve)-1))
    # check curve just once - saves a lot!    
    if check_curve:  # may omit when trying to find a line through islands
        (max_angle,flipped) = check_curve_segments(curve)
        if (max_angle > pi/3):
            raise ValueError,str(
                'curve rotates too quickly (%2.2f deg)'
                ' or changes rotation sense' % (degrees(max_angle)))

    sumsq = 0
    for pt in points:
        s_dist = dist_to_curve(pt, curve,check_curve=False)
        if verbose>2: # debugging info
            if sumsq == 0:  # only plot curve once
                pl.plot(transpose(curve)[0],transpose(curve)[1],'c')
                pl.gca().set_aspect('equal')
            pl.scatter([pt[0]],[pt[1]])    
            pl.plot(pt[0]+abs(s_dist)*cos(linspace(0,2*pi,30)),
                    pt[1]+abs(s_dist)*sin(linspace(0,2*pi,30)))
        sumsq += s_dist**2

    if verbose>2: pl.show()    
    return sqrt(sumsq/len(points))
    

if __name__ == "__main__":
    import doctest
    doctest.testmod()
     
