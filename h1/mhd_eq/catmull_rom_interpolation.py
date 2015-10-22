from numpy import zeros, size, sin, cos, exp, array, dot, sum, diff, mgrid

matcatmull =  0.5*array([[-1,  3, -3,  1],  # 28us
                         [ 2, -5,  4, -1],
                         [-1,  0,  1,  0],
                         [0,   2,  0,  0]])
count=0

power_arr = array([3,2,1,0])

def catmull_int(t, arr, offset=None, mat=matcatmull):
    """ Catmull-Rom cubic spline interp in 1D 
    If t is a vector, assume 0<t<1, and use offset into arr
    If scalar - assume t is a fractional index and calc offset from int(t)
        
    Warning: In both cases, the first element is assumed to have fractional index -1 !!
        This is illustrated in the next line, where the values ARE the indices
        >>> catmull_int(0.5,[-1,0,1,2])
        0.5
        >>> catmull_int(1.5,[1,2,3,2,1])
        2.625
    """

    global count
    count+=1
    #if not hasattr(t, '__len__'):   #  this could be faste - avoid using size
    if size(t) == 1: 
        offset = int(t)
        tt = t - offset
        tvec = tt**power_arr
    elif size(t)==4: tvec = t
    else: raise Exception, 't must have 1 or 4 elements'

    """ Catmull-Rom cubic spline uses finite diff for gradient - from 
    http://en.wikipedia.org/wiki/Cubic_Hermite_spline
    """
    #print type(tvec), type(mat)
    xx = dot(tvec,mat)
    if offset<0:  raise LookupError(' negative index [%d] to arr %d %.4f' %
                                    (offset, t, offset))
    return((xx * arr[offset:offset+4]).sum())

def catmull_int2d(xy, mesh, grids=None):
    """ return the value of a scalar obtained by cubic interpolation into 
    uniform 2d mesh, 
            where xy are fractional indices, (grids=None)
            or if grids=[x2d grid, y2dgrid] x and y values

    Two D case:
    Method is to start with the x value, and interpolate the scalar onto
    the actual x value for four neighbouring mesh y values (two below, two abov).
    Then interpolate in y to get the scalar at the desired x,y
    The first subscript is the x subscript, so that
    makearr2d(5,func=testfnx)[:,0]
    Out[21]: array([-1.,  0.,  1.,  2.,  3.])

    >>> catmull_int2d([1,3],makearr2d(10)) - testfn(1,3)
    0.0
   
    The following creates a grid of coordinate values - shold be easy to get exactly.
    Note that the enumeration starts 1 element inside the array, to make sure all 4 values required are in the array, so i and j need to be adjusted by 1
    >>> gr,gz=mgrid[3:12:10j, -3:2:6j]
    >>> sum(abs(array([[gr[i+1,j+1] - catmull_int2d([r,z],gr,grids=[gr, gz]) for (i,r) in enumerate(gr[1:-3,0])] for (j,z) in enumerate(gz[0,1:-3])])))
    0.0

    # handy for debug: can trace back to x and y  x.0y fn=4.02 -> x=4,y=2
    >>> def fn(x,y): return(x+y/100)  
    >>> gx,gy=mgrid[0:10:21j, 0:5:11j]
    >>> mesh = array([[fn(x,y) for y in gy[0,:]] for x in gx[:,0]])
    >>> sum(abs(array([[fn(x,y) - catmull_int2d([x,y],mesh,grids=[gx, gy]) for x in gx[1:-3,0]] for y in gy[0,1:-3]])))
    0.0

    # test off mesh points - should be true to machine accuracy as fn is linear
    >>> sum(abs(array([[fn(x+0.1,y+0.2) - catmull_int2d([x+0.1,y+0.2],mesh,grids=[gx, gy]) for x in gx[1:-3,0]] for y in gy[0,1:-3]])))  < 1e-13
    True

    Note 1: Needs 5 evaluations of catmull_int (not 4 as stated in wikipedia)   
    Note 2:  Attempts to access data too close to the mesh edge will be detected
    by normal subscript checks for subscripts too large - but we need to check
    explicitly for too small, as a negative subscript in python is legal, and means relative
    to the end of the array.

    """
    inxy = array(xy).copy()  # save for error message if needed
    if grids != None:
        xy = [(xy[0]-grids[0][0,0])/diff(grids[0][0:2,0]) - 1,
              (xy[1]-grids[1][0,0])/diff(grids[1][0,0:2]) - 1]
    xbase = int(xy[0]) ; xt = xy[0] - xbase
    ybase = int(xy[1]) ; yt = xy[1] - ybase

    #for 4 y mesh points, compute scalar at the correct x value
    yint = zeros(4)
    for j in range(0,4):
        if (ybase+j) < 0: 
            raise LookupError,(
                '2D mesh indices - yind= %d for x=%.2g, y=%.2g' % 
                (ybase+j, inxy[0], inxy[1]))
        yint[j] = catmull_int(xy[0], mesh[:,ybase+j])
    val = catmull_int(yt, yint) # yt because the messh is pre-selected
#    x=1/0 for debug
    return(val)

def catmull_int3d(xyz, mesh, grids=None):
    """ return the value of a scalar obtained by cubic interpolation into 
    uniform 3d mesh, where xyz are fractional indices

    Two D case: (above - this is just for explanation)
    Method is to start with the x value, and interpolate the scalar onto
    the actual x value for four neighbouring mesh y values (two below, two abov).
    Then interpolate in y to get the scalar at the desired x,y

    Three D case - use the above to get values at the desired x,y for the 
    4 neighbouring z values, and interpolate to get value at desired x,y,z

    >>> catmull_int3d([1,3,2],makearr3d(10)) - testfn(1,3,2)
    0.0
    >>> abs(catmull_int3d([1.07,3,2],makearr3d(10)) - testfn(1.07,3,2)) < 1e-5
    True
    >>> abs(catmull_int3d([1.07,3,2],makearr3d(10,func=cubicfn)) - cubicfn(1.07,3,2)) < 2e-6
    True

    >>> def fn(x,y,z): return(x+y/100+z/10000)  
    >>> gx,gy,gz=mgrid[0:10:21j, 0:5:11j, -3:3:7j]
    >>> mesh = array([[[fn(x,y,z) for z in gz[0,0,:]] for y in gy[0,:,0]] for x in gx[:,0,0]])
    >>> sum(abs(array([[[fn(x,y,z) - catmull_int3d([x,y,z],mesh,grids=[gx, gy, gz]) for x in gx[1:-3,0,0]] for y in gy[0,1:-3,0]] for z in gz[0,0,1:-3]])))
    0.0
    >>> sum(abs(array([[[fn(x+.1,y+.2,z+.3) - catmull_int3d([x+.1,y+.2,z+.3],mesh,grids=[gx, gy, gz]) for x in gx[1:-3,0,0]] for y in gy[0,1:-3,0]] for z in gz[0,0,1:-3]]))) < 1e-12
    True

    Speed: - i5/760 610us per 3D
       arr=makearr3d(30)
       time for x in arange(0,1,.001): catmull_int3d([x,1,0],arr)

    Note 1: Needs 21 evaluations of catmull_int (not 16 as stated in wikipedia)   
    Note 2:  Attempts to access data too close to the mesh edge will be detected
    by normal subscript checks for subscripts too large - but we need to check
    explicitly for too small, as a negative subscript in python is legal, and means relative
    to the end of the array.

    """
    inxyz = array(xyz).copy()  # save for error message if needed
    if grids != None:
        xyz = [(xyz[0]-grids[0][0,0,0])/diff(grids[0][0:2,0,0]) - 1,
               (xyz[1]-grids[1][0,0,0])/diff(grids[1][0,0:2,0]) - 1,
               (xyz[2]-grids[2][0,0,0])/diff(grids[2][0,0,0:2]) - 1]

    xbase = int(xyz[0]) ; xt = xyz[0] - xbase
    ybase = int(xyz[1]) ; yt = xyz[1] - ybase
    zbase = int(xyz[2]) ; zt = xyz[2] - zbase


    #for 4 z values, compute scalar at the correct y value, each of which
    # were obtained using 4 y values obtained at the correct x value
    zint = zeros(4)
    for k in range(0,4):
        yint = zeros(4)
        #for 4 y values, compute scalar at the correct x value
        for j in range(0,4):
            if ((zbase+k) < 0) or ((ybase+j) < 0): 
                raise LookupError,('negative mesh indices - yind= %d, zind=%d' % 
                                   (ybase+j, zbase+k))
            yint[j] = catmull_int(xyz[0], mesh[:,ybase+j,zbase+k])
        zint[k] = catmull_int(yt, yint)  # yt because the mesh is pre-selected
    val = catmull_int(zt, zint)
#    x=1/0 for debug
    return(val)

def testfn(x,y,z=None,gridsize=None):
    """ A function that cubic interpolation can't do exactly
    will serve for 2D and threeD
    """
    if gridsize==None: gridsize = 0.1
    if z != None:
        return(sin(x*gridsize)*cos(y*gridsize/2)*exp(-z*gridsize/3.))
    else:
        return(sin(x*gridsize)*cos(y*gridsize/2))

def testfnx(x,y,z=None,gridsize=None):
    """ A trivial function that cubic interpolation for debugging
    will serve for 2D and threeD
    """
    if gridsize==None: gridsize = 0.1
    return(x)

def testfny(x,y,z=None,gridsize=None):
    """ A trivial function that cubic interpolation for debugging
    will serve for 2D and threeD
    """
    if gridsize==None: gridsize = 0.1
    return(y)

def cubicfn(x,y,z,gridsize=None):
    """ A function that cubic interpolation should do exactly, but
    catmull-rom is not exact - it approximates the derivative.
    """
    if gridsize==None: gridsize = 0.1
    return(((x*gridsize)**3)*((y*gridsize)**2)*(-(z*gridsize)))

def makearr2d(dim=10, func=testfn, gridsize=None):
    """ prepare test array for catmull_int3d, using function func(x,y)
    - dim is a scalar or a list of three indices
    To simplify testing, the array has extra elements below 0, so that the
    first valid fractional index is 0, and it refers to the second point in 
    the array.  So arr[1,1] is set to fn(0,0)
    """
    if gridsize==None: gridsize = 0.1
    if size(dim) == 1:
        nx = dim ; ny = dim+1
        dim2 = [nx, ny]
    else:
        dim2 = dim
        nx = dim2[0] ;  ny = dim2[1]
    arr = zeros(dim2)

    for i in range(nx):
        for j in range(ny):
            arr[i,j] = func(i-1,j-1,gridsize=gridsize)

    return(arr)

        
def makearr3d(dim=10, func=testfn, gridsize=None):
    """ prepare test array for catmull_int3d, using function func(x,y,z)
    - dim is a scalar or a list of three indices
    To simplify testing, the array has extra elements below 0, so that the
    first valid fractional index is 0, and it refers to the second point in 
    the array.  So arr[1,1,1] is set to fn(0,0,0)
    """
    if gridsize==None: gridsize = 0.1
    if size(dim) == 1:
        nx = dim ; ny = dim+1  ; nz = dim+2
        dim3 = [nx, ny, nz]
    else:
        dim3 = dim
        nx = dim3[0] ;  ny = dim3[1] ; nz = dim3[2] 
    arr = zeros(dim3)

    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                arr[i,j,k] = func(i-1,j-1,k-1,gridsize=gridsize)

    return(arr)

def dotest(n):
    for i in range(n):
        catmull_int3d([20.22,20,20],ss)
        
if __name__ == "__main__":
    import doctest
    doctest.testmod()
    """
    import cProfile
    ss=makearr3d(30,func=cubicfn,gridsize=0.05)
    cProfile.run('dotest(1000)', 'catstats')
    #cProfile.run('dotest(1)', 'catstats')
    import pstats
    p = pstats.Stats('catstats')
    p.sort_stats('time').print_stats()
    """    
"""
This shows order 3 convergence
catmull_int3d([10.11,10,10],makearr3d(15,func=cubicfn,gridsize=0.1))-cubicfn(10.11,10,10,gridsize=.1)
#  -7.6362000000163022e-05
catmull_int3d([20.22,20,20],makearr3d(30,func=cubicfn,gridsize=0.05))-cubicfn(20.22,20,20,gridsize=.05)
#  -1.2012000000227729e-05
catmull_int3d([40.44,40,40],makearr3d(60,func=cubicfn,gridsize=0.025))-cubicfn(40.44,40,40,gridsize=.025)
#  -4.6200000047846856e-07
catmull_int3d([80.88,80,80],makearr3d(120,func=cubicfn,gridsize=0.0125))-cubicfn(80.88,80,80,gridsize=.0125)
# 1.5675000009096607e-07
log(7.6e-5/1.5e-7)/log(8)
# 2.99  (same result log(4.4e-6/8.7e-9)/log(8) for testfn)

It also shows how clumsy the coding is for convergence tests!  The main reason for this was that the 
test function is scaled so that the inputs can be the array indices.

# timing in ipython
ss=makearr3d(30,func=cubicfn,gridsize=0.05)
%timeit catmull_int3d([20.22,20,20],ss)
# h1svr: 560us /point - c.f. 0.7us for Bline ~1000 times slower!

# this way reports a slower time.
from timeit import Timer 
T=Timer('catmull_int3d([20.22,20,20],ss)',setup='from catmull_rom_interpolation import catmull_int3d, makearr3d, cubicfn\nss=makearr3d(30,func=cubicfn,gridsize=0.05)')
T.timeit(1000)
# h1svr: 0.85ms

simple fun: def fn(x,y): x+=y+1+2+3+4;return(x*y)   347ns
def fn(x): return(sin(x))   3.7us
def fn(x): return(sin(linspace(0,x,1))) : 13u
def fn(x): return(sin(linspace(0,x,1000))): 74us
"""
