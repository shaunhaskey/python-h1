
__doc__ = \
"""
for history in netCDF4Tools:

get full user name:
pwd.getpwuid(os.getuid()).pw_gecos
get login name:
getpass.getuser()


Class for 3D grid and it's inverse grid.


History:

08.Aug.2012 bhs
* use H1ErrorClasses

21.Jun.2012 bhs
* start with this script

"""

__version__ = "0.1"
__versionTime__ = "21 June 2012 21:30"
__author__ = "Bernhard Seiwald <bernhard.seiwald@gmail.com>"


import time
import numpy
import random
from scipy import interpolate
import netCDF4Tools
from H1ErrorClasses import *
import miscmath


#WRITE_PROGRESS=True
WRITE_PROGRESS=False


class GRID3D():
    """
    """

    def __init__(self, 
                 nc_filename=None, 
                 compute_grid=True,
                 load_grid=True, 
                 save_grid=False,
                 gridfunction=None, 
                 gnr=None, gnp=None, gnz=None,
                 testmode = False):
        """
                 nc_filename=None, 
                 compute_grid=True,
                 load_grid=True, 
                 save_grid=False,
                 gridfunction=None, 
                 gnr=None, gnp=None, gnz=None,
                 testmode = False

        """
        self.isAvail = False

        if(load_grid and (nc_filename != None)):
            try:
                self.load(nc_filename)
                self.isAvail = True
                return
            except:
                pass
            # end   try:
        # end   if(load_grid and (nc_filename != None)):

        if(testmode):
            gridfunction = gridTestFunc
        # end   if(testmode):

        if(gridfunction == None):
            raise Grid3DError('No function for computing gridvalues provided!')
        # end   if(gridfunction == None):

        if(compute_grid):
            self.computeGrid(gridfunction, nr=gnr, np=gnp, nz=gnz)
            self.computeInverseGrid()
            self.isAvail = True
        # end   if(compute_grid):

        if(save_grid):
            self.save(nc_filename)
        # end   if(save_grid):

        return


    def computeGrid(self, gridfunction, nr=None, np=None, nz=None):
        """
        Compute a 3D grid and the inverse 3D grid of the torus.

        Loop over (s,theta_b,phi_b) and compute
        R(s,theta_b,phi_b), phi(s,theta_b,phi_b), Z(s,theta_b,phi_b)
        in order to compute the 3D grid.

        In the second step invert the grid: compute indices for R,phi,Z
        and store s(R,phi,Z), theta_b(R,phi,Z), phi_b(R,phi,Z). This
        inverse grid (igrid) might be usefull to determine a guess for 
        inverting coordinates.


        ATTENTION: As many evaluations of the gridfunction 
           (nr * np * nz) are necessary to - hopefully - fill
           the grid without holes, the computation might be 
           extremly time consuming!!!

        Call e.g.:
        boo = BOOZER.BOOZER(ncfile)
        compute(boo.fct_RPZ, nr=20, np=90, nz=20)

        INPUT:
          gridfunction ... function to compute values for inverse grid
          nr ..... number of r values; cylindrical coordinates
          np ..... number of phi values; cylindrical coordinates
          nz ..... number of z values; cylindrical coordinates

        RETURN:
          None
        """
        # vector with nan
        nan_vec = numpy.array([numpy.nan,numpy.nan,numpy.nan])

        # some 'magic' numbers for fallback
        if(nr == None):   nr=10  #8
        if(np == None):   np=60  #12
        if(nz == None):   nz=12  #8

        # use same grid dimensions in Boozer coordinates
        ns_b = 1*nr
        nt_b = 1*nz
        np_b = 1*np

        self.nr = nr
        self.np = np
        self.nz = nz
        self.ns_b = ns_b
        self.nt_b = nt_b
        self.np_b = np_b

        print "".join(["computing 3D grid with \n",
                       "   ns_b={0:d}, nt_b={1:d}, np_b={2:d} "
                       .format(self.ns_b, self.nt_b, self.np_b),
                       "(={0:d}) gridpoints on Boozer grid...\n"
                       .format(self.ns_b*self.nt_b*self.np_b),
                       "may take a while until finished..."])

        self.s_b = numpy.linspace(0., 1., self.ns_b)
        self.t_b = numpy.linspace(0., 2*numpy.pi, self.nt_b)
        self.p_b = numpy.linspace(0., 2*numpy.pi, self.np_b)
        self.grid3D = numpy.zeros([self.ns_b,self.nt_b,self.np_b,3])
        self.grid3D[:,:,:] = nan_vec

        if(WRITE_PROGRESS):   print "computing real space grid"
        t1 = time.time()
        # this loop i shall parallelize!!!
        # import multiprocessing
        # see mp_sample3.py
        # also VMECFBScan.py
        # http://www.ibm.com/developerworks/aix/library/au-multiprocessing/
        # BUT: some work around needed according:
        #   http://stackoverflow.com/questions/1816958/cant-pickle-type-instancemethod-when-using-pythons-multiprocessing-pool-ma
        for (i,s) in enumerate(self.s_b):
            for (j,t) in enumerate(self.t_b):
                for (k,p) in enumerate(self.p_b):
                    if(WRITE_PROGRESS):  print "i(s), j(t), k(p) ", i, j, k
                    self.grid3D[i,j,k] = numpy.array(gridfunction(s,t,p,calc_derivs=False))

        t2 = time.time()
        print "finished after {0:.1f}s".format(t2-t1)

        self.r_min = self.grid3D[:,:,:,0].min()
        self.r_max = self.grid3D[:,:,:,0].max()
        self.dr = (self.r_max-self.r_min)/max((self.nr-1),1)
        self.p_min = 0.            #grid3D[:,:,:,1].min()
        self.p_max = 2*numpy.pi    #grid3D[:,:,:,1].max()
        self.dp = (self.p_max-self.p_min)/max((self.np-1),1)
        self.z_min = self.grid3D[:,:,:,2].min()
        self.z_max = self.grid3D[:,:,:,2].max()
        self.dz = (self.z_max-self.z_min)/max((self.nz-1),1)

        return


    def computeInverseGrid(self):
        """
        """
        # vector with nan
        nan_vec = numpy.array([numpy.nan,numpy.nan,numpy.nan])

        self.igrid3D = numpy.zeros([self.nr,self.np,self.nz,3])-1
        self.igrid3D[:,:,:] = nan_vec

        if(WRITE_PROGRESS):
            fd = open('tg1.txt', 'w')

        ig = numpy.zeros(self.igrid3D.shape, int)
        for (i,s) in enumerate(self.s_b):
            for (j,t) in enumerate(self.t_b):
                for (k,p) in enumerate(self.p_b):
                    ir = int((self.grid3D[i,j,k,0]-self.r_min)/self.dr)
                    ip = int((self.grid3D[i,j,k,1]-self.p_min)/self.dp)
                    iz = int((self.grid3D[i,j,k,2]-self.z_min)/self.dz)
                    if(WRITE_PROGRESS):
                        if(self.igrid3D[ir,ip,iz,0] == -1):
                            print "grid3D[",i,j,k,"] unused..."

                    self.igrid3D[ir,ip,iz] = numpy.array([s,t,p])
                    #ig[i,j,k,0] = ir
                    #ig[i,j,k,1] = ip
                    #ig[i,j,k,2] = iz
                    ig[ir,ip,iz,0] += 1
                    ig[ir,ip,iz,1] += 1
                    ig[ir,ip,iz,2] += 1
                    if(WRITE_PROGRESS):
                        if(ig[ir,ip,iz,0] > 1):
                            tstr = "".join(["cell {0:d},{1:d},{2:d}".format(ir,ip,iz), 
                                            " R: {0:.3f},{1:.3f},{2:.3f}".format(self.grid3D[i,j,k,0],self.grid3D[i,j,k,1],self.grid3D[i,j,k,2]), 
                                            " S: {0:.3f},{1:.3f},{2:.3f}".format(s,t,p)])
                            print tstr
                            fd.write("".join([tstr,'\n']))

        if(WRITE_PROGRESS):
            fd.close()
        print "...ready. Keep fingers crossed that there are no holes..."
        d3shape = self.igrid3D.shape
        nel = numpy.prod(d3shape[0:3])
        ind = numpy.isnan(self.igrid3D.reshape([nel,3])[:,0]).reshape(d3shape[0:3])
        holes_ind = numpy.array(numpy.where(ind==True)).T
        try:
            n = holes_ind.size
            holesFound = True
        except:
            holesFound = False

        if(holesFound):
            print "{0:d} holes".format(n)
            if(WRITE_PROGRESS):   print "found at: ", holes_ind
            # try to fill holes
            ret = self.fixHoles(holes_ind)
        else:
            print "no holes found"

        self.is_computed = True
        if(WRITE_PROGRESS):
            return(ig)
        else:
            return

    def fixHoles(self, h_ind):
        """
        This routine tries to fix/fill all the holes in the 
        inverse grid igrid3D.

        Calls fixOneHole() to do the work.

        Method:
        Use ordinary qubic spline for splining data along an axes of
        igrid3D. Use spline evaluation for approximating a value.

        INPUT:
          h_ind ... array with indices of the holes

        RETURN:
          ret ..... success: True  otherwise False
        """
        ret = True
        igrid3Dshape = self.igrid3D.shape
        for indx in h_ind[:]:
            try:
                self.fixOneHole(indx)
            except:
                pass
        return(ret)

    def fixOneHole(self, indx):
        """
        This routine tries to fix/fill one hole in the 
        inverse grid igrid3D.

        Method:
        Use ordinary qubic spline for splining data along an axes of
        igrid3D. Use spline evaluation for approximating a value.


        INPUT:
          indx ... triplet defining the index of the hole in igrid3D

        RETURN:
          None
        """
        for j in range(3):  # for all three co-ordinates (R,phi,Z)
            try:
                x = numpy.where(~numpy.isnan(self.igrid3D[:,indx[1],indx[2],j]))
                y = self.igrid3D[indx[1],indx[2],x,j]
                tck = interpolate.splrep(x,y)
                xin = h_ind[2]
                self.igrid3D[indx[0],indx[1],indx[2],j] = interpolate.splev(xin, tck, der=0, ext=0)
            except:
                if(WRITE_PROGRESS):  
                    estr = "".join(["fixHoles(): hole for ",
                                    "igrid3D{0:s} [{1:d}] ".format(indx, j),
                                    "can't be fixed."])
                    raise Grid3DError(estr)
                raise

        if(WRITE_PROGRESS):  
            print "".join(["fixHoles(): hole for ",
                           "igrid3D{0:s} [{1:d}] ".format(indx, j),
                           "fixed."])
        return


    def stGuessRPZ(self, R, phi, Z):
        """
        Determine a proper start value for the 3d Real2Boozer()

        Method: The torus is mapped to a 3d grid. 

        INPUT:
          R, phi, Z ... cylindrical coordinates
          s, theta_b, phi_b ... Boozer coordinates
        RETURN:
          
        """
        if(not self.is_computed):
            print "3d grid not computed. Computing now. This takes a while..."
            self.compute()

        ir = miscmath.nint([(R   - self.r_min) / self.dr])[0]
        ip = miscmath.nint([(phi - self.p_min) / self.dp])[0]
        iz = miscmath.nint([(Z   - self.z_min) / self.dz])[0]

        (ir,ip,iz, success) = self.nextNeighbour(ir,ip,iz)

        (s, theta_b, phi_b) = self.igrid3D[ir,ip,iz]

        return(numpy.array([s, theta_b, phi_b]))


    def nextNeighbour(self, ir, ip, iz):
        """
        Find next cell in the inverse grid, igrid3D, containing 
        not a nan-tripple.

        Ok, it does not really search the next neighbour, because 
        this routine is somewhat simplified.

        Method:
        Starting at the given point, walk along the axis until a
        good cell is found.

        INPUT:
           ir, ip, iz ... index where to start search for good neighbour
           
        RETURN:
           ir, ip, iz ... index with good neighbour
           success ....... True for a good cell otherwise False 
        """

        # check whether there is at least one good cell in igrid3D
        if(numpy.array(numpy.where(~numpy.isnan(self.igrid3D))).size == 0):
            if(WRITE_PROGRESS):   print 'nextNeighbour() 1: ir, ip, iz ', ir, ip, iz
            return(ir, ip, iz, False)

        # check for boundaries
        # for ir and iz stay within the limits
        """
        if(ir not in range(self.igrid3D.shape[0])):
            if(WRITE_PROGRESS):   
                print 'nextNeighbour() 2: ir, ip, iz ', ir, ip, iz
            return(ir, ip, iz, False)
        if(iz not in range(self.igrid3D.shape[2])):
            if(WRITE_PROGRESS):   
                print 'nextNeighbour() 3: ir, ip, iz ', ir, ip, iz
            return(ir, ip, iz, False)
        """
        if(ir < 0):
            ir = numpy.max([ir, 0])
        if(ir >= self.igrid3D.shape[0]):
            ir = numpy.min([ir, self.igrid3D.shape[0]-1])
        if(iz < 0):
            iz = numpy.max([iz, 0])
        if(ir >= self.igrid3D.shape[2]):
            iz = numpy.min([iz, self.igrid3D.shape[2]-1])

        # for phi we allow for walking around the torus
        if(ip < 0):   
            ip += self.igrid3D.shape[1]
        if(ip >= self.igrid3D.shape[1]):   
            ip -= self.igrid3D.shape[1]

        # let's see whether we found a good cell
        if(not self.igrid3DIsNan(ir, ip, iz)):
            if(WRITE_PROGRESS):   print 'nextNeighbour() 4: ir, ip, iz ', ir, ip, iz
            return(ir, ip, iz, True)

        # In order to reduce biased search direction, change the
        # order of directions for walking through the grid randomly.
        # In other words, do a random walk.
        # number of main directions (axes) for walk
        na = 3
        ser = range(2*na)
        random.shuffle(ser)
        d = numpy.zeros([2*na,na], int)
        for i in range(0, 2*na, 2):
            i0 = ser[i]
            i1 = ser[i+1]
            d[i0,i/2] = 1
            d[i1,i/2] = -1

        success = False

        for i in range(2*na):
            (ir1, ip1, iz1) = (ir+d[i,0], ip+d[i,1], iz+d[i,2])
            (ir2, ip2, iz2, success) = self.nextNeighbour(ir1, ip1, iz1)
            if(success):
                if(WRITE_PROGRESS):
                    print "nextNeighbour()  A{0:d}: ir1, ip1, iz1: {1:d},{2:d},{3:d}".format(i, ir1, ip1, iz1)
                return(ir2, ip2, iz2, success)

        if(WRITE_PROGRESS):   print 'nextNeighbour()  11: ir, ip, iz ', ir, ip, iz
        return(ir2, ip2, iz2, success)

    def igrid3DIsNan(self, ir, ip, iz):
        """
        """
        if(numpy.isnan(self.igrid3D[ir,ip,iz]).tolist().count(True) != 0):
            return(True)
        else:
            return(False)


    def load(self, nc_filename):
        """
        load 3D grid

        INPUT:
          nc_filename ... name of netCDF-4 file

        RETURN:
          None
        """
        try:
            print "try loading ...grid3D: {0:s}".format(nc_filename)
            gr = netCDF4Tools.readNetCDF4(nc_filename)
                #self.is_computed = True
            print "...grid3D loaded"
        except:
            raise Grid3DError("can't load grid3D")
            
        # import all
        for (key, value) in gr.iteritems():
            cmdstr = "self.{0:s} = value".format(key)
            if(WRITE_PROGRESS):  print cmdstr.replace('key', key)
            exec(cmdstr) 

        return


    def save(self, nc_filename):
        """
        save the grid

        INPUT:
          nc_filename ... name of netCDF-4 output file

        RETURN:
          None
        """
        grid3Dstruct = {}
        grid3Dstruct.update({"nr":self.nr})
        grid3Dstruct.update({"np":self.np})
        grid3Dstruct.update({"nz":self.nz})
        grid3Dstruct.update({"ns_b":self.ns_b})
        grid3Dstruct.update({"nt_b":self.nt_b})
        grid3Dstruct.update({"np_b":self.np_b})
        grid3Dstruct.update({"s_b":self.s_b})
        grid3Dstruct.update({"t_b":self.t_b})
        grid3Dstruct.update({"p_b":self.p_b})
        grid3Dstruct.update({"is_computed":self.is_computed})
        grid3Dstruct.update({"r_min":self.r_min})
        grid3Dstruct.update({"r_max":self.r_max})
        grid3Dstruct.update({"dr":self.dr})
        grid3Dstruct.update({"p_min":self.p_min})
        grid3Dstruct.update({"p_max":self.p_max})
        grid3Dstruct.update({"dp":self.dp})
        grid3Dstruct.update({"z_min":self.z_min})
        grid3Dstruct.update({"z_max":self.z_max})
        grid3Dstruct.update({"dz":self.dz})
        grid3Dstruct.update({"grid3D":self.grid3D})
        grid3Dstruct.update({"igrid3D":self.igrid3D})

        netCDF4Tools.writeNetCDF4(grid3Dstruct, nc_filename)
        return




def gridTestFunc(s, t, p, calc_derivs=False):
    """
    This is a very simple function to create a torus.
    Used e.g. for testing GRID3D.

ns=4
nt=10
np=15
sarr=numpy.linspace(0,1,ns)
tarr=numpy.linspace(0, 2*numpy.pi, nt)
parr=numpy.linspace(0, 2*numpy.pi, np)
x=[]
y=[]
z=[]
for s in sarr:
  for t in tarr:
    for p in parr:
       (R,phi,Z)=GRID3D.gridTestFunc(s, t, p)
       x.append(R*numpy.cos(phi))
       y.append(R*numpy.sin(phi))
       z.append(Z)

mi = min(min(x), min(y), min(z))
ma = max(max(x), max(y), max(z))
fig = pylab.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.set_xlim(mi,ma)
ax.set_ylim(mi,ma)
ax.set_zlim(mi,ma)
ax.plot3D(x,y,z)
fig.show()
    """
    R0 = 1.1
    R1 = numpy.sqrt(s) * 0.2
    Z1 = numpy.sqrt(s) * 0.3
    a = (2.*numpy.pi)/3
    m=1.
    n=1.
    arg = m*t + n*p
    R = R0 + R1*numpy.cos(arg)
    phi = p #- a*numpy.sin(arg)
    Z = Z1*numpy.sin(arg)
    return(R,phi,Z)
