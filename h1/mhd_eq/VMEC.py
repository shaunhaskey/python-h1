#! /usr/bin/env python

__doc__ = \
"""
x=numpy.linspace(0,5,10)
y=numpy.linspace(-2,2,9)
gx,gy=numpy.meshgrid(x,y)
z=gy*numpy.sin(gx)
grida=griddata((gx.flatten(),gy.flatten()), z.flatten(),(grid_r,grid_z))
pylab.contourf(grid_r, grid_z, grida, hold=True)
pylab.contour(x,y,z,50)
pylab.show()

A class for VMEC equilibria with some useful methods.

History:

18.Jan.2013 bhs
* removed: calc_rz_coefs2(), calc_rzl_coefs2(), getCrossSectDataOLD()

19.Dec.2012 bhs
* mayavi may be used for plotting fluxsurfaces
* mod B may be plotted on the fluxsurface

13.Dec.2012 bhs
* getFluxsurfData()  meshgrid output is now available (optional)

21.Nov.2012 bhs
* introduced method fct_RPZ(...)
* self.grid3D = GRID3D.GRID3D(...)
both, in order to be compatible with the BOOZER class

08.Aug.2012 bhs
* make use of H1ErrorClasses
* make use of SPLINE
* make use of spline_spec
* implement getCrossSectData similar as in BOOZER

16 Mar.2012
* Make importAll() an own method.

06 Sep.2011 bhs
* start with this script
"""

"""
Summing up Fourier series is done in:
fct_RZ(), fct_RPZ(), fct_RZL(), getFluxSurfData(), getCrossSectData()

"""


__version__ = "0.4"
__versionTime__ = "19 Jan 2013"
__author__ = "Bernhard Seiwald <bernhard.seiwald@gmail.com>, Boyd D. Blackwell <boyd.blackwell@anu.edu.au>"



import os
import sys
import time
import numpy
import scipy
import scipy.interpolate
import pylab
import warnings
from copy import deepcopy

#try:
import mpl_toolkits.mplot3d.axes3d as p3
#except:
#    from mpl_toolkits.mplot3d import Axes3D

try:  # try for scipy 0.9 or greater - does cubic, but sadly, only for 1 and 2d
    from scipy.interpolate import griddata
    interp_method = 'cubic'
    version=0.9# fudge - can't string compare version 0.10.0
except: # if not, settle for matplotlib version which only does 'linear' or 'nn'
    from matplotlib.mlab import griddata
    interp_method = 'linear'
    version=0.7 # fudge

from matplotlib.backends.backend_pdf import PdfPages

try:
    from mayavi import mlab
    ISOK_MAYAVI = True
except:
    ISOK_MAYAVI = False

try:
    from scipy.io import netcdf
    net_cdf_module = 'netcdf'
except ImportError:
    try:
        import netCDF4 as netcdf
        netcdf_module = 'netCDF4'
    except: raise ImportError('Unable to load netcdf module/library')


from H1ErrorClasses import *
import convenience
import netCDF4Tools
import SPLINE
from spline_spec import leadingfunction_mod as leadingfunction
from spline_spec import spline_spec_mod as spline_spec
import GRID3D
import miscmath
from catmull_rom_interpolation import catmull_int2d as interp2d

PI = numpy.pi
PI2 = 2.*PI


debug=0



__all__ = ['VMEC']


__all_splines__ = -1

class VMEC():
    """
    A class for VMEC models.
    This class deals with NetCDF output (wout*.nc) of VMEC.

    Author: Bernhard Seiwald
    Date:   

    """
    from time import clock as secs

    msgs = {}

    def __init__(self, 
                 nc_filename, 
                 import_all=True, 
                 show_vars=False,
                 maxm=None,
                 compute_spline_type=None, 
                 load_spline=True, 
                 save_spline=False,
                 use_spline_type=SPLINE.CONV_SPLINE,
                 compute_grid=False,
                 load_grid=True, 
                 save_grid=False, 
                 gridfunction=None,
                 gnr=None, gnp=None, gnz=None):
        """
        INPUT:
          nc_filename, 
          import_all=True, 
          show_vars=False,
          maxm=None,
          compute_spline_type=__all_splines__, 
          load_spline=True, 
          save_spline=False,
          use_spline_type=SPLINE.CONV_SPLINE,
          compute_grid=False,
          load_grid=True, 
          save_grid=False, 
          gridfunction=None,
          gnr=None, gnp=None, gnz=None):

        """
        print "loading VMEC file " + nc_filename + "..."
        self.nc_filename = nc_filename
        self.ds = netcdf.netcdf_file(self.nc_filename, 'r')

        # import all variables
        self.import_all = import_all
        if import_all:
            self.importAll()
 
        self.DMU0 = 4.0*numpy.pi*1.0e-7        # 1.256637d-06

        self.ncpath,ncfname = os.path.split(self.nc_filename)
        self.config_name, scfilename_ext = os.path.splitext(ncfname)

        self.es = numpy.array(self.ds.variables['phi'][:]/max(self.ds.variables['phipf']))

        self.compute_spline_type = 0
        if(compute_spline_type == None):
            compute_spline_type = 0

        if(compute_spline_type == __all_splines__):
            self.compute_spline_type = self.allSplMeths()
        else:
            self.compute_spline_type = compute_spline_type

        self.splkind = use_spline_type

        # compute coefficients for splines
        self.spline_coefs_done = False
        spline_name = "".join([self.config_name, '.splcoef'])
        if(load_spline):
            try:
                print "try loading spline coefficients..."
                print "".join([spline_name, ".rmnc"])
                self.rSpline = SPLINE.SPLINE(nc_name="".join([spline_name, 
                                                              ".rmnc"]))
                self.lSpline = SPLINE.SPLINE(nc_name="".join([spline_name, 
                                                              ".lmns"]))
                self.zSpline = SPLINE.SPLINE(nc_name="".join([spline_name, 
                                                              ".zmns"]))
                self.bSpline = SPLINE.SPLINE(nc_name="".join([spline_name, 
                                                              ".bmnc"]))
                print "...spline coefficients loaded"
                self.spline_coefs_done = True
            except:
                print "can't load spline coefficients"

        if((compute_spline_type != 0) and not self.spline_coefs_done):
            self.calc_spline_coefs()
            if(save_spline):
                self.rSpline.save("".join([spline_name, ".rmnc"]))
                self.lSpline.save("".join([spline_name, ".lmns"]))
                self.zSpline.save("".join([spline_name, ".zmns"]))
                self.bSpline.save("".join([spline_name, ".bmnc"]))
        # compute 3D grid if asked for

        if(gridfunction == None):
            gridfunction=self.fct_RPZ

        self.maxm = maxm
        grid_name = "".join([self.config_name, '.grid3D.nc'])
        #try:
        self.grid3D = GRID3D.GRID3D(nc_filename=grid_name, 
                                    compute_grid=compute_grid,
                                    load_grid=load_grid, 
                                    save_grid=save_grid,
                                    gridfunction=gridfunction, 
                                    gnr=gnr, gnp=gnp, gnz=gnz)
        #except:
        #    print "grid3D can't be provided..."
        print "3d grid is {0:s}available".format("" if self.grid3D.isAvail else "not ")

        if(show_vars): self.show_variables(self)

        print "...VMEC file " + nc_filename + " loaded"

        return


    def warn_but_not_too_often(self,message, wait=3):
        """ manage a message list to group similar messages together,
        if they come too often.  will return true if it is time to print the
        message again, and print the message
        """
        if self.msgs.has_key(message):
            that_msg = self.msgs[message]
            if  self.secs() > (that_msg['tim']+wait):
                if that_msg['count']>0: 
                    print("%d occurences of " % that_msg['count']),

                print(message)
                tim =self.secs()
                count=0

            else: 
                count = that_msg['count']+1
                tim =  that_msg['tim']
        else:
            count=0
            tim = self.secs()
        self.msgs.update({message:{'count':count, 'tim': tim}})
        return(count==0)
    

    def importAll(self):
        """

        Import all variables.
        Set self.import_all = True

        """

        # dimensions
        dimensions = self.ds.dimensions
        for var in dimensions:
            destvar = convenience.clearFilename(var)
            destvar = destvar.replace('-', '_')
            exestr = "self."+destvar+"=self.ds.dimensions['"+var+"']"
            exec(exestr)

        # variables
        variables = self.ds.variables
        for var in variables:
            if(len(self.ds.variables[var].shape) == 0):
                exestr = "self."+var+"=self.ds.variables['"+var+"'].getValue()"
                exec(exestr)
            else:
                exestr = "self."+var+"=self.ds.variables['"+var+"'][:]"
                exec(exestr)

        self.import_all = True
        return


    def showVariables(self, max_width=80, indent='  '):
        """
        """

        line = indent
        for k in self.ds.variables.keys(): 
            if len(self.ds.variables[k].shape) == 0:
                op = str('%s=%s' % (k,self.ds.variables[k].getValue()))
            else:
                op = str('%s arr%s' % (k,self.ds.variables[k].shape))
            if len(line+op) > max_width:
                print line
                line = indent
            else: line += op + ', '
        if line != indent: print line
    
        return(0)


    ##### some stuff for interpolation start #####

    def allSplMeths(self):
        """
        Determine all methods, which are provided by SPLINE. 


        INPUT:
          None

        RETURN:
          allmeths ... integer number representing all methods
        """
        spl = SPLINE.SPLINE({'x':numpy.arange(5), 'y':numpy.arange(5)})
        allmeths = sum(spl.SPLINE_METHODS_l)
        del(spl)

        return(allmeths)


    def calc_spline_coefs(self):
        """
        Calculate spline coefficients.

        Have a good look at the documentation of the calss SPLINE!

        * for computation of coefficients the flag compute_spline_type
          is taken into account: 
         if SQRT_SPLINE:
            divide odd m terms by sqrt(s), keep even m terms as they are
            use 3rd order qubic spline

         if CONV_SPLINE:
            use just 3rd order qubic spline

         if SPECIAL_SPLINE:
            use special spline equiped with leading function (takes into
            account sqare root characteristics of Rmnc, Zmns).
            make sure that spline_spec is imported.

        * for interpolation/evaluation: 
          use e.g. rSpline(s, deriv=False, kindspec=None)


        INPUT:
          None

        RETURN:
          None
        """

        if(not self.import_all):
            self.importAll()   # we need to import the quantities

        self.rSpline = SPLINE.SPLINE({'x':self.es, 
                                      'y':self.rmnc, 
                                      'kindspec':self.compute_spline_type, 
                                      'ixm':self.xm})
        self.lSpline = SPLINE.SPLINE({'x':self.es, 
                                      'y':self.lmns, 
                                      'kindspec':self.compute_spline_type, 
                                      'ixm':self.xm})
        self.zSpline = SPLINE.SPLINE({'x':self.es, 
                                      'y':self.zmns, 
                                      'kindspec':self.compute_spline_type, 
                                      'ixm':self.xm})
        self.bSpline = SPLINE.SPLINE({'x':self.es, 
                                      'y':self.bmnc, 
                                      'kindspec':self.compute_spline_type, 
                                      'ixm':self.xm})
        self.spline_coefs_done = True
        print "...done"
        return


    def fct_RZ(self, s, theta, phi, calc_derivs=False):

        """
        (r, z, rs, rt, rp, zs, zt, zp) = 
             fct_RZ(self, s, theta, phi, calc_derivs=True)

        (r, z) = 
             fct_RZ(self, s, theta, phi, calc_derivs=False)

        Compute 
        R(s,theta,phi) = SUM(rmnc_b*cos(ixm*theta - ixn*phi))
        and 
        Z(s,theta,phi) = SUM(zmns_b*sin(ixm*theta - ixn*phi))
        and the derivatives with respect to s and theta.
        rmnc, zmns ... VMEC 'spectrum'

        Returns R, Z and derivatives.

        INPUT:
          s ...... flux surface label
          theta, 
          phi .... angle like VMEC variables 
                   (phi equals the reals space cylindrical phi)
          calc_derivs . compute derivatives (logical)

        OUTPUT:
          R, Z
          Rs, Rt, Rp ... derivatives of R with respect to s, theta, phi
          Zs, Zt, Zp ... derivatives of Z with respect to s, theta, phi

        """

        maxm = None

        if(s != 0.0):
            les = s
        else:
            les = 1.0e-30

        phi = numpy.mod(phi,PI2)
        if(calc_derivs):
            (Rmn, Zmn, Rmns, Zmns) = self.calc_rz_coefs(les, calc_derivs=calc_derivs)
        else:
            (Rmn, Zmn) = self.calc_rz_coefs(les, calc_derivs=calc_derivs)

        #if(maxm == None): maxm = self.maxm

        #if(maxm != None):
        #    elts_used = numpy.where(self.xm<maxm)[0]
        #else:
        #    elts_used = numpy.arange(len(self.rmnc[0]))

        ixm = self.xm #[elts_used]
        ixn = self.xn #[elts_used]
        Rmn = Rmn #[elts_used]
        Zmn = Zmn #[elts_used]

        if(calc_derivs):
            Rmns = Rmns #[elts_used]
            Zmns = Zmns #[elts_used]

        argv = ixm*theta - ixn*phi
        cosv = numpy.cos(argv)
        sinv = numpy.sin(argv)
  
        R = sum(Rmn  * cosv)
        Z = sum(Zmn  * sinv)

        if(calc_derivs):
            Rs =  sum(rmns * cosv)
            Rt = -sum(rmn * sinv * ixm)
            Rp =  sum(rmn * sinv * ixn)
            Zs =  sum(zmns * sinv)
            Zt =  sum(zmn * cosv * ixm)
            Zp = -sum(zmn * cosv * ixn)
            return(R, Z, Rs, Rt, Rp, Zs, Zt, Zp)
        else:
            return(R, Z)


    def fct_RPZ(self, s, theta, phi_rs, calc_derivs=False):

        """
        (R, phi, Z, Rs, Rt, Rp, phis, phit, phip, Zs, Zt, Zp) =
             fct_RPZ(self, s, theta, phi_rs, calc_derivs=True)

        (R, phi, Z) = 
             fct_RPZ(self, s, theta, phi_rs, calc_derivs=False)

        Compute 
        R(s,theta,phi) = SUM(rmnc_b*cos(ixm*theta - ixn*phi_rs))
        and 
        Z(s,theta,phi) = SUM(zmns_b*sin(ixm*theta - ixn*phi_rs))
        and the derivatives with respect to s and theta.
        rmnc, zmns ... VMEC 'spectrum'

        Returns R, phi=phi_rs, Z and derivatives.

        ATTENTION: This method is introduced to be compatible with 
        the corresponding method BOOZER.fct_RPZ(...).

        INPUT:
          s ...... flux surface label
          theta, 
          phi_rs . angle like VMEC variables 
                   (phi equals the reals space cylindrical phi)
          calc_derivs . compute derivatives (logical)

        OUTPUT:
          R, phi, Z
          Rs, Rt, Rp ... derivatives of R with respect to s, theta, phi
          phis, phit, phip ... zero
          Zs, Zt, Zp ... derivatives of Z with respect to s, theta, phi

        """

        maxm = None

        if(s != 0.0):
            les = s
        else:
            les = 1.0e-30

        phi_rs = numpy.mod(phi_rs,PI2)
        if(calc_derivs):
            (Rmn, Zmn, Rmns, Zmns) = self.calc_rz_coefs(les, calc_derivs=calc_derivs)
        else:
            (Rmn, Zmn) = self.calc_rz_coefs(les, calc_derivs=calc_derivs)

        #if(maxm == None): maxm = self.maxm

        #if(maxm != None):
        #    elts_used = numpy.where(self.xm<maxm)[0]
        #else:
        #    elts_used = numpy.arange(len(self.rmnc[0]))

        ixm = self.xm #[elts_used]
        ixn = self.xn #[elts_used]
        Rmn = Rmn #[elts_used]
        Zmn = Zmn #[elts_used]

        if(calc_derivs):
            Rmns = Rmns #[elts_used]
            Zmns = Zmns #[elts_used]

        argv = ixm*theta - ixn*phi_rs
        cosv = numpy.cos(argv)
        sinv = numpy.sin(argv)
  
        R = sum(Rmn  * cosv)
        Z = sum(Zmn  * sinv)
        phi = numpy.empty(R.shape)
        phi.fill(phi_rs)

        if(calc_derivs):
            Rs =  sum(rmns * cosv)
            Rt = -sum(rmn * sinv * ixm)
            Rp =  sum(rmn * sinv * ixn)
            Zs =  sum(zmns * sinv)
            Zt =  sum(zmn * cosv * ixm)
            Zp = -sum(zmn * cosv * ixn)
            #return(R, Z, Rs, Rt, Rp, Zs, Zt, Zp)
            phis = numpy.zeros(Rs)
            phit = numpy.zeros(Rs)
            phip = numpy.ones(Rs)
            return(R, phi, Z, Rs, Rt, Rp, phis, phit, phip, Zs, Zt, Zp)
        else:
            #return(R, Z)
            return(R, phi, Z)



    def fct_RZL(self, s, theta, phi, calc_derivs=False):

        """
        (R, Z, l, Rs, Rt, Rp, Zs, Zt, Zp, ls, lt, lp) = 
                  fct_RZL(self, s, theta, phi, calc_derivs=True)

        (R, Z, l) = 
                  fct_RZL(self, s, theta, phi, calc_derivs=False)

        Compute 
        R(s,theta_b,phi_b) = SUM(rmnc*cos(ixm*theta + ixn*phi))
        Z(s,theta_b,phi_b) = SUM(zmns*sin(ixm*theta + ixn*phi))
        and thye streamfunction lambda
        l(s,theta_b,phi_b) = SUM(lmns*sin(ixm*theta + ixn*phi))
        and the derivatives with respect to s, theta_b and phi_b.
        Rmnc, Zmns, lmns ... VMEC 'spectra'

        Returns R, Z, l and derivatives.

        INPUT:
          s ..... flux surface label
          theta, 
          phi ... angle like VMEC variables
          calc_derivs . compute derivatives (logical)

        OUTPUT:
          R, phi, Z
          Rs, Rt, Rp ... derivatives of R with respect to s, theta_b, phi_b
          Zs, Zt, Zp ... derivatives of Z with respect to s, theta_b, phi_b
          ls, lt, lp ... derivatives of Z with respect to s, theta_b, phi_b

        """

        maxm = None

        if(s != 0.0):
            les = s
        else:
            les = 1.0e-30

        phi = numpy.mod(phi,PI2)

        if(calc_derivs):
            (Rmn, Zmn, lmn, Rmns, Zmns, lmns) = self.calc_rzl_coefs(les, calc_derivs=calc_derivs)
        else:
            (Rmn, Zmn, lmn) = self.calc_rzl_coefs(les, calc_derivs=calc_derivs)

        if(maxm == None): maxm = self.maxm

        if(maxm != None):
            elts_used = numpy.where(self.ixm<maxm)[0]
        else:
            elts_used = numpy.arange(len(self.rmnc[0]))

        ixm = self.xm #[elts_used]
        ixn = self.xn #[elts_used]

        Rmn = Rmn #[elts_used]
        Zmn = Zmn #[elts_used]
        lmn = lmn #[elts_used]

        if(calc_derivs):
            Rmns = Rmns #[elts_used]
            Zmns = Zmns #[elts_used]
            lmns = lmns #[elts_used]

        argv = ixm*theta - ixn*phi
        cosv = numpy.cos(argv)
        sinv = numpy.sin(argv)
  
        R = sum(Rmn * cosv)
        Z = sum(Zmn * sinv)
        l = sum(lmn * sinv)

        if(calc_derivs):
            Rs =  sum(Rmns * cosv)
            Rt = -sum(Rmn * sinv * ixm)
            Rp =  sum(Rmn * sinv * ixn)
            Zs =  sum(Zmns * sinv)
            Zt =  sum(Zmn * cosv * ixm)
            Zp = -sum(Zmn * cosv * ixn)
            ls =  sum(lmns * sinv)
            lt =  sum(lmn * cosv * ixm)
            lp = -sum(lmn * cosv * ixn)
            return(R, Z, l, Rs, Rt, Rp, Zs, Zt, Zp, ls, lt, lp)
        else:
            return(R, Z, l)


    def calc_rz_coefs(self, s, calc_derivs=True):
        """
        (rmnc, zmns, rmncs, zmnss) = calc_rz_coefs(self, s, calc_derivs=True)

        Compute coefficients rmnc_b(s), zmns_b(s).
        The derivatives with respect to s are computed 
        if set: calc_derivs=True. 
        All modes are taken into account.
        Number of modes for rmnc_b(s), zmns_b(s) are assumed to 
        be the same and not 'cross checked'!

        RETURN:
          rmnc, zmns, rmncs, zmnss ... if set calc_derivs=True
          rmnc, zmns ................. if set calc_derivs=False
        """

        if (not self.spline_coefs_done):
            print 'calculating spline coefs'
            self.calc_spline_coefs()

        RR = self.rSpline.eval(s, deriv=calc_derivs, kindspec=self.splkind)
        ZZ = self.zSpline.eval(s, deriv=calc_derivs, kindspec=self.splkind)

        if calc_derivs:
            return(RR[0,0,:], ZZ[0,0,:], 
                   RR[0,1,:], ZZ[0,1,:])
        else:
            return(RR[0,:], ZZ[0,:])



    def calc_rzl_coefs(self, s, calc_derivs=True):
        """
        (rmnc, zmns, lmns, rmncs, zmnss, lmnss) = calc_rzl_coefs(self, s, calc_derivs=True)

        Compute coefficients rmnc(s), zmns(s), lmns(s).
        The derivatives with respect to s are computed 
        if set: calc_derivs=True. 
        All modes are taken into account.
        Number of modes for rmnc(s), zmns(s), lmns(s) are assumed to 
        be the same and not 'cross checked'!

        RETURN:
          rmnc, zmns, lmns, rmncs, zmnss, lmnss ... if set calc_derivs=True
          rmnc, zmns, lmns ........................ if set calc_derivs=False
       """

        if (not self.spline_coefs_done):
            print 'calculating spline coefs'
            self.calc_spline_coefs()

        RR = self.rSpline.eval(s, deriv=calc_derivs, kindspec=self.splkind)
        ZZ = self.zSpline.eval(s, deriv=calc_derivs, kindspec=self.splkind)
        ll = self.lSpline.eval(s, deriv=calc_derivs, kindspec=self.splkind)

        if calc_derivs:
            return(RR[0,0,:], ZZ[0,0,:], ll[0,0,:], 
                   RR[0,1,:], ZZ[0,1,:], ll[0,1,:])
        else:
            return(RR[0,:], ll[0,:], ll[0,:])


    ##### some stuff for interpolation end #####


    ##### plotting start #####

    def plotFluxTube(self, s, s_ind=None, no_theta=30, no_phi=180, 
                     adjust_axes=False, figno=None, USE_MAYAVI=True,
                     withB=True, **plkwargs):
        """
        ANSCHAUEN!!!
        http://matplotlib.sourceforge.net/examples/mplot3d/surface3d_demo3.html
        http://matplotlib.sourceforge.net/examples/mplot3d/wire3d_demo.html

        3d plot (geometry) of a flux tube/surface specified by s or s_ind.

        Usage: e.g.:
        boozmn.plotFluxTube(s_ind=-1)

        ATTENTION: Input for s can't be handled until now!
                   Reason: Interpolation (e.g. splines) needed in order
                   to compute Fourier coefficients, e.g. Rmnc(s)...
                   copy/paste from other routines...

        INPUT:
          s .......... flux surface label; array like
          s_ind ...... chosen flux surface index 
          no_theta ... no points on phi=const plane
          no_phi ..... no points in phi direction
          adjust_axes  True/False adjust axes by plotCrossSect()
                       Useful, if equal scaling for the axes is required.
          figno ...... number of figure
          USE_MAYAVI . use mayavi if available
          withB ...... use mod B for color information
          **plkwargs . arguments directly parsed to plot command

        RETURN:
          err ......... 0 for normal operation

        """

        if([s,s_ind].count(None) == 2):
            print "chose either s or s_ind!"
            return(1)
        # end   if([s,s_ind].count(None) == 2):

        if(not(ISOK_MAYAVI & USE_MAYAVI)):
            meshgrid = False
        else:
            meshgrid = True
        # end   if(not(ISOK_MAYAVI & USE_MAYAVI)):

        points = self.getFluxsurfData(s, s_ind=s_ind, 
                                      no_theta=no_theta, no_phi=no_phi,
                                      retB=withB,
                                      coordsys="cart", meshgrid=meshgrid)
        fig = self.doPlot3D(points, phi_b=None, s_b=s, 
                            adjust_axes=adjust_axes, 
                            figno=figno, USE_MAYAVI=USE_MAYAVI,
                            **plkwargs)

        return(fig)


    def doPlot3D(self, points, phi_b=None, s_b=None, 
                 adjust_axes=False, figno=None, USE_MAYAVI=True,
                 **plkwargs):
        """
        Do the 3d plot of a dataset.
        
        INPUT:
          points ..... array of coordinates as returned by getFluxsurfData()
                       ATTENTION: if MAYAVI is used, getFluxsurfData() has to
                       be called with meshrgid=True! In that case
                       no_theta*no_phi -> no_theta,no_phi
                       points.shape:   [nsurf,x,no_theta*no_phi]
                       x might be x=3 for pure coordinates or
                       for x=4 the 4th 'column' contains color information
                       e.g. mod B
          phi_b ...... Boozer angle like variable, defines cross section
          s_b ........ flux surface labels; array like
          adjust_axes  True/False adjust axes by plotCrossSect()
                       Useful, if equal scaling for the axes is required.
          figno ...... number of figure
          USE_MAYAVI . use mayavi if available
          **plkwargs . arguments directly parsed to plot command
        
        RETURN:
          fig
        """
        cols = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
        ncols = len(cols)

        nsurf = points.shape[0]
        if(s_b == None):
            nlegs = 0
        else:
            try:
                nlegs = len(s_b)
            except:
                nlegs = 1  # when we are her, s_b should be a scalar

        if((nlegs != nsurf) & (nlegs != 0)):
            print "Numers of fluxsurfaces and fluxsurface labes do not match."
            print "Don't show fluxfurface labels."
            s_b = None
            nlegs = 0

        # adjust axes limits in order to plot 'equal'
        if(adjust_axes):
            dd = 0.05
            xmin = points[:,0,:].min()
            xmax = points[:,0,:].max()
            ymin = points[:,1,:].min()
            ymax = points[:,1,:].max()
            zmin = points[:,2,:].min()
            zmax = points[:,2,:].max()
            dx = (xmax-xmin)*dd
            dy = (ymax-ymin)*dd
            dz = (zmax-zmin)*dd
            xmin -= dx
            xmax += dx
            ymin -= dy
            ymax += dy
            zmin -= dz
            zmax += dz

            dx = xmax - xmin
            dy = ymax - ymin
            dz = zmax - zmin
            if(dx <= dz):
                xmin = (xmin+xmax)/2 - dz/2.
                xmax = (xmin+xmax)/2 + dz/2.
            else:
                zmin = (zmin+zmax)/2 - dx/2.
                zmax = (zmin+zmax)/2 + dx/2.
            # end   if(dx <= dz):
        # end   if(adjust_axes):

        if(not(ISOK_MAYAVI & USE_MAYAVI)):
            fig = pylab.figure()
            try:
                # try new way. for matplotlib> v1.0.0 
                ax = fig.add_subplot(111, projection='3d')
            except:
                # ok, we need the old way
                ax = p3.Axes3D(fig)
            # end   try:
    
            pylab.rcParams['legend.fontsize']='small'

            s_b_str = ""
            if(s_b != None):
                try:
                    ll = len(s_b)
                    s_b = numpy.array(s_b)
                except:
                    s_b = numpy.array([s_b])
                # end   try:
    
                mv = numpy.abs(s_b).max()
                # stringlength slen:  e.g. a=0.45687
                # we like 3 digs behind komma
                # l = int(numpy.log10(mv)) -> l=0  ;   digits before komma: l+1
                # slen = digits before komma + komma + digs behind komma
                if(mv > 0.0):
                    slen = int(numpy.log10(mv)) + 1 + 1 + 3
                else:
                    slen = 3
                # end   if(mv > 0.0):
    
                acc = 1
            # end   if(s_b != None):
    
            if(phi_b != None):
                title_str = "".join([self.nc_filename.strip(), '\n',
                                     r"$\varphi={0:.2f}^\circ$".format(numpy.rad2deg(phi_b))])
            else:
                 title_str = self.nc_filename.strip()
            # end   if(phi_b != None):

            ax.set_title(title_str)
            ax.set_aspect('equal')
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')
            ax.grid(True)

            surf = []
            legl = []
            legt = []

            if(plkwargs.has_key('marker')):
                predefmarker = plkwargs['marker']
            else:
                predefmarker = ''
            # end   if(plkwargs.has_key('marker')):

            llen = -1
            for i in range(nsurf):
                if(s_b != None):
                    # set marker, if not set, for axis to '+'
                    try:
                        sb = s_b[i]
                    except:
                        sb = s_b
                    if(sb == 0.0):
                        if(not plkwargs.has_key('marker')):
                            plkwargs.update({'marker':'+'})
                    else:
                        if(plkwargs.has_key('marker')):
                            plkwargs.update({'marker':predefmarker})

                nc = i%ncols
                ll, = ax.plot3D(points[i,0,:], points[i,1,:], points[i,2,:],
                                **plkwargs)
                if(nlegs != 0):
                    legl.append(ll)
                    legt.append("s={0:s}".format(convenience.formatNumber(s_b[i], slen, acc)))
                    llen = max(llen, len(legt[-1]))

                ax.hold(True)


            # Put a legend below current axis
            ncols = 40 / llen
            if(nlegs != 0):
                leg = ax.legend(legl, legt, loc='upper center',
                                bbox_to_anchor=(0.5, 0.01), numpoints=1, ncol=ncols)


            ax.set_aspect('equal')
            ax.grid(True)
            if(adjust_axes):
                ax.set_xlim(xmin, xmax)
                ax.set_ylim(zmin, zmax)
            fig.show()
        else: # if(not(ISOK_MAYAVI & USE_MAYAVI)):
            fig = None
            alpha = 1.0
            fig = mlab.figure(fgcolor=(0, 0, 0), bgcolor=(1, 1, 1))
            plotOK = True
            if(points.shape[1] == 3):
                for i in range(nsurf):
                    fluxsurf = mlab.mesh(points[i,0,:,:], 
                                         points[i,1,:,:], 
                                         points[i,2,:,:], 
                                         opacity=alpha) #colormap='jet'
                # end   for i in range(nsurf):
            elif(points.shape[1] == 4):
                for i in range(nsurf):
                    fluxsurf = mlab.mesh(points[i,0,:,:], 
                                         points[i,1,:,:], 
                                         points[i,2,:,:], 
                                         scalars=points[i,3,:,:], 
                                         opacity=alpha, colormap='jet')
                # end   for i in range(nsurf):
                dB = points[:,3,:,:].max() - points[:,3,:,:].min()
                if(dB < 0.5):
                    fmt_str = '%.2f'
                else:
                    fmt_str = '%.1f'
                # end   if(dB < 0.5):
                clbh = mlab.colorbar(title='B [T]', 
                                     orientation='vertical', 
                                     label_fmt=fmt_str)
            else:
                tstr = "".join(["doPlot3D(): Wrong number of 'columns!\n",
                                "  points.shape must be (nsurf, 3/4, ...)\n",
                                "  second dimension is ",
                                "{0:d} ".format(points.shape[1]),
                                "instead of either 3 or 4!"])
                print tstr
                plotOK = False
            # end   if(points.shape[1] == 3):
            if(plotOK):
                olh = mlab.outline(figure=fig)
                axh = mlab.axes(olh,  # xlabel=r'$\alpha$'
                                xlabel='X', ylabel='Y', zlabel='Z')
                if(adjust_axes):
                    (zmin, zmax) = (xmin, xmax)
                    mlab.axes(axh,
                              extent=[xmin, xmax, ymin, ymax, zmin, zmax])
                # end   if(adjust_axes):
                #mlab.axes(fluxsurf)
                viewh = mlab.view(45.0, 60.0, 7, figure=fig)
                mlab.show()
            # end   if(plotOK):
        # end   if(ISOK_MAYAVI & USE_MAYAVI):

        return(fig)




    def getFluxsurfData(self, s_b, s_ind=None, no_theta=30, no_phi=180,
                        retB=False, coordsys=None, meshgrid=False):
        """
        Compute points on specified fluxtubes/fluxsurfaces.

        ATTENTION:
          At least one flux surface label has to be provided for s!!!
          If s_b >= 0 and s_ind are present, both are used.
            e.g.: wout.getCrossSectData(0, [0.1,0.9], [0,1])
            returns cross section data for s_b=0.1, s_b=0.9, the magnetic axis 
            and the last closed surface at phi=0.
          In case s_b == None, s_ind rather than s is used.
            e.g.: wout.getCrossSectData(0, -1, [0,1])
            returns cross section data for the magnetic axis and 
            the last closed surface at phi=0.

        INPUT:
          s_b ...... flux surface labels (array like)   nsurf=len(s_b)
          s_ind ...... indices of flux surfaces to plot (array like)
          no_theta . no points for VMEC angle like variable theta
          no_phi ... no points for VMEC angle variable phi
          retB ..... return mod B  (True/False)
          coordsys ... return points in specified coordinate system
                       "cart" (Cartesian, default), "cyl" (cylindrical)
          meshgrid ... output is meshgrid for each fluxsurface
                       points.shape:   [nsurf,x,no_phi,no_theta]

        RETURN:
          points ... array of coordinates
                     when meshgrid=True  no_theta*no_phi -> no_theta,no_phi
                     either [nsurf,3,no_theta*no_phi] for retB=False
                     or [nsurf,4,no_theta*no_phi] for retB=True

        """

        if([s_b,s_ind].count(None) == 2):
            raise InputDataError("no values for s_b and s_ind")

        if(coordsys == None):
            coordsys = "cart"

        coordsys = coordsys.lower()

        # s_b should be an array (or array like)
        if(s_b != None):
            try:
                ll = len(s_b)
                s_b = numpy.array(s_b)
            except:
                s_b = numpy.array([s_b])
            # end   try:
        # end   if(s_b != None):
        print "s_ind, s_b ", s_ind, s_b

        pureInd = False
        if(s_ind != None):
            s_ind = numpy.array(s_ind)
            if(s_b == None):
                s_b = s_ind  #s_b = self.es[s_ind]
                pureInd = True
            else:
                s_b = numpy.append(s_b, self.es[s_ind])
            # end   if(s_b == None):
        # end   if(s_ind != None):

        # in case s_b is a scalar, make a real array
        if(len(s_b.shape) == 0):
            sb = s_b
            if(pureInd):  
                s_b = numpy.empty([1], dtype=numpy.int)
            else:
                s_b = numpy.empty([1])
            # end   if(pureInd): 
            s_b[0] = sb
        # end   if(len(s_b.shape) == 0):

        points = []

        t = numpy.linspace(0, PI2, no_theta, endpoint=True)
        p = numpy.linspace(0, PI2, no_phi,   endpoint=True)

        (theta,phi_b) = numpy.meshgrid(t, p)
        zerom = numpy.zeros(theta.shape)

        for s in s_b:
            if(pureInd):
                rmnc = self.rmnc[s]
                zmns = self.zmns[s]
                if(retB):
                    bmnc = self.bmnc[s]
                # end   if(retB):
            else:
                rmnc = self.rSpline.eval(s, kindspec=self.splkind)
                zmns = self.zSpline.eval(s, kindspec=self.splkind)
                if(retB):
                    bmnc = self.bSpline.eval(s, kindspec=self.splkind)
                # end   if(retB):
            # end   if(pureInd):


            r_c = numpy.zeros(no_theta*no_phi)
            z_c = numpy.zeros(no_theta*no_phi)
            p_c = numpy.zeros(no_theta*no_phi)
            if(retB):
                B = numpy.zeros(no_theta*no_phi)

            k = 0
            for theta in t:
                for phi_b in p:
                    argv = self.xm*theta - self.xn*phi_b
                    cosv = numpy.cos(argv)
                    sinv = numpy.sin(argv)
                    p_c[k] = phi_b
                    r_c[k] = numpy.sum(rmnc*cosv)
                    z_c[k] = numpy.sum(zmns*sinv)
                    if(retB):
                        B[k] = numpy.sum(bmnc*cosv)
                    # end   if(retB):
                    k += 1
                # end   for phi_b in p:
            # end   for theta in t:

            # make meshgrids   reshape([no_theta,no_phi]).transpose()
            if(meshgrid):
                r_c = numpy.reshape(r_c,[no_theta,no_phi]).transpose()
                p_c = numpy.reshape(p_c,[no_theta,no_phi]).transpose()
                z_c = numpy.reshape(z_c,[no_theta,no_phi]).transpose()
                if(retB):
                    B = numpy.reshape(B,[no_theta,no_phi]).transpose()
            # end   if(meshgrid):

            if(coordsys == "cart"):
                X = r_c * numpy.cos(p_c)
                Y = r_c * numpy.sin(p_c)
                Z = z_c
                if(retB):
                    points.append([X, Y, Z, B])
                else:
                    points.append([X, Y, Z])
                # end   if(retB):
            else:
                if(retB):
                    points.append([r_c, p_c, z_c, B])
                else:
                    points.append([r_c, p_c, z_c])
                # end   if(retB):
            # end   if(coordsys == "cart"):
        # end   for s in s_b:

        return(numpy.array(points))



    def plotCrossSect(self, phi, s, s_ind=None, no_theta=100, 
                      theta_s=0, theta_e=PI2,
                      adjust_axes=False, **plkwargs):
        """
        Plot no_theta points on a cross section for each flux surface 
        specified by s or s_ind.

        ATTENTION:
          At least one flux surface label has to be provided for s!!!
          If s >= 0 and s_ind are present, both are used.
            e.g.: wout.getCrossSectData(0, [0.1,0.9], [0,1])
            returns cross section data for s=0.1, s=0.9,  the magnetic axis 
            and the last closed surface at phi=0.
          In case s == None, s_ind rather than s is used.
            e.g.: wout.getCrossSectData(0, -1, [0,1])
            returns cross section data for the magnetic axis and 
            the last closed surface at phi=0.

        Usage: e.g.:
        VMECwout.plotCrossSect(numpy.deg2rad(30),s_ind=-1))


        INPUT:
          phi ........ real space toroidal angle in rad
          s .......... flux surface label (array like)
          s_ind ...... chosen flux surface index (array like)
          no_theta ... no points on phi=const plane
          theta_s .... start value for theta in rad
          theta_e .... end value for theta in rad
          adjust_axes  True/False adjust axes by plotCrossSect()
                       Useful, if equal scaling for the axes is required.
          **plkwargs . arguments parsed dirctly to plot command

        RETURN:
          err ......... 0 for normal operation

        """

        if([s,s_ind].count(None) == 2):
            raise InputDataError("no values for s and s_ind")
        # end   if([s,s_ind].count(None) == 2):

        points = self.getCrossSectData(phi, s, s_ind=s_ind, 
                                       no_theta=no_theta, 
                                       theta_s=theta_s, 
                                       theta_e=theta_e, 
                                       coordsys="cyl")

        cols = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
        ncols = len(cols)

        nsurf = points.shape[0]
        if(s == None):
            nlegs = 0
        else:
            try:
                nlegs = len(s_b)
            except:
                nlegs = 1  # when we are her, s_b should be a scalar
            # end   Try:
        # end   if(s == None):

        if((nlegs != nsurf) & (nlegs != 0)):
            print "Numers of fluxsurfaces and fluxsurface labes do not match."
            print "Don't show fluxfurface labels."
            s = None
            nlegs = 0
        # end   f((nlegs != nsurf) & (nlegs != 0)):

        # adjust axes limits in order to plot 'equal'
        if(adjust_axes):
            dd = 0.05
            rmin = points[:,0,:].min()
            rmax = points[:,0,:].max()
            zmin = points[:,2,:].min()
            zmax = points[:,2,:].max()
            dr = (rmax-rmin)*dd
            dz = (zmax-zmin)*dd
            rmin -= dr
            rmax += dr
            zmin -= dz
            zmax += dz

            dr = rmax - rmin
            dz = zmax - zmin
            if(dr <= dz):
                rmin = (rmin+rmax)/2 - dz/2.
                rmax = (rmin+rmax)/2 + dz/2.
            else:
                zmin = (zmin+zmax)/2 - dr/2.
                zmax = (zmin+zmax)/2 + dr/2.
            # end   if(dr <= dz):
        # end   if(adjust_axes):


        fig = pylab.figure()
        ax = fig.add_subplot(111)

        pylab.rcParams['legend.fontsize']='small'

        s_b_str = ""
        if(s != None):
            try:
                ll = len(s)
                s = numpy.array(s)
            except:
                s = numpy.array([s])
            # end   try:

            mv = numpy.abs(s).max()
            # stringlength slen:  e.g. a=0.45687
            # we like 3 digs behind komma
            # l = int(numpy.log10(mv)) -> l=0  ;   digits before komma: l+1
            # slen = digits before komma + komma + digs behind komma
            if(mv > 0.0):
                slen = int(numpy.log10(mv)) + 1 + 1 + 3
            else:
                slen = 3
            # end   if(mv > 0.0):

            acc = 1
        # end   if(s != None):

        if(phi != None):
            title_str = "".join([self.nc_filename.strip(), '\n',
                                 r"$\varphi={0:.2f}^\circ$".format(numpy.rad2deg(phi))])
        else:
             title_str = self.nc_filename.strip()
        # end   if(phi != None):

        ax.set_title(title_str)
        ax.set_aspect('equal')

        ax.set_xlabel('R')
        ax.set_ylabel('Z')
        ax.grid(True)

        surf = []
        legl = []
        legt = []

        if(plkwargs.has_key('marker')):
            predefmarker = plkwargs['marker']
        else:
            predefmarker = ''
        # end   if(plkwargs.has_key('marker')):

        llen = -1
        for i in range(nsurf):
            if(s != None):
                # set marker, if not set, for axis to '+'
                try:
                    sb = s[i]
                except:
                    sb = s
                # end   Try:
                if(sb == 0.0):
                    if(not plkwargs.has_key('marker')):
                        plkwargs.update({'marker':'+'})
                else:
                    if(plkwargs.has_key('marker')):
                        plkwargs.update({'marker':predefmarker})
                    # end   if(plkwargs.has_key('marker')):
                # end   if(sb == 0.0):
            # end   if(s != None):

            nc = i%ncols
            ll, = ax.plot(points[i,0,:], points[i,2,:], **plkwargs)
            if(nlegs != 0):
                legl.append(ll)
                legt.append("s={0:s}".format(convenience.formatNumber(s[i], slen, acc)))
                llen = max(llen, len(legt[-1]))
            # end   if(nlegs != 0):

            ax.hold(True)
        # end   for i in range(nsurf):


        # Put a legend below current axis
        ncols = 40 / llen
        if(nlegs != 0):
            leg = ax.legend(legl, legt, loc='upper center',
                            bbox_to_anchor=(0.5, 0.01), numpoints=1, ncol=ncols)
        # end   if(nlegs != 0):

        ax.set_aspect('equal')
        ax.grid(True)
        if(adjust_axes):
            ax.set_xlim(xmin, xmax)
            ax.set_ylim(zmin, zmax)
        # end   if(adjust_axes):
        fig.show()

        return(fig)

            


        """
        # set marker, if not set, for axis to '+'
        if((s == 0.0) or (s_ind == 0)):
            if(not plkwargs.has_key('marker')):
                plkwargs.update({'marker':'+'})
        else:
            if(plkwargs.has_key('marker')):
                plkwargs.update({'marker':predefmarker})

        p = pylab.plot(r_c, z_c, 
                           label="".join([self.nc_filename, '\n'
                                          r"$\varphi$=%.3f" %(phi)]),
                           **plkwargs)

        pylab.rcParams['legend.fontsize']='small'

        #pylab.gca().set_title(self.nc_filename)
        title_str = "".join([self.nc_filename, '\n',
                             r"$\varphi=%.1f^\circ$" %(numpy.rad2deg(phi))])
        pylab.legend()
        pylab.gca().set_title(title_str)
        pylab.gca().set_aspect('equal')
        pylab.gca().set_xlabel('R')
        pylab.gca().set_ylabel('Z')
        pylab.gca().set_aspect('equal')
        pylab.gca().grid(True)
        if(adjust_axes):
            pylab.gca().set_xlim(rmin, rmax)
            pylab.gca().set_ylim(zmin, zmax)

        pylab.show()

        return(0)
        """

    def getCrossSectData(self, phi, s, s_ind=None, theta_arr=None,
                         no_theta=30, theta_s=0, theta_e=PI2, 
                         coordsys=None):
        """
        Calculate points on cross section specified by toroidal 
        VMEC (= cylindrical) angle phi.

        ATTENTION:
          At least one flux surface label has to be provided for s!!!
          If s >= 0 and s_ind are present, both are used.
            e.g.: wout.getCrossSectData(0, [0.1,0.9], [0,1])
            returns cross section data for s=0.1, s=0.9,  the magnetic axis 
            and the last closed surface at phi=0.
          In case s == None, s_ind rather than s is used.
            e.g.: wout.getCrossSectData(0, -1, [0,1])
            returns cross section data for the magnetic axis and 
            the last closed surface at phi=0.
        
        INPUT:
          phi ........ VMEC (= cylindrical) toroidal angle in rad
          s .......... flux surface label (array like)   nsurf=len(s_b)
          s_ind ...... indices of flux surfaces to plot (array like)
          theta_arr .. array with values for theta; 
                       if theta_arr is provided, use theta_arr rather
                       than no_theta, theta_s, theta_e
          no_theta ... no points on phi=const plane
          theta_s .... start value for theta in rad
          theta_e .... end value for theta in rad
          coordsys ... return points in specified coordinate system
                       "cart" (Cartesian, default), "cyl" (cylindrical)

        RETURN:
          points ... array of coordinates
                     points.shape:   [nsurf,3,no_theta]

        """

        s_b = s

        # let's check for s:
        if([s_b,s_ind].count(None) == 2):
            raise InputDataError("no values for s and s_ind")
        # end   if([s_b,s_ind].count(None) == 2):

        if(coordsys == None):
            coordsys = "cart"
        # end   if(coordsys == None):

        coordsys = coordsys.lower()
        phi_b = phi

        # s_b should be an array (or array like)
        if(s_b != None):
            try:
                ll = len(s_b)
                s_b = numpy.array(s_b)
            except:
                s_b = numpy.array([s_b])
            # end   try:
        # end   if(s != None):

        pureInd = False
        if(s_ind != None):
            s_ind = numpy.array(s_ind)
            if(s_b == None):
                s_b = s_ind  #self.es[s_ind]
                pureInd = True
            else:
                s_b = numpy.append(s_b, self.es[s_ind])
            # end   if(s_b == None):
        # end   if(s_ind != None):

        # in case s_b is a scalar, make a real array
        if(len(s_b.shape) == 0):
            sb = s_b
            if(pureInd):  
                s_b = numpy.empty([1], dtype=numpy.int)
            else:
                s_b = numpy.empty([1])
            # end   if(pureInd): 
            s_b[0] = sb
        # end   if(len(s_b.shape) == 0):

        points = []

        #phi_b = numpy.mod(phi_b,PI2)
        if(theta_arr != None):
            theta = numpy.array(theta_arr.copy())
            no_theta = theta.shape[0]
        else:
            theta = numpy.linspace(theta_s, theta_e, no_theta, endpoint=True)
        # end   if(theta_arr != None):

        for s in s_b:
            if(pureInd):
                rmnc = self.rmnc[s]
                zmns = self.zmns[s]
            else:
                rmnc = self.rSpline.eval(s, kindspec=self.splkind)
                zmns = self.zSpline.eval(s, kindspec=self.splkind)
            # end   if(pureInd):

            r_c = numpy.zeros(no_theta)
            p_c = numpy.zeros(no_theta)
            z_c = numpy.zeros(no_theta)
            for (i,thet) in enumerate(theta):
                argv = self.xm*thet - self.xn*phi_b
                cosv = numpy.cos(argv)
                sinv = numpy.sin(argv)
                r_c[i] = numpy.sum(rmnc*cosv)
                p_c[i] = phi_b
                z_c[i] = numpy.sum(zmns*sinv)

            if(coordsys == "cart"):
                X = r_c * numpy.cos(p_c)
                Y = r_c * numpy.sin(p_c)
                Z = z_c

                points.append([X,Y,Z])
            elif(coordsys == "cyl"):
                points.append([r_c, p_c, z_c])
            else:
                raise UnknownCoordSysError("VMEC getCrossSectData(): unknown coordinate system")

        return(numpy.array(points))

    ##### plotting end #####

    ##### misc start #####

    def determineRA(self):
        """
        Determine major radius R and plasma radius a.
        These quantities are necessary for creating a Boozer
        input file for NEO.

        Here, Joachim Geigers approach is used in order to be 
        consistent to the Boozer files (*.bc) Joachim created
        for the benchmark paper. Joachim Geigers input files 
        are used for NEO, NEO2 and NEO-MC runs.

        INPUT:
          None

        RETURN:
          R ... major radius
          a ... plasma radius

Joachims email:
Die Berechnung von aeff ist in diesem
Programm-abschnitt (Fortran) enthalten:

     re00 = cp0
     re36 = cp0
     remv = cp0
     do j = 1, ns
       lj = (j-1)*mnmax
       remv(j) = sum(rmnc(lj+1:lj+mnmax)*zmns(lj+1:lj+mnmax)*ixm)
       do mn = 1, mnmax
         l  = mn + lj
         im = ixm(mn)
         in = ixn(mn)
         do mn2 = 1, mnmax
           l2 = mn2 + lj
           im2= ixm(mn2)
           in2= ixn(mn2)
           if(im.eq.im2) then
             re00(j) = re00(j) + rmnc(l)*zmns(l2)*im
             ifac = 1-2*mod(abs(in-in2),2)
!--------------------------------------------
!        -1   for in-in2 odd
! ifac =
!        +1   for in-in2 even
!--------------------------------------------
             re36(j) = re36(j) + rmnc(l)*zmns(l2)*im*ifac
           endif
         enddo
       enddo
     enddo
     remm = sqrt(abs( (re00 +re36) ) /2. )
     re00 = sqrt(abs( re00 ) )
     re36 = sqrt(abs( re36 ) )
     remv = sqrt(abs( remv ) )

     aeff = remm(ns)
     aefv = remv(ns)

Die Fourierkoeffizienten sind die in vmec-Koordinaten, nicht in Boozerkoordinaten!
Die beiden hier gerechneten Werte sind einmal definiert durch die Querschnitts-
flaechen in den beiden Symmetrieebenen (remm bzw aeff) und zum anderen ueber den
toroidalen Mittelwert der Querschnittsflaechen (remv bzw. aefv).
Ich hoffe, Du kommst damit zurecht. Im uebrigen nehme ich fuer den grossen Radius
den Fourierkoeffizienten R(0,0) am Plasmarand (eigentlich ns-1) um Verzerrungen
durch die Achsverschiebung zu vermeiden.
        """

        if(not self.import_all):
            self.importAll()

        # R approx. R00(ns-1)
        R = self.rmnc[-2][0]

        # aeff
        re00 = numpy.zeros(self.ns)
        re36 = numpy.zeros(self.ns)
        remv = numpy.zeros(self.ns)

        """
        # directly from Joachims fortan routine with falttened rmnc, zmns
        rmnc = self.rmnc.flatten()
        zmns = self.zmns.flatten()
        for j in range(self.ns):
            lj = j*self.mn_mode
            remv[j] = numpy.sum(rmnc[lj:lj+self.mn_mode]*zmns[lj:lj+self.mn_mode]*self.xm)
            for mn in range(self.mn_mode):
                l  = mn + lj
                im1 = self.xm[mn]
                in1 = self.xn[mn]
                for mn2 in range(self.mn_mode):
                    l2 = mn2 + lj
                    im2 = self.xm[mn2]
                    in2 = self.xn[mn2]
                    if(im1 == im2):
                        re00[j] = re00[j] + rmnc[l]*zmns[l2]*im1
                        ifac = 1-2*numpy.mod(numpy.abs(in1-in2),2)
                        #!--------------------------------------------
                        #!        -1   for in1-in2 odd
                        #! ifac =
                        #!        +1   for in1-in2 even
                        #!--------------------------------------------
                        re36[j] = re36[j] + rmnc[l]*zmns[l2]*im1*ifac
                    #endif
                #enddo
            #enddo
        #enddo
        remm = numpy.sqrt(numpy.abs( (re00 +re36) ) /2. )
        re00 = numpy.sqrt(numpy.abs( re00 ) )
        re36 = numpy.sqrt(numpy.abs( re36 ) )
        remv = numpy.sqrt(numpy.abs( remv ) )
        """

        """
        # use the 2d arrays
        for j in range(self.ns):
            remv[j] = numpy.sum(self.rmnc[j,:]*self.zmns[j,:]*self.xm)
            for mn in range(self.mn_mode):
                im1 = self.xm[mn]
                in1 = self.xn[mn]
                for mn2 in range(self.mn_mode):
                    im2 = self.xm[mn2]
                    in2 = self.xn[mn2]
                    if(im1 == im2):
                        re00[j] = re00[j] + self.rmnc[j,mn]*self.zmns[j,mn2]*im1
                        ifac = 1-2*numpy.mod(numpy.abs(in1-in2),2)
                        #!--------------------------------------------
                        #!        -1   for in1-in2 odd
                        #! ifac =
                        #!        +1   for in1-in2 even
                        #!--------------------------------------------
                        re36[j] = re36[j] + self.rmnc[j,mn]*self.zmns[j,mn2]*im1*ifac
                    #endif
                #enddo
            #enddo
        #enddo
        remm = numpy.sqrt(numpy.abs( (re00 +re36) ) /2. )
        re00 = numpy.sqrt(numpy.abs( re00 ) )
        re36 = numpy.sqrt(numpy.abs( re36 ) )
        remv = numpy.sqrt(numpy.abs( remv ) )
        """
        """
        aeff = remm[self.ns-1]
        aefv = remv[self.ns-1]
        """

        # we are interested just in aeff 
        j = -1
        re00 = 0.0
        re36 = 0.0
        for mn in range(self.mn_mode):
            im1 = self.xm[mn]
            in1 = self.xn[mn]
            for mn2 in range(self.mn_mode):
                im2 = self.xm[mn2]
                in2 = self.xn[mn2]
                if(im1 == im2):
                    re00 += self.rmnc[j,mn]*self.zmns[j,mn2]*im1
                    ifac = 1-2*numpy.mod(numpy.abs(in1-in2),2)
                    #!--------------------------------------------
                    #!        -1   for in1-in2 odd
                    #! ifac =
                    #!        +1   for in1-in2 even
                    #!--------------------------------------------
                    re36 += self.rmnc[j,mn]*self.zmns[j,mn2]*im1*ifac

        aeff = numpy.sqrt(numpy.abs( (re00 +re36) ) / 2. )

        return(R, aeff)



    def printForVMECRestart(self, eps=1e-2):
        """
        Print modes for magnetic axis and Rmn, Zmn of the last closed 
        flux surface for restarting VMEC.

        This output may be directly copied into VMEC input file

        INPUT:
          eps ... print modes sqrt(rmnc[-1]**2 + zmns[-1]**2)>eps

        RETURN:
          restparams  dictionary with parameters for VMEC restart
                      restparams = {'AxisParameters': string,
                                    'BoundaryParameters': string}
        """

        ax_str = []
        ax_str.append("!----- Axis Parameters -----\n")
        ax_str.append("RAXIS = %12.8f %12.6f %12.6f \n" 
                      % (tuple(self.ds.variables['rmnc'][0,0:3])))
        ax_str.append("ZAXIS = %12.8f %12.6f %12.6f \n" 
                      % (tuple(self.ds.variables['zmns'][0,0:3])))

        rmnc = self.ds.variables['rmnc'][-1][:]
        zmns = self.ds.variables['zmns'][-1][:]
        xm = self.ds.variables['xm'][:]
        xn = self.ds.variables['xn'][:]
        nfp = self.ds.variables['nfp'].getValue()

        w3 = numpy.where(numpy.sqrt(rmnc**2 + zmns**2)>eps)

        bd_str = []
        bd_str.append("!----- Boundary Parameters ----\n")
        for i in w3[0]:
            n = int(xn[i]/nfp)
            m = int(xm[i])
            bd_str.append("RBC(%d,%d) = %10.6f  ZBS(%d,%d) = %10.6f \n" 
                          % (n,m,rmnc[i], n,m,zmns[i]))

        print "".join(ax_str)[:-1]
        print "".join(bd_str)[:-1]

        restparams = {'AxisParameters': "".join(ax_str), 
                      'BoundaryParameters': "".join(bd_str)}

        return(restparams)

    ##### misc end #####



    ##### some stuff for inverting coordinates start #####

    def real2Mag(self, phi_plane, RZ, acc=1.0e-3, interp_params=None, debug=0):
        """
        simple version; deals with a real space (cylindrical) phi-plane.

        INPUT:
          phi_plane . scalar defining phi-plane in cylindrical coordinates
                      in rad
          RZ ........ R,Z coordinates in specified phi-plane
                      dimension: RZ[n,2]
          acc ....... accuracy in units of R,Z (most likely R,Z are in meters)
          interp_params  parameters for interpolation; distionary
                         e.g.:  default_params = {
            'st_res': [35j,100j], # grid in s,theta space to xform to cart
            'rz_res': [130j, 130j], # the uniform grid in cylindrical space.
            'extent': [1.0,1.38,-.25,.25] # in r,z space rmin,rmax,zmin,zmax
            }
            In case interp_params is not provided by the user the values 
            for 'extent' are determined automatically.

        RETURN:
          stp ... array containing s, theta, phi
                  dimension:  stp[n,3]
        """
        if(debug != 0):
            print "real2Mag(): inversion on phi={0:7.2f}deg plane".format(numpy.rad2deg(phi_plane))
            print "interp_params ", interp_params
            print "RZ ", RZ

        RZ = numpy.array(RZ)
        dim = RZ.shape
        if(len(dim) == 1):
            no_points = 1
            if(dim[0] != 2):
                raise InputDataError("".join(["For a single point only ",
                                              "one pair (RZ) is allowed!"]))
            RZ = RZ.reshape([n,dim[0]])
        elif(len(dim) == 2):
            no_points = dim[0]
        else:
            raise InputDataError("".join(["Input coordinates (R,Z) must ",
                                          "have dimension [n,2]!"]))

        default_params = {        # these do error <=1e-4, for cubic interp
            'st_res': [35j,100j], # grid in s,theta space to xform to cart
            'rz_res': [130j, 130j], # the uniform grid in Cartesian space.
            'extent': [1.0,1.38,-.25,.25] # in r,z space rmin,rmax,zmin,zmax
            }
        if (debug>1): 
            print('interp_params (default)\n %s' % (default_params))

        if(interp_params == None):
            int_pms = deepcopy(default_params)
        else:
            int_pms = deepcopy(interp_params)

        if (debug>1): 
            print('int_pms (input)\n %s' % (int_pms))

        # initialise array with results, stp, with nan
        stp = numpy.empty([no_points,3])
        stp.fill(numpy.nan)

        xandz = [0,2]

        # make an approximate grid for the Cartesian -> Bean transform.
        # normally would only do this once, and save somehow
        # make the sgrid go a little above s=1 so that interpolator
        # has enough margin.  Need to be sure VMEC fourier series doesn;t
        # diverge too much.

        # Try with two grids in theta -> pacman interpolation
        smin = 0.0
        smax = 1.08
        tmin = -PI
        tmax = PI
        (grid_s,grid_theta) = numpy.mgrid[smin:smax:int_pms['st_res'][0], 
                                          tmin:tmax:int_pms['st_res'][1]]
        # this improves the mesh shape near origin
        grid_s = (numpy.array((grid_s)))**2


        # first check that grid is OK - values are imaginary
        #  easiest by making sure grid is not ridiculously coarse
        # need to check just one as they are all the same shape
        if (numpy.shape(grid_s)[0] < 10) or (numpy.shape(grid_s)[1] < 10):
            self.warn_but_not_too_often(
                "st grid small %dx%d! don't forget js in st_res" %
                (numpy.shape(grid_s)))
            if numpy.size(grid_s)<10: raise LookupError('grid_s too small')

        # adjust limits for R and Z using values computed with s=1
        # in case user did not provide interp_params
        if(interp_params == None):
            RpZ = self.getCrossSectData(phi_plane, 1, 
                                        theta_arr=grid_theta[0][:],
                                        coordsys="cyl")
            lim = []
            rmin = RpZ[0,0,:].min()
            rmax = RpZ[0,0,:].max()
            zmin = RpZ[0,2,:].min()
            zmax = RpZ[0,2,:].max()
            lim.append(rmin - 0.1*numpy.abs(rmin))
            lim.append(rmax + 0.1*numpy.abs(rmax))
            lim.append(zmin - 0.1*numpy.abs(zmin))
            lim.append(zmax + 0.1*numpy.abs(zmax))
            int_pms.update({'extent':lim})

        if (debug>1): 
            print('int_pms (actual)\n %s' % (int_pms))

        r_on_bean = numpy.zeros(grid_s.shape)
        z_on_bean = numpy.zeros(grid_s.shape)

        for (i,s) in enumerate(grid_s[:,0]):
            RpZ = self.getCrossSectData(phi_plane, s, 
                                        theta_arr=grid_theta[i][:],
                                        coordsys="cyl")
            r_on_bean[i][:] = RpZ[0,0,:]
            z_on_bean[i][:] = RpZ[0,2,:]


        ext = int_pms['extent']
        print 'ext=', ext
        rz_res = int_pms['rz_res']
        (grid_r,grid_z) = numpy.mgrid[ext[0]:ext[1]:rz_res[0], 
                                      ext[2]:ext[3]:rz_res[1]]

        # check grid as for s and r earlier
        if (numpy.shape(grid_r)[0] < 10) or (numpy.shape(grid_r)[1] < 10):
            self.warn_but_not_too_often(
                "rz grid small %dx%d! don't forget 'j's in st_res" %
                (numpy.shape(grid_r)))
            if numpy.size(grid_r)<10: raise LookupError('grid_r too small')

        rsp = []  # list with indices of points surviving the following check
        for i in range(no_points):
            ir = RZ[i,0]
            iz = RZ[i,1]
            goodp = True
            if not miscmath.contained_in(grid_r, ir):
                print "point ", ir, iz
                print "not in R range: ", min(grid_r.flatten()), max(grid_r.flatten())
                goodp = False
                #raise ValueError('R value outside grid ')
                
            if not miscmath.contained_in(grid_z, iz):
                print "point ", ir, iz
                print "not in Z range: ", min(grid_z.flatten()), max(grid_z.flatten())
                goodp = False
                #raise ValueError('Z value outside grid ')
            if(goodp):
                rsp.append(i)

        #if (debug>0): 
        print "{0:d} out of {1:d} points survived".format(len(rsp), no_points)
            #print "interpolation method=", interp_method

        if version >= 0.9:  
            cart_grid_s =     griddata((numpy.array(r_on_bean).flatten(),
                                        numpy.array(z_on_bean).flatten()), 
                                       grid_s.flatten(), 
                                       (grid_r, grid_z),
                                       method=interp_method)#, fill_value=1)
            cart_grid_theta1 = griddata((numpy.array(r_on_bean).flatten(),
                                        numpy.array(z_on_bean).flatten()), 
                                       (numpy.sqrt(grid_s)*numpy.sin(grid_theta)).flatten(), 
                                       (grid_r, grid_z),
                                       method=interp_method)#,fill_value=999)
            cart_grid_theta2 = griddata((numpy.array(r_on_bean).flatten(),
                                        numpy.array(z_on_bean).flatten()), 
                                       (numpy.sqrt(grid_s)*numpy.cos(grid_theta)).flatten(), 
                                       (grid_r, grid_z),
                                       method=interp_method)#,fill_value=999)
        else:
            cart_grid_s =     griddata(numpy.array(r_on_bean).flatten(), 
                                       numpy.array(z_on_bean).flatten(), 
                                       grid_s.flatten(), 
                                       grid_r, grid_z, 
                                       interp='linear').T#,fill_value=1)
            cart_grid_theta1 = griddata(numpy.array(r_on_bean).flatten(), 
                                       numpy.array(z_on_bean).flatten(), 
                                       numpy.sin(grid_theta.flatten()), 
                                       grid_r, grid_z, 
                                       interp='linear').T
            cart_grid_theta2 = griddata(numpy.array(r_on_bean).flatten(), 
                                       numpy.array(z_on_bean).flatten(), 
                                       numpy.sin(grid_theta.flatten()), 
                                       grid_r, grid_z, 
                                       interp='linear').T

        if(debug>0):
            print "cart grids computed"
        # last closed magnetic surface
        lcms = self.getCrossSectData(phi_plane, 1, no_theta=100, coordsys="cyl")
        # magnetic axis
        ma = self.getCrossSectData(phi_plane, 0, no_theta=1, coordsys="cyl")
        
        if(debug>10): # some plots for checking
            dname = "TestInversionPlots/"
            try:
                os.mkdir(dname)
            except:
                pass
            fig = pylab.figure()
            ax = fig.add_subplot(111)
            ax.set_xlabel(r'$R$')
            ax.set_ylabel(r'$Z$')
            titlestr = "".join(["$\\phi={0:7.2f}^\circ;$ ".format(numpy.rad2deg(phi_plane)),
                                "$no\\_points={0:d}$".format(no_points)])
            ax.set_title(titlestr)
            ax.imshow(cart_grid_s, extent=[min(grid_r.flatten()),max(grid_r.flatten()),min(grid_z.flatten()),max(grid_z.flatten())])
            ax.plot(r_on_bean, z_on_bean,',',color='gray')
            ax.plot(grid_r, grid_z,'b',linewidth='.2')
            ax.plot(grid_r, grid_z,'b',linewidth='.2')
            ax.plot(lcms[0,0,:],lcms[0,2,:], 'k-', linewidth=1)
            ax.plot(ma[0,0,:],ma[0,2,:], 'k+', markersize=10, markeredgewidth=1)

            ax.set_aspect('equal')
            fig.show()
            fname="{0:s}points_p{1:+04.0f}np{2:03d}.pdf".format(dname,numpy.rint(numpy.rad2deg(phi_plane)), no_points)

            fig.savefig(fname)


            fig1 = pylab.figure()
            try:
                # try new way. for matplotlib> v1.0.0 
                ax1 = fig1.add_subplot(111, projection='3d')
            except:
                # ok, we need the old way
                ax1 = p3.Axes3D(fig1)

            ax1.set_xlabel(r'$R$')
            ax1.set_ylabel(r'$Z$')
            ax1.set_zlabel(r'$s$')
            titlestr = "".join(["$\\phi={0:7.2f}^\circ;$ ".format(numpy.rad2deg(phi_plane)),
                                "$no\\_points={0:d}$".format(no_points)])
            ax1.set_title(titlestr)
            ax1.plot3D(grid_r.flatten(), grid_z.flatten(), cart_grid_s.flatten(), 'r.', markersize=0.5)
            ax1.contour(grid_r, grid_z, cart_grid_s, 50) 
            ax1.plot3D(lcms[0,0,:],lcms[0,2,:],numpy.ones(lcms.shape[2]), 'k-', linewidth=5)
            ax1.plot3D(ma[0,0,:],ma[0,2,:],numpy.zeros(ma.shape[2]), 'k+', markersize=10, markeredgewidth=2)
            fig1.show()
            fname="{0:s}s_3d_p{1:+04.0f}np{2:03d}.pdf".format(dname,numpy.rint(numpy.rad2deg(phi_plane)), no_points)

            fig1.savefig(fname)

            #""" have to check for 3d - not working this way on h1svr
            fig2 = pylab.figure()
            try:
                # try new way. for matplotlib> v1.0.0 
                ax2 = fig2.add_subplot(111, projection='3d')
            except:
                # ok, we need the old way
                ax2 = p3.Axes3D(fig2)

            ax2.set_xlabel(r'$R$')
            ax2.set_ylabel(r'$Z$')
            ax2.set_zlabel(r'$\theta$')
            titlestr = "".join(["$\\phi={0:7.2f}^\circ;$ ".format(numpy.rad2deg(phi_plane)),
                                "$no\\_points={0:d}$".format(no_points)])
            ax2.set_title(titlestr)
            ax2.plot3D(grid_r.flatten(), grid_z.flatten(), cart_grid_theta1.flatten(), 'r.', markersize=0.5)
            ax2.contourf(grid_r, grid_z, cart_grid_theta1, 50) 
            ax2.plot3D(lcms[0,0,:],lcms[0,2,:],numpy.ones(lcms.shape[2]), 'k-', linewidth=5)
            ax2.plot3D(ma[0,0,:],ma[0,2,:],numpy.ones(ma.shape[2]), 'k+', markersize=10, markeredgewidth=2)
            fig2.show()
            fname="{0:s}t1_3d_p{1:+04.0f}np{2:03d}.pdf".format(dname,numpy.rint(numpy.rad2deg(phi_plane)), no_points)
            fig2.savefig(fname)

            fig2 = pylab.figure()
            try:
                # try new way. for matplotlib> v1.0.0 
                ax2 = fig2.add_subplot(111, projection='3d')
            except:
                # ok, we need the old way
                ax2 = p3.Axes3D(fig2)

            ax2.set_xlabel(r'$R$')
            ax2.set_ylabel(r'$Z$')
            ax2.set_zlabel(r'$\theta$')
            titlestr = "".join(["$\\phi={0:7.2f}^\circ;$ ".format(numpy.rad2deg(phi_plane)),
                                "$no\\_points={0:d}$".format(no_points)])
            ax2.set_title(titlestr)
            ax2.plot3D(grid_r.flatten(), grid_z.flatten(), cart_grid_theta2.flatten(), 'r.', markersize=0.5)
            ax2.contourf(grid_r, grid_z, cart_grid_theta2, 50) 
            ax2.plot3D(lcms[0,0,:],lcms[0,2,:],numpy.ones(lcms.shape[2]), 'k-', linewidth=5)
            ax2.plot3D(ma[0,0,:],ma[0,2,:],numpy.ones(ma.shape[2]), 'k+', markersize=10, markeredgewidth=2)
            fig2.show()
            fname="{0:s}t2_3d_p{1:+04.0f}np{2:03d}.pdf".format(dname,numpy.rint(numpy.rad2deg(phi_plane)), no_points)
            fig2.savefig(fname)
        # end   if(debug>10): # some plots for checking


        for i in range(no_points):
            if(i in rsp):
                use_grid2 = False
                goodp = True
                try:
                    stp[i,0] = interp2d(RZ[i],cart_grid_s, grids=[grid_r,grid_z])
                except:
                    goodp = False
                    print "s: interp2d() causes problem for: i, RZ ", i, RZ[i]
                    #raise
                try:
                    t1 = interp2d(RZ[i],cart_grid_theta1, grids=[grid_r,grid_z])
                except:
                    goodp = False
                    print "t1: interp2d() causes problem for: i, RZ ", i, RZ[i]
                    #raise
                #stp[i,2] = interp2d(RZ[i],cart_grid_phi, grids=[grid_r,grid_z])

                try:
                    t2 = interp2d(RZ[i],cart_grid_theta2, grids=[grid_r,grid_z])
                except:
                    goodp = False
                    print "t2: interp2d() causes problem for: i, RZ ", i, RZ[i]
                    #raise

                if(goodp):
                    stp[i,1] = numpy.arctan2(t1,t2)
                    stp[i,2] = phi_plane

                    # final check for accuracy
                    # if accuracy is not reached, return "nan"
                    # check also range for s. due to spline extrapolation might be s>1
                    RZt = self.fct_RZ(stp[i,0], stp[i,1], phi_plane)
                    RZt = numpy.array(RZt)
                    dist = numpy.sqrt(numpy.sum((RZ[i]-RZt)**2))
                    if((dist > acc) or (stp[i,0] > 1) or 
                       numpy.isnan(stp[i,0]) or numpy.isnan(stp[i,1])):
                        if (debug>0):
                            print "don't keep point no {0:3d} ({1:12.6e},{2:12.6e}): dist={3:12.6e}, s={4:12.6e}".format( i, RZ[i,0], RZ[i,1], dist, stp[i,0])
                        stp[i,:].fill(numpy.nan)
                else:
                    stp[i,:].fill(numpy.nan)

        if(debug == 0):
            return(stp)
        else:
            return(stp, grid_r, grid_z, cart_grid_s, cart_grid_theta1, cart_grid_theta2, r_on_bean, z_on_bean)



    ##### some stuff for inverting coordinates end #####

