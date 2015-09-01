#! /usr/bin/env python

__doc__ = \
"""
A class for 1-D splines as it is usefull for / imported by 
the classes VMEC and BOOZER.
This class deals with NetCDF4 input/output.

History:

14 June 2012 bhs
* Set 'controlflag' ext to default (ext=0) for splines used 
  from scipy.interpolate. For ext=0 some extrapolation is done. For ext=1
  one may have troubles computing values for the upper boundary... (values=0).

05 June 2012 bhs
* specifiers for spline types are global

16 May 2012 bhs
* start with this script
"""

"""

ncpath,ncfname = os.path.split(self.nc_filename)
scfilename_base, scfilename_ext = os.path.splitext(ncfname)

nc_splcoeffilename = "".join([scfilename_base, '.splcoef', scfilename_ext])
not suitable - try the following...
nc_splcoeffilename = "".join([ncpath, '/', nc_splcoeffilename])

"""


__version__ = "0.1"
__versionTime__ = "16 May 2012 09:50"
__author__ = "Bernhard Seiwald <bernhard.seiwald@gmail.com>"


import numpy
from scipy import interpolate
import convenience 
import netCDF4Tools


try:
    import netCDF4 as netCDF4
    netcdf_module = 'netCDF4'
except: raise Exception, 'Unable to load netCDF4 module/library'


#Commented out for now SRH - problem loading them into Python from the .so 
#from spline_spec import leadingfunction_mod as leadingfunction
#from spline_spec import spline_spec_mod as spline_spec

import collections as collections
if (not collections.__dict__.has_key('OrderedDict')):
    # This is some workaround for h1svr in case there is no
    # collections.OrderedDict() available.
    # The package ordereddict may be found in:
    # bhs112@h1svr:~$/home/bhs112/usr/local/
    try:
        import ordereddict as collections
    except:
        print "OrderedDict not available. Bye..."


debug=0

#WRITE_PROGRESS=True
WRITE_PROGRESS=False


__all__ = ['SPLINE']

# define some constants
# Specifiers for spline methods
# Use powers of 2, so several methods may be accessed at once.
# ATTENTION: When new specifiers are added, the names have(!) 
# to end with _SPLINE !!! In some routines we trigger to this!!!
# Nesd to define here, so after importing SPLINE these constants 
# can be accessed.
CONV_SPLINE = 1
SQRT_SPLINE = 2
SPECIAL_SPLINE = 4


class SPLINE():
    """
    A class for 1-D splines as it is usefull for / imported by 
    the classes VMEC and BOOZER.
    This class deals with NetCDF4 input/output.

    Supported kinds of splines:
    
    CONV_SPLINE: Ordinary 3rd order cubic spline.

    SQRT_SPLINE: Divide odd m terms by sqrt(s), keep even m terms as they are.
                 Use 3rd order qubic spline.

    SPECIAL_SPLINE: Use special spline equiped with leading function (takes 
                 into account sqare root characteristics of Rmnc, Zmns).
                 Make sure that spline_spec is imported!


    ATTENTION: The format for the dictionary data must be as follows:
               data = {'x', x[n], 'y':y[n,m], 
                       'ixm':ixm[m], 'kindspec':specifier}
               with: 
               n = no points;  m = no curves
               ixm = Poloidal mode numbers
               kindspec =  specifier being one of: 
                  CONV_SPLINE=1
                  SQRT_SPLINE=2
                  SPECIAL_SPLINE=4

    Example usage:

    import numpy
    import pylab
    import SPLINE
    import netCDF4Tools as nc4T

    # generate some data
    x=numpy.linspace(0, 2*numpy.pi, 7)
    xin=numpy.linspace(0, 2*numpy.pi, 14)
    y1=numpy.sin(x)
    y2=numpy.sin(2*x)
    y3=numpy.sin(3*x)

    ixm=numpy.zeros(3)

    spl=SPLINE.SPLINE({'x':x, 'y':numpy.array([y1,y2,y3]).T,'kindspec':'SQRT_SPLINE', 'ixm':ixm})
    #or
    spl=SPLINE.SPLINE({'x':x,'y':numpy.array([y1,y2,y3]).T,'kindspec':5,'ixm':ixm})
    # or
    spl=SPLINE.SPLINE({'x':x,'y':numpy.array([y1,y2,y3]).T,'kindspec':"SPECIAL_SPLINE CONV_SPLINE",'ixm':ixm})

    yy2=spl2.eval(xin)
    pylab.plot(x,y1,'b+-', xin,yy2[:,0],'r+-'); pylab.show()

    spl2.save('testspline')

    spl3=SPLINE.SPLINE(nc_name='testspline', import_all=True)
    yy3=spl3.eval(xin)
    pylab.plot(x,y1,'b+-', xin,yy3[:,0],'r+-'); pylab.show()

    --- end examples ---

    """

    def __init__(self, data=None, nc_name=None, show_vars=False):
        """
        INPUT:
          data ... distionary containing data x, y(x) to spline and 
                   the kind specifier.
          nc_name  netCDF4 file name for data to load

        RETURN:
          spline class

        ATTENTION: The format for the dictionary data must be as follows:
                   data = {'x', x[n], 'y':y[n,m], 
                           'ixm':ixm[m], 'kindspec':specifier}
                   with: 
                   n = no points;  m = no curves
                   ixm = Poloidal mode numbers
                   kindspec =  specifier being one of: 
                      CONV_SPLINE=1
                      SQRT_SPLINE=2
                      SPECIAL_SPLINE=4
        
        Note: For kindspec one may use either the numbers or directly 
              the names. Also combining several spline methods
              is possible. 
              Examples what is accepted:
              'kindspec': 2                  equal to
              'kindspec': "SQRT_SPLINE"
              'kindspec': 5                  equal to
              'kindspec':"SPECIAL_SPLINE | CONV_SPLINE"
          
        """
        # define some constants
        # Specifiers for spline methods
        # Use powers of 2, so several methods may be accessed at once.
        # ATTENTION: When new specifiers are added, the names have(!) 
        # to end with _SPLINE !!! In some routines we trigger to this!!!
        self.CONV_SPLINE = CONV_SPLINE
        self.SQRT_SPLINE = SQRT_SPLINE
        self.SPECIAL_SPLINE = SPECIAL_SPLINE

        self.__all_imported__ = False
        self.__kind_computed__ = 0
        self.__qNames__ = []

        self.init_SPLINE_METHODS()

        if(data != None):
            if(type(data) != dict): 
                raise Exception, "data must be dictionary"
        
        if(nc_name != None):
            if(type(nc_name) != str): 
                raise Exception, "nc_name must be string"
        

        if([data, nc_name].count(None) != 1):
            raise Exception, "".join(["Either data for calculating spline ",
                                      "coefficients\n ",
                                      "or a filename for loading ",
                                      "netCDF4 data has to be provided!"])

        if(nc_name != None):
            self.load(nc_name, show_vars)
            return
        
        yshape = data['y'].shape
        ydim = len(yshape)
        if(not ydim in (1,2)):
            raise Exception, "".join(["Dimension for y={0:d}".format(ydim),
                                      "not allowed. Must be either 1 or 2"])

        if(not data.has_key('kindspec')):
            data.update({'kindspec':self.CONV_SPLINE})

        if(type(data['kindspec']) not in [str,int]):
             raise Exception, "".join(["kindspec must be type of ",
                                       "either string or integer"])

        if(type(data['kindspec']) == str):
            kspec = self.kindspecStr2Num(data['kindspec'])
            data.update({'kindspec':kspec})

        if((not data.has_key('ixm')) 
           and (data['kindspec'] != self.CONV_SPLINE)):
            raise Exception, "".join(["For any other kind of spline than ",
                                      "CONV_SPLINE (=cubic spline) ",
                                      "ixm has to be provided!"])

        if(not data.has_key('ixm')):
            if(ydim == 1):
                no = 1
            else:
                no = yshape[1]
            data.update({'ixm':numpy.zeros(no)})
            
        self.importQuantity('x', data['x'])
        self.importQuantity('y', data['y'])
        self.importQuantity('ixm', data['ixm'])
        self.importQuantity('kindspec', data['kindspec'])

        self.setHistory()
        self.setSource()
        self.setComment()

        if(data != None):
            self.calc()

        return


    def init_SPLINE_METHODS(self):
        """
        Look for definition of spline methods self.*_SPLINE 
        and add them to the dictionary self.SPLINE_METHODS_d
        and add values to the list self.SPLINE_METHODS_l.

        INPUT:
          None

        RETURN:
          None
        """
        d = dir(self)
        self.SPLINE_METHODS_d = collections.OrderedDict()
        self.SPLINE_METHODS_l = []
        for dd in d:
            if(dd[-7:]=='_SPLINE'):
                exec('val=self.{0:s}'.format(dd))
                self.SPLINE_METHODS_d.update({dd:val})
                self.SPLINE_METHODS_l.append(val)

        self.SPLINE_METHODS_d = collections.OrderedDict(
            sorted(self.SPLINE_METHODS_d.items(), key=lambda x: x[1]))
        self.SPLINE_METHODS_l.sort()
        return


    def kindspecStr2Num(self, kindspec):
        """
        Taking into account a string, kindspec, specifying the 
        kind of spline, this routine return the corresponding number.

        E.g.
        kindspec = "SPECIAL_SPLINE | CONV_SPLINE"
        will lead to
        self.SPECIAL_SPLINE + self.CONV_SPLINE -> 5

        

        INPUT:
          kindspec ... string specifying the kind of spline
        RETURN:
          kspec ...... integer value representing the kind of spline
        """
        kstrl = kindspec.replace(' ', '|')
        kstrl = kstrl.split('|')
        kspec = 0
        for kstr in kstrl:
            ks = kstr.strip()
            if(ks in self.SPLINE_METHODS_d.keys()):
                kspec |= self.SPLINE_METHODS_d[ks]

        return(kspec)


    def setHistory(self, history=None):
        """
        Set the history of the spline.

        INPUT:
          history ... a string telling something about the history

        RETURN:
          None
        """
        self.history = history
        return


    def setSource(self, data_source=None):
        """
        Set the data source of the spline. 
        This is some kind of additional comment.

        The string for data_source might be a configuration specifier, 
        e.g.: VMEC_fb_kh00.000

        INPUT:
          data_source ... a string telling something about the data source

        RETURN:
          None
        """
        self.data_source = data_source
        return


    def setComment(self, comment=None):
        """
        Set some comment.

        INPUT:
          comment ... a string with the comment

        RETURN:
          None
        """
        self.comment = comment
        return



    def importQuantity(self, qName, qData):
        """
        Import a quantity qName into class.
        This quantity is lateron accessed as: self.qName

        INPUT:
          qName ... name of quantity
          qData ... value to be added

        RETURN:
          None
        """
        exestr = "self.{0:s} = qData".format(qName)
        if(WRITE_PROGRESS):   
            print "importQuantity(): exec({0:s})".format(exestr)
        exec(exestr)
        if(qName not in self.__qNames__):
            self.__qNames__.append(qName)
        return


    def importAll(self):
        """

        Import all variables.
        Set self.__all_imported__ = True

        """

        for (skey, svalue) in self.ds.iteritems():
            cmdstr = "self.{0:s} = svalue".format(skey)
            if(WRITE_PROGRESS):   print cmdstr
            exec(cmdstr)
            if(skey not in self.__qNames__):
                self.__qNames__.append(skey)


        # now check for splines from scipy.interpolate 
        # and restore tck-structure
        if("tck__spl_coefs" in self.__qNames__):
            if(WRITE_PROGRESS):   print 'tck__spl_coefs found'
            if("tck__spl_order" in self.__qNames__):
                if(WRITE_PROGRESS):   print 'tck__spl_order found'
                self.tck = self.unpackTCK(self.tck__spl_coefs, 
                                          self.tck__spl_order)
                # remove the unnecessary ...__spl_...
                del self.tck__spl_coefs
                del self.tck__spl_order
            else:
                raise Exception, "".join(["tck-structure for CONV_SPLINE ",
                                          "cannot be restored"])

        
        if("sqrt_tck__spl_coefs" in self.__qNames__):
            if(WRITE_PROGRESS):   print 'sqrt_tck__spl_coefs found'
            if("sqrt_tck__spl_order" in self.__qNames__):
                if(WRITE_PROGRESS):   print 'sqrt_tck__spl_order found'
                self.sqrt_tck = self.unpackTCK(self.sqrt_tck__spl_coefs, 
                                               self.sqrt_tck__spl_order)
                # remove the unnecessary ...__spl_...
                del self.sqrt_tck__spl_coefs
                del self.sqrt_tck__spl_order
            else:
                raise Exception, "".join(["tck-structure for SQRT_SPLINE ",
                                          "cannot be restored"])

        self.__all_imported__ = True
        return


    def showVariables(self, max_width=80, indent='  '):
        """
        """
        __version__ = '1.0'

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


    def calc(self):
        """
        Compute spline coefficients.

        INPUT:
          None

        RETURN:
          None

        ATTENTION: The format for the dictionary data must be as follows:
                   data = {'x', x[n], 'y':y[n,m], 'kind':specifier}
                   with: 
                   * n = no points;  m = no curves
                   * kind specifier being one of: CONV_SPLINE=1,
                     SQRT_SPLINE=2, SPECIAL_SPLINE=4
        """
        if(self.kindspec == None):   
            self.kindspec = self.CONV_SPLINE

        if(type(self.kindspec) == str):
            self.kindspec = self.kindspecStr2Num(kindspec)

        if((self.ixm == None) and (self.kindspec != self.CONV_SPLINE)):
            raise Exception, "".join(["For any other kind of spline than ",
                                      "CONV_SPLINE (=cubic spline) ",
                                      "ixm has to be provided!"])

        calc_flag = False
        if((self.kindspec & self.CONV_SPLINE) == self.CONV_SPLINE):
            print "computing coefficients for conventional 3rd order spline..."
            self.calc_spline_conv_coefs()
            calc_flag = True
        if((self.kindspec & self.SQRT_SPLINE) == self.SQRT_SPLINE):
            print "".join(["computing coefficients for 3rd order spline ",
                           "taking into account for ",
                           "odd modes mode/sqrt(s)..."])
            self.calc_spline_sqrt_coefs()
            calc_flag = True
        if((self.kindspec & self.SPECIAL_SPLINE) == self.SPECIAL_SPLINE):
            print "computing coefficients for special 3rd order spline..."
            self.calc_spline_spec_coefs()
            calc_flag = True

        if(not calc_flag):
            raise Exception, \
                "Unsupported value {0:d} for kind specifier.".format(self.kindspec)
        
        return


    def eval(self, xin, deriv=False, kindspec=None):
        """
        Evaluate spline at the desired point xin.

        INPUT:
          xin ...... evaluate spline at this point.
          deriv .... True/False compute derivatives up to 3rd order.
          kindspec . kind of spline to use.

        RETURN:
          y[,yp,ypp,yppp] value at point xin and, optional, derivatives
                          up to 3rd order.
        """

        if(self.__kind_computed__ == 0):
            raise Exception, "".join(["No spline coefficients ",
                                      "have been computed so far..."])

        # if kindspec is not provided, find the simplest spline...
        if(kindspec == None):
            for k in self.SPLINE_METHODS_l:
                if((self.__kind_computed__ & k) == k):
                    kindspec = k
                    break

        if(not (kindspec in self.SPLINE_METHODS_l)):
            raise Exception, \
                "Unsupported value {0:d} for kind specifier.".format(kindspec)

        for (key,value) in self.SPLINE_METHODS_d.items():
            if(value == kindspec):
                break

        if((self.__kind_computed__ & kindspec) != kindspec):
            raise Exception, "".join(["No spline coefficients for the ",
                                      "requested type ", key, " ", 
                                      "have been computed so far..."])

        xintype = type(xin)
        if(xintype in [list, tuple]):
            xin = numpy.array(xin)
        elif(xintype is not numpy.ndarray):  # assuming it is a scalar...
            xin = numpy.array([xin])
        
        ret = []
        for xi in xin:
            if(kindspec == self.CONV_SPLINE):
                ret.append(self.splev_conv(xi, deriv))
            elif(kindspec == self.SQRT_SPLINE):
                ret.append(self.splev_sqrt(xi, deriv))
            elif(kindspec == self.SPECIAL_SPLINE):
                ret.append(self.splev_spec(xi, deriv))
            else:
                raise Exception, \
                    "".join(["Unsupported value {0:d}".format(kindspec),
                             " for kind specifier."])

        return(numpy.array(ret))



    def load(self, nc_name, show_vars=False):
        """
        Load input data and computed spline coefficients from a netCDF4 file.

        INPUT:
          nc_name ... name of netCDF4 file
          show_vars . show variables

        RETURN:
          None
        """
        if(nc_name[-3:] != '.nc'):
            nc_name = "".join([nc_name, '.nc'])

        print "loading netCDF4 file " + nc_name + "..."
        #self.ds = netCDF4.Dataset(nc_name, mode='r')
        self.ds = netCDF4Tools.readNetCDF4(nc_name)

        # import all variables
        self.importAll()

        if(show_vars): self.showVariables()

        return


    def save(self, nc_name):
        """
        Save input data and computed spline coefficients to a netCDF4 file.

        INPUT:
          nc_name ... name of netCDF4 file

        RETURN:
          None
        """

        # prepare data
        sdata = {}
        sdata.update({'__kind_computed__':self.__kind_computed__})
        sdata.update({'history':self.history})
        sdata.update({'data_source':self.data_source})
        sdata.update({'comment':self.comment})

        # we like to store all data named in self.__qNames__
        # append these data to sdata
        for sitem in self.__qNames__:
            if(sitem in ["tck", "sqrt_tck"]):
                cmdstr = "(spl_coefs, spl_order) = self.packTCK(self.{0:s})".format(sitem)
                if(WRITE_PROGRESS):   print cmdstr
                exec(cmdstr)
                cmdstr = "sdata.update({{'{0:s}__spl_coefs':spl_coefs}})".format(sitem)
                if(WRITE_PROGRESS):   print cmdstr
                exec(cmdstr)
                cmdstr = "sdata.update({{'{0:s}__spl_order':spl_order}})".format(sitem)
                if(WRITE_PROGRESS):   print cmdstr
                exec(cmdstr)
            else:
                cmdstr = "sdata.update({{'{0:s}':self.{0:s}}})".format(sitem)
                if(WRITE_PROGRESS):   print cmdstr
                exec(cmdstr)

        # check the whole stuff in sdata for None's and replace them by "None"
        for key,value in sdata.items():
            if(value is None):
                if(WRITE_PROGRESS):
                    print "".join(["sdata.update({", key, ": 'None'})"])
                sdata.update({key: 'None'})

        # write data to netCDF4 file
        netCDF4Tools.writeNetCDF4(sdata, nc_name)

        return


    def calc_spline_conv_coefs(self):
        """
        Calculate cubic spline coefficients for a set of curves.

        INPUT:
          None

        RETURN:
          None

        For evaluation/interpolation: 
         Make sure that splev is imported from scipy.interpolate.
         Call like:      ynew = splev(xnew,tck,der=0)
         or for derivs   yder = splev(xnew,tck,der=1).
        """

        # do some checks
        if(self.y.shape[0] != self.x.shape[0]):
            raise Exception, "".join(["dimensions do not agree for: y[n,m]",
                                      "and x[n]"])

        self.tck=[]
        if(len(self.y.shape) == 1):
            self.tck.append(interpolate.splrep(self.x, self.y, s=0))
        else:
            n = self.y.shape[1]
            for i in range(n):
                self.tck.append(interpolate.splrep(self.x, self.y[:,i], s=0))

        self.__kind_computed__ |= self.CONV_SPLINE
        self.importQuantity('tck', self.tck)
        self.importQuantity('__kind_computed__', self.__kind_computed__)

        return


    def calc_spline_sqrt_coefs(self):
        """
        Calculate cubic spline coefficients for a set of curves.

        INPUT:
          None

        RETURN:
          None

        For further usage compute spline coefficients for
        rmnc_b(s), zmns_b(s), pmns_b(s).

        Divide odd m terms by sqrt(s), 
        keep even m terms as they are.
        Use 3rd order qubic spline.

        ATTENTION: For interpolation keep in mind to 're-do' sqrt(s)!

        For evaluation/interpolation: 
         Use splev_sqrt(). This routine calls splev() which is
         to be found in from scipy.interpolate.
         Make sure that splev is imported from scipy.interpolate.
         Call like:      ynew = splev_sqrt(xnew,sqrt_tck,der=0)
         or for derivs   yder = splev_sqrt(xnew,sqrt_tck,der=1).
        """

        # as we shall not divide by zero...
        eps = 1.e-30
        self.x = numpy.where(self.x >= eps, self.x, eps)

        # do some checks
        if(self.y.shape[0] != self.x.shape[0]):
            raise Exception, "".join(["dimensions do not agree for: y[n,m]",
                                      "and x[n]"])
        if(self.y.shape[1] != self.ixm.shape[0]):
            raise Exception, "".join(["dimensions do not agree for: y[n,m]",
                                      "and ixm[m]"])

        self.sqrt_tck=[]
        if(len(self.y.shape) == 1):
            self.sqrt_tck.append(interpolate.splrep(self.x, self.y, s=0))
        else:
            sqx = numpy.sqrt(self.x)
            n = self.y.shape[1]
            for i in range(n):
                if(self.ixm[i]%2 == 0):
                    self.sqrt_tck.append(interpolate.splrep(self.x, 
                                                            self.y[:,i], 
                                                            s=0))
                else:
                    self.sqrt_tck.append(interpolate.splrep(self.x, 
                                                            self.y[:,i]/sqx, 
                                                            s=0))

        self.__kind_computed__ |= self.SQRT_SPLINE
        self.importQuantity('sqrt_tck', self.sqrt_tck)
        self.importQuantity('__kind_computed__', self.__kind_computed__)

        return


    def calc_spline_spec_coefs(self):
        """
        Calculate for a set of curves coefficients for a spline
        which is equipped with a leading function (also noted as 
        testfunction or leading term) lf(x).


        INPUT:
          None

        RETURN:
          None

        For further usage compute spline coefficients for
        rmnc_b(s), zmns_b(s), pmns_b(s).

        This spline is equiped with a special leading function lf(s) 
        in order to better fit the characteristics of the coefficients 
        rmnc_b(s) and zmns_b(s) close to the magnetic axis, which is s=0.

        The characteristics of the coefficients is a square root 
        characteristics. To fit this, the qubic spline polynomial
        P^3(s) is modified as follows:  
                P^3(s)   ->    lf(s)*P^3(s)
        with   
                lf(s) = s^(m/2)   
        with m beeing the poloidal modenumber

        For evaluation/interpolation: 
         Make sure that spline_spec is imported.
         Call like:      ynew = splev(xnew,tck,der=0)
         or for derivs   yder = splev(xnew,tck,der=1).
        """

        c1 = 0.0
        cn = 0.0
        sw1 = 2
        sw2 = 4
        lambda1 = numpy.ones(self.x.shape[0])
        indx = numpy.arange(self.x.shape[0])+1
        lf = leadingfunction.lf

        # do some checks
        if(self.y.shape[0] != self.x.shape[0]):
            raise Exception, "".join(["dimensions do not agree for: y[n,m]",
                                      "and x[n]"])
        if(self.y.shape[1] != self.ixm.shape[0]):
            raise Exception, "".join(["dimensions do not agree for: y[n,m]",
                                      "and ixm[m]"])


        self.splcoef=[]
        if(len(self.y.shape) == 1):
            self.splcoef.append(
                spline_spec.splinecof3(self.x, self.y, c1, cn, lambda1, indx, 
                                       sw1, sw2, self.ixm, lf))
        else:
            for i in range(self.y.shape[1]):
                self.splcoef.append(
                    spline_spec.splinecof3(self.x, self.y[:,i], 
                                           c1, cn, lambda1, indx, 
                                           sw1, sw2, self.ixm[i], lf))

        self.splcoef = numpy.array(self.splcoef)
        self.__kind_computed__ |= self.SPECIAL_SPLINE
        self.importQuantity('splcoef', self.splcoef)
        self.importQuantity('__kind_computed__', self.__kind_computed__)

        return




    def splev_conv(self, xin, deriv=False, ext=0):
        """
        Spline evaluation.
        This routine calls splev which has to be importet from:
        scipy.interpolate.

        For description it is refered to the help of splev, which is 
        partly copied here.

        ATTENTION: This routine has two restrictions:
           * Computation of deriatives: only 1. derivative implemented!
           * Values for x <= eps (e.g. eps=1.e-30):  are set to eps!
             Because for odd m we divide by sqrt(x).

        INPUT:
          xin ... do spline interpolation for this point
          deriv . compute first derivative (True/False)

        RETURN:
          ret ... interpolated values y and first derivative if deriv=True
                  y=ret[0,:]  yp=ret[1,:]
          

--- in the man page for splev we may read:

splev(x, tck, der=0, ext=0)
    Evaluate a B-spline or its derivatives.
    
    Given the knots and coefficients of a B-spline representation, evaluate
    the value of the smoothing polynomial and its derivatives.  This is a
    wrapper around the FORTRAN routines splev and splder of FITPACK.
    
    Parameters
    ----------
    x : array_like
        A 1-D array of points at which to return the value of the smoothed
        spline or its derivatives.  If `tck` was returned from `splprep`,
        then the parameter values, u should be given.
    tck : tuple
        A sequence of length 3 returned by `splrep` or `splprep` containing
        the knots, coefficients, and degree of the spline.
    der : int
        The order of derivative of the spline to compute (must be less than
        or equal to k).
    ext : int
        Controls the value returned for elements of ``x`` not in the
        interval defined by the knot sequence.
    
        * if ext=0, return the extrapolated value.
        * if ext=1, return 0
        * if ext=2, raise a ValueError
    
        The default value is 0.
    
    Returns
    -------
    y : ndarray or list of ndarrays
        An array of values representing the spline function evaluated at
        the points in ``x``.  If `tck` was returned from splrep, then this
        is a list of arrays representing the curve in N-dimensional space.

        """
        #ext = 1

        n = len(self.tck) #self.tck.shape[0]
        y = numpy.zeros(n)
        if(deriv):
            yp = numpy.zeros(n)

        for i in range(n):
            y[i] = interpolate.splev(xin, self.tck[i], der=0, ext=ext)
            if(deriv):
                yp[i] = interpolate.splev(xin, self.tck[i], der=1, ext=ext)

        if(deriv):
            return(numpy.array([y, yp]))
        else:
            return(y)


    def splev_sqrt(self, xin, deriv=False, ext=0):
        """
        Spline evaluation.

        For computing spline coefficients sqrt_tck the routine 
        calc_spline_sqrt_coefs() has to be called first. 
        For computing the coefficients the following has bee taken 
        into account:
        Divide odd m terms by sqrt(s), 
        keep even m (m=poloidal modenumber) terms as they are.
        Use 3rd order qubic spline.
 
        For splev_sqrt() the poloidal modenumber is taken into account 
        in order to compute correct values.

        This routine calls splev which has to be importet from:
        scipy.interpolate.

        For description it is refered to the help of splev, which is 
        partly copied here.

        ATTENTION: This routine has two restrictions:
           * Computation of deriatives: only 1. derivative implemented!
           * Values for x <= eps (e.g. eps=1.e-30):  are set to eps!
             Because for odd m we divide by sqrt(x).

        INPUT:
          xin ... do spline interpolation for this point
          deriv . compute first derivative (True/False)

        RETURN:
          ret ... interpolated values y and first derivative if deriv=True
                  y=ret[0,:]  yp=ret[1,:]

--- in the man page for splev we may read:

splev(x, tck, der=0, ext=0)
    Evaluate a B-spline or its derivatives.
    
    Given the knots and coefficients of a B-spline representation, evaluate
    the value of the smoothing polynomial and its derivatives.  This is a
    wrapper around the FORTRAN routines splev and splder of FITPACK.
    
    Parameters
    ----------
    x : array_like
        A 1-D array of points at which to return the value of the smoothed
        spline or its derivatives.  If `tck` was returned from `splprep`,
        then the parameter values, u should be given.
    tck : tuple
        A sequence of length 3 returned by `splrep` or `splprep` containing
        the knots, coefficients, and degree of the spline.
    der : int
        The order of derivative of the spline to compute (must be less than
        or equal to k).
    ext : int
        Controls the value returned for elements of ``x`` not in the
        interval defined by the knot sequence.
    
        * if ext=0, return the extrapolated value.
        * if ext=1, return 0
        * if ext=2, raise a ValueError
    
        The default value is 0.
    
    Returns
    -------
    y : ndarray or list of ndarrays
        An array of values representing the spline function evaluated at
        the points in ``x``.  If `tck` was returned from splrep, then this
        is a list of arrays representing the curve in N-dimensional space.

        """

        # as we shall not divide by zero...
        eps = 1.e-30
        xin = numpy.where(xin >= eps, xin, eps)

        #ext = 0

        n = len(self.sqrt_tck) #self.sqrt_tck.shape[0]
        y = numpy.zeros(n)
        if(deriv):
            yp = numpy.zeros(n)

        for i in range(n):
            sx = numpy.sqrt(xin)
            if(not deriv):
                y[i] = interpolate.splev(xin, self.sqrt_tck[i], der=0, ext=ext)
                if(self.ixm[i]%2 == 1):
                    y[i] *= sx

            if(deriv):
                yp[i] = interpolate.splev(xin, self.sqrt_tck[i], der=1, ext=ext)
                if(self.ixm[i]%2 != 0):
                    yp[i] = yp[i]*sx + y[i]/(2.0*sx)

        if(deriv):
            return(numpy.array([y, yp]))
        else:
            return(y)



    def splev_spec(self, xin, deriv=False):
        """
        Spline evaluation.

        INPUT:
          xin ... do spline interpolation for this point
          deriv . compute first derivative (True/False)

        RETURN:
          ret ... interpolated values y and first derivative if deriv=True
                  y=ret[0,:]  yp=ret[1,:]

        """
        
        if(self.splcoef.shape[0] != self.ixm.shape[0]):
            raise Exception, "".join(["dimensions do not agree for: y[n,m]",
                                      "and ixm[m]"])

        n = len(self.splcoef) #self.splcoef.shape[0]
        y = numpy.zeros(n)
        if(deriv):
            yp = numpy.zeros(n)

        lf = leadingfunction.lf
        res = spline_spec.splint_horner3_driv(self.x, 
                                              self.splcoef[:,0,:].T,
                                              self.splcoef[:,1,:].T,
                                              self.splcoef[:,2,:].T,
                                              self.splcoef[:,3,:].T,
                                              deriv, self.ixm, xin, lf)
        res=numpy.array(res)

        if(deriv):
            return(numpy.array(res))
        else:
            return(numpy.array(res)[0,:])



    def packTCK(self, tck):
        """
        Re-arange tck structure, which is returned by
        scipy.interpolate.splrep() and used by scipy.interpolate.splev(), 
        in order to have simple arrays. The structure tck contains spline 
        coefficients and the order of the splines.

        tck is structured as follows:
         tck: type list, length no_c equals to number of splined curves
          tck[0:no_c-1]: type tuple, length no_t equals 3 for 3rd order spline
           tck[:][0:no_t-2] spline coefficients
           tck[:][no_t-1]   order of spline (this is last element tck[:][-1]

        Attention: If just one curve is splined, the structure of tck
        starts at the level of touple and the top level list is not present!

        See also help for scipy.interpolate.splrep()
        and unpackTCK().

        INPUT:
          tck ... structure returned by scipy.interpolate.splrep()

        RETURN:
          spl_coefs . array with spline coefficients
          spl_order . array with order of the splines
        """

        # the top level of type list exists only, in case the number
        # of splined curves is greater than 1!
        if(type(tck) == list):
            no_curves = len(tck)
            no_tuples = len(tck[0])
            no_coefs = len(tck[0][0])

            spl_order = numpy.zeros(no_curves, dtype=int)
            spl_coefs = numpy.zeros([no_curves, no_tuples-1, no_coefs])
            for i in range(no_curves):
                spl_order[i] = int(tck[i][-1])
                for j in range(no_tuples-1):
                    spl_coefs[i,j,:] = tck[i][j]
        else:
            no_tuples = len(tck)
            no_coefs = len(tck[0])

            spl_order = int(tck[-1])
            spl_coefs = numpy.zeros([no_tuples-1, no_coefs])
            for j in range(no_tuples-1):
                spl_coefs[j,:] = tck[j]

        return(spl_coefs, spl_order)


    def unpackTCK(self, spl_coefs, spl_order):
        """
        Based on the input, which is returned by
        scipy.interpolate.splrep() and used by scipy.interpolate.splev(), 
        in order to have simple arrays. The structure tck contains spline 
        coefficients and the order of the splines.

        tck is structured as follows:
         tck: type list, length no_c equals to number of splined curves
          tck[0:no_c-1]: type tuple, length no_t equals 3 for 3rd order spline
           tck[:][0:no_t-2] spline coefficients
           tck[:][no_t-1]   order of spline (this is last element tck[:][-1]

        Attention: If just one curve is splined, the structure of tck
        starts at the level of touple and the top level list is not present!

        See also help for scipy.interpolate.splrep()
        and packTCK().

        INPUT:
          spl_coefs . array with spline coefficients
          spl_order . array with order of the splines

        RETURN:
          tck ... structure returned by scipy.interpolate.splrep()
        """
        
        dims = spl_coefs.shape
        
        # The top level of type list exists only, in case the number
        # of splined curves is greater than 1! Means here, that spl_coefs
        # has 3 dimensions when more than one cureve has been splined,
        # otherwise we have 2 dimensions.
        if(len(dims) == 3):
            (no_curves, no_tuples, no_coefs) = dims
            tck = []
            for i in range(no_curves):
                tck.append(tuple())
                for j in range(no_tuples):
                    tck[-1] += (spl_coefs[i,j,:],)
                tck[-1] += (int(spl_order[i]),)
        else:
            (no_tuples, no_coefs) = dims
            tck = tuple()
            for j in range(no_tuples):
                tck += (spl_coefs[j,:],)
            tck += (int(spl_order),)


        return(tck)

