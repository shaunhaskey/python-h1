#! /usr/bin/env python

__doc__ = \
"""
Some useful routines for storing data using netCDF-4.

Intention is to write and read data prepared as dictionary 
into a netCDF-4 file.

Here, for arrays, it is made use of compression. 

ATTENTION: Not all datatypes are supported until now!!!
e.g. complex numbers, several numpy datatypes


* Write dictionary data:

Prepare data e.g. like as follows:

import numpy
import pylab
from scipy import interpolate
import SPLINE
import netCDF4Tools
x=numpy.linspace(0, 2*numpy.pi, 7)
xin=numpy.linspace(0, 2*numpy.pi, 14)
y1=numpy.sin(x)
y2=numpy.sin(2*x)
y3=numpy.sin(3*x)
ixm=numpy.zeros(3)
spl=SPLINE.SPLINE({'x':x, 'y':numpy.array([y1,y2,y3]).T,'kindspec':5, 'ixm':ixm})
idat={'x':x, 'y':numpy.array([y1,y2,y3]).T,'kindspec':5, 'ixm':ixm}
idat.update({'tck':spl.tck, 'splcoef':spl.splcoef})
idat2={'A':'some substructure', 'I':123456789, 'F':0.12345678901234567890}
idat3={'AA': 'a subsubstruct', 'II':42, 'FF':4.7}
idat2.update({'subsubstruct':idat3})
idat.update({'substruct':idat2})
idat.update({'history':'some history... just for fun'})

#and for saving (e.g. to file named ttt.nc:

netCDF4Tools.writeNetCDF4(idat, 'ttt.nc')


# * Read data:

mydata = netCDF4Tools.readNetCDF4('ttt.nc')

# do some check:
# extract spline coefficients
tck = mydata['tck']

yy1=interpolate.splev(xin,tck[0])
pylab.plot(x,y1,'b+-', xin,numpy.sin(xin),'g+-', xin,yy1,'r+-')
pylab.show()

The data structure of the dictionary is restored, so 
mydata has the same structure as idat.



ATTENTION/RECOMMENDATION:
  The example with spline coefficients is just to be seen as an example!
  The STRONGLY RECOMMENDET way for writing and reading coefficients of
  splines computed whith help of the calss SPLINE is for sure to 
  use the corresponding methods of the SPLINE class!
  One of the reasons is that the high level function 
  netCDF4Tools.writeNetCDF4() is for this case not efficient!!!


History:

21. June 2012 bhs
* scalar boolean type added 
  writing to netCDF-4 file: 
    * bool is re-defined to integer
    * name of variable is modified in order to recover the original datatype
  reading from netCDF-4 file:
    * re-defined datatypes are recognised by names of variables
    * restore original datatype and original name of variable

19. June 2012 bhs
* slight problems with datatypes fixed in appendToNetCDF4() fixed

13. June 2012 bhs
* problems restoring datastructure fixed
* use exception classes


28. May 2012 bhs
* start with this script

"""



__version__ = "0.5"
__versionTime__ = "13 June 2012 17:00"
__author__ = "Bernhard Seiwald <bernhard.seiwald@gmail.com>"


import numpy
import time

try:
    import netCDF4
    netcdf_module = 'netCDF4'
except: 
    raise Exception, 'Unable to load netCDF4 module/library'



#WRITE_PROGRESS=True
WRITE_PROGRESS=False


# define basic data types for netCDF4
MY_BASICS = [int, float, str]
MY_REDEFS_SCALAR = [bool, numpy.bool, numpy.bool8, numpy.bool_]
MY_ARRAYS = [numpy.ndarray]
MY_TUPLES = [list, tuple]

EXCLKEYTYPELIST = [str, tuple, list]

# marker for re-defined variables
REDEF_STR = "__redef_"
# separator between module and type. 
# needed because variable names do not allow for '.'
# numpy.bool_  ->  numpy_mrt_bool
MRTSEP_STR = "_mrt_"

# use compression for arrays
COMPRESS_ARR = True


class UnKnownKey(Exception):
    """
    The specifier for a key to exclude must be one defined  
    in EXCLKEYTYPELIST.
    """
    pass


class UnhandledVariable(Exception):
    """
    Can't add nameless variable to dictionary
    """
    pass


class UnhandledDataType(Exception):
    """
    addToTarget():
    Add/append data to a given structure and return extended structure.
    Target data type not handled...
    """
    pass




##### read data from netCDF-4

def readNetCDF4(nc_filename):
    """
    Read a netCDF-4 file which was previously written by
    writeNetCDF4() and restore the python data structure.

    INPUT:
      nc_filename ... name of netCDF-4 file

    RETURN:
      retdata .. dictionary with restored python data structure
    """

    retdata = None

    if(nc_filename[-3:] != ".nc"):
        nc_filename = "".join([nc_filename, ".nc"])

    try:
        with netCDF4.Dataset(nc_filename, format='NETCDF4', 
                             mode='r') as rootgrp:
            # walk through the group an add data
            retdata = walkGroup(rootgrp)
    except (UnKnownKey, UnhandledDataType, UnhandledVariable,
            IOError, RuntimeError):
        print "Can't read data from netCDF file {0:s}".format(nc_filename)
        raise

    return(retdata)


def walkGroup(grp, mydata=None):
    """
    Walk through a group of netCDF-4 data and append data
    to a python dictionary.

    Note:
    Here we shall not use mydata={} for initialising! 
    When using def walkGroup(grp, mydata={}) data from first
    call of this routine are still present when calling again!

    INPUT:
      grp ... netCDF-4 group returned either 
      mydata  structure where to store/append data

    RETURN:
      mydata  python dictionay containing data of the rootgroup
              returned by netCDF4.Dataset() or some subgroup
              e.g. found in rootgroup.groups
    """

    if(mydata == None):
        mydata = dict()

    # deal with attributes
    if(WRITE_PROGRESS):   print "walkGroup(): deal with attributes: "
    myattributes = grp.ncattrs()
    if(len(myattributes) != 0):
        for myattr in myattributes:
            exec("myad=grp.{0:s}".format(myattr))
            mydata = addToTarget(mydata, myad, mname=myattr)
    else:
        if(WRITE_PROGRESS):   print "  no attributes stored"


    # deal with variables
    if(WRITE_PROGRESS):   print "walkGroup(): deal with variables: "
    myvariables = grp.variables
    if(len(myvariables) != 0):
        for myvar in myvariables:
            if(myvariables[myvar].dimensions == ()):
                exec("myval=grp.variables['{0:s}'].getValue()[0]".format(myvar))
            else:
                exec("myval=grp.variables['{0:s}'][:]".format(myvar))
            # check for re-defined variable
            pos = myvar.find(REDEF_STR)
            if(pos > -1):
                origtype = myvar[pos+len(REDEF_STR):-2]
                exestr = "myval = {0:s}(myval)".format(origtype)
                exestr = exestr.replace(MRTSEP_STR, '.')
                exestr = exestr.replace("__builtin__.", '')
                myvar = myvar[:pos]
                if(WRITE_PROGRESS):   print exestr
                exec(exestr)
            mydata = addToTarget(mydata, myval, mname=myvar)
    else:
        if(WRITE_PROGRESS):   print "  no variables stored"


    # deal with groups
    if(WRITE_PROGRESS):   print "walkGroup(): deal with groups: "
    mygroups = grp.groups
    if(len(mygroups) != 0):
        for mygrp in mygroups:
            # The group name gives a hint what we have to expect.
            # Add the correct type to the dictionary.
            # First, deal with goups named like tck__de__list.
            # In this example, tck is the name of an entry in a dictionary
            # and the data are, according __de__list, of type list.
            pos = mygrp.find("__de__")
            groupnamestype = False
            itemtype = "dict"  # default type is dictionary
            exestr = "grpdata={0:s}()".format(itemtype)
            if(WRITE_PROGRESS):   print exestr
            exec(exestr)
            if(pos > -1):
                itemtype = mygrp[pos+len('__de__'):]  # check for named entries
                itemname = mygrp[:pos]
                # does the name of the group tell something about the datatype?
                try:
                    exestr = "grpdata={0:s}()".format(itemtype)
                    if(WRITE_PROGRESS):   print exestr
                    exec(exestr)
                    groupnamestype = True
                except:  # no
                    pass
    
            # Deal with groups named like, e.g. tuple_00000.
            # From the name (example!) we know that this is 
            # entry no 00000 of a tuple. Similar for, e.g. list_0034 and so on.
            if(not groupnamestype):
                itemname = None
                ht = mygrp.split('_')
                if(len(ht) == 2):
                    itemtype = ht[0]
                    try:
                        exestr = "grpdata={0:s}()".format(itemtype)
                        if(WRITE_PROGRESS):   print exestr
                        exec(exestr)
                    except:   # fallback to dict
                        groupnamestype = False
                        itemtype = "dict"  # default type is dictionary
                        exestr = "grpdata={0:s}()".format(itemtype)
                        # raise Exception, "BIG PROBLEM 1 unknown/unhandled type >{0:s}<".format(mygrp)
                else:   # fallback to dict
                    groupnamestype = False
                    itemtype = "dict"  # default type is dictionary
                    exestr = "grpdata={0:s}()".format(itemtype)
                    #raise Exception, "BIG PROBLEM 2 unknown/unhandled type >{0:s}<".format(mygrp)
    
            grpdata = walkGroup(grp.groups[mygrp], grpdata)
            if(WRITE_PROGRESS):   print "NNN walkGroup(): appending: ", grpdata
            addToTarget(mydata, grpdata, mname=itemname)

    else:
        if(WRITE_PROGRESS):   print "  no groups stored"

    return(mydata)



def addToTarget(mtarget, mdata, mname=None):
    """
    Add/append data to a given structure and return extended structure.

    INPUT:
      mtarget ... structure (dict, list, tuple) where to append data
      mdata ..... data to append
      mname ..... name of data. necessary for appending data to dictionary!

    RETURN:
      mtarget ... modified structure with appended data

    """
    mtype = type(mtarget)

    if(mtype == dict):
        if(mname == None):
            raise UnhandledVariable("Can't add nameless variable to dictionary")
        exestr = "mtarget.update({{'{0:s}':mdata}})".format(mname)
    elif(mtype == list):
        exestr = "mtarget.append(mdata)"
    elif(mtype == tuple):
        exestr = "mtarget = mtarget + (mdata,)"
    else:
        raise UnhandledDataType("Target data type {0:s} not handled...".format(mtype))

    if(WRITE_PROGRESS):   print exestr
    exec(exestr)

    return(mtarget)


##### write data to netCDF-4


def writeNetCDF4(mydata, nc_filename):
    """
    Write data (dictionary) to netCDF4 file.
    Arrays are always compressed.

    Intention: Store data of dictionary.
    To read data, stored by writeNetCDF4(), use readNetCDF4(). This will
    restore the original data structure.

    The data to be stored are provided as dictionary.
    e.g. mydata = {"A": "some string", "B": array([-1., 3.4]), "I": 2}
    Data types handled: 
      * scalars: integer, float
      * numpy arrays: integer, float
      * strings

    ATTENTION:
      * No sub-dictionaries are handled!
      * Strings are always global attributes in the netCDF4 file!

    

    INPUT:
      mydata ........ dictionary with data.
      nc_filename ... name of netCDF4 file e.g.: test.nc

    RETURN:
      None
    """

    if(nc_filename[-3:] != ".nc"):
        nc_filename = "".join([nc_filename, ".nc"])
    
    # check whether file exists
    try:
        fd = open(nc_filename, mode='r')
        fd.close()
        print "file '{0:s}' exists".format(nc_filename)
        while(1):
            sinp = raw_input('overwrite (y/N) ').strip().lower()
            if(sinp in ['n', '']):
                return
            if(sinp == 'y'):
                break
    except:
        pass


    try:
        rootgrp = netCDF4.Dataset(nc_filename, format='NETCDF4', mode='w')
    except:
        # some own statement
        print "Unable to open {0:s}".format(nc_filename)
        # and pass the exception to next level
        raise 

    # here it becomes interesting - that's the real action
    # walk through the dictionary of input data and store data
    # to netCDF4
    #walkDict(rootgrp, mydata, exclkeys=None)
    try:
        walkDict(rootgrp, mydata, exclkeys=['history'])
    except:
        rootgrp.close()
        raise
        #raise Exception, "Can't write netCDF4 file {0:s}".format(nc_filename)

    # now add some history
    # We write history here, because if the netCDF4 file can't be
    # written/modified we do not like some entries for history.
    hstr = "".join(['Created: ', time.ctime(time.time())])

    if(mydata.has_key('history')):
        if(mydata['history'] != None):
            hstr = "".join([hstr,
                            '\n ', mydata['history']
                            ])

    rootgrp.history = hstr

    rootgrp.close()

    return


def walkDict(grp, mydata, exclkeys=None):
    """
    Deal with a dictionary. 
    We have to walk through the structure, create groups, until we reach basic
    data types for storage.

    Basic data types and arrays are written directly into the
    corresponding group.
    When an entry of the dictionary is a list, tuple, dictionary,...
    the entry name is appended modified:
    entryname -> entryname__de__<type.__name__>
    with __de__ means 'dictionary entry'.

    INPUT:
      grp ........ group instance as returned e.g. by netCDF4.createGroup()
      mydata ..... data the tuple consists of
      exclkeys ... list of strings; keywords to exclude for storage

    RETURN:
      None
    """

    # prepare list with keys which should be exluded
    if(type(exclkeys) not in EXCLKEYTYPELIST):
        raise UnKnownKey("type of exclkeys must be in ".format(EXCLKEYTYPELIST))

    if(exclkeys is None):
        exclkeys = ''

    exclkeys = list(exclkeys)
    for i in range(len(exclkeys)):
        exclkeys[i] = exclkeys[i].strip()
        if(len(exclkeys[i].strip()) == 0):
            exclkeys[i] = ''

    # now check the items in the dictionary
    for mykey,myvalue in mydata.iteritems():
        if(WRITE_PROGRESS):
            print "".join(["*** walkDict(): grp={0:s}  ".format(grp.path),
                           "key={0:s} {1:s}".format(mykey, type(myvalue))])

        if(mykey not in exclkeys):
            mydatatype = type(myvalue)
            try:
                numpydtype = numpy.dtype(mydatatype)
                isnumpydtype = True
            except:
                isnumpydtype = False

            if(WRITE_PROGRESS):
                print "walkDict(): mykey,myvalue=",mykey,myvalue
                print "walkDict(): mydatatype,isnumpydtype=", mydatatype,isnumpydtype

            if(mydatatype is tuple):
                newgrp = grp.createGroup("".join([mykey, '__de__', 
                                                  mydatatype.__name__]))
                walkTuple(newgrp, myvalue, exclkeys=exclkeys)
                newgrp.sync()
            elif(mydatatype is list):
                newgrp = grp.createGroup("".join([mykey, '__de__', 
                                                  mydatatype.__name__]))
                walkList(newgrp, myvalue, exclkeys=exclkeys)
                newgrp.sync()
            elif(mydatatype is dict):
                newgrp = grp.createGroup("".join([mykey, '__de__', 
                                                  mydatatype.__name__]))
                print "+++ newgrp=", "".join([mykey, '__de__', 
                                                  mydatatype.__name__])
                walkDict(newgrp, myvalue, exclkeys=exclkeys)
                newgrp.sync()
            elif(mydatatype in MY_REDEFS_SCALAR):
                appendToNetCDF4(grp, mykey, myvalue, redeftype=True)
            elif((mydatatype in MY_BASICS) or isnumpydtype):
                appendToNetCDF4(grp, mykey, myvalue)
            elif(mydatatype in MY_ARRAYS):
                appendToNetCDF4(grp, mykey, myvalue, compression=COMPRESS_ARR)
            else:
                print "".join(["walkdict(): {0:s} ".format(mydatatype),
                               "of '{0:s}' not handled...".format(mykey)])
        
    grp.sync

    return


def walkTuple(grp, mydata, exclkeys=None):
    """
    Deal with a tuple. 
    We have to walk through the structure, create groups, until we reach basic
    data types for storage.

    INPUT:
      grp ........ group instance as returned e.g. by netCDF4.createGroup()
      mydata ..... data the tuple consists of
      exclkeys ... keywords to exclude for storage; will be passed to
                   walkDict()

    RETURN:
      None
    """
    i_v = 0
    i_a = 0
    i_l = 0
    i_t = 0
    i_d = 0

    mydatatype = type(mydata)
    if(WRITE_PROGRESS):
        print "".join(["walkTuple(): dealing with : ",
                       "{0:s}".format(mydatatype.__name__)])

    for mydat in mydata:
        mydattype = type(mydat)
        try:
            numpydtype = numpy.dtype(mydatatype)
            isnumpydtype = True
        except:
            isnumpydtype = False

        if(mydattype is tuple):
            newgrpname = "{0:s}_{1:05d}".format(mydattype.__name__, i_t)
            newgrp = grp.createGroup(newgrpname)
            walkTuple(newgrp, mydat, exclkeys=exclkeys)
            newgrp.sync()
            i_t += 1
        elif(mydattype is list):
            newgrpname = "{0:s}_{1:05d}".format(mydattype.__name__, i_l)
            newgrp = grp.createGroup(newgrpname)
            walkList(newgrp, mydat, exclkeys=exclkeys)
            newgrp.sync()
            i_l += 1
        elif(mydattype == dict):
            newgrpname = "{0:s}_{1:05d}".format(mydattype.__name__, i_d)
            newgrp = grp.createGroup(newgrpname)
            walkDict(mydat, newgrp, exclkeys=exclkeys)
            newgrp.sync()
            i_d += 1
        elif((mydattype in MY_REDEFS_SCALAR)  or isnumpydtype):
            myname = "v{0:05d}".format(i_v)
            appendToNetCDF4(grp, myname, mydat, redeftype=True)
            i_v += 1
        elif((mydattype in MY_BASICS)  or isnumpydtype):
            myname = "v{0:05d}".format(i_v)
            appendToNetCDF4(grp, myname, mydat)
            i_v += 1
        elif(mydattype in MY_ARRAYS):
            myname = "a{0:05d}".format(i_a)
            appendToNetCDF4(grp, myname, mydat, compression=COMPRESS_ARR)
            i_a += 1
        else:
            print "".join(["walkTuple(): {0:s} ".format(mydattype),
                           "not handled..."])

    grp.sync
    return


def walkList(grp, mydata, exclkeys=None):
    """
    Deal with a list. 
    We have to walk through the structure, create groups, until we reach basic
    data types for storage.

    INPUT:
      grp ........ group instance as returned e.g. by netCDF4.createGroup()
      mydata ..... data the tuple consists of
      exclkeys ... keywords to exclude for storage; will be passed to
                   walkDict()

    RETURN:
      None
    """
    i_v = 0
    i_a = 0
    i_l = 0
    i_t = 0
    i_d = 0

    mydatatype = type(mydata)
    if(WRITE_PROGRESS):
        print "".join(["walkList(): dealing with : ",
                       "{0:s}".format(mydatatype.__name__)])

    for mydat in mydata:
        mydattype = type(mydat)
        try:
            numpydtype = numpy.dtype(mydatatype)
            isnumpydtype = True
        except:
            isnumpydtype = False

        if(mydattype is tuple):
            newgrpname = "{0:s}_{1:05d}".format(mydattype.__name__, i_t)
            newgrp = grp.createGroup(newgrpname)
            walkTuple(newgrp, mydat, exclkeys=exclkeys)
            newgrp.sync()
            i_t += 1
        elif(mydattype is list):
            newgrpname = "{0:s}_{1:05d}".format(mydattype.__name__, i_l)
            newgrp = grp.createGroup(newgrpname)
            walkList(newgrp, mydat, exclkeys=exclkeys)
            newgrp.sync()
            i_l += 1
        elif(mydattype == dict):
            newgrpname = "{0:s}_{1:05d}".format(mydattype.__name__, i_d)
            newgrp = grp.createGroup(newgrpname)
            walkDict(mydat, newgrp, exclkeys=exclkeys)
            newgrp.sync()
            i_d += 1
        elif(mydattype in MY_REDEFS_SCALAR):
            myname = "v{0:05d}".format(i_v)
            appendToNetCDF4(grp, myname, mydat, redeftype=True)
            i_v += 1
        elif((mydattype in MY_BASICS) or isnumpydtype):
            myname = "v{0:05d}".format(i_v)
            appendToNetCDF4(grp, myname, mydat)
            i_v += 1
        elif(mydattype in MY_ARRAYS):
            myname = "a{0:05d}".format(i_a)
            appendToNetCDF4(grp, myname, mydat, compression=COMPRESS_ARR)
            i_a += 1
        else:
            print "".join(["walkList(): {0:s} ".format(mydattype),
                           "not handled..."])

    grp.sync
    return


def appendToNetCDF4(grp, myname, mydata, compression=False, redeftype=False):
    """
    Append/write a variable to a netCDF4 group.

    INPUT:
      grp ......... group instance as returned e.g. by netCDF4.createGroup()
      myname ...... name of the variable
      mydata ...... a 'basic' data type like int, float, string, numpy.ndarray
                    in principle what's defined by MY_BASICS
      compression . use zip compression (True/False)
      redeftype ... re-define datatype (True/False)
                    is used to store basic python/numpy datatypes 
                    which are not supported by netCDF-4. Such a datatye
                    would be e.g. boolean.

    RETURN:
      None
    """

    if(compression):
        zlib_str = ", zlib=True"   # for compression; ',' is needed, because
                                   # this string is simply added to the 
                                   # the command createVariable(...)
    else:
        zlib_str = ""

    if(WRITE_PROGRESS):
        print "appendToNetCDF4(): adding '{0:s}' to netCDF4 file".format(myname)


    # check whether we have to re-define a datatype...
    if(WRITE_PROGRESS):  print "appendToNetCDF4(): redeftype=", redeftype
    if(redeftype):
        mmtype = type(mydata)
        mmtypestr = mmtype.__name__ 
        # deal with boolean
        if(mmtypestr.find('bool') > -1):  # change to int
            mydata = numpy.int(mydata)
            myname = "".join([myname, 
                              "{0:s}{1:s}{2:s}{3:s}__".format(REDEF_STR, 
                                                              mmtype.__module__,
                                                              MRTSEP_STR, 
                                                              mmtypestr)])
            print "mmtypestr"
            print "myname = ", myname
            print "mydata, type(mydata) = ", mydata, type(mydata)

    #data_n = "".join([myname, '_nc'])
    cmdstr2 = "{0:s}[:] = mydata".format(myname)
    # check for array or scalar,...
    try:
        mshapet = mydata.shape
        isarr = True
    except:
        isarr = False

    iswritten = False
    if(not isarr):
        if(type(mydata) == str):
            cmdstr = "grp.{0:s} = mydata".format(myname)
            if(WRITE_PROGRESS):   print 'appendToNetCDF4(): 1', cmdstr
            exec(cmdstr) 
            iswritten = True
        else:
            cdtype = numpy.dtype(type(mydata)).char
            cmdstr = "{0:s} = grp.createVariable(myname, '{1:s}')".format(myname, cdtype)
            if(WRITE_PROGRESS):   print 'appendToNetCDF4(): 2', cmdstr.replace("myname", myname)
            exec(cmdstr)
            if(WRITE_PROGRESS):   print 'appendToNetCDF4(): 3', cmdstr2
            exec(cmdstr2)
            iswritten = True
    else:
        atype = mydata.dtype
        ashape = mydata.shape
        dn = []
        for i,dim in enumerate(ashape):
            dim_n = "{0:s}_dim{1:d}".format(myname, i+1)
            dn.append(dim_n)
            cmdstr = "{0:s} = grp.createDimension(dim_n, dim)".format(dim_n)
            if(WRITE_PROGRESS):   print 'appendToNetCDF4(): 6', cmdstr
            exec(cmdstr)

        #if(atype == '|S1'):
        #    atype='S1'

        cmdstr = "{0:s} = grp.createVariable(myname, '{1:s}',tuple(dn){2:s})".format(myname, atype.char, zlib_str)
        if(WRITE_PROGRESS):   print 'appendToNetCDF4(): 7', cmdstr
        exec(cmdstr)
        if(WRITE_PROGRESS):   print 'appendToNetCDF4(): 8', cmdstr2
        exec(cmdstr2)
        iswritten = True


    if(not iswritten):
        print "Type {0:s} of variable '{1:s}' not handled...".format(mydata.__class__, myname)

    grp.sync()
    return

