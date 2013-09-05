
"""
some convenience functions

"""

import string
import numpy
import time
import pylab


#WRITE_PROGRESS=True
WRITE_PROGRESS=False


try:
    import netCDF4 as netcdf
    netcdf_module = 'netCDF4'
except: 
    raise Exception, 'Unable to load netCDF4 module/library'



def clearFilename(fname_in):
    """
    Remove/change some characters from/in filename, e.g. ' ' -> '_'

    INPUT:
      fname_in .... filename to clear
    RETURN:
      fname_out ... 'new' filename
    """

    # first replace some characters
    fname_out = fname_in.replace(' ', '_')
    fname_out = fname_out.replace('\\', '')
    fname_out = fname_out.replace('?', '')
    fname_out = fname_out.replace('*', '')
    # keep just the allowed characters
    valid_chars = "-_.()%s%s" % (string.ascii_letters, string.digits)
    fname_out = ''.join(c for c in fname_out if c in valid_chars)

    return(fname_out)



def formatNumber(val, slen, acc, dig=3):
    """
    Format number to string with best fitting to length.

    INPUT:
      val .... value to format
      slen ... max length of string
      acc .... accuracy of formated number compared to val
      dig .... digits behind komma

    RETURN:
      nstr ... string representing val
    """
    #m = 13
    #n = 7
    # iterate m,n until fits
    # probably good: ss="% 7.3g" % (4e6/7.) ; print ss, len(ss)
    #fmtstr= "".join(["%", str(m), "g.", str(n)])
    #s = fmtstr % (4e6/7.) 
    #print s, len(s)
    
    # try it the simple way
    fmtstr = "".join(["%.", str(slen).strip(), "f"])
    nstr = fmtstr %(val)
    nstr = nstr.strip()
    if(len(nstr) > slen):
        nstr = nstr[:slen]

    # check for zero
    if(val == 0.0):
        return(nstr)

    # check accuracy
    h = float(nstr)
    if(h != 0.0):
        if(numpy.abs(val-h) < acc):
            return(nstr)
    
    # OK, not simple enough, construct the string
    po = numpy.int(numpy.round(numpy.log10(numpy.abs(val)),0))
    spo = str(po)
    s1 = str(val/10**po)
    pos = s1.find('.')
    l = pos+1+dig
    reqslen = l +len(spo)+1
    if(reqslen > slen):
        errtxt = "ERROR: can't fullfill requirements for %s; required slen=%d" %(str(val), reqslen)
        raise Exception(errtxt)
        return(None)

    if(reqslen < slen):
        l += slen-reqslen

    if(len(s1) < l):
        s1 = "".join([s1, '0' * (l-len(s1))])

    nstr = "".join([s1[:l], 'e', spo])

    return(nstr)


def eformat(f, prec, exp_digits):
    """
    Format numper with specifying the number of digits for the exponent
    e.g.:
    > print eformat(0.870927939438012, 15, 3)
    >  8.709279394380121e-001

    > print eformat(1.870927939438012, 15, 3)
    >  1.870927939438012e+000

    > print eformat(-1.870927939438012, 15, 3)
    > -1.870927939438012e+000

    INPUT:
      f .......... float number (to be formated)
      prec ....... precision
      exp_digits . number of digits for exponent

    RETURN:
      s ... string with formatted number
    """
    # 
    s = "% .*e"%(prec, f)
    mantissa, exp = s.split('e')
    # add 1 to digits as 1 is taken by sign +/-
    return("%se%+0*d"%(mantissa, exp_digits+1, int(exp)))



def list2string2(mylist):
    """
    Return a string with one line per listelement like:
      "0: element0\n1: element1\n2: element2"
    Last string is not terminated by '\n'

    INPUT:
      mylist ... list of strings

    RETURN:
      s ........ formated string
    """

    s = ""
    for i, v in enumerate(mylist):
        s = "".join([s, '' if i==0 else '\n', "{0:2d}: {1:s}".format(i,v)])

    return(s)


def flist2string(mylist):
    """
    Return a string with one line per listelement like:
      "0: element0\n1: element1\n2: element2"
    Last string is not terminated by '\n'

    INPUT:
      mylist ... list of strings

    RETURN:
      s ........ formated string
    """

    s = ""
    for i, v in enumerate(mylist):
        s = "".join([s, '' if i==0 else '\n', "{0:2d}: {1:f}".format(i,v)])

    return(s)


def warn_but_not_too_often(message, wait=3):
    """ 
    manage a message list to group similar messages together,
    if they come too often.  will return true if it is time to print the
    message again, and print the message
    """
    import time
    # check whether warn_msgs exists and declare in case...
    if(not globals().has_key('warn_but_not_too_often_msgs')):
        print "need to declare warn_..."
        global warn_but_not_too_often_msgs
        warn_but_not_too_often_msgs = {}

    if warn_but_not_too_often_msgs.has_key(message):
        that_msg = warn_but_not_too_often_msgs[message]
        if  time.clock() > (that_msg['tim']+wait):
            if that_msg['count']>0: 
                print("%d occurences of " % that_msg['count']),
            # end   if that_msg['count']>0: 
            print(message)
            tim = time.clock()
            count=0
        else: # if  secs() > (that_msg['tim']+wait):
            count = that_msg['count']+1
            tim =  that_msg['tim']
        # end   if  time.clock() > (that_msg['tim']+wait):
    else: # if warn_but_not_too_often_msgs.has_key(message):
        count=0
        tim = time.clock()
        print(message)
    # end   if warn_but_not_too_often_msgs.has_key(message):

    warn_but_not_too_often_msgs.update({message:{'count':count, 'tim': tim}})
    return(count==0)


##### some netCDF4 stuff   start

def writeNetCDF(mydata, nc_filename, imode='a'):
    """
    Write data to netCDF4 file.
    Arrays are always compressed.

    Intention: Store spline coefficients for usage e.g. for BOOZER or VMEC.

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
      imode ......... mode for opening: append 'a', write 'w'

    RETURN:
      None
    """

    if(nc_filename[-3:] != ".nc"):
        nc_filename = "".join([nc_filename, ".nc"])
        
    newfile = False
        
    try:
        spl = netcdf.Dataset(nc_filename, format='NETCDF4', mode=imode)
    except:
        newfile = True
        spl = netcdf.Dataset(nc_filename, format='NETCDF4', mode='w')


    for k in mydata.iterkeys():
        if(k != 'history'):
            appendToNetCDF(spl, k, mydata[k])
            spl.sync()

    # now add some history
    # We write history here, because if the netCDF4 file can't be
    # written/modified we do not like some entries for history.
    if(newfile):
        hstr = "".join(['Created: ', time.ctime(time.time())])
    else:
        try:
            hstr = spl.history
            hstr = "".join([hstr, '\n ',
                            'Modified: ', time.ctime(time.time()) 
                            ])
        except:
            hstr = "".join(['Modified: ', time.ctime(time.time())])

    if(mydata.has_key('history')):
        if(mydata['history'] != None):
            hstr = "".join([hstr,
                            '\n ', mydata['history']
                            ])

    spl.history = hstr

    spl.close()

    return

def appendToNetCDF(grp, myname, mydata, compression=True):
    """
    Append data to netCDF4 file/group.
    
    Data types handled: 
      * scalars: integer, float
      * numpy arrays: integer, float
      * strings

    INPUT:
      grp ...... handle of the netCDF4 group where to append/store data
      myname ... name of data
      mydata ... data
      compression  compress arrays

    RETURN:
      None
    """

    if(compression):
        zlib_str = ", zlib=True"   # for compression; ',' is needed, because
                                   # this string is simply added to the 
                                   # the command createVariable(...)
    else:
        zlib_str = ""

    data_n = "".join([myname, '_nc'])
    cmdstr2 = "{0:s}[:] = mydata".format(data_n)
    if(mydata.__class__ == str):
        cmdstr = "grp.{0:s} = mydata".format(data_n)
        if(WRITE_PROGRESS):   print 'w1', cmdstr
        exec(cmdstr)
    elif(mydata.__class__ == int):
        cmdstr = "{0:s} = grp.createVariable(myname, numpy.dtype('int32').char)".format(data_n)
        if(WRITE_PROGRESS):   print 'w2', cmdstr
        exec(cmdstr)
        if(WRITE_PROGRESS):   print 'w3', cmdstr2
        exec(cmdstr2)
    elif(mydata.__class__ == float):
        cmdstr = "{0:s} = grp.createVariable(myname, numpy.dtype('float32').char)".format(data_n)
        if(WRITE_PROGRESS):   print 'w4', cmdstr
        exec(cmdstr)
        if(WRITE_PROGRESS):   print 'w5', cmdstr2
        exec(cmdstr2)
    elif(mydata.__class__ == numpy.ndarray):
        atype = mydata.dtype
        ashape = mydata.shape
        dn = []
        for i,dim in enumerate(ashape):
            dim_n = "{0:s}_dim{1:d}".format(myname, i+1)
            dn.append(dim_n)
            cmdstr = "{0:s} = grp.createDimension(dim_n, dim)".format(dim_n)
            if(WRITE_PROGRESS):   print 'w6', cmdstr
            exec(cmdstr)

        datatype = mydata.dtype
        #if(mydata.dtype == '|S1'):
        #    datatype='S1'

        cmdstr = "{0:s} = grp.createVariable(myname, '{1:s}',tuple(dn){2:s})".format(data_n, datatype, zlib_str)
        if(WRITE_PROGRESS):   print 'w7', cmdstr
        exec(cmdstr)
        if(WRITE_PROGRESS):   print 'w8', cmdstr2
        exec(cmdstr2)
    else:
        print "Type {0:s} of {1:s} not handled...".format(mydata.__class__, myname)

    return


##### some netCDF4 stuff   end
