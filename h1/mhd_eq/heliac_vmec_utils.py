""" A collection of utilities for working with vmec and heliac files.
bdb 2011-2011
"""
import StringIO, os
import numpy as np
import pylab as pl
from numpy import sin, cos, sum, linspace, pi, array, average, max, argsort, arctan2, cumsum, diff, sqrt, zeros, mod
from curves_and_points import rms_distance_points_curve

debug = 2

def generate_VMEC_input2(input_file, boundary_details, **kwargs):
    '''Generate the VMEC input file based on the descur output and the various coefficients
    SH : 2Apr2013
    '''

    default_settings = {'phi_edge':'0.05','am':'14, -14', 'ai':'',
                        'ftol_array':'1.e-05,1.e-08,1.0e-10,8.0e-12',
                        'ns_array':'10, 23, 65, 100',
                        'niter':'15000','nzeta':'90', 'ntheta':'90',
                        'ntor':'18', 'mpol':'18',}
    for i in kwargs.keys():
        print i, kwargs[i]
        default_settings[i] = kwargs[i]

    default_settings['boundary'] = boundary_details.rstrip('\n&END\n')

    print os.getpid(), 'generating VMEC input'
    fin=open(input_file[0]+'outcurve','r')
    m=[]; n=[]; rbc=[]; zbs=[]; raxis=[]; zaxis=[];
    lines = fin.readlines(); i = 0
    while lines[i].find('MB')<0 and lines[i].find('NB')<0:i+=1
    i+=1
    tmp = [float(j) for j in lines[i].split() if j!='']
    while len(tmp)>1:
        m.append(int(tmp[0])); n.append(int(tmp[1]))
        rbc.append(tmp[2]); zbs.append(tmp[4])
        if len(tmp)==8:raxis.append(tmp[6]); zaxis.append(tmp[7])
        i+=1
        tmp = [float(j) for j in lines[i].split() if j!='']
    remaining_lines = lines[i+1:]

    with open('/home/srh112/code/python/h1_eq_generation/input.template2','r') as template:
        vmec_input_text = template.read()

    ai_success=0
    if default_settings['ai'] == '':
        iotacoeffs=open("/home/srh112/code/python/h1_eq_generation/iotacoeffs.dat",'r')
        for temp in iotacoeffs:
            temps=shlex.split(temp)
            #if float(temps[0])==float(infile[-3:])/100:
            if np.abs(float(temps[0])-float(input_file[2]))<1.e-3:
                fout.write(' ai             = '+temps[1]+', '+temps[2]+', '+temps[3]+'\n')
                ai_success+=1
        iotacoeffs.close()
    else:
        ai_success+=1
        #fout.write(ai_string)
    if ai_success != 1:
        print '!!!!!!!!!!! no ai was written!!!!!!!', float(input_file[2])

    for i in default_settings.keys():
        vmec_input_text = vmec_input_text.replace('<<{}>>'.format(i), default_settings[i])
    print vmec_input_text
    with open(input_file[0]+'input.'+input_file[1],'w') as fout:
        fout.write(vmec_input_text)
    # for temp in template:
    #     if temp[1:3] == 'am':
    #         break
    #     fout.write(temp)
    # import shlex

    # if ai_string == None:
    #     for temp in iotacoeffs:
    #         temps=shlex.split(temp)
    #         #if float(temps[0])==float(infile[-3:])/100:
    #         if np.abs(float(temps[0])-float(input_file[2]))<1.e-3:
    #             fout.write(' ai             = '+temps[1]+', '+temps[2]+', '+temps[3]+'\n')
    #             ai_success+=1
    # else:
    #     fout.write(ai_string)
    #     ai_success+=1
    #presscoeffs=open("/home/srh112/code/python/h1_eq_generation/pressurecoeffs.dat",'r')
    #fout.write(' am             = 14, -14 \n')
    #presscoeffs.close()
    #fout.write(boundary_details)

    # raxisout='RAXIS = '; zaxisout='ZAXIS = '
    # for i in range(len(raxis)):
    #     raxisout = raxisout + repr(raxis[i]) + ', '
    #     zaxisout = zaxisout + repr(zaxis[i]) + ', '

    # fout.write(raxisout+'\n'+zaxisout+'\n')

    # fout.writelines(remaining_lines)

    # fout.write('/\n&END')

    #fin.close()
    #template.close()
    #fout.close()

def read_descur_output(infile,suffix='.outcurve', fmt='%15.10f', maxm=50, maxn=40, mynan=0,debug=1):
    """ read the rmc etc from the descur outcurve file, and plot the cross 
    section
    This might be better as a class too.
    The outcurve file is a human readable text file, with convergence 
    info as well - the plotout file is more condensed
    """
    # Check that the DESCUR file is OK, and good enough
    #
    # Find the axis info from the descur file, and extract the first 3 pairs
    # Read the rbc zbs info

    df = open(infile,'r')
    lines = df.readlines()
    df.close()

    # simple check of file is - is there RBC(0,0)
    rbcptrs = np.array([max([line.find('RBC(0,0)'),line.find('RBC( 0,0)')]) for line in lines])
    rbcptr = np.where(rbcptrs > -1)[0]
    if len(rbcptr) == 0: raise ValueError(infile + ' incomplete - no RBC(0,0)')

    #find the axis info from the descur file, and extract the first 3 pairs
    axisptrs = np.array([line.find(' MB ') for line in lines])
    axisptr = np.where(axisptrs > -1)[0]
    if len(axisptr) == 0: raise ValueError(infile + ' incomplete, no axis lines')
    axisline = lines[axisptr[0]]
    # look at the headings to find where the axis is
    axistoks = axisline.split()
    RAXISindex = np.where(np.array(axistoks) == 'RAXIS')[0]

    raxis =[] ; zaxis=[]
    for j in range(3):
        try:
            raxis.append(float(lines[axisptr[0]+j+1].split()[RAXISindex]))
            zaxis.append(float(lines[axisptr[0]+j+1].split()[RAXISindex+1]))
        except IndexError:
            print('N too small to compute higher order axis terms %d' % (j))

    #print(raxis,zaxis)    

    # Now use exec(string) to read the lines line RBC(0,0) = 1.23
    # This is quite fast - about 10us per value
    #Need to replace RB with ^JRB, and ZB with ^ZB, then remove all space.
    
    dims = (2*maxn+2, maxm)  # leave a safety strip of width 1
    RBC=mynan+np.zeros(dims)
    RBS=np.zeros(dims)
    ZBC=np.zeros(dims)
    ZBS=np.zeros(dims)

    temp = ""
    for ln in lines[rbcptr[0]:]:
        temp+=(ln.replace("RB","\nRB")
                    .replace("ZB","\nZB")
                    .replace(" ","")
                    .replace("(","[")
                    .replace(")","]")
                    )
    exec(temp)

    max_m = max(np.where(np.sum(abs(ZBS),0)>0)[0])
    max_n = max(np.where(np.sum(abs(ZBS)[0:maxn,:],1)>0)[0])
    check_n = len(np.where(np.sum(abs(ZBS)[0:maxn,:],1)>0)[0])
    if check_n>(max_n+1): raise LookupError('inconsistency in ZBS:'
                                            'max ind=%d, % indicies' 
                                            % (max_n, check_n))
    if debug>0: print('array dimensions are %d x [-%d:%d]' % (max_m, max_n, max_n))

    if debug>3:
        import pdb; pdb.set_trace()
        'debugging, c to continue'
    return(RBC, RBS, ZBC, ZBS, max_m, max_n)

def plot_descur(RBC, RBS, ZBC, ZBS, max_m, max_n, phi_plane=0., nu=150, debug=1, **plkwargs):
    """ Companion plotting routine to read_descur_data
    """
    theta_arr = linspace(0,2*pi,nu,endpoint=True)
    r_c = 0*theta_arr
    z_c = 0*theta_arr
    r_c2 = 0*theta_arr
    z_c2 = 0*theta_arr
    for (i,theta) in enumerate(theta_arr):
        for ixm in range(0, max_m):
            if ixm==0: start_n=0
            else: start_n = -max_n
            for ixn in range(start_n, 1+max_n):
                arg = ixn*phi_plane + ixm*theta
                sinarg = sin(arg)
                cosarg = cos(arg)
                r_c[i] += (RBC[ixn,ixm]*cosarg) + (RBS[ixn,ixm]*sinarg)
                z_c[i] += (ZBC[ixn,ixm]*cosarg) + (ZBS[ixn,ixm]*sinarg)
                r_c2[i] += (RBC[ixn,ixm]*cosarg) #Shaun modification to show the stellarator symetric version
                z_c2[i] += (ZBS[ixn,ixm]*sinarg) #Shaun modification to show the stellarator symetric version

    pl.plot(r_c, z_c, **plkwargs)
    pl.plot(r_c2, z_c2, **plkwargs)
    pl.gca().set_aspect('equal')
    pl.show()
    if debug>3:
        import pdb; pdb.set_trace()
        'debugging, c to continue'

def do_plot(infile, debug=1, phi_plane=0, hold=0, **plkwargs):
    """ convenience function to read_descur_data to plot"""
    (RBC, RBS, ZBC, ZBS,max_m, max_n) = read_descur_output(infile,debug=debug)
    if len(np.shape(phi_plane))==0: phi_plane = array([phi_plane])
    for i in range(len(phi_plane)):
        plot_descur(RBC, RBS, ZBC, ZBS,max_m=max_m, max_n=max_n, debug=debug, phi_plane=phi_plane[i], hold=hold)

# open file
# loop until eof
#   read one line to get w, read with loadtxt (or f?)
#   read 10 lines to get ,m n, coeff - coeff seems redundant
#   read 10 long lines to get each of 10 surfaces

def make_vmec_input2(infile, header, delrbs=True, prefix='input.',suffix='.outcurve', fmt='%15.10f'):
    """ Copy header + the axis and rbc/zbc info from the descur outcurve file
    prefix is the prefix of the output file, and suffix refers to the input
    the stem is the common part.
    Looks like this assumes files are in the current directory.
       make_vmec_input('h1ass027v1_1p00.outcurve', 'header.h1ass027v1_pxx_lowb')

    """
    # Check that the DESCUR file is OK, and good enough
    #
    # Make the output file name from the input filename,
    # Create the output file from header, be it a file or a string
    # Find the axis info from the descur file, and extract the first 3 pairs
    # Copy the rbc zbs info, deleteing the rbs/zbc info if:
    #      delrbs is True AND the coeffs are 0

    # simple check of file is - is there axis info?
    df = open(infile,'r')
    lines = df.readlines()
    rbcptrs = np.array([max([line.find('RBC(0,0)'),line.find('RBC( 0,0)'),line.find('RBC(  0,0)')]) for line in lines])
    rbcptr = np.where(rbcptrs > -1)[0]
    if len(rbcptr) == 0: raise ValueError(infile + ' incomplete - no RBC(0,0)')

    # make the output file name from the input filename,
    ptr = infile.find(suffix)
    stem = infile[:ptr]
    outfile = prefix+stem
    # create the output file from header, be it a file or a string
    #hf = open(header,'r')
    # copy header over to output
    #of = open(outfile,'w')
    #of.writelines(hf.readlines())

    #find the axis info from the descur file, and extract the first 3 pairs
    axisptrs = np.array([line.find(' MB ') for line in lines])
    axisptr = np.where(axisptrs > -1)[0]
    if len(axisptr) == 0: raise ValueError(infile + ' incomplete, no axis lines')
    axisline = lines[axisptr[0]]
    axistoks = axisline.split()
    RAXISindex = np.where(np.array(axistoks) == 'RAXIS')[0]

    raxis =[] ; zaxis=[]
    for j in range(3):
        raxis.append(float(lines[axisptr[0]+j+1].split()[RAXISindex]))
        zaxis.append(float(lines[axisptr[0]+j+1].split()[RAXISindex+1]))

    #print(raxis,zaxis)    
    output_string = ''
    output_string +=("RAXIS = "+fmt+fmt+fmt + " \n") % tuple([r for r in raxis])
    output_string +=("ZAXIS = "+fmt+fmt+fmt + " \n") % tuple([z for z in zaxis])
    #of.write(("RAXIS = "+fmt+fmt+fmt + " \n") % tuple([r for r in raxis]))
    #of.write(("ZAXIS = "+fmt+fmt+fmt + " \n") % tuple([z for z in zaxis]))

    # copy the rbc zbs info, deleteing the rbs/zbc info if:
    #      delrbs is True AND the coeffs are 0
    valptr=(rbcptr[0])
    for line in lines[valptr:]:
        if delrbs:
            start = line.find('RBS')
            end = line.find('ZBS')
            line = line[0:start]+line[end:]
            output_string += line
        #of.write(line)
    output_string += '&END\n'
    #print output_string
    return output_string
    #of.write('&END\n')
    #of.close()

# open file
# loop until eof
#   read one line to get w, read with loadtxt (or f?)
#   read 10 lines to get ,m n, coeff - coeff seems redundant
#   read 10 long lines to get each of 10 surfaces




def read_array(filehandle, rows, columns=None, dtype=None, debug=debug):
    """ This is a practice, will probably write customised routines.
    """
    strobj = StringIO.StringIO()
    for r in range(rows):
        strobj.write(filehandle.readline())
        if debug>2: print r,
    #print(strobj.getvalue())
    #strobj.getvalue()
    strobj.seek(0)    # rewind, needed if write used
    return(np.loadtxt(strobj,dtype=dtype))    

class DescurPlotout():
    """ a class to read and plot descur plotout files, and compute RMS error
    plotout files contain the rz array followed by the rmnc arrays etc
    This overlaps in function with the DescurInData class
    """
    npoints = None
    nplanes = None
    nperiods = None

    def __init__(self, filename=None, plot=False, debug=1):
        
        """ Read a .trace file, such as used for descur input, into the object
        """
        if filename ==None:
            filename = "plotout"
        fp = open(filename,'r')
        self.filename=filename
        params = read_array(fp, rows=1, dtype = int)    
        (self.mpol, self.npoints, self.nplanes, self.mpminus1, self.nphi2, self.nperiods, self.mpnt) = params
        self.rz = read_array(fp, rows=self.npoints*self.nplanes)
        if debug>3:
            import pdb; pdb.set_trace()
            'debugging, c to continue'                
        self.rzmn = read_array(fp, rows=self.mpol*(self.nphi2*2+1))
        self.rmnc = zeros([2*self.nphi2+2, self.mpol])  # allow a 1 wide buffer 
        self.rmns = self.rmnc.copy()
        self.zmns = self.rmnc.copy()
        self.zmnc = self.rmnc.copy()

        nphi=1+2*self.nphi2
        for m in range(self.mpol):  # last is mpol-1
            self.rmnc[-self.nphi2-1:-1,m] = self.rzmn[m*nphi:m*nphi+self.nphi2,0]
            self.rmnc[0:self.nphi2+1,m] = self.rzmn[m*nphi+self.nphi2:(m+1)*nphi,0]

            self.rmns[-self.nphi2-1:-1,m] = self.rzmn[m*nphi:m*nphi+self.nphi2,2]
            self.rmns[0:self.nphi2+1,m] = self.rzmn[m*nphi+self.nphi2:(m+1)*nphi,2]

            self.zmns[-self.nphi2-1:-1,m] = self.rzmn[m*nphi:m*nphi+self.nphi2,1]
            self.zmns[0:self.nphi2+1,m] = self.rzmn[m*nphi+self.nphi2:(m+1)*nphi,1]

            self.zmnc[-self.nphi2-1:-1,m] = self.rzmn[m*nphi:m*nphi+self.nphi2,3]
            self.zmnc[0:self.nphi2+1,m] = self.rzmn[m*nphi+self.nphi2:(m+1)*nphi,3]
        if plot: self.plot_surface()

    def make_curve(self, nu=150, phi=0,debug=debug):
        """ return (r,z) at phi for nu points around the closed curve defined 
        by the descur fourier series (representing one toroidal surface)
        """
        theta = linspace(0,2*pi,endpoint=True, num=nu)
        r = 0*theta
        z = 0*theta
        for n in range(-self.nphi2,self.nphi2+1):
            for m in range(0,self.mpol):
                arg = m*theta - n*phi        # be careful of the +/- sign
                cosarg = cos(arg)
                sinarg = sin(arg)
                # adding in the asym args increased time from 8.5ms to 11.5ms (nu=150) 
                r += self.rmnc[n,m] * cosarg + self.rmns[n,m] * sinarg
                z += self.zmns[n,m] * sinarg + self.zmnc[n,m] * cosarg
        return((r,z))        

    def plot_surface(self, debug=debug, nu=150, phi=0, hold=True, linecolor='gray', marker=',',**plkwargs):
        (r,z)=self.make_curve(nu=nu, phi=phi)
        plkwargs.update({'color':linecolor, 'linewidth':0.3})
        if (linecolor != None) and (linecolor != ''): 
            pl.plot(r,z, hold=hold, **plkwargs)
        pl.gca().set_aspect('equal')
        if marker != '': pl.plot(r,z,marker, hold=1)        

    def plot_points(self, plane=0, hold=True):
        """ plot the points from <plane> """
        pts = self.npoints
        pl.plot(self.rz[plane*pts:(plane+1)*pts,0],self.rz[plane*pts:(plane+1)*pts,1],'+', hold=hold)
        if debug>3:
            import pdb; pdb.set_trace()
            'debugging, c to continue'                

    def RMS_error(self, nu=150, plane = 0, debug=debug, plot=False, hold=0):
        """ Calculate the RMS distance between stored points in the given plane
        and the piecewise linear approximation to the curve representing the 
        toroidal surface stored in the class.
        """
        pts = self.npoints
        points = self.rz[plane*pts:(plane+1)*pts]
        curve = self.make_curve(nu=nu, phi = (2*pi*plane)/(self.nplanes))
        if plot:
            pl.plot(points[:,0],points[:,1],'+', hold=hold)
            pl.plot(curve[0],curve[1])
            pl.gca().set_aspect('equal')
        return(rms_distance_points_curve(points, array(curve).T))
    
    def plot_RMS_errors(self, nu=150, planes = -1, plot_planes=False):
        """ For a range of planes, calculate the RMS distance from the points
        to the curve defined by rmnc etc, using a polygon with nu sides.
        Optionally make a collage of all the planes with curves and points plotted
        """
        if planes == -1:  # -1 means all
            planes = range(0, self.nplanes)
        dist = []
        for p in planes:
            if plot_planes: 
                num = len(planes)
                rows = int(num/(sqrt(num)))
                cols = num/rows
                if rows*cols<num: cols += 1
                pl.subplot(rows,cols, p+1) 
            dist.append(self.RMS_error(nu=nu, plane=p, plot=plot_planes))
            if plot_planes: pl.title('%2.2g' % (dist[p]))

        if plot_planes and len(planes) > 1: pl.figure()
        if len(planes)>1: pl.plot(planes, dist)

    def make_fake_trace(self, filename='test.trace', npoints=None, nplanes=None, nperiods=None, iota=1.23, debug=debug):
        """ Make a file suitable for descur input from rmnc etc.
        Useful for testing descur - this creates a trace which in pronciple
        should be exactly reproducible, at because it is a function that
        DESCUR can represent exactly.
        """     
        if npoints == None: npoints = self.npoints
        if nplanes == None: nplanes = self.nplanes
        if nperiods == None: nperiods = self.nperiods
        f = open(filename,'w')
        f.write(' %10d %10d %10d \n' % (npoints, nplanes, nperiods))
        for p in range(nplanes):
            phi = (2*pi*p)/nplanes
            xzarr = self.make_curve(nu=npoints, phi=phi )
            for (i,xz) in enumerate(np.transpose(xzarr)):
                nrot = int(mod(npoints*i*iota, npoints))  # wild guess
                f.write(' %16.10f %16.10f %16.10f \n' % (xz[0],phi+2*pi*nrot,xz[1]))

        f.close()
        if debug>2:
            import pdb; pdb.set_trace()
            'debugging, c to continue'                


### end DescurPlotout 

def read_bline_trace(fileorstring=None, traces=None, retval='r'):
    """ Given a complete file name or the file as a string, read a bline file
    and return:
        retval='a':  return all the values, (a[0]= 1 row of file)
               'x':  return trace xyz        a[0] = [x,y,z]
               'r':  return trace r/z/phi
        traces=1 or traces = [1] or traces = [1,3,4] or traces = None (all)
        Note: for more than one trace, return is a list of 2d arrays
        This is easier to deal with I think.
    """
    if fileorstring == None:
        fileorstring = '/home/bdb112/heliac/vmec_scan/h1ass027v1_p0.10.bline'

    if np.size(fileorstring)==1:
        f=open(fileorstring)
        lines = f.readlines()
    else: 
        lines=fileorstring
        
    target = 'tart trace'

    ptrs = np.array([s.find(target) for s in lines])
    trace_ptrs = np.where(np.array(ptrs)>-1)[0]
    trace_numbers = [int(lines[trace_ptr][1+len(target):]) for trace_ptr in trace_ptrs]
    if traces == None: # return all
        traces = trace_numbers

    if np.size(traces) > 1:
        arr = []
        for t in traces:
            arr.append(read_bline_trace(lines, t, retval))
        return(arr)
    if len(np.shape(traces))>0: trace = traces[0]
    else: trace=traces

    if trace>len(trace_numbers): raise LookupError(str("%d requested, only %d traces in file" % (trace, len(trace_numbers))))
    if trace == len(trace_numbers): sarr = lines[trace_ptrs[-1]+1:]
    else:sarr = lines[trace_ptrs[trace-1]+1:trace_ptrs[trace]]
# now write the array to a string object so we can use loadtxt
    strobj = StringIO.StringIO()
    for line in sarr:
        strobj.write(line)
#    print(strobj.getvalue())
    strobj.seek(0)    # rewind, needed if write used
    try:
        arr = np.loadtxt(strobj)
    except IOError:
        raise IOError,str(('IOError seeking trace %d. Perhaps plane spacing '
                          'is incompatible.  Check earlier messages for earlier traces') %
                          (trace))

    if retval == 'r':
        return(arr[:,0:3])    
    elif retval == 'x':
        return(np.array([arr[:,0]*np.cos(arr[:,2]),arr[:,0]*np.sin(arr[:,2]),arr[:,1]]).T)    
    else:
        return(arr)    

def fix2pi_skips(phase, sign='+'):
    """ ensure that phase monotonically increases (+) or decreases by
    adding units of 2Pi - 
    New method: highly vectorised and much much faster, doesn't need 
    to know the sign.  2011: fixed loss of endpoints.
    """
    fixed_phase=phase
    fixed_phase[1:] = cumsum((diff(phase) < -pi)*2*pi)+phase[1:]
    fixed_phase[1:] = cumsum((diff(fixed_phase) > pi)*(-2*pi))+fixed_phase[1:]
    return(fixed_phase)

class DescurInData():
    """ Object to hold and operate on .trace data used as input by descur.
    Example: rotate trace 7/12 to almost stand config
       x=-7*pi/12
       D=DescurInData("../vmec/testing/descur_out_of_order/x7.trace")
       xyz=D.xyz.copy()
       xyz[:,0] = cos(x)*D.xyz[:,0] - sin(x)* D.xyz[:,2]
       xyz[:,2] = sin(x)*D.xyz[:,0] + cos(x)* D.xyz[:,2]
       D.xyz=xyz
       D.plot_plane()
       D.rewrite("7rot.trace")

    """
    xyz = None
    npoints = None
    nplanes = None
    nperiods = None
    def __init__(self, file=None, debug=1, myplot=None):
        
        """ Read a .trace file, such as used for descur input, into the object
        d=DescurInData(file={'shape':'ellipse','aspect':3})
        d=DescurInData(file={'shape':'circle', 'npoints':100, 'nplanes':24, 'nperiods':3})

        """
        if file ==None:
            file = "h1ass027v1_0p56.trace"
        if type(file) == type({}): self.create_from_dictionary(dict = file)   

        else:            
        # read from file    
            arr=np.loadtxt(file)
            self.filename = file
            self.xyz = arr[1:]
            self.npoints = int(arr[0,0])
            self.nplanes = int(arr[0,1])
            self.nperiods = int(arr[0,2])

        if myplot == None: myplot=pl.plot
        print('%d points, in %d planes, %d periods' % (self.npoints, self.nplanes, self.nperiods))

    def create_from_dictionary(self,  dict):    
        """ create a trace from a dictionary of shapes
        """
        shapedict = {'circle': {'major_radius': 1, 'minor_radius': 0.2, 'aspect': 1},
                     'ellipse': {'major_radius': 1, 'minor_radius': 0.2, 'aspect': 1.5}
                     }
        # get the params
        shape = shapedict[dict['shape']]
        # put default points etc in
        shape.update({'npoints':100, 'nplanes':24, 'nperiods':3})
        # change according to user input
        shape.update(dict)  # any changed entries
        self.npoints = shape['npoints']
        self.nplanes = shape['nplanes']
        self.nperiods = shape['nperiods']
        xyz=[]
        for  p in range(self.nplanes):
            for pt in range(self.npoints):
                phi = pt*2*np.pi/self.nperiods
                theta = 1.23*phi
                if shape['shape'] == 'circle': xyz.append([1+0.2*cos(theta), phi, 0.2*sin(theta)])
                elif shape['shape'] == 'ellipse': xyz.append([1+0.2*cos(theta), phi, 0.2*shape['aspect']*sin(theta)])
                else: 
                    print('allowed shapes are')
                    print(shapedict)
                    raise ValueError(str(shape) + ' not known')
        self.xyz=array(xyz)        
        self.filename=shape

    def reorder(self, rz=None, fract=0.7, plane=0,debug=0):
        """ Reorder the xyz data in the object - see rewrite to write
            it back to file.  The axis for angle measurement is
            partway (fract) between the average R value and the max R
            value
        """
        xyz = self.xyz[plane*self.npoints:(plane+1)*self.npoints]    
        if rz ==None:
            rz=[fract*max(xyz[:,0])+
                (1-fract)*average(xyz[:,0]), 
                average(xyz[:,2])]
        if debug>1:    
            pl.figure()
            pl.subplot(121)
            self.plot_plane(plane)
            pl.scatter(rz[0],rz[1],marker='+')    
            
        ang = arctan2(xyz[:,0] - rz[0], xyz[:,2] - rz[1])     
        angfix = ang #fix2pi_skips(ang)
        angfixoffs=angfix-min(angfix) # start at zero
        ind = argsort(angfixoffs)
        if debug > 1:
            pl.subplot(122)
            pl.plot(angfixoffs)
            pl.plot(angfixoffs[ind])
            pl.gca().set_aspect('auto')
        if sum(abs(ind - range(len(xyz[:,0])))) == 0:
            print('No reordering')

        for plane in range(self.nplanes):
            xyzp = self.xyz[plane*self.npoints:(plane+1)*self.npoints,:]
            xyzpfix = xyzp.copy()

            for i in range(3):
                xyzpfix[:,i] = xyzp[ind,i]

            self.xyz[plane*self.npoints:(plane+1)*self.npoints,:] = xyzpfix
        if debug>3:
            import pdb; pdb.set_trace()
            'debugging, c to continue'                

    def test_gaps(self, debug=debug):  
        """ Return the biggest gap in theta (actually will use arc length scaled to 2pi)
        """
        dxyz = self.xyz.copy()[0:self.npoints]  # copy the first plane
        dxyz[0:self.npoints-1] = self.xyz[1:self.npoints]-self.xyz[0:self.npoints-1]
        dxyz[-1] = self.xyz[0]-self.xyz[self.npoints-1]
        dist = sqrt(dxyz[:,0]**2 + dxyz[:,2]**2)
        if debug>3:
            import pdb; pdb.set_trace()
            'debugging, c to continue'
        return(2*pi*max(dist)/sum(dist))
#        dist2 = diff(xyz[:0])**2 + diff(xyz[:2])**2
        
    def rewrite(self, filename=None):
        """ Make a .trace file from the current data
        """
        if filename == None: filename=self.filename
        f = open(filename,'w')
        f.write(' %10d %10d %10d \n' % (self.npoints, self.nplanes, self.nperiods))
        for xyz in self.xyz[:]:
            f.write(' %16.10f %16.10f %16.10f \n' % tuple(xyz))
    
        f.close()
    
    def plot_plane(self, debug=1,plane=0, linecolor='gray', hold=0, **plkwargs):
        """ Plot R and Z for a given plane in the trace data
        """
        xyz = self.xyz[plane*self.npoints:(plane+1)*self.npoints]
        pl.plot(xyz[:,0],xyz[:,2], '+', hold=hold, **plkwargs)
        ax=pl.gca()
        ax.set_aspect('equal')
        current_col = ax.get_lines()[-1].get_color()
        pl.scatter(xyz[0,0],xyz[0,2], marker='o', color=current_col, facecolor='')
#color=ax._get_lines.colors[ax._get_lines.count % len(ax._get_lines.colors)])
        
        pl.title("%s plane %d/%d" % (self.filename, plane, self.nplanes))
        plkwargs.update({'color':linecolor, 'linewidth':0.1})
        if (linecolor != None) and (linecolor != ''): 
            pl.plot(xyz[:,0],xyz[:,2], **plkwargs)
        pl.draw()

def make_descur_input(fileorstring, trace=1, nplanes=12, npuncs=100, nperiods=3,max_dth=None,outfile=''):
    """ Really need to refactor reformat_to_trace - solve outfile clumsiness
    """
    if max_dth == None: 
        max_dth = 2*pi/nplanes

    while(1):
        print('try trace = %d' %(trace))
        try:
            outfile=reformat_to_trace(fileorstring, trace, nplanes, npuncs, nperiods,outfile='')
        except ValueError as details:
            print('Value error in reformat - try a smaller trace', details)
            trace -= 1
            continue
        except IndexError as details:
            print('Index error in reformat - try a smaller trace', details)
            trace -= 1
            continue
        except LookupError as details:
            print('Lookup error in reformat - try a smaller trace', details)
            trace -= 1
            continue
        print('reformatting %s' % outfile)
        D=DescurInData(outfile)
        D.reorder()
        dth = D.test_gaps()
        if dth<max_dth:
            D.rewrite()
            return()    
        trace -= 1
        if trace < 0: raise ValueError(
            'Theta gap too large (%2.2g/2pi) in %s' %(max_dth, outfile))
        else: 
            print('Theta gap too large (%2.2g) in %s, try a smaller surface' %
                  (max_dth, outfile))
    

def reformat_to_trace(fileorstring, trace=1, nplanes=12, npuncs=100, nperiods=3,outfile=''):
    """ read a trace from the corresponding bline file, check that there 
    are enough points, and make a .trace file derived from the full input file name
    """
    rzphi=read_bline_trace(fileorstring, traces=trace, retval='r')

    ## check that the array size is OK.
    dphi = np.diff(np.sort(rzphi[:,2]))[0]
    fstep = 2*np.pi/(nperiods*nplanes*dphi)
    istep = int(round(fstep))
    if abs(fstep-istep)>.001:  
        raise LookupError('%d planes not possible with spacing of %8.6f, step is %8.6f' 
                          % ( nplanes, dphi, fstep))

    input_is_filename = len(fileorstring)>1
    if outfile == '':
        if not input_is_filename: outfile = 'out.trace'
        else: # make name from input filename, e.g. replace ".bline"
              # with .trace
            outfile=fileorstring
            ptr = outfile.find('bline')
            if ptr >= -1: outfile = outfile[0:ptr]
            outfile += 'trace'
    of = open(outfile,'w')
    of.write(str("%d  %d  %d\n" % (npuncs, nplanes, nperiods)))

# all the punctures for one plane, then the next, etc
# in the order they came out of bline
# perhaps we should reorder so they don't jump around.
    for p in range(nplanes):
        sel_rzphi = rzphi[range(p*istep,istep*nplanes*npuncs,istep*nplanes)]
        for rzp in sel_rzphi:
            of.write("%15.10f %15.10f %15.10f \n" % (rzp[0], rzp[2], rzp[1]))
    of.close()
    return(outfile)

""" read_bline_trace() """

if __name__ == '__main__':
    """Example usage..."""
    pl.rcParams['legend.fontsize']='medium'
    import sys
    args = sys.argv
    if len(args) > 1: 
        print args[1]
        filename = args[1]
    else: filename = None    

#    D=DescurInData(filename)
#    D.plot_plane()

    P=DescurPlotout(plot=True)
    pl.show()
