#!/usr/bin/env python
#from importAgilentBin import importAgilentBin as readAgilentBin
from numpy import double, min, max, angle, pi
import numpy as np
import pylab as pl
from h1.diagnostics.signal_processing import fix2pi_skips, analytic_phase
from scipy.fftpack import hilbert, rfft, irfft
from time import sleep
import MDSplus as MDS
import sys

"""  Process and optionally plot interfeomter signals

 from plot_Agilent
"""
try: 
    from debug_ import debug_
except:
    from warnings import warn
    def debug_(debug, msg='', *args, **kwargs):
        if str(debug) != '0': warn('attempt to debug ' + msg + 
              " need boyd's debug_.py to debug most effectively")

def extract_phase(filename=None, sig=None, ref=None, t_raw=None, method='hilbert', plot=1, channels=[1], hold=1, time_offset=0, f=1e5, loc=1, bw=0.1, baseline=None, label=None, verbose=1, debug=1, **kwargs):
    """ Extract phase from an interferometer signal.  
    Hilbert with bw: 26ms for 50k points.
    """
    if debug is None:
        debug = verbose

    if not filename is None: 
        (t_raw, dat) = readAgilentBin(filename,[channels],debug=verbose)
        tim = t_raw - time_offset
        sig = dat[channels[0]]   
        if len(channels)>1: ref = dat[channels[1]]
    else:
        tim = t_raw - time_offset

    dt = np.average(np.diff(t_raw))
    if f is None:
        if ref is None:
            ref = sig
        f = np.argmax(np.fft.rfft(ref))/(dt*len(sig))

    indf = int(round(f*(dt*len(sig))))
    lowf = int(indf*(1-bw))
    highf = int(indf*(1+bw))

    if method == 'fourier_shift':
        # potentially fast, but it doesn't work???
        ft = rfft(sig)      # see comment about scipy.fftpack below
        ft[0:2*indf] = 0    # we remove what will be the high frequencies for this version
        ft[2*highf:] = 0    # so we need even numbers, and a factor of two
        ft = np.roll(ft, -2*indf) 
        raw_phs = irfft(ft)
        phs = raw_phs
        debug_(debug,1)

    else: # hilbert
        if bw is None:
            anal_sig = -hilbert(hilbert(sig)) + 1j * hilbert(sig)
        else:
            ft = rfft(sig)
            # for scipy.fftpack, the complex data are interleaved as reals
            ft[0:2*lowf] = 0    # so we need even numbers, and a factor of two
            ft[2*highf:] = 0
            filter
            sig = irfft(ft)
            anal_sig = sig + 1j * hilbert(sig)
            raw_phs = fix2pi_skips(angle(anal_sig))
            phs = fix2pi_skips(raw_phs + 2*pi*f*t_raw)
            
    if baseline is None:  # default to signal from 0.01 to 0.02 of total
        baseline = (0.01*np.array([1,2])*len(phs)).astype(int)
    if len(np.shape(baseline))>0:
        phs = phs - np.average(phs[baseline])
    if plot:
        if hold == 0:
            pl.clf()
        pl.plot(tim, phs,label='phase, rad')
        if label is not None: lab = label
        elif filename is not None: lab = filename.split('/')[-1]
        else: lab = ''
        if loc != 0: 
            pl.legend(loc=loc, title=lab, **kwargs)
            pl.rcParams['legend.fontsize']='small'
        elif loc == 0 : pl.title(lab)
        else: pass
            
        pl.show()
    return(phs)

""" for i in range(0,21): subplot(4,6,i+1);extract_phase(filename='/home/bdb112/mydocs/MDF/data/interferometer/4th July 2013/scope_{i}.bin'.format(i=i),loc=0)

# processing the 21 channel interferometer
"""
put=1
pl.rcParams['legend.fontsize']='small'
count = 0
plots=3  # number of plots shown at beginning
shotlist = [int(sys.argv[1])]


def run_demod(shotlist, put = 1, count = 0, plots = 3):
    print(shotlist)
    for shot in shotlist: #(78518,79055): #[77865,77866,77867]:
        try:
            tr = MDS.Tree('electr_dens',shot)
            chnd=tr.getNode('\electr_dens::top.ne_het:chan_list')
            chan_list=chnd.data()
            if count<plots: pl.figure(str(shot))
            for (n,ch) in enumerate(chan_list):  # chan_list[0:1]
                tries = 2000 if len(shotlist)==1 else 2
                for trie in range(tries):
                    try:
                        nd=tr.getNode('\electr_dens::top.camac:'+ch);
                        sig=nd.data() 
                    except Exception,reason:
                        print('failed accessing interf data: \n{r}'
                              .format(r=reason)),
                    else:
                        break
                    if trie < tries:
                        print('- waiting ')
                        sleep(5)
                    else:
                        print('yep')
                        break
                t_raw=nd.dim_of().data()
                (sgn, bw) = (1, .8)
                #if 'A14_22:INPUT_3' in ch: sgn = -1  # why is this different?
                acrms = np.sqrt(sig.var())
                if acrms<0.02: 
                    continue
                elif acrms<0.03: 
                    bw = .1

                # frequency here!!
                phs = sgn*extract_phase(None, sig=sig,t_raw=t_raw,
                                        f=5e3,bw=bw,plot=0)
                if count < plots:
                    pl.plot(t_raw, phs,
                            label=str('{ch}: {rms:.2f}V'.format(ch=ch,rms=acrms)))
                if n==0:
                    count += 1
                    if (count == 5) and (n ==0): print('suppressing further plots')

                if put: 
                    pnode = tr.getNode('\electr_dens::top.ne_het:ne_'+str(n+1))
                    convExpr = MDS.Data.compile("0.35*$VALUE")
                    rawMdsData = MDS.Float32Array(phs)
                    rawMdsData.setUnits("rad")
                    convExpr.setUnits("1e18/m-3")
                    #build the signal object
                    signal = MDS.Signal(convExpr, rawMdsData, nd.dim_of())
                    pnode.putData(signal)
                    if n==0:
                        centrenode = tr.getNode('\electr_dens::top.ne_het:ne_centre')
                        centrenode.putData(signal)
    #    except Exception, reason:
        except  None, reason:
            print('Exception on shot {s}, "{r}"'.format(s=shot,r=reason))

        pl.legend()
#pl.show()    
print('process_interferometer finished')

"""
#dim = MDS.Dimension(MDS.Window(0, 120000, 0), MDS.Range(None, None, 1e-6))
#convExpr = Data.compile("10.* $VALUE/32768."
convExpr = MDS.Data.compile("0.3*$VALUE")
rawMdsData = MDS.Float32Array(sgn*sig)
rawMdsData.setUnits("rad")
convExpr.setUnits("1e18/m-3")
#build the signal object
signal = MDS.Signal(convExpr, rawMdsData, nd.dim_of())
pnode.putData(signal)
"""

