#small edit
#Example Usage (This is outdated, and they probably won't work):
#a=compareFFT(pf.getDevice("PRL").acq.getdata(69819,"H1ToroidalOutward"))
#for i in range(69819,69824):compareFFT(pf.getDevice("PRL").acq.getdata(i,"H1ToroidalOutward"))
#for i in range(69819,69824):a=compareFFT(pf.getDevice("PRL").acq.getdata(i,"H1ToroidalOutward"))
#plotSVD(69862,'H1ToroidalAxial',SVDtime=[0.05,0.0505])
#for i in range(69816,69820):specWidget(pf.getDevice("PRL").acq.getdata(i,"H1ToroidalMirnov_1x"))

#example : specWidget('H1','H1ToroidalMirnov_4x',71000)
#Buttons : 
#csvAv : related to cross spectral density - not sure quite what
#csvAv2 : Related to cross spectral denisty - not sure quite what
#Peak - Doesn't seem to do anything
#PeakS : Plots another spectrogram - not sure what of though
#Fluc : See if there is a fluc struc at the location of the vertical line
#FlucS : Same as above, but it steps through all the time values and outputs the identified fluc strucs for the shot
#SVD : SVD from pyfusion + plot cumulative phase for any identified flucstrucs + filtered time signals, will also print out the main SVs
#FFT : Not working very well...
#Axplot : Axial coils from HMA
#AllCoils : Plot time traces of all coils in HMA
#Pol1 : Plot time traces of all coils in Poloidal array 1
#Pol2 : Plot time traces of all coils in Poloidal array 2

#specWidget('H1','H1ToroidalMirnov_1z',71686) - high frequency shot
import h1.diagnostics.HMA_funcs as HMA_funcs
#import matplotlib.cm as cm
#import scatter_picker
#from scipy import *
import numpy as num
import numpy
import pyfusion as pf
import matplotlib.pyplot as pt
from matplotlib.widgets import Slider, Button, RadioButtons
import MDSplus as MDS
import matplotlib
import copy
font={'family' : 'sans-serif',
      'weight' : 'normal',
      'size' : 12}
matplotlib.rc('font',**font)

#For the text box
#import Tkinter
#import ScrolledText
#import Tkconstants
#import tkFont
import operator
import time
import csv

def checkDigitiserSync(shot):
    mynode1=MDS.Tree("mirnov",shot).getNode('ACQ132_8:INPUT_32')
    mynode2=MDS.Tree("mirnov",shot).getNode('ACQ132_7:INPUT_32')
    mynode3=MDS.Tree("mirnov",shot).getNode('ACQ132_9:INPUT_32')
    voltage1=mynode1.data()
    time1=mynode1.dim_of().data()
    print mynode1.dim_of()
    voltage2=mynode2.data()
    time2=mynode2.dim_of().data()
    print mynode2.dim_of()
    voltage3=mynode3.data()
    time3=mynode3.dim_of().data()
    print mynode3.dim_of()
    kwargs={'marker':'.', 'markersize':12}
    pt.plot(time1,voltage1,'b,',**kwargs)
    pt.plot(time2,voltage2,'k,',**kwargs)
    pt.plot(time3,voltage3,'r,',**kwargs)
    pt.show()

def plot_signals(shot,dataset):
    downsamplefactor=10
    n_columns=3
    fig = pt.figure()
    h1=pf.getDevice("H1")
    input_data=h1.acq.getdata(shot,dataset)

    n_ch = input_data.signal.n_channels()
    n_rows = int(n_ch/n_columns)
    if n_rows<float(n_ch)/n_columns: n_rows+=1
    print str(n_rows) + ' ' + str(n_columns) + ' n_ch ' + str(n_ch)
    axes_list=[]
    count = 0
    for row in range(n_rows):
        for col in range(n_columns):
            print (row)*n_columns+col+1
            count+=1
            if row==0 and col==0:
                axes_list.append(fig.add_subplot(n_rows, n_columns, row*n_columns+col+1))
            else:
                axes_list.append(fig.add_subplot(n_rows, n_columns, row*n_columns+col+1,sharex=axes_list[0],sharey=axes_list[0]))

            if downsamplefactor==1:
                axes_list[-1].plot(input_data.timebase, input_data.signal.get_channel(row*n_columns+col))
                axes_list[-1].axis([-0.01,0.1,-5, 5])
            elif count<=n_ch:
                plotdata=input_data.signal.get_channel(row*n_columns+col)
                timedata=input_data.timebase
                axes_list[-1].plot(timedata[0:len(timedata):downsamplefactor], plotdata[0:len(timedata):downsamplefactor])
                axes_list[-1].axis([-0.01,0.1,-5,5])
                pt.setp(axes_list[-1].get_xticklabels(), visible=False)
                pt.setp(axes_list[-1].get_yticklabels(), visible=False)
    fig.subplots_adjust(hspace=0)#,wspace=0)
    fig.canvas.draw()
    fig.canvas.show()
    #data.plot_signals(downsamplefactor=50,n_columns=3)



def plotSVD(shot,dataset,SVDtime=[0.05,0.0505]):
    try:
        pt.figure()
        h1=pf.getDevice("H1")
        data=h1.acq.getdata(shot,dataset)
        #print data.signal
        data_reduced_time=data.reduce_time(SVDtime).subtract_mean(copy=False).normalise(method='v',separate=True,copy=False)
        #print data_reduced_time
        fs_set=data_reduced_time.flucstruc()
        print fs_set
        title='shot       t0    kHz     Ep  a12     p     H   Strucs\n'
        print title
        for fs in fs_set:
            if fs.p>0.05:
                line= "%d %7.4g %6.3g %6.3f %.2f %.3f %.3f %s\n" % (shot, fs.t0, fs.freq/1000., num.sqrt(fs.E*fs.p), fs.a12, fs.p, fs.H, str(fs.svs()))   
                print line
        svd_data=data_reduced_time.svd()
        svd_data.svdplot()
    except:
        print "Error on Shot "+ str(shot) 





class SpecWidgetData():
    def __init__(self,device,channel,shot, inc_polar_data = 0):
        '''
        Widget for examining a shot - in particular the Mirnov data
        Device : This is the device for pyfusion - probably 'H1'
        channel : Channel as defined in the pyfusion .cfg file. i.e 'H1ToroidalMirnov_1z'
        shot : shot number

        SRH : 1Oct2014
        '''
        print 'in spec widget'
        self.device=device
        self.channel=channel
        self.shot=shot
        self.upward_array_shot=0
        self.outward_array_shot=0
        self.axial_array_shot=0
        #self.MDSTree=MDS.Tree('h1data',self.shot)
        self.inc_polar_data = inc_polar_data
        print 'start data extract'
        self.extract_data()
        print 'finished data extract'
        self.NFFT=4096
        self.fs=1.0/(self.timebase[1]-self.timebase[0])
        self.timeslice=0.01
        self.poloidal1_array_shot=0
        self.axial_array_shot=0
        self.poloidal2_array_shot=0
    def get_polar(self):
        self.polar_signals = []
        self.polar_signals.append(self.pf_device.acq.getdata(self.shot,'H1ToroidalMirnov_1x').signal)
        self.polar_signals.append(self.pf_device.acq.getdata(self.shot,'H1ToroidalMirnov_1y').signal)
        self.polar_signals.append(self.pf_device.acq.getdata(self.shot,'H1ToroidalMirnov_1z').signal)
        print 'finished getting polarisation data'

    def extract_data(self):
        try:
            self.pf_device = pf.getDevice(self.device)
            print 'got device'
            self.pf_data=self.pf_device.acq.getdata(self.shot,self.channel)
            print 'got data'
            self.shot=self.pf_data.meta['shot']
            self.coil_name=self.pf_data.channels.name       
            self.original_data=self.pf_data.signal
            self.timebase=self.pf_data.timebase
            if self.inc_polar_data:
                print 'getting polarisation data'
                self.get_polar()
        except None:
            print('Error extracting new data from %s shot %d, channel %s' %
                  (self.device, self.shot, self.channel))
        count=0
        success=0
        ne_success=0
        kh_success=0
        antenna_success=0
        self.MDSTree = None
        while count<6 and (kh_success==0 or ne_success==0 or antenna_success==0):
            if kh_success==0:
                try:
                    if self.MDSTree == None:
                        self.MDSTree=MDS.Tree('h1data',self.shot)
                    self.main_current=self.MDSTree.getNode('.operations.magnetsupply.lcu.setup_main:I2').data()
                    self.sec_current=self.MDSTree.getNode('.operations.magnetsupply.lcu.setup_sec:I2').data()
                    self.kh=float(self.sec_current)/float(self.main_current)
                    rf_top_node=self.MDSTree.getNode('.rf:i_top')
                    rf_bot_node=self.MDSTree.getNode('.rf:i_bot')
                    rf_bot_node=self.MDSTree.getNode('.rf:i_bot')

                    print 'hello'
                    T2phase_node = self.MDSTree.getNode('.LOG.HEATING.SNMP.T2.OPERATIONAL.LLRF:STALLRFPHD')
                    T1phase_node = self.MDSTree.getNode('.LOG.HEATING.SNMP.T1.OPERATIONAL.LLRF:STALLRFPHD')
                    print 'hello'
                    self.antenna_phase_diff = (T2phase_node.data() - T1phase_node.data())/1000.
                    print self.antenna_phase_diff
                    print 'extracting data'
                    self.rf_top_value = 1
                    self.rf_bot_value = 1
                    self.rf_top_time = 1
                    self.rf_bot_time = 1
                    #self.rf_top_value=rf_top_node.record.data()
                    #self.rf_bot_value=rf_bot_node.record.data()
                    print 'extracting times'
                    #self.rf_top_time=rf_top_node.record.dim_of().data()
                    #self.rf_bot_time=rf_bot_node.record.dim_of().data()
                    #print 'rf_top : ', self.rf_top_value
                    #print 'rf_bot : ', self.rf_bot_value
                    #print 'rf_top time : ', self.rf_top_time
                    #print 'rf_top time: ', self.rf_bot_time


                    #test=self.MDSTree.getNode('.electr_dens:NE_HET:NE_CENTRE')
                    #ne_node=self.MDSTree.getNode('.electr_dens:NE_HET:NE_CENTRE')
                    kh_success=1
                    antenna_success=1
                    print 'got kh value'
                except:
                    print 'Error obtaining kh or antenna values'
                    self.kh=3
                    self.rf_top_time=1
                    self.rf_bot_time=1
                    self.rf_top_value=1
                    self.rf_bot_value=1
                    self.antenna_phase_diff = None
            if ne_success==0:
                try:
                    if self.MDSTree == None:
                        print 'MDSTree need to get again'
                        self.MDSTree = MDS.Tree('h1data',self.shot)
                    print 'getting node'
                    
                    #ne_node=self.MDSTree.getNode('.electr_dens:NE_HET:NE_CENTRE')
                    ne_node=MDS.Tree('electr_dens',self.shot).getNode('.NE_HET:NE_CENTRE')
                    #ne_node=MDS.Tree('electr_dens',self.shot).getNode('.NE_HET:NE_5')
                    print 'got node'
                    self.ne_value=ne_node.record.data()
                    print 'got data'
                    self.ne_time=ne_node.record.dim_of().data()
                    print 'got time'
                    ne_success=1
                except:
                    print 'Error obtaining Ne values'
                    self.ne_value=1
                    self.ne_time=1


            count=count+1
            time.sleep(0.5)

    def extract_axial_array(self):
        if self.axial_array_shot!=self.shot:
            #print 'Extracting Axial Array Data'
            self.axial_data=pf.getDevice("H1").acq.getdata(self.shot,"H1ToroidalAxial")
            self.axial_array_shot=self.shot

    def extract_poloidal1_array(self):
        if self.poloidal1_array_shot!=self.shot:
            #print 'Extracting Poloidal1 Array Data'
            self.poloidal1_data=pf.getDevice("H1").acq.getdata(self.shot,"H1Poloidal1")
            self.poloidal1_array_shot=self.shot

    def extract_poloidal2_array(self):
        if self.poloidal2_array_shot!=self.shot:
            #print 'Extracting Poloidal2 Array Data'
            self.poloidal2_data=pf.getDevice("H1").acq.getdata(self.shot,"H1Poloidal2")
            self.poloidal2_array_shot=self.shot

    def extract_upward_array(self):
        if self.upward_array_shot!=self.shot:
            print 'Extracting Upward Array Data'
            self.upward_data=pf.getDevice("H1").acq.getdata(self.shot,"H1ToroidalUpward")
            self.upward_array_shot=self.shot

    def extract_outward_array(self):
        if self.outward_array_shot!=self.shot:
            print 'Extracting Outward Array Data'
            self.outward_data=pf.getDevice("H1").acq.getdata(self.shot,"H1ToroidalOutward")
            self.outward_array_shot=self.shot

    def average_coherence(self, NFFT,start_time=0.05):
        #self.extract_axial_array()
        #self.average_csd_func3(start_time=start_time)
        signals=self.axial_data.signal.copy()
        timebase=self.axial_data.timebase.copy()

        #!!!!!!!!!!! Need a better way than this
        num_slices=range(0,NFFT*7,100)

        num_freqs=NFFT/2+1
        window_values=numpy.hanning(NFFT)
        window_normalise_factor=(numpy.abs(window_values)**2).sum()
        start_point=numpy.searchsorted(timebase,start_time)
        signal_fft={}
        fft_slices={}
        fft_conj_slices={}
        Pxx={}
        #print 'start Pxx calcs'
        for i in range(0,signals.shape[0]):
            slices=numpy.zeros((len(num_slices),num_freqs),dtype=numpy.complex_)
            Pxx_slices=numpy.zeros((len(num_slices),num_freqs),dtype=float)
            location=0
            for step in num_slices:
                this_slice=(signals[i,start_point+step:start_point+step+NFFT])*window_values
                slices[location,:]=numpy.fft.fft(this_slice,n=NFFT)[:NFFT/2+1]
                Pxx_slices[location,:]=slices[location,:]*numpy.conjugate(slices[location,:])/window_normalise_factor*2/2000000
                location=location+1
            fft_slices[i]=slices
            fft_conj_slices[i]=numpy.conjugate(slices)
            Pxx[i]=numpy.mean(Pxx_slices,axis=0)
            del slices
        #print 'start Pxy calcs'
        Phase={}
        Pxy={}
        Cxy={}
        Cxy_average=numpy.zeros((num_freqs),dtype=float)
        total=0

        for i in range(0,signals.shape[0]):
            for j in range(i+1,signals.shape[0]):
                temp_i=fft_slices[i]
                temp_j=fft_slices[j]
                Pxy_slices=numpy.zeros((len(num_slices),num_freqs),dtype=numpy.complex_)
                location=0
                for step in num_slices:
                    Pxy_slices[location,:] = temp_i[location,:]*numpy.conj(temp_j[location,:])/window_normalise_factor*2/2000000
                    location=location+1
                Pxy[i,j]=numpy.mean(Pxy_slices,axis=0)
                Cxy[i,j]=((numpy.abs(Pxy[i,j]))**2)/(Pxx[i]*Pxx[j])
                Phase[i,j]=numpy.arctan2(Pxy[i,j].imag,Pxy[i,j].real)
                Cxy_average=Cxy_average+Cxy[i,j]
                total=total+1
        Cxy_average=Cxy_average/total
        return Pxx, Pxy, Cxy, Cxy_average, Phase

    def coherence_peak_finding(self,Cxy_average,Phase,NFFT,cut_off=0.7):
        #print 'filter Coherence values'
        Cxy_average_filtered=Cxy_average.copy()
        Phase_filtered=Phase[1,2].copy()
        for i in range(0,Cxy_average_filtered.shape[0]):
            if Cxy_average_filtered[i]>cut_off:
                pass
            else:
                Cxy_average_filtered[i]=0
                Phase_filtered[i]=0
        #Create filtered versions where everything except those above threshold are zero
 
        #Peak identification        
        #print 'start Peak identification'
        freqs=numpy.fft.fftfreq(NFFT,d=1./self.fs)[0:NFFT/2+1]
        freqs[-1]=freqs[-1]*-1 #change sign of last relevant one
        index=2
        identified_peak_list=[]
        identified_number=0
        while index < Cxy_average_filtered.shape[0]:
            if Cxy_average_filtered[index]>0:
                start_point2=index
                index=index+1
                while Cxy_average_filtered[index]> 0 and index < Cxy_average_filtered.shape[0]:
                    index=index+1
                end_point=index-1
                identified_peak_list.append([start_point2,end_point,freqs[start_point2],freqs[end_point]])
            index=index+1
        #print identified_peak_list
        #Finished Generating identified peaks
        return identified_peak_list,Cxy_average_filtered,Phase_filtered


    def coherence_analytic_phase(self,starting_point,ending_point,long_NFFT,signal_fft_values,start_time):
        #expand the window slightly(?? need a better way of doing this)
        print 'hello',starting_point, ending_point
        starting_point=starting_point-15
        ending_point=ending_point+15

        fft_window=numpy.hanning(ending_point-starting_point)
        mask=numpy.zeros(long_NFFT,dtype=complex)
        mask[starting_point:ending_point]=fft_window

        mean_phase_diff=[]
        for i in range(0,len(signal_fft_values)-1):
            signal1_ifft=numpy.fft.ifft(signal_fft_values[i]*mask,n=long_NFFT)
            signal2_ifft=numpy.fft.ifft(signal_fft_values[i+1]*mask,n=long_NFFT)
            temp=signal1_ifft*numpy.conjugate(signal2_ifft)
            if i==0:
                phase_diff=numpy.zeros([14,len(temp)],dtype=float)
            phase_diff[i,:]=numpy.arctan2(temp.imag,temp.real)
            phase_diff_time=numpy.arange(start_time,start_time+long_NFFT*1./self.fs,1./self.fs)
            mean_phase_diff.append(numpy.mean(phase_diff[i,:]))

        #print mean_phase_diff
        phase_plot_csd=[0]
        phase_location_csd=[0]
        for j in range(0,len(mean_phase_diff)):
            phase_location_csd.append(j+1)
            phase_plot_csd.append(mean_phase_diff[j]/(2*numpy.pi)+phase_plot_csd[j])
        return phase_diff, phase_diff_time,mean_phase_diff,phase_location_csd,phase_plot_csd

    def plot_coherent_phase(self,Cxy_average,Phase,Cxy_average_filtered,Phase_filtered,Pxx,freqs):
        coherent_phase_fig=matplotlib.pyplot.figure()
        ax=coherent_phase_fig.add_subplot(311)
        ax2=ax.twinx()
        ax.plot(freqs/1000,Cxy_average,'r',label='CoUnfilt')
        ax2.plot(freqs/1000,Phase[1,2],label='PhaseUnfilt')
        ax.set_title('Average Unfiltered Covariance and Phase')
        ax.set_xlim([0,self.fs/1000/2])

        ax3=coherent_phase_fig.add_subplot(312)
        ax4=ax3.twinx()
        ax3.plot(freqs/1000,Cxy_average_filtered,'r',label='CoherenceFilt')
        ax4.plot(freqs/1000,Phase_filtered,label='PhaseFilt')
        ax3.set_title('Filtered Coherence and Phase')
        ax3.set_xlim([0,self.fs/1000/2])
        ax5=coherent_phase_fig.add_subplot(313)
        ax6=ax5.twinx()

        for i in range(0,len(Pxx)):
            ax5.plot(freqs/1000,Pxx[i],'k')
        ax5.set_xlim([0,self.fs/1000/2])
        ax6.plot(freqs/1000,Cxy_average_filtered,'r',label='CoherenceFilt')
        ax5.set_title('PSD for some channels and CoherenceFilt')
        ax5.set_xlim([0,self.fs/1000/2])
        coherent_phase_fig.canvas.draw()
        coherent_phase_fig.show()

    #Plot the Analytic Phase Results for each Peak        
    def plot_analytic_phase(self,phase_diff, phase_diff_time,freqs2,signal_fft_values,starting_point,ending_point,phase_location_csd,phase_plot_csd):
        fig9=pt.figure()
        #fig10=pt.figure()
        ax50,ax60,ax70=(fig9.add_subplot(221),fig9.add_subplot(223),fig9.add_subplot(122))

        for i in range(0,len(phase_diff)):
            ax60.plot(phase_diff_time,phase_diff[i,:],label='%d,%d,$\phi$=%.2f'%(i,i+1,num.mean(phase_diff)))
        ax50.plot(freqs2/1000,num.absolute(signal_fft_values[1]))
        ax50.plot(freqs2[starting_point]/1000,num.absolute(signal_fft_values[1][starting_point]),'or')
        ax50.plot(freqs2[ending_point]/1000,num.absolute(signal_fft_values[1][ending_point]),'ok')
        ax50.set_title('Start Freq %.1fkHz, End Freq %.1fkHz'%(freqs2[starting_point]/1000,freqs2[ending_point]/1000))
        ax50.set_xlim([0,self.fs/1000/2])

        ax70.plot(phase_location_csd,phase_plot_csd,'o-',label='Phase')
        ax70.set_ylim([-3.5,0])
        ax70.grid(b=True)
        x=num.arange(0,15,0.5)
        #self.plot_mode_phase_lines(ax70,x)
        ax70.legend(ncol=1,loc=3)
        ax70.set_title('Shot %d, Coil Phase, %.1f-%.1fkHz'%(self.shot,freqs2[starting_point]/1000,freqs2[ending_point]/1000))
        ax70.set_xlabel('Coil Number/Location')
        ax70.set_ylabel('Cumulative Phase (2pi rad)')
        ax70.set_ylim([-4,0])
        ax60.set_xlabel('time')
        ax60.set_ylabel('Phase (rad)')
        ax60.set_ylim([-3.14,3.14])
        ax60.set_title('AnalyticPhase')
        ax60.grid(b=True)
        #ax60.legend(ncol=3)
        fig9.canvas.draw()
        fig9.show()


    def average_csd_func2(self,start_time=0.05,coherent_plot=0,an_phase_plot=0):
        #print 'running average coherence'
        NFFT=self.NFFT/2
        freqs=numpy.fft.fftfreq(n=NFFT,d=1./self.fs)[0:NFFT/2+1]
        freqs[-1]=freqs[-1]*-1 #change sign of last relevant one
        #self.extract_axial_array()
        signals=self.axial_data.signal.copy()
        timebase=self.axial_data.timebase.copy()
        start_point=numpy.searchsorted(timebase,start_time)

        #Calculate,filter and plot the average coherent phase
        Pxx, Pxy, Cxy, Cxy_average, Phase = self.average_coherence(NFFT,start_time=start_time)
        identified_peak_list,Cxy_average_filtered,Phase_filtered=self.coherence_peak_finding(Cxy_average,Phase,NFFT)
        if coherent_plot==1:
            self.plot_coherent_phase(Cxy_average,Phase,Cxy_average_filtered,Phase_filtered,Pxx,freqs)

        #Generate the FFT for the Analytic Phase calculation
        long_NFFT=NFFT*8
        window_values2=numpy.hanning(long_NFFT)
        freqs2=numpy.fft.fftfreq(n=long_NFFT,d=1./self.fs)
        signal_fft_values={}
        for i in range(0,signals.shape[0]):
            signal_fft_values[i]=numpy.fft.fft((signals[i,start_point:start_point+long_NFFT])*window_values2,n=long_NFFT)

        #Generate the Analytic Phase for each identified peak
        output_data=[]
        for a in identified_peak_list:
            starting_point=numpy.searchsorted(freqs2,a[2])
            ending_point=numpy.searchsorted(freqs2,a[3])

            #Obtain the Analytic Phase Results for each Peak and then plot the results
            phase_diff, phase_diff_time,mean_phase_diff,phase_location_csd,phase_plot_csd=self.coherence_analytic_phase(starting_point,ending_point,long_NFFT,signal_fft_values,start_time)
            if an_phase_plot==1:
                self.plot_analytic_phase(phase_diff, phase_diff_time,freqs2,signal_fft_values,starting_point,ending_point,phase_location_csd,phase_plot_csd)
            output_current=[a[2],a[3]]
            for i in range(0,len(mean_phase_diff)):
                output_current.append(mean_phase_diff[i])
            output_data.append(output_current)
            #output_data.append([(a[2]+a[3])/2.,numpy.mean(phase_diff),numpy.std(phase_diff)])
        return output_data

    def average_csd_func3(self,start_time=0.05):
        self.extract_axial_array()
        signal_array=self.axial_data.signal.copy()
        signal_timebase=self.axial_data.timebase.copy()
        time=[start_time,start_time+2*self.NFFT*1.0/self.fs]
        reduced_signal_array,reduced_signal_timebase=self.reduce_time_shaun(signal_array,signal_timebase,time)
        PxxList=[]

        #create a cache of all PSDs to limit computation time - but use more memory!
        #create all the ffts

        for i in range(0,reduced_signal_array.shape[0]):
            PxxReturn,freqs=matplotlib.mlab.psd(reduced_signal_array[i,:],NFFT=self.NFFT/2,Fs=self.fs)
            PxxList.append(PxxReturn)
        n=0
        for i in range(0,reduced_signal_array.shape[0]):
            for j in range(i+1,reduced_signal_array.shape[0]):
                signal1=reduced_signal_array[i,:]
                signal2=reduced_signal_array[j,:]
                Pxy,freqs=matplotlib.mlab.csd(signal1,signal2,NFFT=self.NFFT/2,Fs=self.fs)
                if i==0 and j==1:
                    result=num.divide(num.absolute(Pxy)**2,(PxxList[i]*PxxList[j]))
                result=num.divide(num.absolute(Pxy)**2,(PxxList[i]*PxxList[j]))+result
                n=n+1
        result=result/n
        fig=pt.figure()
        ax=fig.add_subplot(111)
        ax.plot(freqs,result)
        ax.set_xlim([0,400000])
        ax.set_ylim([0,1])
        ax.set_title('Average Spectral Coherence : Shot %d,Time : %.2fms'%(self.shot,start_time*1000))
        fig.canvas.draw()
        fig.show()

    def average_csd_func(self,start_time=0.05):
        fig=pt.figure()
        ax=fig.add_subplot(311)
        self.extract_axial_array()
        signal_array=self.axial_data.signal
        signal_timebase=self.axial_data.timebase
        time=[start_time,start_time+self.NFFT*1.0/self.fs]
        reduced_signal_array,reduced_signal_timebase=self.reduce_time_shaun(signal_array,signal_timebase,time)
        runs=0
        for i in range(0,reduced_signal_array.shape[0],3):
            for j in range(i,reduced_signal_array.shape[0],3):
                csd_power,csd_freqs=ax.csd(reduced_signal_array[i,:],reduced_signal_array[j,:],NFFT=self.NFFT,Fs=self.fs)
                #csd_power,csd_freqs=matplotlib.mlab.csd(reduced_signal_array[i,:],reduced_signal_array[j,:],NFFT=self.NFFT,Fs=self.fs)
                if i==0 and j==0:
                    csd_average=num.zeros(csd_power.shape[0],dtype=float)
                csd_total=csd_average+csd_power
                runs=runs+1
        csd_average=csd_total/runs
        ax2=fig.add_subplot(312)
        ax2.plot(csd_freqs,csd_average)
        ax3=fig.add_subplot(313)
        for j in range(0,reduced_signal_array.shape[0]):
            psd_single_power,psd_single_freq=ax3.psd(reduced_signal_array[j,:],NFFT=self.NFFT,Fs=self.fs)
            if j==0:
                psd_array=num.zeros([reduced_signal_array.shape[0],psd_single_power.shape[0]],dtype=float)
            psd_array[j,:]=psd_single_power
        ax3.cla()
        for j in range(0,reduced_signal_array.shape[0]):
            ax3.plot(psd_single_freq,psd_array[j,:])
        ax3.set_xlim([0,300000])
        ax2.set_xlim([0,300000])
        ax.set_xlim([0,300000])
        #fig.canvas.draw()
        return csd_freqs,csd_average,fig

    def reduce_time_shaun(self,signals,timebase, bounds):
        output_signals=signals[:,numpy.searchsorted(timebase,bounds[0]):numpy.searchsorted(timebase,bounds[1])]
        output_timebase=timebase[numpy.searchsorted(timebase,bounds[0]):numpy.searchsorted(timebase,bounds[1])]
        return output_signals,output_timebase

    def coherence_step_through(self,shots):
        f=csv.writer(open('/home/srh112/outputlog.txt','a'),delimiter=' ')
        for j in shots:
            self.shot=j

            print '-------------------Shot %d------------------' %(j)
            tries=0
            success=0
            while tries<10 and success==0:
                try:
                    self.extract_data()
                    self.extract_axial_array()
                    success=1
                except:
                    print 'Error extracting data on Shot : %d'%(j)
                    tries=tries+1
                    time.sleep(0.5)
            for i in num.arange(0,0.07,0.002):
                ne=self.ne_value[num.searchsorted(self.ne_time,i)]
                output_data=self.average_csd_func2(start_time=i,coherent_plot=0,an_phase_plot=0)
                for a in output_data:
                    mean_phase_array=num.array([a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],a[10],a[11],a[13],a[14],a[15]])
                    string= 'Shot = %d, Time : %d ms , kh = %.2f, f = %.2f, Phase-mean = %.2f, Phase-std = %.2f, ne = %.2f'%(self.shot,i*1000,self.kh,(a[0]+a[1])/2.,num.mean(mean_phase_array),num.std(mean_phase_array),ne)
                    print string
                    f.writerow([self.shot,i*1000,self.kh,ne,a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],a[10],a[11],a[13],a[14],a[15]])
                    #f.write(string+'\n')
        del f
                
    def readback_text_file(self,filename='/home/srh112/outputlog.txt',ratio_cutoff=0.5):
        csv_file=csv.reader(open(filename,'rb'),delimiter=' ')
        file_contents=[]
        for a in csv_file:
            file_contents.append(a)

        fig=pt.figure()
        ax=fig.add_subplot(411)
        ax1=fig.add_subplot(412)
        ax2=fig.add_subplot(413)
        ax3=fig.add_subplot(414)
        fig2=pt.figure()
        phase_plot=fig2.add_subplot(111)
        total=0
        for i in file_contents:
            current_i=num.array(i,dtype=float)
            mean_phase_array=current_i[6:20]
            ne=float(current_i[3])
            mean_phase=num.abs(num.mean(mean_phase_array))
            std_phase=num.abs(num.std(mean_phase_array))
            frequency=(current_i[4]+current_i[5])/2.0
            time=current_i[1]
            kh=current_i[2]
            #print 'total :%d, mean: %.2f, std: %.2f' %(total,mean_phase,std_phase)
            frequency_normalised=frequency*(num.abs(ne)**0.5)
            if time>10:
                ax.plot(kh,frequency,'k.')
                ax1.plot(kh,ne,'k.')
                ax2.plot(kh,frequency_normalised,'k.')
                #a=num.abs(float(current_i[4]))
                #b=float(current_i[5])
                if (mean_phase/std_phase)>ratio_cutoff:
                    print 'total :%d, mean: %.2f, std: %.2f, ne: %.2f' %(total,mean_phase,std_phase,ne)
                    total=total+1
                    ax3.plot(kh,frequency_normalised,'k.')
                    phase_plot.plot(kh,mean_phase,'k.')
        #ax.set_xlabel('kh')
        ax.set_ylabel('freq')
        ax.set_title('For time > 10ms')
#        ax.set_ylim([0,160000])
#        ax1.set_ylim([0,2.5])
#        ax2.set_ylim([0,160000*(2**0.5)])
#        ax3.set_ylim([0,160000*(2**0.5)])

        ax.set_xlim([0,1.1])
        ax1.set_xlim([0,1.1])
        ax2.set_xlim([0,1.1])
        ax3.set_xlim([0,1.1])

        #ax1.set_xlabel('kh')
        ax1.set_ylabel('freq')
        ax1.set_title('For time > 10ms and mean phase/std phase > %.2f'%(ratio_cutoff))
        ax3.set_xlabel('kh')
        ax3.set_ylabel('freq')
        ax3.set_title('For time > 10ms and mean phase/std phase > %.2f'%(ratio_cutoff/2))
        ax.grid(b=True)
        ax1.grid(b=True)
        ax2.grid(b=True)
        ax3.grid(b=True)
        fig.canvas.draw()
        fig.show()
        del csv_file


    def plot_rationals(self,axis):
        def rotational_transform(radius,kh):
            coefficients=num.array([[1.24403098, 0.29927867,-0.04178176, -0.0113835 , 0.01371373],
                                   [-0.06438457, 0.17743677,-0.00568132, 0.11426079, -0.0981305],
                                   [0.16757832, -0.41083898,0.00136293,-0.23903926,  0.22891545],
                                   [-0.21602304,  0.16208048,0.05840499,0.1875845 , -0.21617175],
                                   [0.12705246, -0.00544844, -0.03210589,-0.05116255,  0.07173953]])
            c = num.array([-0.25938664,  0.34786773])
            total=num.zeros(len(radius),dtype=float)
            for i in range(0,5):
                for j in range(0,5):
                    component=coefficients[i,j]*((radius-c[0])**i)*((kh-c[1])**j)
                    total=total+component            
            return total

        #Generate the values
        radius=num.arange(0, 1, 0.0005)
        rotational_transforms=num.zeros((len(num.arange(0,1.2,0.001)),len(radius)),dtype=float)
        kh_values_record=num.zeros((len(num.arange(0,1.2,0.001)),1),dtype=float)
        location=0
        for kh_input in num.arange(0,1.2,0.001):
            answer=rotational_transform(radius,kh_input)
            rotational_transforms[location,:]=answer
            kh_values_record[location]=kh_input
            location = location + 1

        #plot the values
        ax=axis.twinx()
        #axis.imshow(rotational_transforms.transpose(),cmap=cm.gray,origin='lower',extent=[0,1.2,0,1])

        label_dictionary={4./3:'4/3',5./4:'5/4',6./5: '6/5',7./5:'7/5',7./6:'7/6'}
        
        #ax.contour(rotational_transforms.transpose(),[4./3,5./4,6./5,7./5,7./6],extent=(0,1.2,0,1),labels=True,label_fmt=label_dictionary)
        #contour_levels=[4./3,5./4,6./5,7./5,7./6]
        contour_levels=[]
        for temp in label_dictionary:
            contour_levels.append(temp)
        CS=ax.contour(rotational_transforms.transpose(),contour_levels,colors='k',extent=(0,1.2,0,1))
        pt.clabel(CS,contour_levels,inline=0,fmt=label_dictionary,fontsize=11,inline_spacing=0)

        ax.set_xlim([0,1.2])
        ax.set_ylim([0,1])
        #axis.set_ylim([0,1])


        

    def readback_flucstrucs_text_file(self,filename='/home/srh112/outputlog_flucstruc.txt',ratio_cutoff=0.1,time_cutoff=0):
        csv_file=csv.reader(open(filename,'rb'),delimiter=' ')
        file_contents=[]
        for a in csv_file:
            file_contents.append(a)

        fig=pt.figure()
        ax=fig.add_subplot(411)
        ax1=fig.add_subplot(412)
        ax2=fig.add_subplot(413)
        ax3=fig.add_subplot(414)
        fig2=pt.figure()
#        phase_plot=fig2.add_subplot(111)
#        iota_plot=phase_plot.twinx()

        phase_plot=fig2.add_subplot(311)
        iota_plot=fig2.add_subplot(312)
        iota_plot_ne_scaled=fig2.add_subplot(313)

        iota_kh_plot=ax.twinx()
        iota_kh_scaled_plot=ax3.twinx()

        freq_values_unfilt=[]
        kh_values_unfilt=[]
        mean_phase_values_unfilt=[]
        freq_values_unfilt_normalised=[]
        ne_values_unfilt=[]


        def filter_array(input_array,index,value,inequality='equal'):
            first_time=0
            for i in range(0,input_array.shape[0]):
                if inequality=='equal':
                    condition_check=input_array[i,index]==value
                elif inequality=='greater':
                    condition_check=input_array[i,index]>=value
                elif inequality=='less':
                    condition_check=input_array[i,index]<=value
                if condition_check==True:
                    if first_time==0:
                        output_array=input_array[i,:].copy()
                        first_time=1
                        print 'first time'
                    elif first_time==1:
                        output_array=num.append([output_array],[input_array[i,:]],axis=0)
                        first_time=2
                    else:
                        output_array=num.append(output_array,[input_array[i,:]],axis=0)
            return output_array


        first_time=0

        #Association of mode numbers with colours
        scale_value=2.*num.pi/15.
        rat_num_dict={'8,7': [4.333*scale_value,'b','o','--',0],
                            '9,7': [4*scale_value,'k','o','--',1],
                            '7,6': [3.667*scale_value,'y','o','--',2],
                            '6,5': [3*scale_value,'k','o','-',3],
                            '7,5': [2.667*scale_value,'b','o','-',4],
                            '5,4': [2.333*scale_value,'b','^','-',5],
                            '4,3': [1.667*scale_value,'y','o','-',6],
                            'unclassified':[0.5,'k','.','-',7]}

        dictionary={'shot':0,'time':1,'kh':2,'ne':3,'fsP':4,'a12':5,'fsH':6,'freq':7,'RMS':8,
                    'fsE':9,'num_SV':10,'num_phases':11,'mean_phase':12,'mode':13}

        for i in file_contents:
            current_i=num.array(i,dtype=float)
            mean_phase_array=current_i[12:12+current_i[dictionary['num_phases']]]

            #manual fix to move "strange" +ve values back to negative - is this valid?:
            for temp_locator in range(0,len(mean_phase_array)):
                if mean_phase_array[temp_locator]>num.pi*0.5:
                    mean_phase_array[temp_locator]=mean_phase_array[temp_locator]-2*num.pi

            mean_phase=num.mean(mean_phase_array)*-1.
            std_phase=num.abs(num.std(mean_phase_array))

            frequency_normalised=current_i[dictionary['freq']]*(num.abs(current_i[dictionary['ne']])**0.5)

            if current_i[dictionary['time']]>time_cutoff and current_i[dictionary['fsP']]>ratio_cutoff:
                if current_i[dictionary['num_SV']]>=2:
                    #find the closest mode
                    overall_min=1000
                    for testing_tmp in rat_num_dict:
                        current_min=num.abs(num.abs(mean_phase)-num.abs(rat_num_dict[testing_tmp][0]))
                        if overall_min>current_min:
                            overall_min=current_min
                            best_fit_mode=rat_num_dict[testing_tmp][4]
                    new_line=num.append(current_i[0:12],mean_phase)
                    new_line=num.append(new_line,best_fit_mode)
                    new_line=num.append(new_line,mean_phase_array)
                    #build the large_array
                    if first_time==0:
                        large_array=new_line.copy()
                        first_time=1
                    elif first_time==1:
                        large_array=num.append([large_array],[new_line],axis=0)
                        first_time=2
                    else:
                        large_array=num.append(large_array,[new_line],axis=0)

                #non SV answers
                freq_values_unfilt.append(current_i[dictionary['freq']])
                kh_values_unfilt.append(current_i[dictionary['kh']])
                freq_values_unfilt_normalised.append(frequency_normalised)
                mean_phase_values_unfilt.append(mean_phase)
                ne_values_unfilt.append(current_i[dictionary['ne']])

        print 'finished, the large array info is'
        print large_array.shape
        

        #Plots
        ax.plot(large_array[:,dictionary['kh']],large_array[:,dictionary['freq']],'bo')

        #unfiltered sv values
        ax1.plot(kh_values_unfilt,ne_values_unfilt,'k.')
        ax2.plot(kh_values_unfilt,freq_values_unfilt_normalised,'k.')

        #can do this one
        ax2.plot(large_array[:,dictionary['kh']],large_array[:,dictionary['freq']]*(num.abs(large_array[:,dictionary['ne']])**0.5),'r.')        
        ax3.plot(large_array[:,dictionary['kh']],large_array[:,dictionary['freq']]*(num.abs(large_array[:,dictionary['ne']])**0.5),'r.')
        phase_plot.scatter(large_array[:,dictionary['kh']],large_array[:,dictionary['mean_phase']],
                           c='k',marker='+',s=large_array[:,dictionary['fsP']]/1.*25,label='all')
        iota_plot.scatter(large_array[:,dictionary['kh']],large_array[:,dictionary['freq']],
                           c='k',marker='+',s=large_array[:,dictionary['fsP']]/1.*25,label='all')
        iota_plot_ne_scaled.scatter(large_array[:,dictionary['kh']],(large_array[:,dictionary['freq']]*(large_array[:,dictionary['ne']]**0.5)),
                           c='k',marker='+',s=large_array[:,dictionary['fsP']]/1.*25,label='all')

        #plot the positive and negative horizontal lines
        for i in rat_num_dict:
            if i != 'unclassified':
                phase_plot.axhline(y=rat_num_dict[i][0],color=rat_num_dict[i][1],linestyle=rat_num_dict[i][3])
                phase_plot.axhline(y=rat_num_dict[i][0]*-1,color=rat_num_dict[i][1],linestyle=rat_num_dict[i][3])
                phase_plot.text(0.01,rat_num_dict[i][0],i,fontsize=9)
                pass

        #plot the mode coloured dots
        for placeholder in ['6,5','5,4','4,3']:
            mode_filtered=filter_array(large_array,dictionary['mode'],rat_num_dict[placeholder][4],inequality='equal')
            mkr_size=mode_filtered[:,dictionary['fsP']]/1.*25
            phase_plot.scatter(mode_filtered[:,dictionary['kh']],mode_filtered[:,dictionary['mean_phase']],
                               s=mkr_size, c=rat_num_dict[placeholder][1], marker=rat_num_dict[i][2],label=placeholder)
            iota_plot.scatter(mode_filtered[:,dictionary['kh']],mode_filtered[:,dictionary['freq']],
                              s=mkr_size, c=rat_num_dict[placeholder][1], marker=rat_num_dict[i][2],label=placeholder)
            iota_plot_ne_scaled.scatter(mode_filtered[:,dictionary['kh']],(mode_filtered[:,dictionary['freq']]*(mode_filtered[:,dictionary['ne']]**0.5)),s=mkr_size, c=rat_num_dict[placeholder][1], marker=rat_num_dict[i][2],label=placeholder)


    
        self.plot_rationals(iota_plot)
        self.plot_rationals(iota_plot_ne_scaled)
        self.plot_rationals(phase_plot)


        iota_kh_plot.set_xlim([0,1.3])
        iota_kh_scaled_plot.set_xlim([0,1.3])
        iota_plot.set_xlim([0,1.2])
        iota_plot_ne_scaled.set_xlim([0,1.2])
        phase_plot.set_xlim([0,1.2])
        ax.set_xlim([0,1.1]);ax1.set_xlim([0,1.1]);ax2.set_xlim([0,1.1]);ax3.set_xlim([0,1.1])

        iota_kh_plot.grid(b=True)
        iota_kh_scaled_plot.grid(b=True)
        iota_plot.grid(b=True)
        phase_plot.grid(b=True)
        iota_plot_ne_scaled.grid(b=True)

        ax.grid(b=True);ax1.grid(b=True);ax2.grid(b=True);ax3.grid(b=True)

        phase_plot.set_ylim([0,2])
        ax.set_ylim([0,250000])
        ax1.set_ylim([0,2.5])
        iota_plot_ne_scaled.set_ylim([0,95000])
        iota_plot.set_ylim([0,160000])

        #iota_plot.set_ylim([0,300000])


        iota_plot.legend()
        iota_plot_ne_scaled.legend()

        phase_plot.legend()

        #iota_plot.set_ylabel('Normalised Radius')

        phase_plot.set_xlabel('kh')
        ax3.set_xlabel('kh')

        phase_plot.set_title('Mode Identification fsP> %.2f, t> %.2fms'%(ratio_cutoff, time_cutoff))

        phase_plot.set_ylabel('Average Phase Diff (rad)')
        ax.set_ylabel('freq (Hz)')
        ax1.set_ylabel('ne (10^18)')
        ax3.set_ylabel('freq * ne^0.5')
        ax2.set_ylabel('freq * ne^0.5')
        iota_plot_ne_scaled.set_ylabel('Hz x sqrt(ne)')
        iota_plot.set_ylabel('Hz')


        ax.set_title('Time > %d ms, p> %.2f'%(time_cutoff,ratio_cutoff))
        ax1.set_title('Time > 10ms and p>%.2f' %(ratio_cutoff))


        
#ax1.set_xlabel('kh')
        def cluster_plot_clicked(event):
            if event.button==1:
                print event.button
                print 'click right button, left button is for zooming'
            else:
                print event.button
                print 'correct button clicked'
                if event.inaxes==phase_plot:
                    print 'hello phase plot'
                    kh_click_value=num.round(event.xdata,2)
                    print kh_click_value
                    print event.ydata
                if event.inaxes==iota_plot:
                    print 'hello phase plot'
                    kh_click_value=num.round(event.xdata,2)
                    kh_filtered=filter_array(large_array,dictionary['kh'],kh_click_value,inequality='equal')
                    print kh_filtered.shape
                    min_location=(num.abs(kh_filtered[:,dictionary['freq']]-event.ydata)).argmin()

                    print event.xdata
                    print event.ydata

                    print '---------- nearest match found-----------'
                    print kh_filtered[min_location,:]
                    print 'freq'
                    print kh_filtered[min_location,dictionary['freq']]
                    print 'kh'
                    print kh_filtered[min_location,dictionary['kh']]

                    relevant_phases=kh_filtered[min_location,14:]
                    print relevant_phases
                    new_fig=pt.figure()
                    phase_axis=new_fig.add_subplot(111)
                    cumul_phases=[0]
                    for i in range(0,len(relevant_phases)):
                        cumul_phases.append(cumul_phases[-1]+relevant_phases[i]/(2*num.pi))
                    print cumul_phases
                    phase_axis.plot(range(0,len(cumul_phases)),cumul_phases,'bo-')
                    self.plot_mode_phase_lines(phase_axis,num.arange(0,15,0.1))
                    phase_axis.set_xlim([0,15])
                    phase_axis.set_ylim([-4,0.5])
                    title_string='Shot:%d, kh:%.2f, t:%dms, nSV:%d, freq:%.2fkHz,\n a12:%.2f, E:%.2f, ne:%.2f'%(kh_filtered[min_location,dictionary['shot']], kh_filtered[min_location,dictionary['kh']],kh_filtered[min_location,dictionary['time']], kh_filtered[min_location,dictionary['num_SV']], kh_filtered[min_location,dictionary['freq']]/1000, kh_filtered[min_location,dictionary['a12']], kh_filtered[min_location,dictionary['fsP']], kh_filtered[min_location,dictionary['ne']])
                    phase_axis.set_title(title_string)
                    phase_axis.grid(b=True)
                    phase_axis.legend()
                    new_fig.canvas.draw()

                if event.inaxes==iota_plot_ne_scaled:
                    print 'hello iota plot ne scaled'
                    kh_click_value=num.round(event.xdata,2)

                    print event.xdata
                    print event.ydata
                    iota_plot_ne_scaled.plot(kh_click_value,event.ydata,'yo')

                fig2.canvas.draw()
        #fig2.add_axes(iota_plot_ne_scaled)
        fig2.add_axes(iota_plot)
        

        fig.canvas.draw()
        fig.show()
        fig2.canvas.mpl_connect('button_press_event', cluster_plot_clicked)

        del csv_file


    def coherence_step_through(self,shots,step_value=0.01):
        f=csv.writer(open('/home/srh112/outputlog.txt','a'),delimiter=' ')
        for j in shots:
            self.shot=j

            print '-------------------Shot %d------------------' %(j)
            tries=0
            success=0
            while tries<10 and success==0:
                try:
                    self.extract_data()
                    self.extract_axial_array()
                    success=1
                except:
                    print 'Error extracting data on Shot : %d'%(j)
                    tries=tries+1
                    time.sleep(0.5)
            for i in num.arange(0,0.07,step_value):
                ne=self.ne_value[num.searchsorted(self.ne_time,i)]
                output_data=self.average_csd_func2(start_time=i,coherent_plot=0,an_phase_plot=0)
                for a in output_data:
                    mean_phase_array=num.array([a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],a[10],a[11],a[13],a[14],a[15]])
                    string= 'Shot = %d, Time : %d ms , kh = %.2f, f = %.2f, Phase-mean = %.2f, Phase-std = %.2f, ne = %.2f'%(self.shot,i*1000,self.kh,(a[0]+a[1])/2.,num.mean(mean_phase_array),num.std(mean_phase_array),ne)
                    print string
                    f.writerow([self.shot,i*1000,self.kh,ne,a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],a[10],a[11],a[13],a[14],a[15]])
                    #f.write(string+'\n')
        del f
        

    def flucstruc_step_through(self,shots,step_value,array='H1ToroidalAxial',file_name='/home/srh112/outputlog_flucstruc_misc.txt',samples=4096,time_bounds=[0,0.08]):
        f=csv.writer(open(file_name,'a'),delimiter=' ')
        for j in shots:
            start_shot_time=time.time()
            self.shot=j
            #print '-------------------Shot %d------------------' %(j)
            count = 0
            success = 0
            ne_success = 0
            kh_success = 0
            while count<6 and kh_success==0 and ne_success==0:
                if kh_success==0:
                    try:
                        self.MDSTree=MDS.Tree('h1data',self.shot)
                        self.main_current=self.MDSTree.getNode('.operations.magnetsupply.lcu.setup_main:I2').data()
                        self.sec_current=self.MDSTree.getNode('.operations.magnetsupply.lcu.setup_sec:I2').data()
                        self.kh=float(self.sec_current)/float(self.main_current)
                        kh_success=1
                    except:
                        print 'Error obtaining kh value'
                        self.kh=3
                if ne_success==0:
                    try:
                        print 'getting node'
                        ne_node=self.MDSTree.getNode('.electr_dens:NE_HET:NE_CENTRE')
                        print 'got node'
                        self.ne_value=ne_node.record.data()
                        print 'got data'
                        self.ne_time=ne_node.record.dim_of().data()
                        print 'got time'
                        ne_success=1
                    except:
                        print 'Error obtaining Ne values'
                        self.ne_value=1
                        self.ne_time=1
                count=count+1
                time.sleep(0.5)

            tries=0
            success=0
            while tries<10 and success==0:
                print array
                try:
                    data=pf.getDevice('H1').acq.getdata(self.shot,array)
                    '''
                    if array=='poloidal1':
                        self.extract_poloidal1_array() #!!!!!!!!!!!
                    elif array=='poloidal2':
                        self.extract_poloidal2_array()
                    else:
                        self.extract_axial_array()
                    '''
                    success=1
                except:
                    print 'Error extracting data on Shot : %d'%(j)
                    tries=tries+1
                    time.sleep(0.5)
            data_extract_time=time.time()-start_shot_time
            '''
            if array=='poloidal1':
                data=copy.deepcopy(self.poloidal1_data)
            elif array=='poloidal2':
                data=copy.deepcopy(self.poloidal2_data)                    
            else:
                data=copy.deepcopy(self.axial_data)
            '''
            time_bounds=[0.005,0.08]
            data_reduced_time=data.reduce_time(time_bounds,copy=True).subtract_mean(copy=False).normalise(method='v',separate=True,copy=False)

            for t_seg in data_reduced_time.segment(samples,overlap=4):
                #print 'hello'
                time_seg_start_time=t_seg.timebase[0]
                time_seg_end_time=t_seg.timebase[-1]
                time_seg_average_time=num.mean([time_seg_start_time,time_seg_end_time])

                times=[]
                times_label=[]

                times.append(time.time())
                times_label.append(' start :')
                previous=time.time()

                times.append(time.time()-previous)
                times_label.append(' copy+reduce :')
                previous=time.time()

                fs_set=t_seg.flucstruc()
                #fs_set=data_reduced_time.flucstruc()

                times.append(time.time()-previous)
                times_label.append(' fluccalc :')
                previous=time.time()

                #######!!!!!
                ne=self.ne_value[num.searchsorted(self.ne_time,time_seg_average_time)]
                for fs in fs_set:
                    if fs.p>0.0:
                        phases=num.zeros(len(fs.dphase),dtype=float)
                        csv_row=[]
                        csv_row.append(self.shot)
                        csv_row.append(time_seg_average_time*1000)
                        csv_row.append(self.kh)
                        csv_row.append(ne)
                        csv_row.append(fs.p)
                        csv_row.append(fs.a12)
                        csv_row.append(fs.H)
                        csv_row.append(fs.freq)
                        csv_row.append((num.mean(data_reduced_time.scales**2))**0.5)
                        csv_row.append(fs.E)
                        csv_row.append(len(fs.svs()))
                        csv_row.append(len(fs.dphase))
                        for fs_phase in range(0,len(fs.dphase)):
                            if num.abs(fs.dphase[fs_phase].delta)<0.001:
                                phases[fs_phase]=0
                            else:
                                phases[fs_phase] = fs.dphase[fs_phase].delta
                            csv_row.append(phases[fs_phase])
                        f.writerow(csv_row)
                        del csv_row
                times.append(time.time()-previous)
                times_label.append(' extract fs :')
                previous=time.time()

                #print 'Shot Times : %.4f to %.4f; comp times:' %(time_seg_start_time,time_seg_end_time),
                for time_spot in range(1,len(times_label)):
                    pass
                    #print ' %s %.3fs ' %(times_label[time_spot],times[time_spot]),
                #print '\n',
                #del times_label
                #del times
            print 'Finish shot : %d , Data Extract: %.2f, Computation: %.2f'%(self.shot,data_extract_time,time.time()-start_shot_time-data_extract_time)
        del f

        
    def plot_raw_signals(self, fig_time_phase, filtered_data, data_reduced_time_axial, fs, freq_fft, phase_diff, include_phase_plot=1, pub_fig = 0, decimate=1):
        cumul_phase_4_3 = num.array([0.00,-0.14,-0.26,-0.35,-0.44,-0.53,-0.62,-0.72,-0.82,-0.94,-1.07,-1.21,-1.36,-1.51,-1.66])
        counter = 0
        if pub_fig:
            cm_to_inch=0.393701
            fig_time_phase.set_figwidth(8.48*cm_to_inch)
            fig_time_phase.set_figheight(8.48*1.3*cm_to_inch*1.2*1./decimate)
            import matplotlib as mpl
            old_rc_Params = mpl.rcParams
            mpl.rcParams['font.size']=8.0
            mpl.rcParams['axes.titlesize']=8.0#'medium'
            mpl.rcParams['xtick.labelsize']=8.0
            mpl.rcParams['ytick.labelsize']=8.0
            mpl.rcParams['lines.markersize']=1.0
            mpl.rcParams['savefig.dpi']=300

        axes_list = []
        if include_phase_plot==0:
            ncols = 1
        else:
            ncols =2
        plot_channels = range(0,len(filtered_data.signal), decimate)
        fluc_struc_phases = [i.delta for i in fs.dphase]
        print len(fluc_struc_phases), len(plot_channels), len(phase_diff)
        for plot_num, i in enumerate(plot_channels):
            if i==0:
                axes_list.append(fig_time_phase.add_subplot(len(plot_channels),ncols, plot_num *ncols+ ncols))
                #axes_list[-1].set_title('Raw(red) filt(black), SVD(blue) FFT(yellow)\n Zoom is synchronised')
            else:
                axes_list.append(fig_time_phase.add_subplot(len(plot_channels), ncols, plot_num*ncols+ncols, sharex=axes_list[0],sharey=axes_list[0]))

            print i, plot_num *ncols+ ncols, phase_diff[i], fluc_struc_phases[i-1]
            axes_list[-1].plot(filtered_data.timebase,filtered_data.signal[i],'k')
            #axes_list[-1].plot(filtered_data_upward.timebase,filtered_data_upward.signal[i],'y')

            axes_list[-1].plot(data_reduced_time_axial.timebase,data_reduced_time_axial.signal[i],'r')

            axes_list[-1].set_ylabel('%d'%(i+1))
            if i==0:
                vert_line_SVD=num.mean([data_reduced_time_axial.timebase[0],data_reduced_time_axial.timebase[-1]])
                vert_line_FFT=num.mean([data_reduced_time_axial.timebase[0],data_reduced_time_axial.timebase[-1]])
            else:
                #vert_line_SVD = vert_line_SVD + num.sum([fs.dphase[i-(j+1)].delta/(2*num.pi)/fs.freq for j in range(decimate)])
                #print len(fs.dphase), len(plot_channels)
                vert_line_SVD -= fluc_struc_phases[i-1]/(2*num.pi)/fs.freq
                vert_line_FFT -= (phase_diff[i])/(2*num.pi)/freq_fft
                print vert_line_SVD - vert_line_FFT, freq_fft, fs.freq
            if counter == 0:
                simul_mode_initial = vert_line_SVD
                simul_mode_phase = simul_mode_initial
                counter+=decimate

            else:
                simul_mode_phase = simul_mode_initial + cumul_phase_4_3[counter]/freq_fft
                counter+=decimate
            #print '############'
            #print vert_line_SVD, vert_line_FFT
            #print fs.freq, freq_fft
            #vert_lines_SVD = vert_line_SVD + num.arange(-10,10)*1./fs.freq
            vert_lines_SVD = vert_line_SVD + num.arange(-10,10)*1./fs.freq
            vert_lines_FFT = vert_line_FFT + num.arange(-10,10)*1./freq_fft
            vert_lines_simul = simul_mode_phase + num.arange(-10,10)*1./freq_fft
            #print vert_lines_SVD, vert_lines_FFT
            #print '############'
            for tmp_locator in range(0,len(vert_lines_SVD)):
                axes_list[-1].axvline(x=vert_lines_SVD[tmp_locator],color='b')
                axes_list[-1].axvline(x=vert_lines_simul[tmp_locator],color='k',linestyle='--')
                axes_list[-1].axvline(x=vert_lines_FFT[tmp_locator],color='y')
            #remove ticks (should include ticks on the bottom one...
            pt.setp(axes_list[-1],yticklabels=[])
            if i<(len(filtered_data.signal)-1):
                pt.setp(axes_list[-1].get_xticklabels(), visible=False)
                #pt.setp(axes_list[0],xticklabels=[])
            else:
                axes_list[-1].set_xlabel('Time(s)')
            axes_list[-1].set_xlim([data_reduced_time_axial.timebase[0],data_reduced_time_axial.timebase[-1]])
        return axes_list




    def plot_raw_signals_single_ax(self, fig_time_phase, filtered_data, data_reduced_time_axial, fs, freq_fft, phase_diff, include_phase_plot=1, pub_fig = 0, decimate=1):
        cumul_phase_4_3 = num.array([0.00,-0.14,-0.26,-0.35,-0.44,-0.53,-0.62,-0.72,-0.82,-0.94,-1.07,-1.21,-1.36,-1.51,-1.66])
        cumul_phase_5_4 = num.array([0.00,-0.19,-0.36,-0.50,-0.62,-0.74,-0.88,-1.02,-1.16,-1.32,-1.51,-1.70,-1.91,-2.12,-2.32])
        cumul_phase_6_5 = num.array([0.00,-0.25,-0.46,-0.64,-0.80,-0.96,-1.13,-1.31,-1.50,-1.71,-1.94,-2.20,-2.46,-2.73,-2.99])
        counter = 0
        if pub_fig:
            cm_to_inch=0.393701
            fig_time_phase.set_figwidth(8.48*cm_to_inch)
            fig_time_phase.set_figheight(8.48*1.3*cm_to_inch*1.2*1./2)
            import matplotlib as mpl
            old_rc_Params = mpl.rcParams
            mpl.rcParams['font.size']=8.0
            mpl.rcParams['axes.titlesize']=8.0#'medium'
            mpl.rcParams['xtick.labelsize']=8.0
            mpl.rcParams['ytick.labelsize']=8.0
            mpl.rcParams['lines.markersize']=1.0
            mpl.rcParams['savefig.dpi']=300

        ax = fig_time_phase.add_subplot(111)
        plot_channels = range(0,len(filtered_data.signal), decimate)
        linear_line_SVD = []
        linear_line_y = []
        simul_line_4_3 = []
        simul_line_5_4 = []
        simul_line_6_5 = []
        for plot_num, i in enumerate(plot_channels):
            plot_sig1 = (filtered_data.signal[i] - num.mean(filtered_data.signal[i]))
            plot_sig1_rms = num.sqrt(num.mean(plot_sig1**2))
            plot_sig1 = plot_sig1 / plot_sig1_rms/3 + i + 1
            ax.plot(filtered_data.timebase*1000,plot_sig1,'k')
            #axes_list[-1].plot(filtered_data_upward.timebase,filtered_data_upward.signal[i],'y')
            plot_sig1 = (data_reduced_time_axial.signal[i] - num.mean(data_reduced_time_axial.signal[i]))
            plot_sig1_rms = num.sqrt(num.mean(plot_sig1**2))
            plot_sig1 = (plot_sig1)/plot_sig1_rms/3 + i+1
            ax.plot(data_reduced_time_axial.timebase*1000,plot_sig1, 'r')
            #axes_list[-1].set_ylabel('%d'%(i+1))
            if i==0:
                vert_line_SVD=num.mean([data_reduced_time_axial.timebase[0],data_reduced_time_axial.timebase[-1]])
                vert_line_FFT=num.mean([data_reduced_time_axial.timebase[0],data_reduced_time_axial.timebase[-1]])
            else:
                vert_line_SVD = vert_line_SVD + num.sum([fs.dphase[i-(j+1)].delta/(2*num.pi)/fs.freq for j in range(decimate)])
                #vert_line_SVD=vert_line_SVD+fs.dphase[i-1].delta/(2*num.pi)/fs.freq
                vert_line_FFT=vert_line_FFT+(phase_diff[i])/(2*num.pi)/freq_fft
            if counter == 0:
                simul_mode_initial = vert_line_SVD
                simul_mode_phase_4_3 = simul_mode_initial
                simul_mode_phase_5_4 = simul_mode_initial
                simul_mode_phase_6_5 = simul_mode_initial
                counter+=decimate

            else:
                simul_mode_phase_4_3 = simul_mode_initial + cumul_phase_4_3[counter]/freq_fft
                simul_mode_phase_5_4 = simul_mode_initial + cumul_phase_5_4[counter]/freq_fft
                simul_mode_phase_6_5 = simul_mode_initial + cumul_phase_6_5[counter]/freq_fft
                counter+=decimate
            #print '############'
            #print vert_line_SVD, vert_line_FFT
            #print fs.freq, freq_fft
            linear_line_SVD.append(vert_line_SVD)
            linear_line_y.append(i+1)
            simul_line_4_3.append(simul_mode_phase_4_3)
            simul_line_5_4.append(simul_mode_phase_5_4)
            simul_line_6_5.append(simul_mode_phase_6_5)

        for i in num.arange(-10,10):
            if i==-10:
                ax.plot((num.array(linear_line_SVD) + i*1./fs.freq)*1000, linear_line_y,'b-',label='expt')
                ax.plot((num.array(simul_line_4_3) + i*1./fs.freq)*1000, linear_line_y,'k--', label = '-4,3')
                #ax.plot((num.array(simul_line_5_4) + i*1./fs.freq)*1000, linear_line_y,'y-', label = '-5,4')
            else:
                ax.plot((num.array(linear_line_SVD) + i*1./fs.freq)*1000, linear_line_y,'b-')
                ax.plot((num.array(simul_line_4_3) + i*1./fs.freq)*1000, linear_line_y,'k--')
                #ax.plot((num.array(simul_line_5_4) + i*1./fs.freq)*1000, linear_line_y,'y-')
        ax.legend(loc='upper right')
            #ax.plot((num.array(simul_line_6_5) + i*1./fs.freq)*1000, linear_line_y,'k--')
            #vert_lines_SVD = vert_line_SVD + num.arange(-10,10)*1./fs.freq
            #vert_lines_FFT = vert_line_FFT + num.arange(-10,10)*1./freq_fft
            #vert_lines_simul = simul_mode_phase + num.arange(-10,10)*1./freq_fft
            #print vert_lines_SVD, vert_lines_FFT
            #print '############'
            #for tmp_locator in range(0,len(vert_lines_SVD)):
            #    axes_list[-1].axvline(x=vert_lines_SVD[tmp_locator],color='b')
            #    axes_list[-1].axvline(x=vert_lines_simul[tmp_locator],color='k',linestyle='--')
                #axes_list[-1].axvline(x=vert_lines_FFT[tmp_locator],color='y')
            #remove ticks (should include ticks on the bottom one...
            #pt.setp(axes_list[-1],yticklabels=[])
            #if i<(len(filtered_data.signal)-1):
            #    pt.setp(axes_list[-1].get_xticklabels(), visible=False)
                #pt.setp(axes_list[0],xticklabels=[])
            #else:
            #    axes_list[-1].set_xlabel('Time(s)')
            #axes_list[-1].set_xlim([data_reduced_time_axial.timebase[0],data_reduced_time_axial.timebase[-1]])
        ax.set_xlim([0.01292*1000,0.01312*1000])
        ax.grid()
        return ax



    def output_fluc_strucs(self,shot,dataset,svd_time=[0.05,0.0505],plot='YES',amplitude_plot='YES',sep_normalisation=True):
        data_axial=copy.deepcopy(self.axial_data)
        print 'copied data'
        #data_reduced_time_axial=data_axial.reduce_time(svd_time,copy=True).subtract_mean(copy=False)
        data_reduced_time_axial=data_axial.reduce_time(svd_time,copy=True).subtract_mean(copy=False).normalise(method='v',separate=sep_normalisation,copy=False)
        print 'reduced time'
        fs_set=data_reduced_time_axial.flucstruc()
        print 'got fs'
        print fs_set
        title='%-7s %-7s %-7s %-7s %-7s %-7s %-7s %-s' %('shot','t0','kHz','Ep','a12','p', 'H', 'Strucs')
        #self.stext.insert(Tkconstants.END, title + '\n')
        fluc_struc_list=[]
        for fs in fs_set:
            if fs.p>0.03:
                phases = '|'.join(["%3.2f" % fs.dphase[j].delta for j in range(0,len(fs.dphase))])

                fluc_struc_list.append([shot, fs.t0, fs.freq/1000., num.sqrt(fs.E*fs.p), fs.a12, fs.p, fs.H, str(fs.svs()),phases])
                #if the fluc struc has 2 SVs-need a better indication...
                if len(str(fs.svs()))>4:
                    #fig10=pt.figure()

#Generate the time - phase plot----------------------------------------------------------
#                    self.extract_upward_array()

#                    data_upward=copy.deepcopy(self.upward_data)
#                    data_reduced_time_upward=data_upward.reduce_time(svd_time,copy=True).subtract_mean(copy=False).normalise(method='v',separate=True,copy=False)

                    temp1=copy.deepcopy(data_reduced_time_axial)
                    filtered_data=temp1.sp_filter_butterworth_bandpass([fs.freq-2000,fs.freq+2000],[fs.freq-10000,fs.freq+10000],3.,10.)
#                    filtered_data_upward=data_reduced_time_upward.sp_filter_butterworth_bandpass([fs.freq-2000,fs.freq+2000],[fs.freq-10000,fs.freq+10000],3.,10.)

                    #filtered_data=data_reduced_time_axial.sp_filter_butterworth_bandpass([fs.freq-2000,fs.freq+2000],[fs.freq-10000,fs.freq+10000],3.,10.,copy=True)

                    temp_timebase=data_reduced_time_axial.timebase
                    temp_signals=data_reduced_time_axial.signal

                    freq_fft,phase,phase_diff,cumulative_phase,amplitude=self.plot_phases(temp_signals,temp_timebase,fs.freq)#!!!!!!!!!!
                    del temp_signals,temp_timebase

                    fig_time_phase=pt.figure()
                    #axes_list=[]
                    axes_list = self.plot_raw_signals(fig_time_phase, filtered_data, data_reduced_time_axial, fs, freq_fft, phase_diff)
                    fig_time_phase2 = pt.figure()
                    axes_list2 = self.plot_raw_signals(fig_time_phase2, filtered_data, data_reduced_time_axial, fs, freq_fft, phase_diff,include_phase_plot=0)
                    fig_time_phase2.canvas.draw(); fig_time_phase2.show()
                    self.relative_amplitude_plot=0 #Temporary Measure!!!!!!!!!!!!!!
                    #setup the figure axes dependent on the two cases
                    if self.relative_amplitude_plot==1: #!!!!!!!!!!!!!!
                        phase_plot_axis=fig10.add_subplot(121)
                        phase_plot_axis.grid(b=True)
                        amplitude_plot_axis=fig10.add_subplot(122)
                        amplitude_plot_axis.grid(b=True)
                    else:
                        #phase_plot_axis=fig10.add_subplot(111)
                        phase_plot_axis=fig_time_phase.add_subplot(1,2,1)
                        phase_plot_axis.grid(b=True)
                    phase_plot=[0]
                    phase_location=[0]
                    for j in range(0,len(fs.dphase)):
                        phase_location.append(j+1)
                        phase_plot.append(fs.dphase[j].delta/(2*num.pi)+phase_plot[j])
                    phase_plot_axis.plot(phase_location,phase_plot,'o-',label='SVD')
                    phase_plot_axis.set_title('Shot=%d, kh=%.2f, time=%dms, freq=%.2fkHz\n Strucs=%s,a12=%.2f, E=%.1f' %(self.shot,self.kh,fs.t0*1000,fs.freq/1000.,str(fs.svs()),fs.a12,fs.p))
                    phase_plot_axis.set_xlabel('Coil Position')
                    phase_plot_axis.set_ylabel('Cumulative Phase (2pi rad)')
                    phase_plot_axis.set_ylim([0,3.5])
                    x=num.arange(0,15,0.5)
                    scaling_value=15
                    self.plot_mode_phase_lines(phase_plot_axis,x)#!!!!!!!!!!!!!!!!

                    signals,timebase=self.reduce_time_shaun(self.axial_data.signal,self.axial_data.timebase,svd_time)
                    freq_fft,phase,phase_diff,cumulative_phase,amplitude=self.plot_phases(signals,timebase,fs.freq)#!!!!!!!!!!
                    if self.relative_amplitude_plot==1: #!!!!!!!!!!!!
                        self.extract_upward_array()
                        self.extract_outward_array()
                        signals_upward,timebase_upward=self.reduce_time_shaun(self.upward_data.signal,self.upward_data.timebase,svd_time)
                        signals_outward,timebase_outward=self.reduce_time_shaun(self.outward_data.signal,self.outward_data.timebase,svd_time)
                        freq_fft, phase_upward,phase_diff_upward,cumulative_phase_upward,amplitude_upward=self.plot_phases(signals_upward,timebase_upward,fs.freq)
                        freq_fft,phase_outward,phase_diff_outward,cumulative_phase_outward,amplitude_outward=self.plot_phases(signals_outward,timebase_outward,fs.freq)
                        amplitude_total=(amplitude**2+amplitude_outward**2+amplitude_upward**2)**0.5

                    phase_plot_axis.plot(cumulative_phase/2/num.pi,'ko-',label='rawfft')

                    if self.relative_amplitude_plot==1:
                        rel_amplitude_axis=amplitude_plot_axis.twinx()
                        width=0.25
                        bar_left = range(len(amplitude_total))
                        amplitude_plot_axis.bar(bar_left,amplitude_total,width,color='y',label='Total')
                        amplitude_plot_axis.bar(bar_left,amplitude,width,color='k',label='Axial')
                        rel_amplitude_axis.plot(amplitude/(amplitude_total),'bo-',label='Ratio')
                        amplitude_plot_axis.set_xlim([0,16])
                        rel_amplitude_axis.set_ylim([0,1])
                        amplitude_plot_axis.legend()

                    #amplitude_plot_axis.xticks(in
                    phase_plot_axis.set_xlim([0,16])
                    phase_plot_axis.legend(ncol=2,loc='lower right')
                    fig_time_phase.canvas.draw(); fig_time_phase.show()
        fluc_struc_list.sort(key=operator.itemgetter(5),reverse=True)
        print title
        for i in fluc_struc_list:
            line= "%-7d %-7.4g %-7.3g %-7.3f %-7.2f %-7.3f %-7.3f %-s %s" % (i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7],i[8])   
            print line
        if plot=='YES':
            pt.figure()
            svd_data=data_reduced_time_axial.svd()
            svd_data.svdplot()                    
        return None


    def plot_mode_phase_lines(self,plot_axis,array_of_coils):
        scaling_value=15
        plot_axis.plot(array_of_coils,1.667*array_of_coils/scaling_value,'y',label='4,3 mode')
        plot_axis.plot(array_of_coils,2.333*array_of_coils/scaling_value,'r',label='5,4 mode')
        plot_axis.plot(array_of_coils,2.667*array_of_coils/scaling_value,'b',label='7,5 mode')
        plot_axis.plot(array_of_coils,3*array_of_coils/scaling_value,'k',label='6,5 mode')
        plot_axis.plot(array_of_coils,3.667*array_of_coils/scaling_value,'m',label='7,6 mode')
        plot_axis.plot(array_of_coils,4*array_of_coils/scaling_value,'k--',label='9,7 mode')
        plot_axis.plot(array_of_coils,4.333*array_of_coils/scaling_value,'b-.',label='8,7 mode')

    def plot_phases(self,orig_signals,orig_timebase,freq):
        amplitude=num.zeros(orig_signals.shape[0],dtype=float)
        phase=num.zeros(orig_signals.shape[0],dtype=float)
        for i in range(0,orig_signals.shape[0]):
            fft_transform=num.fft.fft(orig_signals[i,:])
            freq_array=num.fft.fftfreq(len(fft_transform),d=orig_timebase[1]-orig_timebase[0])
            freq_fft=freq_array[num.searchsorted(freq_array[0:len(freq_array)/2],freq)]
            tmp_val = fft_transform[num.searchsorted(freq_array[0:len(freq_array)/2],freq)]
            #real=fft_transform[num.searchsorted(freq_array[0:len(freq_array)/2],freq)].real
            #imag=fft_transform[num.searchsorted(freq_array[0:len(freq_array)/2],freq)].imag
            amplitude[i]=num.abs(tmp_val)#(real**2+imag**2)**0.5
            phase[i]=num.angle(tmp_val)
            #phase[i]=num.arctan2(real,imag)
            #tmp = num.arctan2(imag,real)
            #phase[i]=tmp

            phase_diff=num.zeros(orig_signals.shape[0],dtype=float)
            cumulative_phase=num.zeros(orig_signals.shape[0],dtype=float)
            for i in range(1,len(phase)):
                phase_diff[i]=phase[i]-phase[i-1]
                while phase_diff[i]>num.pi:
                    phase_diff[i]=phase_diff[i]-2*num.pi
                while phase_diff[i]<(-1*num.pi):
                    phase_diff[i]=phase_diff[i]+2*num.pi
                cumulative_phase[i]=cumulative_phase[i-1]+phase_diff[i]
        return freq_fft,phase, phase_diff, cumulative_phase, amplitude



class specWidget:
    def __init__(self,device,channel,shot):
        #self.root=Tkinter.Tk()
        #self.stext=ScrolledText.ScrolledText(self.root,bg='white',height=10,width=300)
        #self.f=tkFont.Font(family='times',size=10)
        #self.stext.tag_config('font', font=self.f)
        #self.stext.pack(fill=Tkconstants.BOTH,side=Tkconstants.LEFT,expand=True)
        print 'initialising'
        self.data_class=SpecWidgetData(device,channel,shot)
        print 'data_class made'
        self.clim=[-90,0]
        self.freq_axis=[0,self.data_class.fs/2]
        self.noverlap=11.0/16
        self.amplitude=10000
        self.relative_amplitude_plot=0
        self.plot_ne_scaling=0
        self.ne_scaling_axis_range=[-2,0]
        self.plot_flucs = 0
        self.plot_polar = 0

        #Create Axes
        self.fig=pt.figure(figsize=(10,8))
        self.spec_axis=self.fig.add_axes([0.25,0.1,0.3,0.45])
        #self.ne_scaling_axis=self.spec_axis.twinx()

        #pt.setp(self.ne_scaling_axis,yticklabels=[])


        self.fft_axis=self.fig.add_axes([0.6,0.1,0.3,0.45],sharey=self.spec_axis)
        self.ne_axis=self.fig.add_axes([0.25,0.65,0.3,0.2],sharex=self.spec_axis)
        self.antenna_axis=self.ne_axis.twinx()


        self.signal_axis=self.fig.add_axes([0.6,0.65,0.3,0.2],sharex=self.spec_axis)
        self.ax_fft_amp=self.fig.add_axes([0.25,0.02,0.3,0.03])
        self.ax_clim_upper=self.fig.add_axes([0.05,0.8,0.1,0.05])
        self.ax_clim_lower=self.fig.add_axes([0.05,0.73,0.1,0.05])
        self.ax_nfft_radio=self.fig.add_axes([0.05,0.2,0.1,0.15])
        self.ax_noverlap_radio=self.fig.add_axes([0.05,0.05,0.1,0.15])
        self.ax_relative_amplitude_radio=self.fig.add_axes([0.05,0.66,0.1,0.05])
        self.ax_frequency_axis=self.fig.add_axes([0.25,0.95,0.3,0.03])
        self.ax_frequency_axis2=self.fig.add_axes([0.6,0.95,0.3,0.03])
        self.ax_shot_plus_one=self.fig.add_axes([0.6,0.02,0.05,0.03])
        self.ax_shot_minus_one=self.fig.add_axes([0.65,0.02,0.05,0.03])
        self.ax_shot_plus_ten=self.fig.add_axes([0.7,0.02,0.05,0.03])
        self.shot_minus_ten=self.fig.add_axes([0.75,0.02,0.05,0.03])
        self.shot_zero=self.fig.add_axes([0.80,0.02,0.05,0.03])
        self.shot_go=self.fig.add_axes([0.85,0.02,0.05,0.03])


        self.ax_svd_spawn=self.fig.add_axes([0.05, 0.49, 0.05, 0.03])
        self.ax_axial_plot_spawn=self.fig.add_axes([0.05, 0.45, 0.05, 0.03])
        self.ax_fft_compare_spawn=self.fig.add_axes([0.11, 0.49, 0.05, 0.03])
        self.ax_all_coils_plot_spawn=self.fig.add_axes([0.11, 0.45, 0.05, 0.03])
        self.ax_poloidal1_plot_spawn=self.fig.add_axes([0.11, 0.41, 0.05, 0.03])
        self.ax_poloidal2_plot_spawn=self.fig.add_axes([0.05, 0.41, 0.05, 0.03])


        self.ax_fluc_calc=self.fig.add_axes([0.05, 0.53, 0.05, 0.03])
        self.ax_fluc_step=self.fig.add_axes([0.11, 0.53, 0.05, 0.03])
        self.ax_peak_finding=self.fig.add_axes([0.05, 0.57, 0.05, 0.03])
        self.ax_peak_finding_step=self.fig.add_axes([0.11, 0.57, 0.05, 0.03])
        self.ax_csd_averaging=self.fig.add_axes([0.05, 0.61, 0.05, 0.03])
        self.ax_csd_averaging2=self.fig.add_axes([0.11, 0.61, 0.05, 0.03])

        self.plot_ne_scaling_ax=self.fig.add_axes([0.05, 0.90, 0.03, 0.08])
        self.inc_flucstrucs=self.fig.add_axes([0.08, 0.90, 0.03, 0.08])
        self.plot_polar=self.fig.add_axes([0.11, 0.90, 0.03, 0.08])


        self.plot_ne_scaling_slider_ax=self.fig.add_axes([0.15, 0.95, 0.1, 0.02])
        self.plot_ne_scaling_slider2_ax=self.fig.add_axes([0.15, 0.92, 0.1, 0.02])
        self.plot_ne_scaling_slider_time_ax=self.fig.add_axes([0.15, 0.89, 0.1, 0.02])

        #Create Widgets
        self.clim_slider1=Slider(self.ax_clim_upper,'',0,200,valinit=141)
        self.clim_slider2=Slider(self.ax_clim_lower,'',-130,0,valinit=-105)
        self.fft_amp=Slider(self.ax_fft_amp,'Freq',0,10000, valinit=500)
        self.frequency_axis_slider=Slider(self.ax_frequency_axis,'',0,self.data_class.fs/2*0.99,valinit=0)
        self.frequency_axis_slider2=Slider(self.ax_frequency_axis2,'',0,self.data_class.fs/2,valinit=self.data_class.fs/2)
        self.shot_plus_one=Button(self.ax_shot_plus_one,'+1',hovercolor='0.975')
        self.shot_minus_one=Button(self.ax_shot_minus_one,'-1',hovercolor='0.975')
        self.shot_plus_ten=Button(self.ax_shot_plus_ten,'+10',hovercolor='0.975')
        self.shot_minus_ten=Button(self.shot_minus_ten,'-10',hovercolor='0.975')
        self.shot_zero=Button(self.shot_zero,'0',hovercolor='0.975')
        self.shot_go=Button(self.shot_go,'go',hovercolor='0.975')



        self.svd_spawn_button=Button(self.ax_svd_spawn,'SVD',hovercolor='0.975')
        self.axial_plot_spawn_button=Button(self.ax_axial_plot_spawn,'AxPlot',hovercolor='0.975')
        self.ax_fft_spawn_button=Button(self.ax_fft_compare_spawn,'FFT',hovercolor='0.975')
        self.ax_all_coils_plot_spawnbutton=Button(self.ax_all_coils_plot_spawn,'AllCoils',hovercolor='0.975')
        self.ax_fluc_calcbutton=Button(self.ax_fluc_calc,'Fluc',hovercolor='0.975')
        self.ax_fluc_step=Button(self.ax_fluc_step,'FlucS',hovercolor='0.975')
        self.ax_csd_averaging_button=Button(self.ax_csd_averaging,'csdAv',hovercolor='0.975')
        self.ax_csd_averaging_button2=Button(self.ax_csd_averaging2,'csdAv2',hovercolor='0.975')
        self.ax_poloidal1_plot_button=Button(self.ax_poloidal1_plot_spawn,'Pol1',hovercolor='0.975')
        self.ax_poloidal2_plot_button=Button(self.ax_poloidal2_plot_spawn,'Pol2',hovercolor='0.975')

        self.peak_finding_step=Button(self.ax_peak_finding_step,'PeakS',hovercolor='0.975')
        self.peak_finding=Button(self.ax_peak_finding,'Peak',hovercolor='0.975')
        self.nfft_radio=RadioButtons(self.ax_nfft_radio,('16384','8192','4096','2048','1024','512'),active=2)
        self.noverlap_radio=RadioButtons(self.ax_noverlap_radio,('15','11','9','5', '0'),active=1)
        self.relative_amplitude_radio=RadioButtons(self.ax_relative_amplitude_radio,('1','0'),active=1)
        self.plot_ne_scaling_radio=RadioButtons(self.plot_ne_scaling_ax,('1','0'),active=1)
        self.inc_flucstrucs_radio=RadioButtons(self.inc_flucstrucs,('1','0'),active=1)
        self.plot_polar_radio=RadioButtons(self.plot_polar,('1','0'),active=1)

        self.plot_ne_scaling_slider=Slider(self.plot_ne_scaling_slider_ax,'',0,5,valinit=1)
        self.plot_ne_scaling_slider2=Slider(self.plot_ne_scaling_slider2_ax,'',1,4,valinit=2.5)
        self.plot_ne_scaling_slider_time=Slider(self.plot_ne_scaling_slider_time_ax,'',-10,10,valinit=0)



        self.fig.canvas.draw()


        #Callback Definitions
        self.clim_slider1.on_changed(self.update_clim)
        self.clim_slider2.on_changed(self.update_clim)
        self.frequency_axis_slider.on_changed(self.update_freq_axis)
        self.frequency_axis_slider2.on_changed(self.update_freq_axis)
        self.nfft_radio.on_clicked(self.update_nfft)
        self.noverlap_radio.on_clicked(self.update_noverlap)
        self.relative_amplitude_radio.on_clicked(self.update_relative_amplitude)
        self.shot_plus_one.on_clicked(self.increment_shot_plus_one)
        self.shot_minus_one.on_clicked(self.increment_shot_minus_one)
        self.shot_plus_ten.on_clicked(self.increment_shot_plus_ten)
        self.shot_minus_ten.on_clicked(self.increment_shot_minus_ten)
        self.shot_zero.on_clicked(self.increment_shot_zero)
        self.shot_go.on_clicked(self.increment_go)


        self.svd_spawn_button.on_clicked(self.svd_spawn)
        self.axial_plot_spawn_button.on_clicked(self.axial_plot_spawn)
        self.ax_fft_spawn_button.on_clicked(self.fft_spawn)
        self.ax_all_coils_plot_spawnbutton.on_clicked(self.all_coils_plot_spawn)
        self.ax_fluc_calcbutton.on_clicked(self.fluc_calc)
        self.ax_fluc_step.on_clicked(self.fluc_step)
        self.fft_amp.on_changed(self.update_fft_amp)
        self.peak_finding.on_clicked(self.peak_finding_single)
        self.peak_finding_step.on_clicked(self.peak_finding_step_func)
        self.ax_csd_averaging_button.on_clicked(self.average_csd_button)
        self.ax_csd_averaging_button2.on_clicked(self.average_csd_button2)
        self.ax_poloidal2_plot_button.on_clicked(self.poloidal2_plot_spawn)
        self.ax_poloidal1_plot_button.on_clicked(self.poloidal1_plot_spawn)

        self.plot_ne_scaling_radio.on_clicked(self.update_ne_scaling_plot)
        self.inc_flucstrucs_radio.on_clicked(self.update_inc_flucstrucs)
        self.plot_polar_radio.on_clicked(self.update_plot_polar_radio)

        self.plot_ne_scaling_slider.on_changed(self.update_ne_scaling_axis_range)
        self.plot_ne_scaling_slider2.on_changed(self.update_ne_scaling_axis_range)
        self.plot_ne_scaling_slider_time.on_changed(self.update_ne_scaling_axis_range)



        #Initialise plots
        #self.update_spec()
        self.update_fft()
        self.shot_text_update=pt.figtext(0.90,0.02,str(self.data_class.shot),figure=self.fig)
        self.redraw_figure(default = 1, update_spec = 1)
        self.fig.canvas.mpl_connect('button_press_event', self.figure_clicked)

    def update_relative_amplitude(self,label):
        self.relative_amplitude_plot=int(label)

    def update_ne_scaling_plot(self,label):
        self.plot_ne_scaling=int(label)
        self.redraw_figure()

    def update_inc_flucstrucs(self,label):
        self.plot_flucs = int(label)
        print 'update inc flucs', self.plot_flucs
        self.redraw_figure(update_spec = 1, get_flucs = int(label))

    def update_plot_polar_radio(self,label):
        self.data_class.inc_polar_data = int(label)
        print 'update plot polar', self.data_class.inc_polar_data
        if self.data_class.inc_polar_data:
            self.data_class.get_polar()
        self.update_fft()
        self.redraw_figure()


    def average_csd_button(self,event):
        self.data_class.average_csd_func(start_time=self.data_class.timeslice)

    def average_csd_button2(self,event):
        self.data_class.average_csd_func2(start_time=self.data_class.timeslice,coherent_plot=1,an_phase_plot=1)

    def peak_finding_step_func(self,event):
        fig=pt.figure()
        ax=fig.add_subplot(111)
        a,b,c,spectrogram_plot=ax.specgram(self.data_class.original_data[num.searchsorted(self.data_class.timebase,0):num.searchsorted(self.data_class.timebase,0.1)],NFFT=self.data_class.NFFT,Fs=self.data_class.fs,noverlap=int(self.noverlap*self.data_class.NFFT),window=num.hanning(self.data_class.NFFT))
        for i in num.arange(0,0.08,0.001):
            print 'time : %.2fms' %(i*1000)
            returned_peaks=self.peak_finding_func(start_time=i)
            for b in returned_peaks:
                ax.plot(i,b[0],'k.')
        ax.set_ylim([0,500000])
        fig.canvas.draw()

    def peak_finding_single(self,event):
        self.peak_finding_func(start_time=self.data_class.timeslice)

    def peak_finding_func(self,start_time=0.05):
        fft_xaxis,fft_yaxis,returned_fig=self.data_class.average_csd_func(start_time=start_time)
        pt.close(returned_fig)
        overall_average=num.average(fft_yaxis)
        min_freq=0
        step=10000
        start_location=num.searchsorted(fft_xaxis,min_freq)
        end_location=num.searchsorted(fft_xaxis,min_freq+step)
        peaks=[]
        print fft_xaxis
        print fft_yaxis
        print len(fft_xaxis)
        print len(fft_yaxis)
        while end_location<len(fft_xaxis):
            current_frequencies=num.array(fft_xaxis[start_location:end_location])
            current_fft=num.array(fft_yaxis[start_location:end_location])
            print 'size of arrays'
            print len(current_frequencies)
            print len(current_fft)
            print current_fft.argmax()
            if (current_fft.max()/num.average(current_fft))>3:
                a=[current_frequencies[current_fft.argmax()],current_fft.max()]
                if a in peaks: 
                    pass
                else:
                    appended=0
                    too_small=0
                    for b in peaks:
                        if a[0]>b[0]-10000 and a[0]<b[0]+10000:
                            if a[1]>b[1]:
                                peaks[peaks.index(b)]=a
                                appended=1
                                #print 'Start: %d, End: %d, Average: %.2f, Max: %.2f,Ratio: %.2f, Frequency: %.2f'%(min_freq,min_freq+step,num.average(current_fft),current_fft.max(),current_fft.max()/num.average(current_fft),current_frequencies[current_fft.argmax()])
                            else:
                                too_small=1
                    if appended!=1 and too_small==0:
                        peaks.append(a)
                        appended=0
                        too_small==0
                        #print 'Start: %d, End: %d, Average: %.2f, Max: %.2f,Ratio: %.2f, Frequency: %.2f'%(fft_xaxis[start_location],fft_xaxis[end_location],num.average(current_fft),current_fft.max(),current_fft.max()/num.average(current_fft),current_frequencies[current_fft.argmax()])

            start_location=start_location+1
            end_location=end_location+1

        filtered_peaks=[]
        for b in peaks:
            if b[1]>overall_average*10 and b[0]>2000:
                filtered_peaks.append(b)
        return filtered_peaks

    def fluc_step(self,event):
        #self.stext.insert(Tkconstants.END, '----------------------------------\n')
        #self.stext.see(Tkconstants.END)
        self.data_class.extract_axial_array()
        for i in range(0,80,5):
            self.data_class.output_fluc_strucs(self.data_class.shot,'H1ToroidalAxial',svd_time=[i/1000.0,i/1000.0+self.data_class.NFFT*1.0/self.data_class.fs],plot="NO")

    def fluc_calc(self,event):
        #self.stext.insert(Tkconstants.END, '----------------------------------\n')
        #self.stext.see(Tkconstants.END)
        self.data_class.extract_axial_array()
        self.data_class.output_fluc_strucs(self.data_class.shot,'H1ToroidalAxial',svd_time=[self.data_class.timeslice,self.data_class.timeslice+self.data_class.NFFT*1.0/self.data_class.fs],plot="NO")

    def axial_plot_spawn(self,event):
        plot_signals(self.data_class.shot,'H1ToroidalAxial')

    def fft_spawn(self,event):
        pass

    def all_coils_plot_spawn(self,event):
        plot_signals(self.data_class.shot,'H1Toroidal')

    def poloidal1_plot_spawn(self,event):
        plot_signals(self.data_class.shot,'H1Poloidal1')

    def poloidal2_plot_spawn(self,event):
        plot_signals(self.data_class.shot,'H1Poloidal2')

    def svd_spawn(self,event):
        self.data_class.extract_axial_array()
        self.data_class.output_fluc_strucs(self.data_class.shot,'H1ToroidalAxial',svd_time=[self.data_class.timeslice,self.data_class.timeslice+self.data_class.NFFT*1.0/self.data_class.fs])

        #tempdata=pf.getDevice(self.data_class.device).acq.getdata(self.data_class.shot,"H1ToroidalAxial")
        
        #print 'performing SVD on ' + str(self.data_class.shot)
        #tempdata_reduced_time=tempdata.reduce_time([0.05, 0.0505]).subtract_mean().normalise(method='v',separate=True)
        #tempdata_reduced_time.svd().svdplot()
        #self.fig2.canvas.draw()
     #   [self.data_class.timeslice, self.data_class.timeslice+0.0005]

    def figure_shot_text(self):
        try:
            self.shot_text_update.remove()
        except:
            pass
        self.shot_text_update.set_text(str(self.data_class.shot))
        self.fig.canvas.draw()

    def increment_shot_plus_one(self,event):
        self.data_class.shot=self.data_class.shot+1
        self.figure_shot_text()
        if event.button==1:
            self.increment()

    def increment_shot_minus_one(self,event):
        self.data_class.shot=self.data_class.shot-1
        self.figure_shot_text()
        if event.button==1:
            self.increment()

    def increment_shot_plus_ten(self,event):
        self.data_class.shot=self.data_class.shot+10
        self.figure_shot_text()
        if event.button==1:
            self.increment()

    def increment_shot_minus_ten(self,event):
        self.data_class.shot=self.data_class.shot-10
        self.figure_shot_text()
        if event.button==1:
            self.increment()

    def increment_shot_zero(self,event):
        h1tree = MDS.Tree('h1data',0,'READONLY')
        h1shotnumber = h1tree.getCurrent('h1data')
        self.data_class.shot=h1shotnumber
        self.figure_shot_text()
        if event.button==1:
            self.increment()

    def increment_go(self,event):
        self.increment()
        pass

    def increment(self):
        self.data_class.extract_data()
        self.update_fft()
        #self.update_spec()
        self.redraw_figure(update_spec = 1, get_flucs = 1)

    def figure_clicked(self,event):
        print 'clicked', event.button
        if event.button==3:
            if (event.inaxes==self.spec_axis):# or (event.inaxes==self.ne_scaling_axis):
                self.data_class.timeslice=event.xdata/1000.
                #self.stext.insert(Tkconstants.END, 'Time = '+ str(self.data_class.timeslice)+'\n')
                #self.stext.see(Tkconstants.END)
                self.old_spec_ylim = self.spec_plot.get_axes().get_ylim()
                self.old_spec_xlim = self.spec_plot.get_axes().get_xlim()

                self.update_fft()
                self.spec_plot.set_clim(self.clim)
                self.redraw_figure()
                if self.data_class.axial_array_shot==self.data_class.shot:
                    self.data_class.output_fluc_strucs(self.data_class.shot,'H1ToroidalAxial',svd_time=[self.data_class.timeslice,self.data_class.timeslice+self.data_class.NFFT*1.0/self.data_class.fs],plot="NO")
                    #self.data_class.average_csd_func2(start_time=self.data_class.timeslice)
        elif event.button == 2:
            self.redraw_figure(default=1)
    def update_fft_amp(self,val):
        self.amplitude=val
        self.redraw_figure()

    def update_freq_axis(self,val):
        self.freq_axis=[self.frequency_axis_slider.val,self.frequency_axis_slider2.val]
        self.redraw_figure()

    def update_ne_scaling_axis_range(self,val):
        self.ne_scaling_axis_range=[self.plot_ne_scaling_slider.val,self.plot_ne_scaling_slider2.val]
        self.redraw_figure()

    def update_noverlap(self,label):
        self.noverlap=1.0*float(label)/16.0
        #self.update_spec()
        self.redraw_figure(update_spec = 1)

    def update_nfft(self,label):
        self.data_class.NFFT=int(label)
        #self.update_spec()
        self.update_fft()
        self.redraw_figure(update_spec = 1)

    def update_clim(self,val):
        self.clim=[self.clim_slider2.val,self.clim_slider1.val+self.clim_slider2.val]
        self.spec_plot.set_clim(self.clim)
        self.redraw_figure()

    def update_timeslice(self,val):
        self.data_class.timeslice=val
        self.update_fft()
        self.redraw_figure()
    
    #Redraw the whole widget with the updated settings when something has changed
    def redraw_figure(self, default = 0, update_spec = 0, get_flucs = 0):
        if default == 0:
            self.old_spec_ylim = self.spec_plot.get_axes().get_ylim()
            self.old_spec_xlim = self.spec_plot.get_axes().get_xlim()
            print 'update spec axis', self.old_spec_ylim, self.old_spec_xlim
            tmp_spec_ylim = self.spec_plot.get_axes().get_ylim()
            tmp_spec_xlim = self.spec_plot.get_axes().get_xlim()
            print 'spec axis values', tmp_spec_ylim, tmp_spec_xlim
        try:
            self.spec_axisvertline.remove()
            self.ne_axisvertline.remove()
            self.signal_axisverline.remove()
        except:
            pass

        #clear axes
        self.ne_axis.cla()
        self.antenna_axis.cla()
        self.signal_axis.cla()
        self.fft_axis.cla()

        #replot lines
        self.fft_axis.plot(self.fft_yaxis,self.fft_xaxis/1000., 'r-')

        if self.data_class.inc_polar_data:
            colors = ['b','k','0.5']
            for i, tmp in enumerate(self.fft_yaxis_naked):
                self.fft_axis.plot(tmp,self.fft_xaxis/1000.,color = colors[i])

        self.ne_axis.plot(self.data_class.ne_time*1000.,self.data_class.ne_value)
        self.antenna_axis.plot(self.data_class.rf_top_time*1000.,self.data_class.rf_top_value,'r-')
        self.antenna_axis.plot(self.data_class.rf_bot_time*1000.,self.data_class.rf_bot_value,'k-')
        self.signal_axis.plot(self.data_class.timebase*1000.,self.data_class.original_data)

        if update_spec:
            self.update_spec()

        if self.plot_flucs == 1:
            if get_flucs:
                print 'getting fluc strucs'
                fs_dictionary, serial_number, success = HMA_funcs.single_shot_fluc_strucs(self.data_class.shot, 'H1ToroidalUpward', [0,0.1], 2048, power_cutoff = 0.0, n_svs = 1, overlap = 2)#fft_length)
                self.flucs_time_list2 = []; self.flucs_freq_list2 = []
                self.flucs_time_list1 = []; self.flucs_freq_list1 = []

                for serial in fs_dictionary.keys():

                    if fs_dictionary[serial]['SVs']>1:

                        self.flucs_time_list2.append(fs_dictionary[serial]['time'])
                        self.flucs_freq_list2.append(fs_dictionary[serial]['freq']/1000.)
                    else:
                        self.flucs_time_list1.append(fs_dictionary[serial]['time'])
                        self.flucs_freq_list1.append(fs_dictionary[serial]['freq']/1000.)
                        
            self.spec_axis.plot(self.flucs_time_list2, self.flucs_freq_list2, 'kx')
            self.spec_axis.plot(self.flucs_time_list1, self.flucs_freq_list1, 'k,')

        #include an ne predicted frequency
        if self.plot_ne_scaling==1:
            try:
                self.ne_scale_plot.pop(0).remove()
                print 'removed ne plot'
            except:
                print 'failed to remove ne plot'
                pass
            #self.ne_scaling_axis_range=[self.plot_ne_scaling_slider.val,self.plot_ne_scaling_slider2.val]
            constants = (0.5/(2*num.pi))/((4*num.pi*1.673*(10**(-16))**0.5))
            #print 'scale factor, ', constants
            #self.ne_scaling_axis.plot(self.data_class.ne_time,self.plot_ne_scaling_slider.val*1000*((self.plot_ne_scaling_slider2.val*self.data_class.ne_value)**(-0.5)),'y')
            #self.plot_ne_scaling_slider_time.on_changed(self.update_ne_scaling_axis_range)
            iota_config = 1.3; n_mode = 4; m_mode = 3; major_radius = 1.22
            k_parallel = (num.abs(n_mode - m_mode * iota_config))/major_radius
            B0 = 0.5
            rho = self.plot_ne_scaling_slider2.val*1.672*(10**(-27))*self.data_class.ne_value*(10**18)
            perm = 4*num.pi*(10**(-7))
            v_a = B0 / ((rho * perm)**(0.5))
            scaling_factor = self.plot_ne_scaling_slider.val
            predicted_frequency = scaling_factor * k_parallel * v_a /(2.0 * num.pi)
            print 'k_parallel,  rho, min(v_a) ',k_parallel,rho, num.min(v_a)
            print 'predicted_frequency : ', predicted_frequency
            #self.ne_scaling_axis.plot(self.data_class.ne_time*1000+(self.plot_ne_scaling_slider_time.val),predicted_frequency,'k')
            self.ne_scale_plot = self.spec_axis.plot(self.data_class.ne_time*1000+(self.plot_ne_scaling_slider_time.val),predicted_frequency/1000.,'k')
            #self.ne_scaling_axis.set_ylim(self.freq_axis)
        else:
            try:
                self.ne_scale_plot.pop(0).remove()
                print 'removed ne plot'
            except:
                print 'failed to remove ne plot'

        #redo vert lines
        self.ne_axisvertline=self.ne_axis.axvline(x=self.data_class.timeslice *1000,color='y')
        self.signal_axisvertline=self.signal_axis.axvline(x=self.data_class.timeslice*1000,color='y')
        self.spec_axisvertline=self.spec_axis.axvline(x=self.data_class.timeslice*1000,color='y')
        #redo the scales depending on default or not

        self.fft_axis.grid(b=True)

        if default == 1:
            self.spec_axis.set_ylim([0,1000])
            self.spec_axis.set_xlim([-15,120])
        else:
            self.spec_axis.set_xlim(tmp_spec_xlim)
            self.spec_axis.set_ylim(self.old_spec_ylim)

        self.spec_plot.set_clim(self.clim)
        self.ne_axis.set_ylim([0,2.5])
        self.antenna_axis.set_ylim([0, 100])
        self.signal_axis.set_ylim([-5.5,5.5])
        self.fft_axis.set_xlim([0,self.amplitude])

        #misc labelling
        self.spec_axis.set_xlabel('Time (ms)')
        self.spec_axis.set_ylabel('Freq (kHz)')
        self.ne_axis.set_xlabel('Time (ms)')
        self.ne_axis.set_ylabel('ne density (10^18/m^3)')
        self.antenna_axis.set_ylabel('Antenna Current (A)')
        self.antenna_axis.yaxis.set_ticks_position('right')
        self.fft_axis.set_xlabel('Amplitude (a.u)')
        self.fft_axis.set_ylabel('Freq (Hz)')
        self.spec_axis.set_title(str(self.data_class.shot) + ' kh=' + str(self.data_class.kh) + ',(T2 - T1):' +str(self.data_class.antenna_phase_diff) + 'deg', fontsize='medium')
        self.signal_axis.set_xlabel('Time (s)')
        self.signal_axis.set_ylabel('Volts')
        self.signal_axis.grid(b=True)
        self.ne_axis.grid(b=True)
        self.fft_axis.set_title(self.data_class.coil_name[-17:].replace('_','-'))
        #pt.setp(self.ne_scaling_axis,yticklabels=[])
        self.fig.canvas.draw()
        #self.fig.savefig(str(self.data_class.shot)+'.png') - edit this so that lots of pictures can be made quickly

    def update_fft(self):
        start_location=num.searchsorted(self.data_class.timebase,self.data_class.timeslice)
        self.fft_data=num.fft.fft(num.array(self.data_class.original_data[start_location:start_location+self.data_class.NFFT])*num.hanning(self.data_class.NFFT))
        k=num.arange(-self.data_class.NFFT/2,self.data_class.NFFT/2)
        self.fft_xaxis=k*self.data_class.fs/self.data_class.NFFT
        self.fft_yaxis=abs(num.concatenate((self.fft_data[self.data_class.NFFT/2:],self.fft_data[:self.data_class.NFFT/2])))
        if self.data_class.inc_polar_data:
            self.fft_yaxis_naked = []
            for tmp_data in self.data_class.polar_signals:
                tmp_fft = num.fft.fft(num.array(tmp_data[start_location:start_location+self.data_class.NFFT])*num.hanning(self.data_class.NFFT))
                self.fft_yaxis_naked.append(abs(num.concatenate((tmp_fft[self.data_class.NFFT/2:],tmp_fft[:self.data_class.NFFT/2]))))

    def update_spec(self):
        try:
            self.spec_axis.cla()
        except:
            pass
        start_loc = 0
        end_loc = -1
        #start_loc = searchsorted(self.data_class.timebase,-0.005)
        #end_loc  = searchsorted(self.data_class.timebase,0.1) -1
        #a,b,c,self.spec_plot=self.spec_axis.specgram(self.data_class.original_data[searchsorted(self.data_class.timebase,-0.005):searchsorted(self.data_class.timebase,0.1)],NFFT=self.data_class.NFFT,Fs=self.data_class.fs,noverlap=int(self.noverlap*self.data_class.NFFT),window=num.hanning(self.data_class.NFFT), cmap='jet', xextent = [self.data_class.timebase[searchsorted(self.data_class.timebase,-0.005)], self.data_class.timebase[searchsorted(self.data_class.timebase,0.1)]])
        print self.data_class.timebase.shape, start_loc, end_loc, self.data_class.timebase[start_loc], self.data_class.timebase[end_loc]
        a,b,c,self.spec_plot=self.spec_axis.specgram(self.data_class.original_data[start_loc: end_loc], NFFT=self.data_class.NFFT, Fs=self.data_class.fs/1000., noverlap=int(self.noverlap*self.data_class.NFFT), window=num.hanning(self.data_class.NFFT), cmap='jet', xextent = [self.data_class.timebase[start_loc]*1000, self.data_class.timebase[end_loc]*1000])
