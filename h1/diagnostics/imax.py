import scipy.signal
import MDSplus as MDS
import numpy as np
import matplotlib.pyplot as pt

class ImaxData():
    def __init__(self,shot_list, plot_data = 0, get_mirnov_triggers = 1, plot_mirnov_triggers = 0):
        self.image_array = np.zeros((16,512,512),dtype=float)
        self.shot_list = shot_list
        if plot_data: fig, ax = pt.subplots(nrows= 4, ncols = 4,sharex = 1, sharey = 1); ax = ax.flatten()
        for i, shot in enumerate(self.shot_list):
            self.image_array[i,:,:] = MDS.Tree('imax',shot).getNode('.images').data()[0,:,:]
            if plot_data: 
                im = ax[i].imshow(self.image_array[i,:,:], interpolation = 'none', aspect = 'auto', origin='lower')
                im.set_clim(247,65535)
        if get_mirnov_triggers:
            self.get_mirnov_triggers(plot_mirnov_triggers)
        else:
            self.phase_average = np.linspace(0, 2.*np.pi, endpoint = False)
            self.phase_std = self.phase_average * 0
            self.amp_average = self.phase_average * 0 + 1
            self.amp_std = self.phase_average * 0
            self.electron_dens = self.phase_average * 0 + 1
        if plot_data:
            ax[-1].set_xlim([0,512])
            ax[-1].set_ylim([0,512])
            fig.subplots_adjust(hspace=0.0, wspace=0.0,left=0., bottom=0.,top=0.95, right=0.95)
            fig.suptitle('Raw Images')
            fig.canvas.draw(); fig.show()

        #self.get_john_analysis(plot_john_data = 0)
        #self.calibrate_data(dark_shot=1103, white_shot=1107, plot_calibration_im = 0, clip = 0.2, plot_cal_images = 1, cal_sum = 1, med_filt = 3, mode_amp_filt = self.amp_average, electron_density = self.electron_dens)
        #self.fourier_decomp(plot_amps = 1, plot_phases = 1)

    def get_mirnov_triggers(self, plot_mirnov_triggers):
        if plot_mirnov_triggers: fig_trig, ax_trig = pt.subplots(nrows= 4, ncols = 4,sharex = 1, sharey = 1); ax_trig = ax_trig.flatten()
        if plot_mirnov_triggers: fig_phase, ax_phase = pt.subplots(nrows = 2); 
        if plot_mirnov_triggers: fig_hilb, ax_hilb = pt.subplots(nrows= 4, ncols = 4,sharex = 1, sharey = 1); ax_hilb = ax_hilb.flatten()
        phase_average_list = []; phase_std_list = []
        amp_average_list = []; amp_std_list = []
        density_list = []
        for i, shot in enumerate(self.shot_list):
            PLL_node = MDS.Tree('h1data',shot).getNode('.operations.mirnov.a14_14:input_6')
            mirnov_node = MDS.Tree('h1data',shot).getNode('.operations.mirnov:a14_15:input_1')
            electron_density = MDS.Tree('h1data',shot).getNode('.electr_dens.NE_HET:NE_CENTRE')

            #provide 4ms padding infront of and after the PLL signal goes high
            a = PLL_node.data() > 1.5
            padding = int(0.004/(PLL_node.dim_of().data()[1] - PLL_node.dim_of().data()[0]))
            start_ind = np.argmax(a) - padding
            end_ind = len(a) - np.argmax(a[::-1]) + padding

            mirnov_data, mirnov_time, PLL_data, PLL_time, electr_data, hilb_sig, maxima = self.hilb_transform_check(PLL_node, mirnov_node, electron_density, start_ind, end_ind)
            phase_complex = np.exp(1j*np.angle(hilb_sig[maxima]))
            phase_average_list.append(np.angle(np.mean(phase_complex)))
            phase_std_list.append(np.sqrt(np.log(1./np.abs(np.mean(phase_complex)))))
            amp_std_list.append(np.std(np.abs(hilb_sig[maxima])))
            amp_average_list.append(np.mean(np.abs(hilb_sig[maxima])))
            density_list.append(np.mean(electr_data[maxima]))
            if plot_mirnov_triggers:
                ax_trig[i].plot(mirnov_node.dim_of().data(),mirnov_node.data())
                ax_trig[i].plot(PLL_node.dim_of().data(),PLL_node.data())
                ax_hilb[i].plot(mirnov_time[maxima], np.angle(hilb_sig)[maxima],'o')

        self.phase_average = np.array(phase_average_list)
        self.phase_std = np.array(phase_std_list)
        self.amp_average = np.array(amp_average_list)
        self.amp_std = np.array(amp_std_list)
        self.electron_dens = np.array(density_list)
        if plot_mirnov_triggers:
            fig_trig.canvas.draw(); fig_trig.show()
            fig_hilb.canvas.draw(); fig_hilb.show()
            #ax_phase[0].plot(range(16), phase_average_list, '-o')
            ax_phase[0].errorbar(range(16), phase_average_list, yerr= phase_std_list,fmt='-o')
            ax_phase[1].errorbar(range(16), amp_average_list, yerr= amp_std_list,fmt='-o')
            ax_phase[1].set_xlabel('shot number relative to first shot')
            ax_phase[1].set_ylabel('amplitude')
            ax_phase[0].set_ylabel('phase')
            ax_phase[1].set_ylim([0,np.max(np.array(amp_average_list)+np.array(amp_std_list)) + 0.1])
            fig_phase.suptitle('{} - {} hilbert amps and phases at trigger times'.format(np.min(self.shot_list), np.max(self.shot_list)))
            fig_phase.canvas.draw(); fig_phase.show()

    def hilb_transform_check(self, PLL_node, mirnov_node, electron_node, start_ind, end_ind):#start_time = 0.005, end_time = 0.04):
        '''Hilbert transform of the Mirnov signal

        SH : 9July2013
        '''
        PLL_data = PLL_node.data()
        PLL_time = PLL_node.dim_of().data()
        electr_data = electron_node.data()
        electr_time = electron_node.dim_of().data()
        #start_loc = np.argmin(np.abs(PLL_time - start_time))
        #end_loc = np.argmin(np.abs(PLL_time - end_time))
        print start_ind, end_ind, end_ind - start_ind
        PLL_data = PLL_data[start_ind:end_ind]
        PLL_time = PLL_time[start_ind:end_ind]
        mirnov_data = mirnov_node.data()
        mirnov_time = mirnov_node.dim_of().data()
        mirnov_data = np.interp(PLL_time, mirnov_time, mirnov_data)
        electr_data = np.interp(PLL_time, electr_time, electr_data)
        mirnov_data = mirnov_data - np.mean(mirnov_data)
        mirnov_time = PLL_time
        hilb_sig = scipy.signal.hilbert(mirnov_data)
        maxima = np.zeros(PLL_data.shape,dtype=bool)
        maxima[1:-1] = (PLL_data>4.5)[1:-1] * (PLL_data[1:-1]>PLL_data[0:-2])* (PLL_data[1:-1]>PLL_data[2:])
        return mirnov_data, mirnov_time, PLL_data, PLL_time, electr_data, hilb_sig, maxima

    def fourier_decomp(self,plot_amps = 0, plot_phases = 0):
        '''Perform an FFT on the image and plot the results if asked to
        '''
        self.fft_values = np.fft.fft(self.image_array_cal, axis = 0)
        if plot_amps: fig, ax = pt.subplots(nrows= 3, ncols = 3,sharex = 1, sharey = 1); ax = ax.flatten()
        if plot_phases: fig2, ax2 = pt.subplots(nrows= 3, ncols = 3,sharex = 1, sharey = 1); ax2 = ax2.flatten()
        im_list = []
        for i in range(8):
            if plot_amps:
                im_list.append(ax[i].imshow(np.abs(self.fft_values[i,:,:])*2, interpolation = 'nearest', aspect = 'auto',origin='lower'))
            if plot_phases:
                im2 = ax2[i].imshow(np.angle(self.fft_values[i,:,:]), interpolation = 'nearest', aspect = 'auto', origin='lower')
                im2.set_clim([-np.pi,np.pi])
        if plot_amps:
            for i in range(2,len(im_list)):im_list[i].set_clim(im_list[1].get_clim())
            ax[-1].set_xlim([0,512])
            ax[-1].set_ylim([0,512])

            fig.suptitle('Fourier amplitude (color range 0, 0.5)')
            fig.subplots_adjust(hspace=0.0, wspace=0.0,left=0., bottom=0.,top=0.95, right=0.95)
            fig.canvas.draw(); fig.show()
            #fig.savefig('imax_Fourier_amplitude.png')
        if plot_phases:
            fig2.subplots_adjust(hspace=0.0, wspace=0.0,left=0., bottom=0.,top=0.95, right=0.95)
            fig2.suptitle('Fourier phase (color range -pi,pi)')
            #fig2.savefig('imax_Fourier_phase.png')
            fig2.canvas.draw(); fig2.show()

    def calibrate_data(self, dark_shot=1103, white_shot=1107, plot_calibration_im = 0, clip = 0.2, plot_cal_images=0, cal_sum = 1, med_filt = 3, mode_amp_filt = None):
        '''Calibrate the imax image using a dark and white shot

        need to provide a dark_shot number, white_shot number calibrate by
        sum of all pixels (cal_sum) apply a median filter with side length
        med_filt amply a scaling based on the mode amplitudes from the
        Mirnov signal (mode_amp_filt is an array)

        SH: 3Jul2013
        '''
        self.image_array_cal = np.zeros(self.image_array.shape,dtype=float)
        self.dark_image = MDS.Tree('imax',dark_shot).getNode('.images').data()[0,:,:]
        self.white_image = MDS.Tree('imax',white_shot).getNode('.images').data()[0,:,:]
        if plot_calibration_im:
            fig, ax = pt.subplots(nrows= 2,sharex = 1, sharey = 1); ax = ax.flatten()
            im = ax[0].imshow(self.white_image, interpolation = 'none', aspect = 'auto',origin='lower')
            pt.colorbar(im,ax=ax[0])
            im = ax[1].imshow(self.dark_image, interpolation = 'none', aspect = 'auto', origin='lower')
            pt.colorbar(im,ax=ax[1])
            fig.canvas.draw(); fig.show()
        full_scale = (self.white_image - self.dark_image)
        max_value = np.max(full_scale)
        if plot_cal_images: 
            fig, ax = pt.subplots(nrows= 4, ncols = 4,sharex = 1, sharey = 1); ax = ax.flatten()
            im_list = []
        for i in range(16):
            self.image_array_cal[i,:,:] = (self.image_array[i,:,:] - self.dark_image)/np.clip(self.white_image - self.dark_image, max_value*clip,65000)
            print i,
            if cal_sum:
                print 'cal_sum,', np.sum(self.image_array_cal[i,:,:]), self.electron_dens[i], self.amp_average[i]
                self.image_array_cal[i,:,:] = self.image_array_cal[i,:,:]/np.sum(self.image_array_cal[i,:,:])* 65000.
            if med_filt!=0:
                print 'med_filt,',
                self.image_array_cal[i,:,:] = scipy.signal.medfilt(self.image_array_cal[i,:,:], med_filt)
            if mode_amp_filt:
                print 'mode_amp_filt',
                self.image_array_cal[i,:,:] = self.image_array_cal[i,:,:] / self.amp_average[i]
            if plot_cal_images:
                im_list.append(ax[i].imshow(self.image_array_cal[i,:,:], interpolation = 'none', aspect = 'auto', origin='lower'))
            print ''
        if plot_cal_images: 
            for i in range(0,len(im_list)):im_list[i].set_clim(im_list[0].get_clim())
            ax[-1].set_xlim([0,512])
            ax[-1].set_ylim([0,512])
            fig.subplots_adjust(hspace=0.0, wspace=0.0,left=0., bottom=0.,top=0.95, right=0.95)
            fig.suptitle('imax_corrected Images')
            fig.canvas.draw(); fig.show()
            fig.savefig('imax_corrected_images.png')

    def get_john_analysis(self, plot_john_data = 0):
        '''Extracts Johns results from the MDSplus tree
        SH : 9July2013
        '''
        imax_tree = MDS.Tree('imax',self.shot_list[0])
        self.john_amp_data = imax_tree.getNode('.processed.amplitude').data()
        self.john_phase_data = imax_tree.getNode('.processed.phase').data()
        if plot_john_data:
            fig, ax = pt.subplots(nrows = 3,ncols = 2, sharex = 1, sharey = 1)
            im_list = []
            for i in range(john_amp_data.shape[0]):
                im_list.append(ax[i,0].imshow(self.john_amp_data[i,:,:],interpolation = 'nearest', aspect = 'auto',origin='lower'))
                im_ang = ax[i,1].imshow(self.john_phase_data[i,:,:],interpolation = 'nearest', aspect = 'auto',origin='lower')
                im_ang.set_clim([-np.pi,np.pi])
            im_list[-1].set_clim(im_list[1].get_clim())
            fig.canvas.draw(); fig.show()

    def plot_fft_values(self, list_of_slices):
        fig,ax = pt.subplots(nrows=2)
        for i in range(16):
            for j in list_of_slices:
                ax[0].plot(np.abs(self.fft_values[i,:,j]))
                ax[1].plot(np.angle(self.fft_values[i,:,j]))
        fig.canvas.draw(); fig.show()

def make_animation(start_shot_list, titles, harmonic = 1, base_directory = '', prefix = 'hello'):
    '''
    Creates an animation of the L and R side of the 4/3 and 5/4 whale tail
    Need to provide the starting shot list and the titles

    SRH: 12July2013
    '''
    phases = np.linspace(0,2.*np.pi,20,endpoint=True)
    fft_list = []
    for start_shot, title in zip(start_shot_list, titles):
        tmp1 = ImaxData(range(start_shot,start_shot + 16), plot_data = 0, plot_mirnov_triggers = 0,get_mirnov_triggers = 0)
        tmp1.calibrate_data(dark_shot=1103, white_shot=1107, plot_calibration_im = 0, clip = 0.2, plot_cal_images = 0, cal_sum = 1, med_filt = 3, mode_amp_filt = None)
        tmp1.fourier_decomp(plot_amps = 0, plot_phases = 0)
        fft_list.append(tmp1.fft_values)
    fig, ax = pt.subplots(nrows = 2, ncols = 2,sharex = 1, sharey = 1); ax = ax.flatten()
    clim_list = []
    for i, phase in enumerate(phases):
        for j, fft_values in enumerate(fft_list):
            tmp = np.real(fft_values[harmonic,:,:] * np.exp(1j*phase))
            ax[j].cla()
            im = ax[j].imshow(tmp, interpolation = 'none', aspect = 'auto', origin='lower',cmap='hot')
            if i==0:
                clim=im.get_clim()
                clim_val = np.max(np.abs(clim))
                clim_list.append([-clim_val, clim_val])
                print clim
                #ax[j].title('{}'.format(titles[j]))
            im.set_clim(clim_list[-1])
            fig.canvas.draw(); fig.show()
            fig.savefig('{}{}_{:02d}.png'.format(base_directory, prefix, i))
    
