import scipy.signal
import MDSplus as MDS
import numpy as np
import matplotlib.pyplot as pt
import h1.diagnostics.winspec as winspec
import h1.diagnostics.SPE_reader as SPE_reader
import copy
class ImaxData():
    def __init__(self,shot_list=None, SPE_file = None, single_mdsplus=0, plot_data = 0, get_mirnov_triggers = 0, plot_mirnov_triggers = 0, mdsplus_tree = 'imax', mdsplus_node='.pimax.images', flipud=0, fliplr=0, rot90=0):
        self.shot_list = shot_list
        if single_mdsplus!=0:
            #self.image_array = np.zeros((len(shot_list),512,512),dtype=float)
            self.image_array = MDS.Tree(mdsplus_tree,single_mdsplus).getNode(mdsplus_node).data()
            if flipud: self.image_array = np.flipud(self.image_array)
            if fliplr: self.image_array = np.fliplr(self.image_array)
            if rot90: self.image_array = np.rot90(self.image_array)
            self.shot_list = [single_mdsplus]
        elif shot_list!=None:
            self.image_array = np.zeros((len(shot_list),512,512),dtype=float)
            for i, shot in enumerate(self.shot_list):
                tmp_data = MDS.Tree(mdsplus_tree,shot).getNode(mdsplus_node).data()[0,:,:]
                if flipud: tmp_data = np.flipud(tmp_data)
                if fliplr: tmp_data = np.fliplr(tmp_data)
                if rot90: tmp_data = np.rot90(tmp_data)

                if i==0: self.image_array = np.zeros((len(shot_list),tmp_data.shape[0],tmp_data.shape[1]),dtype=float)
                self.image_array[i,:,:] = tmp_data.copy()
        elif SPE_file!=None:
            if SPE_file.__class__==list:
                n_files = len(SPE_file)
                self.image_array = np.zeros((n_files,512,512),dtype=float)
                for i in range(n_files):
                    tmp_data = winspec.SpeFile(SPE_file[i]).data[0,:,:]
                    if flipud: tmp_data = np.flipud(tmp_data)
                    if fliplr: tmp_data = np.fliplr(tmp_data)
                    if rot90: tmp_data = np.rot90(tmp_data)
                    self.image_array[i,:,:] = tmp_data.copy()
            else:
                self.image_array = np.zeros((1,512,512),dtype=float)
                cal_file = winspec.SpeFile(SPE_file)
                #self.image_array[0,:,:] = np.rot90(cal_file.data[0,:,:])
                self.image_array = cal_file.data
                for i in range(self.image_array.shape[0]):
                    if flipud: self.image_array = np.flipud(self.image_array)
                    if fliplr: self.image_array = np.fliplr(self.image_array)
                    if rot90: self.image_array = np.rot90(self.image_array)
                    self.image_array[i,:,:] = self.image_array[i, :, :].copy()
        else:
            raise ValueError("either shot_list, or SPE_file must not be None")
        if get_mirnov_triggers:
            self.get_mirnov_triggers(plot_mirnov_triggers)
        else:
            self.phase_average = np.linspace(0, 2.*np.pi, endpoint = False)
            self.phase_std = self.phase_average * 0
            self.amp_average = self.phase_average * 0 + 1
            self.amp_std = self.phase_average * 0
            self.electron_dens = self.phase_average * 0 + 1
        if plot_data: self.plot_images()

    def get_phase_std(self,n_phases = 1, plot = False, ax = None):
        phase_ave_list = []
        phase_std_list = []
        for loc, i in enumerate(self.shot_list):
            print i
            if loc==1: 
                tmp_plot = True
            else:
                tmp_plot = False
            tmp1, tmp2 = PLL_performance2(i, n_phases = n_phases, plot=tmp_plot)
            phase_ave_list.extend(tmp1)
            phase_std_list.extend(tmp2)
        if plot:
            if ax == None:
                fig, ax = pt.subplots()
            phases = np.array(phase_ave_list)
            stds = np.array(phase_std_list)
            ax.errorbar(np.arange(phases.shape[0]) + 1, phases, yerr = stds)
            if ax ==None:
                fig.canvas.draw(); fig.show()


    def plot_fourier_harmonics_pub(self,pub_fig = 0, clim=None, savefig= None, fig = None, ax = None, draw_show = 1, fliplr = 1, save_fig=None,n_harmonics =3):
        fig,ax = pt.subplots(nrows=n_harmonics,ncols=2)
        if n_harmonics == 1: ax = ax[np.newaxis,:]
        if pub_fig:
            cm_to_inch=0.393701
            import matplotlib as mpl
            old_rc_Params = mpl.rcParams
            mpl.rcParams['font.size']=8.0
            mpl.rcParams['axes.titlesize']=8.0#'medium'
            mpl.rcParams['xtick.labelsize']=8.0
            mpl.rcParams['ytick.labelsize']=8.0
            mpl.rcParams['lines.markersize']=5.0
            mpl.rcParams['savefig.dpi']=300
            fig.set_figwidth(8.48*cm_to_inch)
            fig.set_figheight((8.48*(0.5*n_harmonics)+0.1)*cm_to_inch)
        amp_im_list = []
        phase_im_list = []
        for i in range(n_harmonics):
            if fliplr:
                tmp = np.fliplr(self.fft_values[i+1,:,:])
            else:
                tmp = self.fft_values[i,:,:]
            amp_im_list.append(ax[i,0].imshow(np.abs(tmp)*2, interpolation = 'nearest', aspect = 'auto',origin='lower',cmap='spectral'))  
            ax[i,0].set_ylabel('harmonic {}'.format(i+1))
            phase_im_list.append(ax[i,1].imshow(np.angle(tmp), interpolation = 'nearest', aspect = 'auto',origin='lower',cmap='RdBu'))
            phase_im_list[-1].set_clim([-np.pi,np.pi])
            if i==0:
                clim = np.array(amp_im_list[-1].get_clim())
            else:
                tmp_clim = np.array(amp_im_list[-1].get_clim())
                if tmp_clim[0]<clim[0]: clim[0]=tmp_clim[0]
                if tmp_clim[1]>clim[1]: clim[1]=tmp_clim[1]

        for i in amp_im_list:i.set_clim(clim)
        fig.subplots_adjust(hspace=0.05, wspace=0.05,left=0.05, bottom=0.05,top=0.95, right=0.95)
        tmp1 = [ax[i,0] for i in range(n_harmonics)]
        tmp2 = [ax[i,1] for i in range(n_harmonics)]
        amp_cbar = pt.colorbar(amp_im_list[-1],ax = tmp1, orientation='horizontal',pad=0.01)
        phase_cbar = pt.colorbar(phase_im_list[-1],ax = tmp2, orientation='horizontal',pad=0.01)
        amp_cbar.set_label('Amplitude (a.u)')
        phase_cbar.set_label('Phase (rad)')
        clim_cur = amp_im_list[-1].get_clim()
        amp_cbar.set_ticks(np.round(np.linspace(clim_cur[0],clim_cur[1],5),1)[0:-1])
        clim_cur = [-np.pi,np.pi]
        phase_cbar.set_ticks(np.round(np.linspace(clim_cur[0],clim_cur[1],5),1)[0:])
        phase_cbar.ax.tick_params(labelsize=7) 
        amp_cbar.ax.tick_params(labelsize=7) 
        for j in ax.flatten():
            xticks = j.xaxis.get_major_ticks()
            yticks = j.yaxis.get_major_ticks()
            for jj in xticks: jj.label1.set_visible(False)
            for jj in yticks: jj.label1.set_visible(False)
        # if plot_amps:
        #     for i in range(2,len(im_list)):im_list[i].set_clim(im_list[1].get_clim())
        #     ax[-1].set_xlim([0,self.fft_values.shape[1]])
        #     ax[-1].set_ylim([0,self.fft_values.shape[2]])
        #     fig.subplots_adjust(hspace=0.0, wspace=0.0,left=0., bottom=0.,top=0.95, right=0.95)
        #     if draw:
        #         fig.suptitle('Fourier amplitude (color range 0, 0.5)')
        #         fig.canvas.draw(); fig.show()
        #     #fig.savefig('imax_Fourier_amplitude.png')
        # if plot_phases:
        #     fig2.subplots_adjust(hspace=0.0, wspace=0.0,left=0., bottom=0.,top=0.95, right=0.95)
        #     if draw:
        #         fig2.suptitle('Fourier phase (color range -pi,pi)')
        #         fig2.canvas.draw(); fig2.show()

        if save_fig!=None:
            fig.savefig(save_fig, bbox_inches='tight',pad_inches=0.05)
        fig.canvas.draw(); fig.show()


    def plot_images_subset(self,cal=0, pub_fig = 0, image_indices = None, clim=None, savefig= None, fig = None, ax = None, draw_show = 1, fliplr = 1, save_fig=None, inc_DC_diff = False):
        fig = pt.figure()
        if image_indices==None:
            image_indices = [0,2,4,6,8,10,12,14]
        if pub_fig:
            cm_to_inch=0.393701
            import matplotlib as mpl
            old_rc_Params = mpl.rcParams
            mpl.rcParams['font.size']=8.0
            mpl.rcParams['axes.titlesize']=8.0#'medium'
            mpl.rcParams['xtick.labelsize']=8.0
            mpl.rcParams['ytick.labelsize']=8.0
            mpl.rcParams['lines.markersize']=5.0
            mpl.rcParams['savefig.dpi']=300
            fig.set_figwidth(8.48*cm_to_inch)
            fig.set_figheight(8.48*1.3*cm_to_inch)

        if cal==0: 
            print 'Uncalibrated images'
            plot_array = self.image_array
        else:
            print 'Calibrated images'
            plot_array = self.image_array_cal
        ax = []; im_list = []
        rows = len(image_indices)/2+1
        if inc_DC_diff: rows+=1
        row = 0; col = 0
        for i,index in enumerate(image_indices):
            print row, col
            if i==0:
                ax.append(pt.subplot2grid((rows,2), (row,col)))
            else:
                ax.append(pt.subplot2grid((rows,2), (row,col), sharex=ax[0],sharey=ax[0]))
            if col==0:
                col+=1
            else:
                col = 0
                row += 1
        #ax.append(pt.subplot2grid((rows,2), (0,1), sharex=ax[0],sharey=ax[0]))
        #ax.append(pt.subplot2grid((rows,2), (1,0), sharex=ax[0],sharey=ax[0]))
        #ax.append(pt.subplot2grid((rows,2), (1,1), sharex=ax[0],sharey=ax[0]))
        if inc_DC_diff:
            ax.append(pt.subplot2grid((rows,2), (rows-2,0), sharex=ax[0],sharey=ax[0]))
            ax.append(pt.subplot2grid((rows,2), (rows-2,1), sharex=ax[0],sharey=ax[0]))
        ax2 = pt.subplot2grid((rows,2), (rows-1,0), colspan=2)
        #ax2 = pt.subplot2grid((rows,2), (rows-1,0), colspan=1)
        
        for i, index in enumerate(image_indices):
            if fliplr:
                tmp = np.fliplr(plot_array[index,:,:])
            else:
                tmp = plot_array[index,:,:]
            im_list.append(ax[i].imshow(tmp, interpolation = 'none', aspect = 'auto', origin='lower', cmap = 'spectral'))
            if i==0:
                clim = np.array(im_list[-1].get_clim())
            else:
                tmp_clim = np.array(im_list[-1].get_clim())
                if tmp_clim[0]<clim[0]: clim[0]=tmp_clim[0]
                if tmp_clim[1]>clim[1]: clim[1]=tmp_clim[1]
            ax[i].set_ylabel('frame {}'.format(index+1))
        for i in im_list: i.set_clim(clim)
        if inc_DC_diff:
            if fliplr:
                tmp = np.fliplr(np.abs(self.fft_values[0,:,:]))
                tmp2 = np.fliplr(plot_array[image_indices[1],:,:] - plot_array[image_indices[0],:,:])
            else:
                tmp = np.abs(self.fft_values[0,:,:])
                tmp2 = plot_array[image_indices[1],:,:] - plot_array[image_indices[0],:,:]
            im_DC = ax[-2].imshow(tmp, interpolation = 'none', aspect = 'auto', origin='lower', cmap = 'spectral')
            im_diff = ax[-1].imshow(tmp2, interpolation = 'none', aspect = 'auto', origin='lower', cmap = 'RdBu')
            val = np.min(np.abs(im_diff.get_clim()))*0.8
            #print tmp_min, tmp_max
            im_diff.set_clim([-val,val])
            ax[-2].set_ylabel('DC')
            ax[-1].set_ylabel('frame {} - frame {}'.format(image_indices[1]+1, image_indices[0]+1))
        ax[0].text(30,140,'TFC',color='white')
        ax[0].text(177,140,'TFC',color='white')
        #for i in im_list: i.set_clim([3.7,4.4])
        plot_locs = [[120,150], [140,60]]
        periods = 3
        for i in plot_locs:
            single_line = []
            for j in range(plot_array.shape[0]):
                if fliplr:
                    single_line.append(np.fliplr(plot_array[j,:,:])[i[1],i[0]])
                else:
                    single_line.append(plot_array[j,i[1],i[0]])
            orig_single_line = copy.deepcopy(single_line)
            for jjj in range(periods-1):
                single_line.extend(orig_single_line)
            ax2.plot(range(16*periods),single_line,'.-')
            #ax2.plot(np.arange(16)+16,single_line,'b.-')
            #ax2.plot(range(16),plot_array[:,i[0],i[1]],'o-')
            #ax2.plot(np.arange(16)+16,plot_array[:,i[0],i[1]],'o-')
            for j in ax:
                j.plot(i[0],i[1], 'kx')
        ax2.set_xlim([0,16*periods-1])
        ymin, ymax = ax2.get_ylim()
        ax2.vlines(np.arange(0,periods*16,16),ymin,ymax)
        for j in ax:
            xticks = j.xaxis.get_major_ticks()
            yticks = j.yaxis.get_major_ticks()
            for jj in xticks: jj.label1.set_visible(False)
            for jj in yticks: jj.label1.set_visible(False)
            j.set_xlim([0,256])
            j.set_ylim([0,256])
        ax2.set_xlabel('frame number (or time)')
        ax2.set_ylabel('pixel intensity (a.u)')
        ax2.grid()
        if save_fig!=None:
            fig.subplots_adjust(hspace=0.05, wspace=0.05,left=0.05, bottom=0.05,top=0.95, right=0.95)
            fig.tight_layout()
            fig.savefig(save_fig)
        fig.canvas.draw(); fig.show()


    def get_mirnov_triggers_new(self,):
        if plot_mirnov_triggers: fig_trig, ax_trig = pt.subplots(nrows= 4, ncols = 4,sharex = True, sharey = True); ax_trig = ax_trig.flatten()
        if plot_mirnov_triggers: fig_phase, ax_phase = pt.subplots(nrows = 2); 
        if plot_mirnov_triggers: fig_hilb, ax_hilb = pt.subplots(nrows= 4, ncols = 4,sharex = True, sharey = True); ax_hilb = ax_hilb.flatten()
        #PLL_node = MDS.Tree('imax',shot).getNode('.PLL
        mirnov_node = MDS.Tree('imax',shot).getNode('.operations.mirnov:a14_15:input_1')

        
    def get_mirnov_triggers(self, plot_mirnov_triggers, old_way = 1):
        if plot_mirnov_triggers: fig_trig, ax_trig = pt.subplots(nrows= 4, ncols = 4,sharex = True, sharey = True); ax_trig = ax_trig.flatten()
        if plot_mirnov_triggers: fig_phase, ax_phase = pt.subplots(nrows = 2); 
        if plot_mirnov_triggers: fig_hilb, ax_hilb = pt.subplots(nrows= 4, ncols = 4,sharex = True, sharey = True); ax_hilb = ax_hilb.flatten()
        phase_average_list = []; phase_std_list = []
        amp_average_list = []; amp_std_list = []
        density_list = []
        if old_way:
            PLL_node_path = '.operations.mirnov.a14_14:input_6'
            PLL_tree = 'h1data'
            mirnov_node_path = '.operations.mirnov:a14_15:input_1'
            mirnov_tree = 'h1data'
        else:
            PLL_node_path = '.waveforms.PLL'
            PLL_tree = 'imax'
            mirnov_node_path = '.waveforms.lock_signal'
            mirnov_tree = 'imax'

        for i, shot in enumerate(self.shot_list):
            PLL_node = MDS.Tree(PLL_tree,shot).getNode(PLL_node_path)
            mirnov_node = MDS.Tree(mirnov_tree,shot).getNode(mirnov_node_path)
            electron_density = MDS.Tree('h1data',shot).getNode('.electr_dens.NE_HET:NE_CENTRE')

            #provide 4ms padding infront of and after the PLL signal goes high
            a = PLL_node.data() > 1.5
            padding = int(0.004/(PLL_node.dim_of().data()[1] - PLL_node.dim_of().data()[0]))
            start_ind = np.argmax(a) - padding
            end_ind = len(a) - np.argmax(a[::-1]) + padding

            mirnov_data, mirnov_time, PLL_data, PLL_time, electr_data, hilb_sig, maxima = self.hilb_transform_check(PLL_node, mirnov_node, electron_density, start_ind, end_ind)
            fig_tmp,ax_tmp = pt.subplots(nrows = 3)
            ax_tmp[0].plot(mirnov_time, mirnov_data)
            ax_tmp[0].plot(PLL_time, PLL_data)
            ax_tmp[1].plot(mirnov_time, np.abs(hilb_sig))
            ax_tmp[2].plot(mirnov_time, np.angle(hilb_sig))
            fig_tmp.canvas.draw(); fig_tmp.show()
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

        if plot_mirnov_triggers:
            fig_trig.canvas.draw(); fig_trig.show()
            fig_hilb.canvas.draw(); fig_hilb.show()
        self.phase_average = np.array(phase_average_list)
        self.phase_std = np.array(phase_std_list)
        self.amp_average = np.array(amp_average_list)
        self.amp_std = np.array(amp_std_list)
        self.electron_dens = np.array(density_list)
        if plot_mirnov_triggers:
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

    def decimate_data(self,decimate_pixel=1):
        tmp_count = 0
        for tmp1 in range(decimate_pixel):
            for tmp2 in range(decimate_pixel):
                if tmp_count==0:
                    tmp_decimate = self.image_array_cal[:,tmp1::decimate_pixel, tmp2::decimate_pixel]
                else:
                    tmp_decimate += self.image_array_cal[:,tmp1::decimate_pixel, tmp2::decimate_pixel]
                tmp_count += 1
        self.image_array_cal = tmp_decimate/tmp_count
        

    def fourier_decomp(self,plot_amps = 0, plot_phases = 0, amp_fig=None, amp_ax = None, phase_fig = None, phase_ax = None, draw = 1, amp_clims = None, fliplr = 1):
        '''Perform an FFT on the image and plot the results if asked to
        '''
        self.fft_values = np.fft.fft(self.image_array_cal, axis = 0)

        n_subplots = self.fft_values.shape[0]/2
        n_cols = int(np.ceil(n_subplots**0.5))
        if n_subplots/float(n_cols)>n_subplots/n_cols:
            n_rows = n_subplots/n_cols + 1
        else:
            n_rows = n_subplots/n_cols

        if plot_amps:
            if amp_fig==None or amp_ax ==None:
                fig, ax = pt.subplots(nrows=n_rows, ncols = n_cols,sharex = True, sharey = True); ax = ax.flatten()
            else:
                fig = amp_fig
                ax = amp_ax

        if plot_phases: 
            if phase_fig==None or phase_ax ==None:
                fig2, ax2 = pt.subplots(nrows=n_rows, ncols = n_cols,sharex = True, sharey = True); ax2 = ax2.flatten()
            else:
                fig2 = phase_fig
                ax2 = phase_ax
        im_list = []
        for i in range(self.fft_values.shape[0]/2):
            if fliplr:
                tmp = np.fliplr(self.fft_values[i,:,:])
            else:
                tmp = self.fft_values[i,:,:]
            if plot_amps:

                im_list.append(ax[i].imshow(np.abs(tmp)*2, interpolation = 'nearest', aspect = 'auto',origin='lower'))                
                if amp_clims!=None:
                    im_list[-1].set_clim(amp_clims[i])

            if plot_phases:
                im2 = ax2[i].imshow(np.angle(tmp), interpolation = 'nearest', aspect = 'auto', origin='lower')
                im2.set_clim([-np.pi,np.pi])
        if plot_amps:
            for i in range(2,len(im_list)):im_list[i].set_clim(im_list[1].get_clim())
            ax[-1].set_xlim([0,self.fft_values.shape[1]])
            ax[-1].set_ylim([0,self.fft_values.shape[2]])
            fig.subplots_adjust(hspace=0.0, wspace=0.0,left=0., bottom=0.,top=0.95, right=0.95)
            if draw:
                fig.suptitle('Fourier amplitude (color range 0, 0.5)')
                fig.canvas.draw(); fig.show()
            #fig.savefig('imax_Fourier_amplitude.png')
        if plot_phases:
            fig2.subplots_adjust(hspace=0.0, wspace=0.0,left=0., bottom=0.,top=0.95, right=0.95)
            if draw:
                fig2.suptitle('Fourier phase (color range -pi,pi)')
                fig2.canvas.draw(); fig2.show()

    def plot_fft_overlay(self,harmonic=1):
        fig, ax = pt.subplots(ncols = 3)
        base_im = ax[0].imshow(np.abs(self.fft_values[0,:,:])*2, cmap = 'bone', interpolation = 'nearest', aspect = 'auto',origin='lower')
        top_im = ax[1].imshow(np.abs(self.fft_values[harmonic,:,:])*2, interpolation = 'nearest', aspect = 'auto',origin='lower', alpha = 0.5)

        top_im = ax[2].imshow(np.abs(self.fft_values[harmonic,:,:])*2, interpolation = 'nearest', aspect = 'auto',origin='lower')
        base_im = ax[2].imshow(np.abs(self.fft_values[0,:,:])*2, cmap = 'bone', interpolation = 'nearest', aspect = 'auto',origin='lower', alpha = 0.9)
        fig.subplots_adjust(hspace=0.0, wspace=0.0,left=0., bottom=0.,top=0.95, right=0.95)
        fig.canvas.draw(); fig.show()

    def calibrate_data(self, dark_shot=1103, white_shot=1107, dark_SPE = None, white_SPE= None, plot_calibration_im = 0, clip = 0.2, plot_cal_images=0, cal_sum = 1, med_filt = 3, mode_amp_filt = None, old_imax = 1, clip_image = 20000):
        '''Calibrate the imax image using a dark and white shot

        need to provide a dark_shot number, white_shot number calibrate by
        sum of all pixels (cal_sum) apply a median filter with side length
        med_filt amply a scaling based on the mode amplitudes from the
        Mirnov signal (mode_amp_filt is an array)

        SH: 3Jul2013
        '''
        self.image_array_cal = np.zeros(self.image_array.shape,dtype=float)
        if old_imax:
            node = '.images'
        else:
            node = '.pimax.images'

        if dark_SPE==None:
            self.dark_image = MDS.Tree('imax',dark_shot).getNode(node).data()[0,:,:]
        else:
            tmp0, tmp1 = SPE_reader.extract_SPE_data(dark_SPE)
            self.dark_image = +tmp0[0,:,:]

        if white_SPE==None:
            self.white_image = MDS.Tree('imax',white_shot).getNode(node).data()[0,:,:]
        else:
            tmp0, tmp1 = SPE_reader.extract_SPE_data(white_SPE)
            self.white_image = +tmp0[0,:,:]

        if plot_calibration_im:
            fig, ax = pt.subplots(nrows= 2,sharex = True, sharey = True); ax = ax.flatten()
            im = ax[0].imshow(self.white_image, interpolation = 'none', aspect = 'auto',origin='lower')
            pt.colorbar(im,ax=ax[0])
            im = ax[1].imshow(self.dark_image, interpolation = 'none', aspect = 'auto', origin='lower')
            pt.colorbar(im,ax=ax[1])
            fig.canvas.draw(); fig.show()
        full_scale = (self.white_image - self.dark_image)
        max_value = np.max(full_scale)
        print 'max value : ', max_value
        # if plot_cal_images: 
        #     fig, ax = pt.subplots(nrows= 4, ncols = 4,sharex = True, sharey = True); ax = ax.flatten()
        #     im_list = []

        if med_filt!=0:
            self.white_image = scipy.signal.medfilt(self.white_image, med_filt)
            self.dark_image = scipy.signal.medfilt(self.dark_image, med_filt)

        for i in range(self.image_array.shape[0]):
            #med filt, subtract dark, divide by clip(white-dark)
            if med_filt!=0:
                print 'med_filt,',
                self.image_array_cal[i,:,:] = scipy.signal.medfilt(self.image_array[i,:,:], med_filt)
            else:
                self.image_array_cal[i,:,:] = self.image_array[i,:,:].copy()

            #subtract the dark image
            self.image_array_cal[i,:,:] = np.clip(self.image_array_cal[i,:,:] - self.dark_image,0,65536)
            tmp = np.clip(self.white_image - self.dark_image,0.2*np.max(self.white_image - self.dark_image), 65536) 
            self.image_array_cal[i,:,:] = self.image_array_cal[i,:,:] / tmp
            
            #print 'clipped image min,max,mean', np.min(tmp), np.max(tmp), np.mean(tmp)
            #np.clip(self.white_image - self.dark_image, max_value*clip,65000)
            #self.image_array_cal[i,:,:] = np.clip(self.image_array[i,:,:] - self.dark_image, 0,50000)/np.clip(self.white_image - self.dark_image, max_value*clip,65000)
            #print np.min(np.clip(self.white_image - self.dark_image, max_value*clip,65000)), np.max(np.clip(self.white_image - self.dark_image, max_value*clip,65000))
            print i,
            if cal_sum:
                print 'cal_sum,', np.sum(self.image_array_cal[i,:,:]), self.electron_dens[i], self.amp_average[i]
                self.image_array_cal[i,:,:] = self.image_array_cal[i,:,:]/np.sum(self.image_array_cal[i,:,:])* 65000.
                print 'cal_sum,', np.sum(self.image_array_cal[i,:,:]), self.electron_dens[i], self.amp_average[i]
            if mode_amp_filt:
                print 'mode_amp_filt',
                self.image_array_cal[i,:,:] = self.image_array_cal[i,:,:] / self.amp_average[i]
            print 'cal_sum,', np.sum(self.image_array_cal[i,:,:]), np.min(self.image_array_cal[i,:,:]),np.max(self.image_array_cal[i,:,:]), np.mean(self.image_array_cal[i,:,:]),self.electron_dens[i], self.amp_average[i]
            # if plot_cal_images:
            #     im_list.append(ax[i].imshow(self.image_array_cal[i,:,:], interpolation = 'none', aspect = 'auto', origin='lower'))
            # print ''
        if plot_cal_images: 
            self.plot_images(cal=1, clim=0)
            # for i in range(0,len(im_list)):im_list[i].set_clim(im_list[0].get_clim())
            # ax[-1].set_xlim([0,512])
            # ax[-1].set_ylim([0,512])
            # fig.subplots_adjust(hspace=0.0, wspace=0.0,left=0., bottom=0.,top=0.95, right=0.95)
            # fig.suptitle('imax_corrected Images')
            # fig.canvas.draw(); fig.show()
            # fig.savefig('imax_corrected_images.png')

    def get_john_analysis(self, plot_john_data = 0):
        '''Extracts Johns results from the MDSplus tree
        SH : 9July2013
        '''
        imax_tree = MDS.Tree('imax',self.shot_list[0])
        self.john_amp_data = imax_tree.getNode('.processed.amplitude').data()
        self.john_phase_data = imax_tree.getNode('.processed.phase').data()
        if plot_john_data:
            fig, ax = pt.subplots(nrows = 3,ncols = 2, sharex = True, sharey = True)
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


def extract_PLL_data(shot_number):
    PLL_node_path = '.mirnov.acq132_8:input_32'
    mirnov_node_path = '.mirnov.acq132_7:input_24'
    read_node_path = '.mirnov.acq132_9:input_32'
    monitor_node_path = '.mirnov.acq132_7:input_32'

    #PLL_node_path = '.waveforms.PLL'
    PLL_tree = 'h1data'
    #mirnov_node_path = '.waveforms.lock_signal'
    mirnov_tree = 'h1data'

    Tree = MDS.Tree(PLL_tree,shot_number)
    PLL_node = Tree.getNode(PLL_node_path)
    mirnov_node = Tree.getNode(mirnov_node_path)
    read_node = Tree.getNode(read_node_path)
    monitor_node = Tree.getNode(monitor_node_path)
    monitor_data = monitor_node.data()
    read_data = read_node.data()
    time_ax = PLL_node.dim_of().data()
    mirnov_data = mirnov_node.data()
    mirnov_data = mirnov_data - np.mean(mirnov_data)
    PLL_data = PLL_node.data()
    return mirnov_data, PLL_data, time_ax, read_data, monitor_data


def PLL_performance(shot_number, n_phases = 4):
    mirnov_data, PLL_data, time_ax, read_data, monitor_data = extract_PLL_data(shot_number)
    time_ax2 = np.arange(time_ax[4],time_ax[-4],1./500000)
    read_data_interp = np.interp(time_ax2, time_ax, read_data)
    tmp = read_data_interp[1:] - read_data_interp[0:-1]
    tmp_argsorted = np.argsort(tmp)
    pos_edges = []; pos_edges_amps = []
    neg_edges = []; neg_edges_amps = []
    for i in tmp_argsorted[0:n_phases]:
        neg_edges.append(time_ax2[i])
        neg_edges_amps.append(0)
    for i in tmp_argsorted[-n_phases:]:
        pos_edges.append(time_ax2[i])
        pos_edges_amps.append(0)

    hilb_sig = scipy.signal.hilbert(mirnov_data)
    #peak = PLL_data>1.5
    maxima = np.zeros(PLL_data.shape,dtype=bool)
    pos_edge = np.zeros(PLL_data.shape,dtype=bool)
    pos_edge[1:] = (PLL_data[1:] - PLL_data[0:-1])>1.5

    #maxima[1:-1] = (PLL_data>2.5)[1:-1] * (PLL_data[1:-1]>PLL_data[0:-2])* (PLL_data[1:-1]>PLL_data[2:])
    #peak = peak * maxima

    peak = pos_edge

    peak_mon = monitor_data>0.5
    maxima_mon = np.zeros(peak.shape,dtype=bool)
    maxima_mon[1:-1] = (monitor_data>0.5)[1:-1] * (monitor_data[1:-1]>monitor_data[0:-2])* (monitor_data[1:-1]>monitor_data[2:])
    peak_mon = peak_mon * maxima_mon
    peak_mon_phases = np.angle(hilb_sig)[peak_mon]
    start = 0
    print neg_edges

    neg_edges.sort()
    pos_edges.sort()
    print neg_edges
    print pos_edges
    phase_average_list = []; phase_std_list = []; phase_time_list = []
    for i in range(len(neg_edges)):
        end = neg_edges[i]
        end = pos_edges[i]
        print start, end
        tmp_truth = peak_mon * (time_ax>start) * (time_ax<end)
        #angles = np.angle(hilb_sig)[tmp_truth]

        phase_complex = np.exp(1j*np.angle(hilb_sig[tmp_truth]))
        phase_average_list.append(np.angle(np.mean(phase_complex)))
        phase_std_list.append(np.sqrt(np.log(1./np.abs(np.mean(phase_complex)))))
        phase_time_list.append(np.mean(time_ax[tmp_truth]))
        #amp_std_list.append(np.std(np.abs(hilb_sig[maxima])))

        print 'n {}, mean {:.2f}deg, std {:.2f}deg'.format(np.sum(tmp_truth),np.rad2deg(phase_average_list[-1]), np.rad2deg(phase_std_list[-1]))
        #print angles, len(angles)
        start = neg_edges[i]
        start = pos_edges[i]
    return phase_average_list, phase_std_list


def PLL_performance2(shot_number, n_phases=4, plot = False):
    mirnov_data, PLL_data, time_ax, read_data, monitor_data = extract_PLL_data(shot_number)
    time_ax2 = np.arange(time_ax[4],time_ax[-4],1./500000)
    #self.extract_PLL_data()
    #time_ax2 = np.arange(self.time_ax[4],self.time_ax[-4],1./500000)
    read_data_interp = np.interp(time_ax2, time_ax, read_data)
    tmp = read_data_interp[1:] - read_data_interp[0:-1]
    tmp_argsorted = np.argsort(tmp)
    pos_edges = []; pos_edges_amps = []
    neg_edges = []; neg_edges_amps = []
    for i in tmp_argsorted[0:n_phases]:
        neg_edges.append(time_ax2[i])
        neg_edges_amps.append(0)
    for i in tmp_argsorted[-n_phases:]:
        pos_edges.append(time_ax2[i])
        pos_edges_amps.append(0)

    hilb_sig = scipy.signal.hilbert(mirnov_data)
    #peak = self.PLL_data>1.5
    phase_data = PLL_data
    phase_data = monitor_data
    maxima = np.zeros(phase_data.shape,dtype=bool)
    pos_edge = np.zeros(phase_data.shape,dtype=bool)
    pos_edge[1:] = (phase_data[1:] - phase_data[0:-1])>1.5

    #maxima[1:-1] = (self.PLL_data>2.5)[1:-1] * (self.PLL_data[1:-1]>self.PLL_data[0:-2])* (self.PLL_data[1:-1]>self.PLL_data[2:])
    #peak = peak * maxima

    peak = pos_edge
    print 'hello!!'
    peak_mon = monitor_data>2.5
    logic_high = monitor_data>0.2
    logic_low = monitor_data<0.1
    logic_high_edges = logic_low[0:-1]*logic_high[1:]
    #logic_high_edges[logic_high_edges[1:]*logic_high_edges[0:-1]]==False
    #logic_high_edges[logic_high_edges[2:]*logic_high_edges[0:-2]]==False
    #logic_high_edges[logic_high_edges[3:]*logic_high_edges[0:-3]]==False
    print 'logic high edges', np.sum(logic_high_edges)


    maxima_mon = np.zeros(peak.shape,dtype=bool)
    maxima_mon[1:-1] = (monitor_data>0.5)[1:-1] * (monitor_data[1:-1]>monitor_data[0:-2])* (monitor_data[1:-1]>monitor_data[2:])
    peak_mon = peak_mon * maxima_mon
    peak_mon_phases = np.angle(hilb_sig)[peak_mon]
    print np.sum(monitor_data[peak_mon]<2.5)
    start = 0
    print neg_edges

    neg_edges.sort()
    pos_edges.sort()
    print 'neg edges: ', neg_edges
    print 'pos edges: ', pos_edges
    phase_average_list = []; phase_std_list = []; phase_time_list = []
    for i in range(len(neg_edges)):
        end = neg_edges[i]
        end = pos_edges[i]
        print start, end
        tmp_truth = peak_mon * (time_ax>start) * (time_ax<end)
        tmp_truth = logic_high * (time_ax>start) * (time_ax<end)
        #angles = np.angle(hilb_sig)[tmp_truth]
        phase_complex = np.exp(1j*np.angle(hilb_sig[tmp_truth]))
        phase_average_list.append(np.angle(np.mean(phase_complex)))
        phase_std_list.append(np.sqrt(np.log(1./np.abs(np.mean(phase_complex)))))
        phase_time_list.append(np.mean(time_ax[tmp_truth]))
        #amp_std_list.append(np.std(np.abs(hilb_sig[maxima])))

        print 'n {}, mean {:.2f}deg, std {:.2f}deg'.format(np.sum(tmp_truth),np.rad2deg(phase_average_list[-1]), np.rad2deg(phase_std_list[-1]))
        #print angles, len(angles)
        start = neg_edges[i]
        start = pos_edges[i]
    if plot:
        PLL_fig, PLL_ax = pt.subplots(nrows=3, sharex = True)
        pub_fig = 1
        if pub_fig:
            cm_to_inch=0.393701
            import matplotlib as mpl
            old_rc_Params = mpl.rcParams
            mpl.rcParams['font.size']=8.0
            mpl.rcParams['axes.titlesize']=8.0#'medium'
            mpl.rcParams['xtick.labelsize']=8.0
            mpl.rcParams['ytick.labelsize']=8.0
            mpl.rcParams['lines.markersize']=5.0
            mpl.rcParams['savefig.dpi']=150
            PLL_fig.set_figwidth(8.48*cm_to_inch)
            PLL_fig.set_figheight(8.48*1.5*cm_to_inch)

        im = PLL_ax[0].specgram(mirnov_data, NFFT = 4096, Fs = 2000, noverlap = 1000,cmap='jet',rasterized = True)
        print im
        print im[-1].get_clim()
        im[-1].set_clim([-80,10])
        #PLL_ax[0].plot(time_ax[::2], PLL_data[::2], label = 'PLL')
        PLL_ax[1].plot(time_ax[::2]*1000, monitor_data[::2], label = 'Camera Gate', rasterized=True)
        PLL_ax[1].plot(time_ax[::2]*1000, read_data[::2], label = 'Camera Readout', rasterized=True)
        PLL_ax[1].plot(time_ax[logic_high_edges]*1000, monitor_data[logic_high_edges], 'o',label = 'Camera Gates')
        #PLL_ax[0].plot(time_ax2, read_data_interp, label = 'read2')
        #PLL_ax[0].plot(pos_edges, pos_edges_amps, 'x')
        #PLL_ax[0].plot(neg_edges, neg_edges_amps, 'o')
        PLL_ax[1].legend(prop={'size':7}, loc='upper left')
        #PLL_ax[1].plot(time_ax[::4], mirnov_data[::4], label = 'Mirnov')
        #PLL_ax[1].plot(time_ax[peak], mirnov_data[peak], 'o')
        #PLL_ax[1].legend(prop={'size':8}, loc='upper left')

        #ax_phase[0].errorbar(range(16), phase_average_list, yerr= phase_std_list,fmt='-o')
        #ax_phase[1].errorbar(range(16), amp_average_list, yerr= amp_std_list,fmt='-o')
        #PLL_ax[2].errorbar(phase_time_list, phase_average_list, yerr = phase_std_list,fmt='-o', label='phase w std')
        #PLL_ax[2].legend(prop={'size':8}, loc='upper left')
        #PLL_ax[2].set_ylim([-np.pi,np.pi])

        #PLL_ax[2].plot(time_ax[::2], np.abs(hilb_sig)[::2], label = 'Mirnov hilb amp')
        #PLL_ax[2].legend()

        #PLL_ax[1].plot(time_ax[::2], np.angle(hilb_sig)[::2], label = 'Mirnov hilb phase')
        PLL_ax[2].plot(time_ax[peak]*1000, np.angle(hilb_sig)[peak], 'o',label = 'Camera Gates')
        #PLL_ax[3].plot(PLL_node.dim_of().data(), PLL_node.data())
        PLL_ax[2].legend(prop={'size':7}, loc='best')
        #PLL_fig.subplots_adjust(hspace=0.0, wspace=0.0,left=0.05, bottom=0.,top=0.95, right=0.95)
        PLL_ax[2].set_xlim([0.03*1000,0.10*1000])
        PLL_ax[0].set_ylim([0,120])
        PLL_ax[1].set_ylim([0,5])
        PLL_ax[2].grid()
        PLL_ax[-1].set_xlabel('Time (ms)')
        PLL_ax[2].set_ylabel('Phase (rad)')
        PLL_ax[1].set_ylabel('Logic Signals (V)')
        PLL_ax[0].set_ylabel('Frequency (kHz)')
        PLL_ax[1].grid()
        PLL_fig.savefig('PLL_multiple_phases.pdf',bbox_inches='tight', pad_inches=0.1)
        PLL_fig.savefig('PLL_multiple_phases.eps',bbox_inches='tight', pad_inches=0.1)
        PLL_fig.canvas.draw(); PLL_fig.show()


    return phase_average_list, phase_std_list


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
    fig, ax = pt.subplots(nrows = 2, ncols = 2,sharex = True, sharey = True); ax = ax.flatten()
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
    



def database_of_shots():
    dict_of_shots = {}

    ##### Broadband, Nandi and John #########
    ##This has a different orientation to the 514nm data..... ##
    dict_of_shots['broadband'] = {}
    dict_of_shots['broadband'] = {'0.33':{},'0.44':{},'0.63':{},'0.83':{}}
    for i in dict_of_shots['broadband']: i = {}
    dict_of_shots['broadband']['0.33']['center'] = {'shot_list': range(71235,71235+16), 'comment':'LHS', 'n':-5, 'm':4}
    dict_of_shots['broadband']['0.44']['center'] = {'shot_list': range(71651,71651+16), 'comment':'RHS', 'n':-5, 'm':4}
    dict_of_shots['broadband']['0.63']['center'] = {'shot_list': range(71102,71102+16), 'comment':'LHS', 'n':-4, 'm':3}
    dict_of_shots['broadband']['0.83']['center'] = {'shot_list': range(70912,70912+16), 'comment':'RHS', 'n':-4, 'm':3}
    for i in dict_of_shots['broadband'].keys():
        dict_of_shots['broadband'][i]['center']['mdsplus_tree'] = 'imax'
        dict_of_shots['broadband'][i]['center']['mdsplus_node'] = '.images'
        dict_of_shots['broadband'][i]['center']['mdsplus_tree_path'] = '/home/srh112/IMAX_DATA/'
        
        dict_of_shots['broadband'][i]['center']['dark_SPE'] = None
        dict_of_shots['broadband'][i]['center']['white_SPE'] = None
        dict_of_shots['broadband'][i]['center']['dark_shot'] = 1103
        dict_of_shots['broadband'][i]['center']['white_shot'] = 1107
        dict_of_shots['broadband'][i]['center']['CCD_side_length'] = 0.01575
        dict_of_shots['broadband'][i]['center']['camera'] = 'pimax'
        dict_of_shots['broadband'][i]['center']['fliplr'] = 0
        dict_of_shots['broadband'][i]['center']['flipud'] = 1
        dict_of_shots['broadband'][i]['center']['rot90'] = 0

    ######## 514nm ############# 
    dict_of_shots['514nm'] = {}
    dict_of_shots['514nm']['0.83'] = {}
    dict_of_shots['514nm']['0.44'] = {}
    dict_of_shots['514nm']['0.33'] = {}
    dict_of_shots['514nm']['0.63'] = {}
    dict_of_shots['514nm']['0.56'] = {}
    dict_of_shots['514nm']['0.69'] = {}
    dict_of_shots['514nm']['0.58'] = {}
    dict_of_shots['514nm']['0.37'] = {}
    dict_of_shots['514nm']['0.28'] = {}

    shot_list = [79887, 79888, 79880, 79889, 79881, 79890, 79882, 79891, 79883, 79892, 79884, 79893, 79885, 79894, 79886, 79895]
    dict_of_shots['514nm']['0.44']['top'] = {'shot_list':shot_list, 'comment':'RHS 514nm', 'n':-5,'m':4}
    shot_list = [79911, 79896, 79904, 79897, 79905, 79898, 79906, 79899, 79907, 79900, 79908, 79901, 79909, 79902, 79910, 79903]
    dict_of_shots['514nm']['0.44']['bottom'] = {'shot_list':shot_list, 'comment':'RHS 514nm', 'n':-5,'m':4}
    shot_list = [79927, 79912, 79920, 79913, 79921, 79914, 79922, 79915, 79923, 79916, 79924, 79917, 79925, 79918, 79926, 79919]
    dict_of_shots['514nm']['0.44']['center'] = {'shot_list':shot_list, 'comment':'RHS 514nm', 'n':-5,'m':4}

    shot_list = [79680, 79712, 79681, 79713, 79682, 79714, 79683, 79715, 79684, 79716, 79685, 79717, 79686, 79718, 79687, 79719]
    dict_of_shots['514nm']['0.83']['top'] = {'shot_list':shot_list, 'comment':'RHS 514nm', 'n':-4,'m':3}
    shot_list = [79688, 79696, 79689, 79697, 79690, 79698, 79691, 79699, 79692, 79700, 79693, 79701, 79694, 79702, 79695, 79703]
    dict_of_shots['514nm']['0.83']['bottom'] = {'shot_list':shot_list, 'comment':'RHS 514nm', 'n':-4,'m':3}
    shot_list = [79677, 79704, 79678, 79705, 79671, 79706, 79672, 79707, 79673, 79708, 79674, 79709, 79675, 79710, 79676, 79711]
    dict_of_shots['514nm']['0.83']['center'] = {'shot_list':shot_list, 'comment':'RHS 514nm', 'n':-4,'m':3}

    shot_list = [80028, 80036, 80029, 80037, 80030, 80038, 80031, 80039, 80032, 80040, 80033, 80041, 80034, 80042, 80035, 80043]
    dict_of_shots['514nm']['0.33']['top'] = {'shot_list':shot_list, 'comment':'RHS 514nm', 'n':-5,'m':4}
    shot_list = [80012, 80020, 80013, 80021, 80014, 80022, 80015, 80023, 80016, 80024, 80017, 80025, 80018, 80026, 80019, 80027]
    dict_of_shots['514nm']['0.33']['bottom'] = {'shot_list':shot_list, 'comment':'RHS 514nm', 'n':-5,'m':4}
    shot_list = [79996, 80004, 79997, 80005, 79998, 80006, 79999, 80007, 80000, 80008, 80001, 80009, 80002, 80010, 80003, 80011]
    dict_of_shots['514nm']['0.33']['center'] = {'shot_list':shot_list, 'comment':'RHS 514nm', 'n':-5,'m':4}

    shot_list = [80056, 80064, 80057, 80065, 80058, 80066, 80059, 80067, 80060, 80068, 80061, 80069, 80062, 80070, 80063, 80071]
    dict_of_shots['514nm']['0.63']['top'] = {'shot_list':shot_list, 'comment':'RHS 514nm', 'n':-4,'m':3}
    shot_list = [80072, 80080, 80073, 80081, 80074, 80082, 80075, 80083, 80076, 80084, 80077, 80085, 80078, 80086, 80079, 80087]
    dict_of_shots['514nm']['0.63']['center'] = {'shot_list':shot_list, 'comment':'RHS 514nm', 'n':-4,'m':3}
    shot_list = [80088, 80096, 80089, 80097, 80090, 80098, 80091, 80099, 80092, 80100, 80093, 80101, 80094, 80102, 80095, 80103]
    dict_of_shots['514nm']['0.63']['bottom'] = {'shot_list':shot_list, 'comment':'RHS 514nm', 'n':-4,'m':3}

    shot_list = [80134, 80142, 80118, 80135, 80119, 80136, 80120, 80137, 80121, 80138, 80122, 80139, 80123, 80140, 80133, 80141]
    dict_of_shots['514nm']['0.56']['bottom'] = {'shot_list':shot_list, 'comment':'RHS 514nm', 'n':-4,'m':3}
    shot_list = [80150, 80158, 80143, 80151, 80144, 80152, 80145, 80153, 80146, 80154, 80147, 80155, 80148, 80156, 80149, 80157]
    dict_of_shots['514nm']['0.56']['center'] = {'shot_list':shot_list, 'comment':'RHS 514nm', 'n':-4,'m':3}
    shot_list = [80166, 80176, 80159, 80167, 80160, 80168, 80161, 80169, 80162, 80170, 80163, 80171, 80164, 80172, 80165, 80173]
    dict_of_shots['514nm']['0.56']['top'] = {'shot_list':shot_list, 'comment':'RHS 514nm', 'n':-4,'m':3}

    shot_list = [80444, 80437, 80445, 80438, 80446, 80439, 80447, 80440, 80448, 80441, 80449, 80442, 80450, 80435, 80443, 80436]
    dict_of_shots['514nm']['0.69']['bottom'] = {'shot_list':shot_list, 'comment':'RHS 514nm', 'n':-4,'m':3}
    shot_list = [80428, 80421, 80429, 80422, 80430, 80423, 80431, 80424, 80432, 80425, 80433, 80426, 80434, 80419, 80427, 80420]
    dict_of_shots['514nm']['0.69']['center'] = {'shot_list':shot_list, 'comment':'RHS 514nm', 'n':-4,'m':3}
    shot_list = [80461, 80454, 80462, 80455, 80463, 80456, 80464, 80457, 80465, 80458, 80466, 80459, 80467, 80452, 80460, 80453]
    dict_of_shots['514nm']['0.69']['top'] = {'shot_list':shot_list, 'comment':'RHS 514nm', 'n':-4,'m':3}

    shot_list = [80522, 80507, 80515, 80508, 80516, 80509, 80517, 80510, 80518, 80511, 80523, 80512, 80520, 80513, 80521, 80514]
    dict_of_shots['514nm']['0.58']['bottom'] = {'shot_list':shot_list, 'comment':'RHS 514nm', 'n':-4,'m':3}
    shot_list = [80491, 80499, 80492, 80500, 80493, 80501, 80494, 80502, 80495, 80503, 80496, 80504, 80497, 80505, 80498, 80506]
    dict_of_shots['514nm']['0.58']['center'] = {'shot_list':shot_list, 'comment':'RHS 514nm', 'n':-4,'m':3}
    shot_list = [80490, 80475, 80483, 80476, 80484, 80477, 80485, 80478, 80486, 80479, 80487, 80480, 80488, 80481, 80489, 80482]
    dict_of_shots['514nm']['0.58']['top'] = {'shot_list':shot_list, 'comment':'RHS 514nm', 'n':-4,'m':3}

    shot_list = [80557, 80565, 80558, 80551, 80559, 80552, 80560, 80553, 80561, 80554, 80562, 80555, 80563, 80556, 80564, 80549]
    dict_of_shots['514nm']['0.37']['bottom'] = {'shot_list':shot_list, 'comment':'RHS 514nm', 'n':-5,'m':4}
    shot_list = [80574, 80567, 80575, 80568, 80576, 80569, 80577, 80570, 80578, 80571, 80579, 80572, 80580, 80573, 80581, 80566]
    dict_of_shots['514nm']['0.37']['center'] = {'shot_list':shot_list, 'comment':'RHS 514nm', 'n':-5,'m':4}
    shot_list = [80590, 80583, 80591, 80584, 80592, 80585, 80593, 80586, 80594, 80587, 80595, 80588, 80596, 80589, 80597, 80582]
    dict_of_shots['514nm']['0.37']['top'] = {'shot_list':shot_list, 'comment':'RHS 514nm', 'n':-5,'m':4}

    shot_list = [80640, 80649, 80641, 80650, 80642, 80643, 80651, 80652, 80644, 80653, 80645, 80654, 80646, 80655, 80656, 80647]
    dict_of_shots['514nm']['0.28']['bottom'] = {'shot_list':shot_list, 'comment':'RHS 514nm', 'n':-5,'m':4}
    shot_list = [80637, 80622, 80630, 80623, 80631, 80624, 80632, 80625, 80633, 80626, 80634, 80627, 80635, 80628, 80636, 80629]
    dict_of_shots['514nm']['0.28']['center'] = {'shot_list':shot_list, 'comment':'RHS 514nm', 'n':-5,'m':4}
    shot_list = [80620, 80612, 80605, 80621, 80613, 80614, 80607, 80608, 80615, 80616, 80609, 80617, 80610, 80618, 80611, 80619]
    dict_of_shots['514nm']['0.28']['top'] = {'shot_list':shot_list, 'comment':'RHS 514nm', 'n':-5,'m':4}

    for i in dict_of_shots['514nm'].keys():
        for j in dict_of_shots['514nm'][i].keys():
            dict_of_shots['514nm'][i][j]['mdsplus_tree'] = 'h1data'
            dict_of_shots['514nm'][i][j]['mdsplus_node'] = '.mirnov.pimax.pimax.images'
            dict_of_shots['514nm'][i][j]['dark_SPE'] = '/home/srh112/code/python/imax/imax_calibrations/dark_514_x1_400.spe'
            dict_of_shots['514nm'][i][j]['white_SPE'] = '/home/srh112/code/python/imax/imax_calibrations/white_514_1x.spe'
            #dict_of_shots['514nm'][i][j]['dark_SPE'] = 'pimax4_1000.spe'
            #dict_of_shots['514nm'][i][j]['white_SPE'] = 'pimax4_1004.spe'
            dict_of_shots['514nm'][i][j]['dark_shot'] = None
            dict_of_shots['514nm'][i][j]['white_shot'] = None
            dict_of_shots['514nm'][i][j]['mdsplus_tree_path'] = r'h1data::/data/h1/recent;h1data::/data/h1/sorted/~f~e~d~c'
            dict_of_shots['514nm'][i][j]['CCD_side_length'] = 0.0131
            dict_of_shots['514nm'][i][j]['camera'] = 'pimax4'
            dict_of_shots['514nm'][i][j]['fliplr'] = 0
            dict_of_shots['514nm'][i][j]['flipud'] = 0
            dict_of_shots['514nm'][i][j]['rot90'] = 0


    base_dir = '/home/srh112/Desktop/tomo_stuff/IMAX_H_1_Shots/706_728_imaging_Dec_2010/'
    dict_of_shots['707nm'] = {}
    dict_of_shots['707nm']['0.83'] = {}

    #CENTRE VIEW 707nm
    file_list = ['{}/{}.SPE'.format(base_dir,i) for i in range(69009,69009+16)]
    dict_of_shots['707nm']['0.83']['center'] = {'shot_list': None, 'comment':'RHS', 'n':-4, 'm':3, 'SPE_files' : file_list}

    #TOP VIEW
    shot_list = range(69026, 69043); shot_list.remove(69039)
    file_list = ['{}/{}.SPE'.format(base_dir,i) for i in shot_list]
    dict_of_shots['707nm']['0.83']['top'] = {'shot_list': None, 'comment':'RHS', 'n':-4, 'm':3, 'SPE_files' : file_list}
    

    #BOTTOM VIEW
    shot_list = [69043, 69044, 69045, 69046, 69047, 69048, 69049,69050,69051,69052,69053,69054,69055,69056,69057,69058]
    file_list = ['{}/{}.SPE'.format(base_dir,i) for i in shot_list]
    file_list = ['{}/{}.SPE'.format(base_dir,i) for i in range(69043,69043+16)]
    dict_of_shots['707nm']['0.83']['bottom'] = {'shot_list': None, 'comment':'RHS', 'n':-4, 'm':3, 'SPE_files' : file_list}
    return dict_of_shots


        
