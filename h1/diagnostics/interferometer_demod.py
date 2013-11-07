import scipy.special as spec
import MDSplus as MDS
import numpy as np
import matplotlib.pyplot as pt
import scipy.signal as sig
import scipy.optimize as opt


class interf_demod(object):
    def __init__(self,shot_number, f_m, dig_mult, start_samples = 10000, force_symmetric = 1, individual_eps_calc = 1):
        '''Class for demodulating the sinusoidally modulated interferometer
        SRH : 5Nov2013
        '''
        self.shot_number = shot_number
        self.f_m = f_m
        self.dig_mult = dig_mult
        self.start_samples = start_samples
        self.f_s = self.f_m * self.dig_mult
        self.get_up_to_carriers()

    def get_up_to_carriers(self,force_symmetric = 1, individual_eps_calc = 1):
        '''Convenience function to do several things at once - it gets you the carriers as a function of time
        SRH: 7Nov2013
        '''
        self.extract_mdsplus_data()
        self.calc_carriers_time(force_symmetric = force_symmetric, individual_eps_calc = individual_eps_calc)
        self.calc_carriers_time_rfft()

    def extract_mdsplus_data(self,):
        '''Get the data out of MDSplus
        SRH: 5Nov2013
        '''
        T = MDS.Tree('h1data',self.shot_number)
        count = 0
        self.chan_list = range(21)
        self.rms_pre_shot = []; self.dc_pre_shot = []
        self.rms_full_shot = []; self.dc_full_shot = []
        for i, chan_num in enumerate(self.chan_list):
            digi, chan = self.chan_num_to_digi(chan_num)
            n = T.getNode('.electr_dens.camac.A14_{}.input_{}'.format(digi,chan))
            signal = n.data()
            if count == 0:
                keep_length = len(signal) - (len(signal)%self.dig_mult)
                self.raw_signals = np.zeros((len(self.chan_list), keep_length),dtype = float)
            self.raw_signals[i,:] = +signal[:keep_length]
            #signal = signal[:keep_length]
            dc, ac_rms = self.calculate_ac_dc(signal, full = 0)
            self.rms_pre_shot.append(ac_rms)
            self.dc_pre_shot.append(dc)
            dc, ac_rms = self.calculate_ac_dc(signal, full = 1)
            self.rms_full_shot.append(ac_rms)
            self.dc_full_shot.append(dc)

            #self.rms_pre_shot.append(np.sqrt(1./self.start_samples*np.sum((signal[:self.start_samples]-np.mean(signal[:self.start_samples]))**2)))
            count += 1
        self.t = np.arange(-self.start_samples, keep_length - self.start_samples) * 1./self.f_s

    def calculate_ac_dc(self,signal, full=1):
        '''Calculate the ac (rms) and dc components of the raw
        interferometer channels 
        SRH: 5Nov2013 
        '''
        if full:
            stop_loc = len(signal)
        else:
            stop_loc = self.start_samples
        ac_rms = np.sqrt(1./self.start_samples*np.sum((signal[:stop_loc]-np.mean(signal[:stop_loc]))**2))
        dc = np.mean(signal[:stop_loc])
        return dc, ac_rms

    def plot_raw_signals(self, subtract_dc = 1):
        '''Plot the raw digitised signals
        SRH: 5Nov2013
        '''
        n_chans = self.raw_signals.shape[0]
        fig, ax = pt.subplots(nrows = n_chans,sharex =True, sharey=True)
        for i in range(self.raw_signals.shape[0]):
            if subtract_dc:
                ax[i].plot(self.t, self.raw_signals[i,:] - self.dc_pre_shot[i])
            else:
                ax[i].plot(self.t, self.raw_signals[i,:])
        fig.subplots_adjust(hspace=0, wspace=0,left=0.05, bottom=0.05,top=0.95, right=0.95)
        ax[0].set_ylim([-0.4,0.4])
        ax[-1].set_xlabel('Sample Number')
        fig.canvas.draw(); fig.show()

    def chan_num_to_digi(self,chan_num):
        '''Maps the logical channel number with 0 being plasma top and
        20 plasma bottom to the digitiser channels based on
        A14_XX:input_X 
        SRH: 5Nov2013
        '''
        chan_num_dict = {'0':[22,5],'1':[22,4], '2':[22,3],'3':[22,2],
                         '4':[22,1], '5':[21,6], '6':[21,5],
                         '7':[21,4], '8':[21,3], '9':[21,2],
                         '10':[21,1], '11':[22,6], '12':[23,1],
                         '13':[23,2], '14':[23,3], '15':[23,4],
                         '16':[23,5], '17':[23,6], '18':[24,1],
                         '19':[24,2], '20':[24,3]}
        return chan_num_dict[str(chan_num)]

    def calc_epsilon(self,):
        '''calculate epsilon, the timing error
        SRH: 5Nov2013
        '''
        self.e_recon = np.angle(self.signal_fft[:,self.f_m_loc]/1j)
        self.e_recon2 = np.angle(self.signal_fft[:,self.f_m_loc*2]/1j)
        self.e_recon3 = np.angle(self.signal_fft[:,self.f_m_loc*3]/1j)

    def calc_carriers_time(self, n_carriers = 6, force_symmetric = 0, individual_eps_calc = 0, tophat_width_prop = 1):
        '''Calculate the carriers as a function of time
        SRH: 5Nov2013
        '''
        self.signal_fft = np.fft.fft(self.raw_signals)
        self.fft_ax = np.fft.fftfreq(self.signal_fft.shape[1], 1./self.f_s)

        #Find the modulation frequency
        self.f_m_loc = np.argmin(np.abs(self.fft_ax - self.f_m))
        self.f_m2_loc = np.argmin(np.abs(self.fft_ax - 2*self.f_m))

        #Get epsilon for the first harmonic
        self.e_recon = np.angle(self.signal_fft[:,self.f_m_loc]/1j)

        self.carriers = np.zeros((n_carriers,self.signal_fft.shape[0], self.signal_fft.shape[1]),dtype = complex)
        #Width of the tophat
        self.width = int(self.f_m_loc) * tophat_width_prop
        if (self.width%2)!=0: self.width = self.width - 1
        for n in range(n_carriers):
            if (n%2)==0: 
                mult = 1
            else:
                mult = -1j
            n_freq_loc = np.argmin(np.abs(self.fft_ax - n*self.f_m))
            tmp = self.signal_fft *0
            if n==0:
                tmp[:,0:self.width/2] = self.signal_fft[:, n_freq_loc:n_freq_loc+self.width/2]
                #print len(np.conj(signal_fft[n_freq_loc:n_freq_loc+self.width/2]))
                tmp[:,-self.width/2+1:] = self.signal_fft[:,-self.width/2+1:]
            else:
                tmp[:,0:self.width/2] = self.signal_fft[:,n_freq_loc:n_freq_loc+self.width/2]
                #tmp[0]= signal_fft[n_freq_loc]
                #tmp[1:self.width/2] = signal_fft[n_freq_loc+1:n_freq_loc+self.width/2]
                #tmp[-self.width/2+1:] = signal_fft[n_freq_loc-self.width/2+1:n_freq_loc]
                tmp[:,-self.width/2+1:] = self.signal_fft[:,n_freq_loc-self.width/2+1:n_freq_loc]
            if n==0:
                e_recon_tmp = self.e_recon
            elif (n%2)==0:
                e_recon_tmp = np.angle(self.signal_fft[:,self.f_m_loc*n])/n
            else:
                e_recon_tmp = np.angle(self.signal_fft[:,self.f_m_loc*n]/-1j)/n
            #e_recon_tmp = e_recon
            a = tmp* np.exp(-1j*n*e_recon_tmp[:,np.newaxis])*mult
            #print n, a[0,:]
            if force_symmetric:
                a[:,-self.width/2+1:] = np.conj(a[:,1:self.width/2])[:,::-1]
            #h_n.append(np.fft.ifft(a))
            self.carriers[n,:,:] = np.fft.ifft(a)

    def calc_carriers_time_rfft(self, n_carriers = 6, force_symmetric = 0, individual_eps_calc = 0, tophat_width_prop = 1):
        '''Same as calc_carriers_time except using the real fft and
        real ifft which forces symmetry, and removes some of the
        problems with getting the frequency shifting wrong 
        SRH:5Nov2013
        '''
        self.e_recon = np.angle(self.signal_fft[:,self.f_m_loc]/1j)
        self.carriers_rfft = np.zeros((n_carriers,self.signal_fft.shape[0], self.signal_fft.shape[1]),dtype = complex)
        self.signal_rfft = np.fft.rfft(self.raw_signals)
        #Width of the tophat
        self.width = int(self.f_m_loc) * tophat_width_prop
        if (self.width%2)!=0: self.width = self.width - 1
        for n in range(n_carriers):
            if (n%2)==0: 
                mult = 1
            else:
                mult = -1j
            n_freq_loc = np.argmin(np.abs(self.fft_ax - n*self.f_m))
            tmp = self.signal_rfft *0
            tmp[:,0:self.width/2] = self.signal_fft[:, n_freq_loc:n_freq_loc+self.width/2]
            if n==0:
                e_recon_tmp = self.e_recon
            elif (n%2)==0:
                e_recon_tmp = np.angle(self.signal_fft[:,self.f_m_loc*n])/n
            else:
                e_recon_tmp = np.angle(self.signal_fft[:,self.f_m_loc*n]/-1j)/n
            a = tmp* np.exp(-1j*n*self.e_recon[:,np.newaxis])*mult
            #h_n.append(np.fft.ifft(a))
            self.carriers_rfft[n,:,:] = np.fft.irfft(a)

    def plot_carriers(self,):
        pass

    def phi1_min(self,phi_1_test,ratio, bessel_numerator, bessel_denomenator):
        '''Minimisation function for optimizer to calculate phi1
        SRH: 5Nov2013
        '''
        return np.abs(ratio - spec.jn(bessel_numerator,phi_1_test)/spec.jn(bessel_denomenator,phi_1_test))

    def calc_phi1_means_std(self, harm1, harm2, start_point = 0, end_point = None, plot_fig = True):
        '''Calculates phi1 - modulation depth for each of the channels individually using the mean value from start_point to end_point
        Produces plots of the mean values, with error bars

        SRH : 7Nov2013
        '''
        if end_point == None: end_point = self.start_samples
        self.ratio = np.abs(self.carriers[harm1,:,start_point:end_point]) / np.abs(self.carriers[harm2,:,start_point:end_point])
        self.ratio_means = np.mean(self.ratio, axis = -1)
        self.ratio_stds = np.std(self.ratio, axis = -1)
        if plot_fig:
            fig, ax = pt.subplots(nrows = 2, sharex = True)
            ax[0].errorbar(range(len(self.ratio_means)), self.ratio_means, yerr=self.ratio_stds)
            ax[0].set_ylim([-3,3])
        self.phi1_mean = []
        self.phi1_std_upper = []
        self.phi1_std_lower = []
        for i in range(self.ratio_means.shape[0]):
            print opt.fmin(self.phi1_min, np.deg2rad(135), args=(self.ratio_means[i],harm1, harm2),disp=0, xtol=0.00001)
            self.phi1_mean.append(opt.fmin(self.phi1_min, np.deg2rad(135), args=(self.ratio_means[i],harm1, harm2),disp=0, xtol=0.00001)[0])
            self.phi1_std_upper.append(opt.fmin(self.phi1_min, np.deg2rad(135), args=(self.ratio_means[i]+self.ratio_stds[i], harm1, harm2),disp=0, xtol=0.00001)[0])
            self.phi1_std_lower.append(opt.fmin(self.phi1_min, np.deg2rad(135), args=(self.ratio_means[i]-self.ratio_stds[i], harm1, harm2),disp=0, xtol=0.00001)[0])
        if plot_fig:
            tmp1 = np.abs(np.array(self.phi1_mean) - np.array(self.phi1_std_lower))
            tmp2 = np.abs(np.array(self.phi1_mean) - np.array(self.phi1_std_upper))
            err_vals = np.maximum(tmp1, tmp2)
            ax[1].errorbar(range(len(self.ratio_means)), np.rad2deg(self.phi1_mean), yerr=np.rad2deg(err_vals))
            ax[1].set_ylim([0,180])
            fig.canvas.draw(); fig.show()

    def extract_harmonic_amps(self,start_point = 0, end_point = None, show_fig = True):
        '''Get the amplitude of the harmonics to check how well the
        sinusoidal modulation is working 
        SRH: 5Nov2013 
        '''
        if end_point == None: end_point = self.start_samples
        end_point += (start_point - end_point)%self.dig_mult
        tmp_fft = np.fft.fft(self.raw_signals[:,start_point:end_point])
        tmp_fft_ax = np.fft.fftfreq(tmp_fft.shape[1], 1./self.f_s)
        f_m_loc_tmp = np.argmin(np.abs(tmp_fft_ax-self.f_m))
        n_harmonics = int((self.f_s/2.)/self.f_m)
        if (self.f_m*n_harmonics)>=(self.f_s/2.): n_harmonics -= 1
        self.harmonic_list = np.arange(n_harmonics+1)*f_m_loc_tmp
        self.harmonics = tmp_fft[:,self.harmonic_list]
        if show_fig:
            fig, ax = pt.subplots(nrows=5, ncols = 5, sharex = True, sharey = True); ax = ax.flatten()
            for i in range(tmp_fft.shape[0]):
                ax[i].plot(tmp_fft_ax, np.abs(tmp_fft[i,:]))
                ax[i].plot(tmp_fft_ax[self.harmonic_list], np.abs(self.harmonics[i,:]), 'o')
            ax[-1].set_xlim([0,np.max(tmp_fft_ax)])
            ax[-1].set_ylim([0,np.max(np.abs(self.harmonics))])
            fig.subplots_adjust(hspace=0.015, wspace=0.015,left=0.10, bottom=0.10,top=0.95, right=0.95)
            fig.canvas.draw(); fig.show()

    def compare_harmonics(self,phi1 = None, show_fig = True):
        '''Compare the harmonic strengths with the amplitudes they
        should be. Can force a phi1 value to use, or use the already
        calculated one

        SRH: 5Nov2013
        '''
        chans, max_harmonic = self.harmonics.shape
        if phi1==None: 
            phi1 = np.array(self.phi1_mean)
        else:
            phi1 = np.array([phi1] * chans)
        bessel_amps = []
        for i in phi1:
            bessel_amps.append(spec.jn(range(max_harmonic), i))
        bessel_amps = np.array(bessel_amps)
        print bessel_amps.shape, self.harmonics.shape

        if show_fig:
            fig, ax = pt.subplots(nrows=5, ncols = 5, sharex = True, sharey = True); ax = ax.flatten()
            for i in range(chans):
                ax[i].plot(range(1,max_harmonic),np.abs(self.harmonics[i,1:])*np.abs(bessel_amps[i,1])/ np.abs(self.harmonics[i,1]),'o-')
                #ax[i].plot(range(1,max_harmonic),np.sign(np.real(self.harmonics[i,1:]))*np.abs(self.harmonics[i,1:])*np.abs(bessel_amps[i,1])/ np.abs(self.harmonics[i,1]),'o-')
                ax[i].plot(range(1,max_harmonic), bessel_amps[i,1:],'x-')
            ax[-1].set_xlim([0,max_harmonic])
            ax[-1].set_ylim([0,np.max(bessel_amps)])
        
            fig.subplots_adjust(hspace=0.015, wspace=0.015,left=0.10, bottom=0.10,top=0.95, right=0.95)
            fig.canvas.draw(); fig.show()

    def calculate_phase(self, phi1, odd, even, use_rfft = 0, clim = [0,5], show_fig = True, save_fig_name = None):
        '''Calculate the actual phase shift caused by the plasma. Can use the rfft
        SRH: 5Nov2013
        '''
        if use_rfft: 
            carriers = self.carriers_rfft
        else:
            carriers = self.carriers
        odd_part = np.real(carriers[odd,:,:]/spec.jn(odd,phi1))
        even_part = np.real(carriers[even,:,:]/spec.jn(even,phi1))
        output = np.unwrap(np.arctan2(odd_part, even_part), axis = -1)
        self.initial_phase = np.mean(output[:,30:1000], axis = -1)
        output = (output - self.initial_phase[:,np.newaxis])
        #if np.mean(output, axis = -1)<0: output*=-1
        self.output = +output
        if show_fig:
            self.im_fig, self.im_ax = pt.subplots(ncols = 2, sharey=True)
            self.im =self.im_ax[0].imshow(np.abs(self.output),interpolation = 'nearest', aspect='auto',cmap = 'hot', extent=[self.t[0],self.t[-1],self.output.shape[0],0], origin = 'upper')
        times = [0.02,0.04,0.06]
        for t in times:
            t_loc = np.argmin(np.abs(self.t - t))
            plot_vals = np.abs(self.output[:,t_loc])<10
            if show_fig:
                self.im_ax[1].plot(np.abs(self.output[:,t_loc][plot_vals]),np.arange(self.output.shape[0])[plot_vals],'x-')
                self.im_ax[0].vlines(t, self.output.shape[0],0, colors = 'b')
        if show_fig:
            self.im_ax[1].grid()
            self.im_ax[1].set_xlim(clim)
            #self.im =self.im_ax.imshow(np.abs(np.flipud(self.output)), aspect='auto',cmap = 'hot', extent=[self.t[0],self.t[-1],self.output.shape[0],0])
            self.im.set_clim(clim)
            self.im_ax[0].set_xlim([-0.005,0.12])
            self.im_ax[1].set_ylim([0.,self.output.shape[0]])
            pt.colorbar(self.im,ax = self.im_ax[0])
            if save_fig_name!=None:
                self.im_fig.suptitle(save_fig_name.rstrip('.png'))
                self.im_fig.savefig(save_fig_name)
            self.im_fig.canvas.draw(); self.im_fig.show()


    def calculate_phase_simul(self, phi1, use_rfft = 0, clim = [0,5], show_fig = True):
        '''Calculate the actual phase shift caused by the plasma using all available carriers
        SRH: 5Nov2013
        '''
        if use_rfft: 
            carriers = self.carriers_rfft
        else:
            carriers = self.carriers

        harms, chans, Ns = carriers.shape
        self.J = np.zeros((harms, 3), dtype = float)
        self.J[::2,2] = spec.jn(np.arange(0, harms,2), phi1)
        self.J[1::2,1] = spec.jn(np.arange(1, harms,2), phi1)
        self.J[0,0] = 1
        self.J_inv = np.linalg.pinv(self.J)
        self.h = np.zeros((harms, Ns),dtype = float)
        self.result_list = []
        if show_fig:
            fig, ax = pt.subplots(nrows = 1); ax = [ax]
        self.output_phase_simul = np.zeros((chans, Ns), dtype = float)
        for i in range(chans):
            tmp = np.dot(self.J_inv, np.real(carriers[:,i,:]))
            self.result_list.append(tmp)
            output = np.unwrap(np.arctan2(tmp[2,:], tmp[1,:]))
            initial_phase = np.mean(output[30:1000])
            print initial_phase,
            self.output_phase_simul[i,:] = (output - initial_phase)
            tmp_mean = np.mean(self.output_phase_simul[i,:])
            print np.cos(tmp_mean)*np.sin(tmp_mean), tmp_mean, np.sign(np.cos(tmp_mean)*np.sin(tmp_mean)) == np.sign(tmp_mean)
            
            if show_fig:
                ax[0].plot(self.t, np.abs(self.output_phase_simul[i,:]))
        print ''

        if show_fig:
            self.im_fig2, self.im_ax2 = pt.subplots(ncols = 2, sharey=True)
            self.im2 =self.im_ax2[0].imshow(np.abs(self.output_phase_simul),interpolation = 'nearest', aspect='auto',cmap = 'hot', extent=[self.t[0],self.t[-1],self.output.shape[0],0], origin = 'upper')
            self.im2.set_clim(clim)
        times = [0.02,0.04,0.06]
        for t in times:
            t_loc = np.argmin(np.abs(self.t - t))
            plot_vals = np.abs(self.output_phase_simul[:,t_loc])<10
            if show_fig:
                self.im_ax2[1].plot(np.abs(self.output_phase_simul[:,t_loc][plot_vals]),np.arange(self.output_phase_simul.shape[0])[plot_vals],'x-')
                self.im_ax2[0].vlines(t, self.output_phase_simul.shape[0],0, colors = 'b')
        if show_fig:
            self.im_ax2[1].grid()
            self.im_ax2[1].set_xlim(clim)
            self.im_ax2[0].set_xlim([-0.005,0.12])
            self.im_ax2[1].set_ylim([0.,self.output_phase_simul.shape[0]])
            pt.colorbar(self.im2,ax = self.im_ax2[0])
            self.im_fig2.canvas.draw(); self.im_fig2.show()
            fig.canvas.draw(); fig.show()
        

    def min_func_phi1_simul(self, phi1, harms, chans, h):
        '''Function to mimise the difference between h and J dot S

        Outstanding questions: How is the affected by bad channels?

        SRH: 7Nov2013
        '''
        J = np.zeros((harms, 3), dtype = float)
        J[::2,2] = spec.jn(np.arange(0, harms,2), phi1)
        J[1::2,1] = spec.jn(np.arange(1, harms,2), phi1)
        J[0,0] = +1
        J_tile = np.tile(J,(chans,1))
        J_inv = np.linalg.pinv(J_tile)
        S = np.dot(J_inv, h)
        return np.sum((np.dot(J_tile, S) - h)**2)


    def calculate_phase_phi_simul(self, phi1, use_rfft = 0, clim = [0,5], show_fig = True, channel_mask = None):
        '''Calculate phi1 (modulation depth) simulatneously using all interferometer channels

        SRH: 5Nov2013
        '''
        if use_rfft: 
            carriers = self.carriers_rfft
        else:
            carriers = self.carriers
        n_samples = 2000
        harms, chans, Ns = carriers.shape

        self.J = np.zeros((harms, 3), dtype = float)

        #self.h = np.zeros((harms*chans, Ns),dtype = float)
        self.h = np.zeros((harms*chans, n_samples),dtype = float)
        for i in range(chans):
            self.h[i*harms:(i+1)*harms,:] = np.real(carriers[:,i,:n_samples])
        phi1_list = np.linspace(0.3,1.3*np.pi,50)
        error_list = []
        for phi1 in phi1_list:
            self.J[::2,2] = spec.jn(np.arange(0, harms,2), phi1)
            self.J[1::2,1] = spec.jn(np.arange(1, harms,2), phi1)
            self.J[0,0] = 1
            self.J_tile = np.tile(self.J,(chans,1))
            self.J_inv = np.linalg.pinv(self.J_tile)
            self.S = np.dot(self.J_inv, self.h)
            tmp = np.dot(self.J_tile, self.S)
            error_list.append(np.sum((tmp - self.h)**2))
            print phi1, error_list[-1]
        phi_useful = phi1_list[np.argmin(error_list)]
        self.J[::2,2] = spec.jn(np.arange(0, harms,2), phi_useful)
        self.J[1::2,1] = spec.jn(np.arange(1, harms,2), phi_useful)
        self.J[0,0] = 1
        self.J_tile = np.tile(self.J,(chans,1))
        self.J_inv = np.linalg.pinv(self.J_tile)
        self.S = np.dot(self.J_inv, self.h)
        tmp = np.dot(self.J_tile, self.S)
        
        a = opt.fmin(self.min_func_phi1_simul, np.deg2rad(135), args=(harms, chans, self.h),disp=1, xtol=0.00001)
        a2 = opt.minimize(self.min_func_phi1_simul, np.deg2rad(135), args=(harms, chans, self.h),method = 'L-BFGS-B', bounds = [(0.05,np.pi)])
        a3 = opt.minimize(self.phase_phi1_phi2_error, [np.deg2rad(135), np.deg2rad(140),0.5,0.5], args=(harms, chans, self.h),method = 'L-BFGS-B', bounds = [(0.05,np.pi), (0.05,np.pi),(0.001,20),(0.001,20)])
        print a
        print a2
        print a3
        if show_fig:
            fig, ax = pt.subplots()
            ax.plot(np.rad2deg(phi1_list), error_list)
            ax.plot(np.rad2deg(a[0]), self.min_func_phi1_simul(a[0],harms, chans, self.h) ,'o')
            ax.plot(np.rad2deg(a2['x'][0]), self.min_func_phi1_simul(a2['x'][0],harms, chans, self.h) ,'o')
            fig.canvas.draw(); fig.show()


    def phase_phi1_phi2_error(self, guess, harms, chans, h):
        '''Attempt at using two different modulation depths

        SRH: 7Nov2013
        '''
        J = np.zeros((harms, 3), dtype = float)
        phi1, phi2, amp1, amp2 = guess
        J[::2,2] = amp1*spec.jn(np.arange(0, harms,2), phi1) + amp2*spec.jn(np.arange(0, harms,2), phi1)
        J[1::2,1] = amp1*spec.jn(np.arange(1, harms,2), phi1) + amp2*spec.jn(np.arange(1, harms,2), phi2)
        J[0,0] = 1
        J_tile = np.tile(J,(chans,1))
        J_inv = np.linalg.pinv(J_tile)
        S = np.dot(J_inv, h)
        tmp = np.dot(J_tile, S)
        return np.sum((tmp - self.h)**2)


    def calculate_phase_phi1_phi2(self, use_rfft = 0, clim = [0,5]):
        '''Calculate the actual phase shift caused by the plasma using multiple carriers
        SRH: 5Nov2013
        '''
        if use_rfft: 
            carriers = self.carriers_rfft
        else:
            carriers = self.carriers
        n_samples = 2000
        harms, chans, Ns = carriers.shape

        self.J = np.zeros((harms, 3), dtype = float)

        #self.h = np.zeros((harms*chans, Ns),dtype = float)
        self.h = np.zeros((harms*chans, n_samples),dtype = float)
        for i in range(chans):
            self.h[i*harms:(i+1)*harms,:] = np.real(carriers[:,i,:n_samples])
        min_rad = 0.3
        max_rad = 1.3*np.pi
        
        phi1_list = np.linspace(min_rad,max_rad,50)
        error_list = []
        error_matrix = np.zeros((phi1_list.shape[0], phi1_list.shape[0]),dtype=float)

        for i, phi1 in enumerate(phi1_list):
            for j, phi2 in enumerate(phi1_list):
                self.J[::2,2] = spec.jn(np.arange(0, harms,2), phi1) + spec.jn(np.arange(0, harms,2), phi1)
                self.J[1::2,1] = spec.jn(np.arange(1, harms,2), phi1) + spec.jn(np.arange(1, harms,2), phi2)
                self.J[0,0] = 1
                self.J_tile = np.tile(self.J,(chans,1))
                self.J_inv = np.linalg.pinv(self.J_tile)
                self.S = np.dot(self.J_inv, self.h)
                tmp = np.dot(self.J_tile, self.S)
                error_matrix[i,j] = np.sum((tmp - self.h)**2)
            print phi1, np.min(error_matrix[i,:])
        fig, ax = pt.subplots()
        self.im_tmp = ax.imshow(error_matrix,extent = [np.rad2deg(min_rad), np.rad2deg(max_rad),np.rad2deg(min_rad), np.rad2deg(max_rad)])
        fig.canvas.draw(); fig.show()

    def get_C_Q_initial_vals(self,):
        pass
def compare_linear_sinusoidal(sinusoidal_shot, sinusoidal_freq, sinusoidal_mult, linear_shot):
    inter_obj = interf_demod(sinusoidal_shot, sinusoidal_freq, sinusoidal_mult)
    inter_obj.extract_mdsplus_data()
    #inter_obj.fft_signal()
    inter_obj.calc_epsilon()
    inter_obj.calc_carriers_time(force_symmetric = 1, individual_eps_calc = 1)
    inter_obj.calc_carriers_time_rfft()
    inter_obj.calc_phi1_means_std(3, 1, plot_fig = False)

    clim = [0,8]
    inter_obj.calculate_phase(np.deg2rad(125),1,2, clim=clim, show_fig = False)
    inter_obj.calculate_phase_simul(np.deg2rad(125), clim=clim, show_fig = False)

    #inter_obj.extract_harmonic_amps(start_point = 0, end_point = None, show_fig = False)
    #inter_obj.calculate_phase_phi_simul(np.deg2rad(125), use_rfft = 0, clim = [0,5], show_fig = show_phi_min)
    #inter_obj.compare_harmonics(show_fig = show_compare_harms)

    import MDSplus as MDS
    import matplotlib.pyplot as pt
    fig,ax = pt.subplots()
    fig2,ax2 = pt.subplots()
    #for ch in range(1,22):
    colors = ['b','k','r','y','m']
    colors.extend(colors)
    colors.extend(colors)

    fig_big, ax_big = pt.subplots(nrows=5, ncols = 5, sharex = True, sharey = True); ax_big = ax_big.flatten()

    #for i, ch in enumerate(range(12,22,1)):
    for i, ch in enumerate(range(21)):
        if ch<=10:
            fenton_ch = 11 - ch 
        else:
            fenton_ch = ch + 1
        T = MDS.Tree('h1data',linear_shot)
        n = T.getNode('.electr_dens.ne_het:ne_{}'.format(fenton_ch))
        try:
            ax_big[i].plot(inter_obj.t, np.abs(inter_obj.output[ch,:]),'--')
            tmp_data = n.raw_of().data()
            t2 = np.arange(-10000, len(tmp_data) - 10000) * 1./1000000.
            ax.plot(t2, tmp_data)
            ax_big[i].plot(t2, tmp_data)
            #if ch==1:
            print i, ch, fenton_ch
            if True:
                ax2.plot(t2, n.raw_of().data(),'{}-'.format(colors[i]))
                ax2.plot(inter_obj.t, np.abs(inter_obj.output[ch,:]),'--{}'.format(colors[i]))
                #ax2.plot(inter_obj.output_phase_simul[10,6000::],'--')
        except Exception,e:
            print 'problem',ch
            print e
    ax_big[-1].set_ylim([0,2.*np.pi])
    ax_big[-1].set_xlim([0,0.1])
    ax.set_ylim([0,2.*np.pi])
    ax2.set_ylim([0,2.*np.pi])
    fig.canvas.draw(); fig.show()
    fig2.canvas.draw(); fig2.show()

    fig_big.subplots_adjust(hspace=0.015, wspace=0.015,left=0.10, bottom=0.10,top=0.95, right=0.95)
    fig_big.canvas.draw(); fig.show()

def compare_sinusoidal_sinusoidal(sinusoidal1_shot, sinusoidal1_freq, sinusoidal1_mult, sinusoidal2_shot, sinusoidal2_freq, sinusoidal2_mult):
    inter_obj = interf_demod(sinusoidal1_shot, sinusoidal1_freq, sinusoidal1_mult)
    inter_obj.extract_mdsplus_data()
    #inter_obj.fft_signal()
    #inter_obj.calc_epsilon()
    inter_obj.calc_carriers_time(force_symmetric = 1, individual_eps_calc = 1)
    inter_obj.calc_carriers_time_rfft()

    inter_obj.calc_phi1_means_std(3, 1, plot_fig = False)

    clim = [0,8]
    inter_obj.calculate_phase(np.deg2rad(125),1,2, clim=clim, show_fig = False)
    inter_obj.calculate_phase_simul(np.deg2rad(125), clim=clim, show_fig = False)

    inter_obj2 = interf_demod(sinusoidal2_shot, sinusoidal2_freq, sinusoidal2_mult)
    inter_obj2.extract_mdsplus_data()
    #inter_obj2.fft_signal()
    #inter_obj2.calc_epsilon()
    inter_obj2.calc_carriers_time(force_symmetric = 1, individual_eps_calc = 1)
    inter_obj2.calc_carriers_time_rfft()
    inter_obj2.calc_phi1_means_std(3, 1, plot_fig = False)

    clim = [0,8]
    inter_obj2.calculate_phase(np.deg2rad(125),1,2, clim=clim, show_fig = False)
    inter_obj2.calculate_phase_simul(np.deg2rad(125), clim=clim, show_fig = False)

    #inter_obj.extract_harmonic_amps(start_point = 0, end_point = None, show_fig = False)
    #inter_obj.calculate_phase_phi_simul(np.deg2rad(125), use_rfft = 0, clim = [0,5], show_fig = show_phi_min)
    #inter_obj.compare_harmonics(show_fig = show_compare_harms)

    fig_big, ax_big = pt.subplots(nrows=5, ncols = 5, sharex = True, sharey = True); ax_big = ax_big.flatten()

    #for i, ch in enumerate(range(12,22,1)):
    for i, ch in enumerate(range(21)):
        try:
            ax_big[i].plot(inter_obj.t, np.abs(inter_obj.output[ch,:]),'-')
            ax_big[i].plot(inter_obj2.t, np.abs(inter_obj2.output[ch,:]),'--')
        except Exception,e:
            print 'problem',ch
            print e
    ax_big[-1].set_ylim([0,2.*np.pi])
    ax_big[-1].set_xlim([0,0.1])
    fig_big.subplots_adjust(hspace=0.015, wspace=0.015,left=0.10, bottom=0.10,top=0.95, right=0.95)
    fig_big.canvas.draw(); fig_big.show()
