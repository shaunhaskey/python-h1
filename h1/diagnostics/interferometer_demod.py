import scipy.special as spec
import MDSplus as MDS
import numpy as np
import matplotlib.pyplot as pt
import scipy.signal as sig
import scipy.optimize as opt


class interf_demod(object):
    def __init__(self,shot_number, f_m, dig_mult, start_samples = 10000, force_symmetric = 1, individual_eps_calc = 1, eps_method = 'using_1st_carrier', tophat_width_prop = 1, chan_list_dig = None, ch_z = None):
        '''Class for demodulating the sinusoidally modulated interferometer
        SRH : 5Nov2013
        '''
        self.shot_number = shot_number
        self.f_m = f_m
        self.dig_mult = dig_mult
        self.start_samples = start_samples
        print self.start_samples,
        self.start_samples = self.start_samples - self.start_samples%(dig_mult)
        print self.start_samples
        self.f_s = self.f_m * self.dig_mult

        tr = MDS.Tree('electr_dens',shot_number)
        if chan_list_dig==None:
            self.chan_list_dig = tr.getNode('\electr_dens::top.ne_het:chan_list').data().tolist()
        else:
            self.chan_list_dig = chan_list_dig
        if chan_list_dig==None:
            self.ch_z = tr.getNode('\electr_dens::top.ne_het:chan_z').data()
        else:
            self.ch_z = ch_z
        if len(self.chan_list_dig)!=len(self.ch_z):
            raise Exception("chan_list_dig and ch_z are different lengths, each channel needs a z coordinate")
        top_down_order = (np.argsort(self.ch_z)[::-1]).tolist()
        self.chan_list_dig_top_down = np.array(self.chan_list_dig)[top_down_order]

        self.get_up_to_carriers(force_symmetric = force_symmetric, individual_eps_calc = individual_eps_calc, eps_method = eps_method, tophat_width_prop = tophat_width_prop)

    def get_up_to_carriers(self,force_symmetric = 1, individual_eps_calc = 1, eps_method = 'using_1st_carrier', tophat_width_prop = 1):
        '''Convenience function to do several things at once - it gets you the carriers as a function of time
        SRH: 7Nov2013
        '''
        self.extract_mdsplus_data()
        self.calc_carriers_time(force_symmetric = force_symmetric, eps_method = eps_method, tophat_width_prop = tophat_width_prop)
        self.calc_carriers_time_rfft(force_symmetric = force_symmetric, individual_eps_calc = individual_eps_calc, tophat_width_prop = tophat_width_prop)

    def extract_mdsplus_data(self,):
        '''Get the data out of MDSplus
        SRH: 5Nov2013
        '''
        T = MDS.Tree('h1data',self.shot_number)
        count = 0
        #self.chan_list = range(21)
        self.rms_pre_shot = []; self.dc_pre_shot = []
        self.rms_full_shot = []; self.dc_full_shot = []
        for i, chan_loc in enumerate(self.chan_list_dig_top_down):
            n = T.getNode('.electr_dens.camac.{}'.format(chan_loc))
            #digi, chan = self.chan_num_to_digi(chan_num)
            #n = T.getNode('.electr_dens.camac.A14_{}.input_{}'.format(digi,chan))
            signal = n.data()
            self.orig_signal_length = len(signal)
            if count == 0:
                self.keep_length = self.orig_signal_length - (self.orig_signal_length%self.dig_mult)
                self.raw_signals = np.zeros((len(self.chan_list_dig_top_down), self.keep_length),dtype = float)
            self.raw_signals[i,:] = +signal[:self.keep_length]
            #signal = signal[:self.keep_length]
            dc, ac_rms = self.calculate_ac_dc(signal, full = 0)
            self.rms_pre_shot.append(ac_rms)
            self.dc_pre_shot.append(dc)
            dc, ac_rms = self.calculate_ac_dc(signal, full = 1)
            self.rms_full_shot.append(ac_rms)
            self.dc_full_shot.append(dc)
            #self.rms_pre_shot.append(np.sqrt(1./self.start_samples*np.sum((signal[:self.start_samples]-np.mean(signal[:self.start_samples]))**2)))
            count += 1
        self.t = np.arange(-self.start_samples, self.keep_length - self.start_samples) * 1./self.f_s

    def extract_epsilon_sync_signal(self,plot=False):
        '''Get epsilon from the digitised fm modulation signal 
        '''
        T = MDS.Tree('h1data',self.shot_number)
        n = T.getNode('.electr_dens.camac.A14_24.input_5')
        #-ve is because of a polarity problem on the BNC-> lemo connector
        sync_signal = -n.data()[:self.raw_signals.shape[1]]
        sync_signal_fft = np.fft.fft(sync_signal)
        fft_ax = np.fft.fftfreq(len(sync_signal_fft),1./self.f_s)
        f_m_loc = np.argmin(np.abs(fft_ax - self.f_m))
        #+np.pi/2 because the modulation is represented by a sine wave, not a cos, 
        #+np.pi,  mod 2pi,  - np.pi, to be within -np.pi, and np.pi
        self.sync_epsilon = (np.angle(sync_signal_fft[f_m_loc])+np.pi/2+np.pi)%(2.*np.pi) - np.pi
        if plot:
            fig, ax = pt.subplots()
            ax.plot(sync_signal)
            fig.canvas.draw(); fig.show()
        return self.sync_epsilon

    def epsilon_analysis(self,ch = 7, use_full = 0):
        self.extract_epsilon_sync_signal()
        if use_full:
            signal = self.signal_fft
            fm_loc = self.f_m_loc
        else:
            signal = self.signal_start_fft
            fft_ax_tmp = np.fft.fftfreq(len(signal),1./self.f_s)
            fm_loc = np.argmin(np.abs(fft_ax_tmp - self.f_m))
        self.e_recon = np.angle(signal[:,fm_loc]/1j)
        self.e_reconb = np.angle(-signal[:,fm_loc]/1j)
        self.e_overall = self.e_reconb * 0
        self.e_overall2 = np.angle(signal[:,fm_loc]/signal[:,3*fm_loc])/(-2.)

        self.e_overall[self.e_recon>0] = self.e_recon[self.e_recon>0]
        self.e_overall[self.e_reconb>0] = self.e_reconb[self.e_reconb>0]
        self.e_recon2 = np.angle(signal[:,fm_loc*2])/2
        self.e_recon3 = np.angle(signal[:,fm_loc*3]/1j)/3
        self.e_recon2b = np.angle(-signal[:,fm_loc*2])/2
        self.e_recon3b = np.angle(-signal[:,fm_loc*3]/1j)/3
        print('sync:{:.2f}, -1sync:{:.2f}, overall:{:.2f},e_recon1:{:.2f}, e_recon1b:{:.2f}, e_recon2:{:.2f}, e_recon2b:{:.2f}'.format(self.sync_epsilon, self.sync_epsilon - np.pi, self.e_overall[ch], self.e_recon[ch], self.e_reconb[ch],self.e_recon2[ch], self.e_recon2b[ch]))
        print('{},{}'.format((self.sync_epsilon - self.e_recon[ch])/np.pi, (self.sync_epsilon - self.e_reconb[ch])/np.pi))

    def plot_raw_signals(self, subtract_dc = 1):
        '''Plot the raw digitised signals
        SRH: 5Nov2013
        '''
        n_chans = self.raw_signals.shape[0]
        fig, ax = pt.subplots(nrows = 5, ncols = 5,sharex =True, sharey=False); ax = ax.flatten()
        for i in range(self.raw_signals.shape[0]):
            if subtract_dc:
                ax[i].plot(self.t, self.raw_signals[i,:] - self.dc_pre_shot[i])
            else:
                ax[i].plot(self.t, self.raw_signals[i,:])
        fig.subplots_adjust(hspace=0, wspace=0,left=0.05, bottom=0.05,top=0.95, right=0.95)
        ax[0].set_ylim([-0.4,0.4])
        ax[0].set_xlim([-0.005,0.005])
        ax[-1].set_xlabel('Sample Number')
        fig.canvas.draw(); fig.show()

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
        self.e_reconb = np.angle(-self.signal_fft[:,self.f_m_loc]/1j)
        self.e_overall = self.e_reconb * 0
        self.e_overall2 = np.angle(self.signal_fft[:,self.f_m_loc]/self.signal_fft[:,3*self.f_m_loc])/(-2.)

        self.e_overall[self.e_recon>0] = self.e_recon[self.e_recon>0]
        self.e_overall[self.e_reconb>0] = self.e_reconb[self.e_reconb>0]
        self.e_recon2 = np.angle(self.signal_fft[:,self.f_m_loc*2])
        self.e_recon3 = np.angle(self.signal_fft[:,self.f_m_loc*3]/1j)
        self.e_recon2b = np.angle(-self.signal_fft[:,self.f_m_loc*2])
        self.e_recon3b = np.angle(-self.signal_fft[:,self.f_m_loc*3]/1j)


    def get_e_recon_matrix(self,n_carriers):
        '''Calculate the epsilon values for each of the carriers
        SRH: 11Nov2013
        '''
        self.e_recon_matrix = np.zeros((self.signal_fft.shape[0], n_carriers),dtype = float)
        self.e_recon_matrix_start = np.zeros((self.signal_fft.shape[0], n_carriers),dtype = float)
        for n in range(1,n_carriers):
            mult = 1
            if (n%2)!=0: mult = 1j
            self.e_recon_matrix[:,n] = np.angle(self.signal_fft[:,n*self.f_m_loc]/mult)/n
            self.e_recon_matrix_start[:,n] = np.angle(self.signal_start_fft[:,n*self.f_m_loc_start]/mult)/n
        harmonic_strength = np.abs(self.signal_fft[:,self.f_m_loc])
        harmonic_strength_sort = np.argsort(harmonic_strength)
        important_harms = self.e_recon_matrix[harmonic_strength_sort[-3:],1]
        important_harms[important_harms<0]+=np.pi
        print '########################'
        print important_harms, np.mean(important_harms), np.std(important_harms)

        harmonic_strength = np.abs(self.signal_start_fft[:,self.f_m_loc_start])
        harmonic_strength_sort = np.argsort(harmonic_strength)
        important_harms = self.e_recon_matrix_start[harmonic_strength_sort[-3:],1]
        important_harms[important_harms<0]+=np.pi
        print('##### Epsilon ########')
        print('3 values {}, mean:{:.3f}, std:{:.3f}, std/mean:{:.3f}%'.format(important_harms, np.mean(important_harms), np.std(important_harms),np.std(important_harms)/np.mean(important_harms)*100))
        self.e_ave = np.mean(important_harms)

    def calc_carriers_time(self, n_carriers = 6, force_symmetric = 0, tophat_width_prop = 1, eps_method = 'using_1st_carrier'):
        '''Calculate the carriers as a function of time
        SRH: 5Nov2013
        '''
        #FFT's
        self.signal_fft = np.fft.fft(self.raw_signals)
        self.signal_start_fft = np.fft.fft(self.raw_signals[:,:self.start_samples])

        #FFT frequencies
        self.fft_ax = np.fft.fftfreq(self.signal_fft.shape[1], 1./self.f_s)
        self.fft_ax_start = np.fft.fftfreq(self.signal_start_fft.shape[1], 1./self.f_s)

        self.carriers = np.zeros((n_carriers,self.signal_fft.shape[0], self.signal_fft.shape[1]),dtype = complex)

        # Find the modulation frequency
        self.f_m_loc = self.signal_fft.shape[1]/self.dig_mult
        self.f_m_loc_start = self.signal_start_fft.shape[1]/self.dig_mult
        print('check f_m loc is right :', (float(self.signal_fft.shape[1])/self.dig_mult) == self.f_m_loc, (float(self.signal_start_fft.shape[1])/self.dig_mult) == self.f_m_loc_start)
        #np.argmin(np.abs(self.fft_ax - self.f_m))
        
        # Get epsilon for the first harmonic
        self.get_e_recon_matrix(n_carriers)
        if eps_method == 'best_channels_fixed':
            print('using best_channels_fixed')
            self.use_eps = self.e_recon_matrix_start * 0 + self.e_ave
        elif eps_method == 'individual':
            print('using individual')
            self.use_eps = +self.e_recon_matrix_start
            self.use_eps = +self.e_recon_matrix
        elif eps_method == 'using_1st_carrier':
            print('using 1st carrier')
            tmp = self.e_recon_matrix_start[:,1]
            self.use_eps = self.e_recon_matrix_start * 0 + tmp[:,np.newaxis]

        #self.e_recon = np.angle(self.signal_fft[:,self.f_m_loc]/1j)

        # self.e_recon = np.angle(self.signal_fft[:,self.f_m_loc]/self.signal_fft[:,3*self.f_m_loc])/(-2.)
        # enforce the use of a single epsilon...
        
        #self.e_recon = self.e_recon *0 + self.e_recon[10]

        #Width of the tophat and make it an odd number
        self.width = int(self.f_m_loc * tophat_width_prop)
        if (self.width%2)!=0: self.width = self.width - 1

        for n in range(n_carriers):
            mult = 1
            if (n%2)!=0: mult = -1j

            #n_freq_loc = np.argmin(np.abs(self.fft_ax - n*self.f_m))
            n_freq_loc = self.f_m_loc * n
            fft_shifted = self.signal_fft *0

            #Shift the frequencies
            if n==0:
                fft_shifted[:,0:self.width/2] = self.signal_fft[:, n_freq_loc:n_freq_loc+self.width/2]
                fft_shifted[:,-self.width/2+1:] = self.signal_fft[:,-self.width/2+1:]
            else:
                fft_shifted[:,0:self.width/2] = self.signal_fft[:,n_freq_loc:n_freq_loc+self.width/2]
                fft_shifted[:,-self.width/2+1:] = self.signal_fft[:,n_freq_loc-self.width/2+1:n_freq_loc]

            #Apply the timing correction
            #Need to check the best epsilon for each case.... 
            e_recon_tmp = self.use_eps[:,n]
            # if n==0:
            #     e_recon_tmp = self.e_recon
            # elif (n%2)==0:
            #     if individual_eps_calc:
            #         e_recon_tmp = np.angle(self.signal_fft[:,self.f_m_loc*n])/n
            #     else:
            #         e_recon_tmp = +self.e_recon
            # else:
            #     if individual_eps_calc:
            #         e_recon_tmp = np.angle(self.signal_fft[:,self.f_m_loc*n]/-1j)/n
            #         print 'using individual carrier epsilon calc'
            #     else:
            #         print 'using same epsilon for each carrier'
            #         e_recon_tmp = +self.e_recon
            self.carriers_epsilon_used = +e_recon_tmp
            fft_shifted_corrected = fft_shifted * np.exp(-1j*n*e_recon_tmp[:,np.newaxis])*mult
            if force_symmetric:
                fft_shifted_corrected[:, -self.width/2+1:] = np.conj(fft_shifted_corrected[:, 1:self.width/2])[:,::-1]
            self.carriers[n,:,:] = np.fft.ifft(fft_shifted_corrected)

    def access_carriers(self,):
        fig, ax = pt.subplots(nrows = 5, ncols = 5, sharex = True, sharey = True); ax = ax.flatten()
        for i in range(self.carriers.shape[1]):
            real_part = np.sum(np.abs(np.real(self.carriers[:,i,:])), axis = 1)
            imag_part = np.sum(np.abs(np.imag(self.carriers[:,i,:])), axis = 1)
            ax[i].bar(range(self.carriers.shape[0]), imag_part*0+1, color='b')
            ax[i].bar(range(self.carriers.shape[0]), real_part/(real_part + imag_part), color='r')
        ax[0].set_xlim([0,self.carriers.shape[0]])
        ax[0].set_ylim([0,1])
        fig.canvas.draw(); fig.show()

    def compare_all_epsilons(self,):
        fig, ax = pt.subplots(nrows = 5, ncols = 5, sharex = True, sharey = True); ax = ax.flatten()
        for i in range(self.carriers.shape[1]):
            real_part = np.sum(np.abs(np.real(self.carriers[:,i,:])), axis = 1)
            imag_part = np.sum(np.abs(np.imag(self.carriers[:,i,:])), axis = 1)
            ax[i].bar(range(self.carriers.shape[0]), imag_part*0+1, color='b')
            ax[i].bar(range(self.carriers.shape[0]), real_part/(real_part + imag_part), color='r')
        ax[0].set_xlim([0,self.carriers.shape[0]])
        ax[0].set_ylim([0,1])
        fig.canvas.draw(); fig.show()


    def calc_carriers_time_rfft(self, n_carriers = 6, force_symmetric = 0, individual_eps_calc = 0, tophat_width_prop = 1):
        '''Same as calc_carriers_time except using the real fft and
        real ifft which forces symmetry, and removes some of the
        problems with getting the frequency shifting wrong 
        SRH:5Nov2013
        '''
        self.e_recon = np.angle(self.signal_fft[:,self.f_m_loc]/1j)
        #self.e_recon = np.angle(self.signal_fft[:,self.f_m_loc]/self.signal_fft[:,3*self.f_m_loc])/(-2.)
        self.e_recon = self.e_recon *0 + self.e_recon[10]
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
                e_recon_tmp = +self.e_recon
            elif (n%2)==0:
                if individual_eps_calc:
                    e_recon_tmp = np.angle(self.signal_fft[:,self.f_m_loc*n])/n
                else:
                    e_recon_tmp = +self.e_recon
            else:
                if individual_eps_calc:
                    e_recon_tmp = np.angle(self.signal_fft[:,self.f_m_loc*n]/-1j)/n
                    print 'using individual carrier epsilon calc'
                else:
                    print 'using same epsilon for each carrier'
                    e_recon_tmp = +self.e_recon
            self.carriers_rfft_epsilon_used = +e_recon_tmp
            #if n==0:
            #    e_recon_tmp = self.e_recon
            #elif (n%2)==0:
            #    e_recon_tmp = np.angle(self.signal_fft[:,self.f_m_loc*n])/n
            #else:
            #    e_recon_tmp = np.angle(self.signal_fft[:,self.f_m_loc*n]/-1j)/n
            a = tmp* np.exp(-1j*n*self.e_recon[:,np.newaxis])*mult
            #h_n.append(np.fft.ifft(a))
            self.carriers_rfft[n,:,:] = np.fft.irfft(a)

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
        end_point -= (start_point - end_point)%self.dig_mult
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
        if phi1 == None: 
            phi1 = np.array(self.phi1_mean)
            title_text = 'Using individual modulation depths for each channel'
        else:
            phi1 = np.array([phi1] * chans)
            title_text = 'Using fixed modulation depth for all channels'
        bessel_amps = []
        for i in phi1:
            bessel_amps.append(spec.jn(range(max_harmonic), i))
        bessel_amps = np.array(bessel_amps)
        print bessel_amps.shape, self.harmonics.shape
        
        #eps_to_use = np.angle(self.harmonics[:,1]/1j)[8]
        eps_to_use = self.e_ave
        if show_fig:
            fig, ax = pt.subplots(nrows=5, ncols = 5, sharex = True, sharey = True); ax = ax.flatten()
            for i in range(chans):
                strengths = self.harmonics[0,:]*0.

                Q_0 = self.harmonics[i,1]/bessel_amps[i,1]/1j/np.exp(1j*eps_to_use)
                C_0 = self.harmonics[i,2]/bessel_amps[i,2]/np.exp(1j*2*eps_to_use)
                Q_0_prop_real = np.abs(np.real(Q_0)/np.abs(Q_0))*100
                C_0_prop_real = np.abs(np.real(C_0)/np.abs(C_0))*100
                print('ch:{}, Q_O_prop:{:.2f}%, C_0_prop:{:.2f}%, Initial Phase:{:.2f}deg'.format(i, Q_0_prop_real, C_0_prop_real, np.rad2deg(np.angle(1j*np.real(Q_0) + np.real(C_0)))))
                #strengths[1::2] = self.harmonics[i,1::2]*bessel_amps[i,1]/ self.harmonics[i,1]
                n_vals = np.arange(self.harmonics.shape[1])
                strengths[1::2] = self.harmonics[i,1::2]/np.exp(1j*n_vals[1::2]*eps_to_use)/np.real(Q_0)/1j
                #strengths[2::2] = self.harmonics[i,2::2]*bessel_amps[i,2]/ self.harmonics[i,2]
                strengths[2::2] = self.harmonics[i,2::2]/np.exp(1j*n_vals[2::2]*eps_to_use)/np.real(C_0)
                strengths = 1*strengths[1:]
                print np.real(strengths)/np.abs(strengths)*100
                #ax[i].plot(range(1,max_harmonic),np.abs(self.harmonics[i,1:])*np.abs(bessel_amps[i,1])/ np.abs(self.harmonics[i,1]),'o-')
                ax[i].plot(range(1,max_harmonic),strengths,'o-')
                #ax[i].plot(range(1,max_harmonic),np.sign(np.real(self.harmonics[i,1:]))*np.abs(self.harmonics[i,1:])*np.abs(bessel_amps[i,1])/ np.abs(self.harmonics[i,1]),'o-')
                ax[i].plot(range(1,max_harmonic), bessel_amps[i,1:],'x-')
            ax[-1].set_xlim([0,max_harmonic])
            ax[-1].set_ylim([0,np.max(bessel_amps)])
        
            fig.subplots_adjust(hspace=0.015, wspace=0.015,left=0.10, bottom=0.10,top=0.95, right=0.95)
            
            fig.suptitle('predicted harmonic amps from bessel funcs (x), measured (o) - dividing by Q i exp(i n eps) for odd and C exp(i n eps)for even\nNote first 2 harmonics must agree perfectly because they are used to calc Q(0) and C(0). {}'.format(title_text))
            fig.canvas.draw(); fig.show()

    def calculate_phase(self, phi1, odd, even, use_rfft = 0, clim = [0,5], show_fig = True, save_fig_name = None, times = None):
        '''Calculate the actual phase shift caused by the plasma. Can use the rfft
        SRH: 5Nov2013
        '''
        if use_rfft: 
            carriers = self.carriers_rfft
            epsilon_used = +self.carriers_rfft_epsilon_used
        else:
            carriers = self.carriers
            epsilon_used = +self.carriers_epsilon_used
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
            self.im_ax[0].set_xlabel('time (s)')
            self.im_ax[0].set_ylabel('Interferometer Channel')
            fig, ax = pt.subplots(ncols = 2); 
            for i in range(self.output.shape[0]):
                ax[0].plot(self.t, self.output[i,:])
            ax[0].set_ylim([-2.*np.pi, 2.*np.pi])
            ax[1].plot(epsilon_used, np.mean(self.output,axis = -1), 'o')
            for i in range(len(epsilon_used)):
                ax[1].text(epsilon_used[i], np.mean(self.output,axis = -1)[i], str(i))
            fig.canvas.draw(); fig.show()
        
        if times==None: times = [0.02,0.04,0.063]
        colours = ['b','g','k']
        for col, t in zip(colours, times):
            t_loc = np.argmin(np.abs(self.t - t))
            plot_vals = np.abs(self.output[:,t_loc])<10
            if show_fig:
                self.im_ax[1].plot(np.abs(self.output[:,t_loc][plot_vals]),np.arange(self.output.shape[0])[plot_vals],'{}x-'.format(col))
                self.im_ax[0].vlines(t, self.output.shape[0],0, colors = col)
        if show_fig:
            self.im_ax[1].grid()
            self.im_ax[1].set_xlim(clim)
            #self.im =self.im_ax.imshow(np.abs(np.flipud(self.output)), aspect='auto',cmap = 'hot', extent=[self.t[0],self.t[-1],self.output.shape[0],0])
            self.im.set_clim(clim)
            self.im_ax[0].set_xlim([-0.005,0.12])
            self.im_ax[1].set_ylim([0.,self.output.shape[0]])
            cbar = pt.colorbar(self.im,ax = self.im_ax[0])
            cbar.set_label('plasma phase shift (rad)')
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
            epsilon_used = +self.carriers_rfft_epsilon_used
        else:
            carriers = self.carriers
            epsilon_used = +self.carriers_epsilon_used

        harms, chans, Ns = carriers.shape
        self.J = np.zeros((harms, 3), dtype = float)
        self.J[::2,2] = spec.jn(np.arange(0, harms,2), phi1)
        self.J[1::2,1] = spec.jn(np.arange(1, harms,2), phi1)
        self.J[0,0] = 1
        self.J_inv = np.linalg.pinv(self.J)
        self.h = np.zeros((harms, Ns),dtype = float)
        self.result_list = []
        if show_fig:
            fig, ax = pt.subplots(ncols = 2); 
        self.output_phase_simul = np.zeros((chans, Ns), dtype = float)
        for i in range(chans):
            tmp = np.dot(self.J_inv, np.real(carriers[:,i,:]))
            self.result_list.append(tmp)
            output = np.unwrap(np.arctan2(tmp[2,:], tmp[1,:]))
            initial_phase = np.mean(output[30:1000])
            #print initial_phase,
            self.output_phase_simul[i,:] = (output - initial_phase)
            tmp_mean = np.mean(self.output_phase_simul[i,:])
            #print np.cos(tmp_mean)*np.sin(tmp_mean), tmp_mean, np.sign(np.cos(tmp_mean)*np.sin(tmp_mean)) == np.sign(tmp_mean)
            
            if show_fig:
                ax[0].plot(self.t, self.output_phase_simul[i,:])

        
        print ''

        if show_fig:
            ax[0].set_ylim([-2.*np.pi, 2.*np.pi])
            ax[1].plot(epsilon_used, np.mean(self.output_phase_simul,axis = -1), 'o')
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


    def calculate_phase_phi_simul(self, use_rfft = 0, clim = [0,5], show_fig = True, channel_mask = None, max_harm = None):
        '''Calculate phi1 (modulation depth) simulatneously using all interferometer channels

        max_harm is the maximum harmonic that is used for the phi simultaneous fit
        SRH: 5Nov2013
        '''
        if use_rfft: 
            carriers = self.carriers_rfft
        else:
            carriers = self.carriers
        n_samples = 2000

        harms, chans, Ns = carriers.shape
        if max_harm!=None: harms = max_harm
        self.J = np.zeros((harms, 3), dtype = float)

        #self.h = np.zeros((harms*chans, Ns),dtype = float)
        self.h = np.zeros((harms*chans, n_samples),dtype = float)
        for i in range(chans):
            self.h[i*harms:(i+1)*harms,:] = np.real(carriers[:max_harm,i,:n_samples])
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
            #print phi1, error_list[-1]
        self.phi1_min_simul = phi1_list[np.argmin(error_list)]
        self.J[::2,2] = spec.jn(np.arange(0, harms,2), self.phi1_min_simul)
        self.J[1::2,1] = spec.jn(np.arange(1, harms,2), self.phi1_min_simul)
        self.J[0, 0] = 1
        self.J_tile = np.tile(self.J,(chans,1))
        self.J_inv = np.linalg.pinv(self.J_tile)
        self.S = np.dot(self.J_inv, self.h)
        tmp = np.dot(self.J_tile, self.S)
        
        a = opt.fmin(self.min_func_phi1_simul, np.deg2rad(135), args=(harms, chans, self.h),disp=1, xtol=0.00001)
        a2 = opt.minimize(self.min_func_phi1_simul, np.deg2rad(135), args=(harms, chans, self.h),method = 'L-BFGS-B', bounds = [(0.05,1.3*np.pi)])
        a3 = opt.minimize(self.phase_phi1_phi2_error, [np.deg2rad(135), np.deg2rad(140),0.5,0.5], args=(harms, chans, self.h),method = 'L-BFGS-B', bounds = [(0.05,1.3*np.pi), (0.05,1.3*np.pi),(0.001,20),(0.001,20)])
        print '######## fmin min_func_phi1_simul ########'
        print a, self.phi1_min_simul
        print '######## minimize min_func_phi1_simul #########'
        print a2, self.phi1_min_simul
        self.phi1_min_simul = a2['x']
        print '######## minimize min_func_phi1_phi2_error ########'
        print a3, self.phi1_min_simul
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
    inter_obj.calc_carriers_time(force_symmetric = 1, eps_method = 'using_1st_carrier')
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
    inter_obj.calc_carriers_time(force_symmetric = 1, eps_method = 'using_1st_carrier')
    inter_obj.calc_carriers_time_rfft()

    inter_obj.calc_phi1_means_std(3, 1, plot_fig = False)

    clim = [0,8]
    inter_obj.calculate_phase(np.deg2rad(125),1,2, clim=clim, show_fig = False)
    inter_obj.calculate_phase_simul(np.deg2rad(125), clim=clim, show_fig = False)

    inter_obj2 = interf_demod(sinusoidal2_shot, sinusoidal2_freq, sinusoidal2_mult)
    inter_obj2.extract_mdsplus_data()
    #inter_obj2.fft_signal()
    #inter_obj2.calc_epsilon()
    inter_obj2.calc_carriers_time(force_symmetric = 1, eps_method = 'using_1st_carrier')
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


def plot_from_MDSplus(shot):
    fig,ax = pt.subplots()
    #fig2,ax2 = pt.subplots()
    #for ch in range(1,22):
    #colors = ['b','k','r','y','m']
    #colors.extend(colors)
    #colors.extend(colors)

    fig_big, ax_big = pt.subplots(nrows=5, ncols = 5, sharex = True, sharey = True); ax_big = ax_big.flatten()
    tr = MDS.Tree('electr_dens',shot)
    chan_list_dig = tr.getNode('\electr_dens::top.ne_het:chan_list').data().tolist()
    ch_z = tr.getNode('\electr_dens::top.ne_het:chan_z').data()
    top_down_order = (np.argsort(ch_z)[::-1]).tolist()
    #for i, ch in enumerate(range(12,22,1)):
    for i, ch in enumerate(top_down_order):
        print i,ch
        n = tr.getNode('\electr_dens::top.ne_het:ne_{}'.format(ch+1))
        sig = n.data()
        if i==0:
            output_data = np.zeros((len(top_down_order),len(sig)),dtype=float)
            t = n.dim_of().data()
            #t = np.arange(len(sig))
        output_data[i,:] = +sig
    im =ax.imshow(output_data,interpolation = 'nearest', aspect='auto',cmap = 'hot', extent=[t[0],t[-1],output_data.shape[0],0], origin = 'upper')
    im.set_clim([0,2.5])
    fig.canvas.draw(); fig.show()
    #fig_big.subplots_adjust(hspace=0.015, wspace=0.015,left=0.10, bottom=0.10,top=0.95, right=0.95)
    #fig_big.canvas.draw(); fig.show()


def put_data_MDSplus(shot, f_m=None, dig_mult=None):
    count = 0
    plots=0  # number of plots shown at beginning
    print('hello world')
    print f_m, dig_mult
    try:
            
        tr = MDS.Tree('electr_dens',shot)
        chnd = tr.getNode('\electr_dens::top.ne_het:chan_list')
        chan_list = chnd.data()
        ch_z_nd = tr.getNode('\electr_dens::top.ne_het:chan_z')
        ch_z = ch_z_nd.data()
        if dig_mult==None:
            dig_mult = tr.getNode('\electr_dens::top.ne_het:digi_mult').data()
            print('Got dig multiplier out of the tree : {}'.format(dig_mult))
        if f_m==None:
            f_m = tr.getNode('\electr_dens::top.ne_het:fm_freq').data()
            print('Got freq modulation out of the tree : {}'.format(f_m))
            
        #Get the data
        print shot, f_m, dig_mult
        inter_obj = interf_demod(shot, f_m, dig_mult, start_samples = 10000, force_symmetric = 1, individual_eps_calc = 0, eps_method = 'using_1st_carrier', tophat_width_prop = 1)

        clim = [0,5]
        inter_obj.calculate_phase_phi_simul(use_rfft = 0, clim = [0,5], show_fig = 0)
        inter_obj.calculate_phase(inter_obj.phi1_min_simul, 1, 2, clim=clim, show_fig = 0)
        #inter_obj.calculate_phase_simul(inter_obj.phi1_min_simul, clim=clim, show_fig = 0)
        chan_list = inter_obj.chan_list_dig
        chan_list_dig_top_down = inter_obj.chan_list_dig_top_down.tolist()
        print 'hello'
        for (n,ch) in enumerate(chan_list): 
            if ch[0:3]=='A14':
                digitiser_A14 = True
            else:
                digitiser_A14 = False
            print 'hello', n
            nd = tr.getNode('\electr_dens::top.camac:'+ch);
            old_dim_of = nd.dim_of()
            phs = np.zeros((inter_obj.orig_signal_length,), dtype=float)
            #sig=nd.data() 
            #t_raw=nd.dim_of().data()
            #(sgn, bw) = (1, .8)
            #acrms = np.sqrt(sig.var())
            #if acrms<0.02: 
            #    continue
            #elif acrms<0.03: 
            #    bw = .1
            #phs = sgn*extract_phase(None, sig=sig,t_raw=t_raw,
            #                        f=5e3,bw=bw,plot=0)
            index = chan_list_dig_top_down.index(ch)
            #print n, ch, index, chan_list_dig_top_down[index]
            
            phs[:inter_obj.keep_length] = inter_obj.output[index,:]
            if np.mean(phs)<0: 
                print 'correcting sign'
                phs = -1.*phs
            #dim = Dimension(Window(startIdx[chan], endIdx[chan], trigTime), Range(None, None, clockPeriod))
            startIdx=tr.getNode('\electr_dens::top.camac:{}:startidx'.format(ch)).data()
            if digitiser_A14:
                trigTime = tr.getNode('\electr_dens::top.camac:{}:stop_trig'.format(ch.split(':')[0])).data()
            else:
                trigTime = tr.getNode('\electr_dens::top.camac:{}:trigger'.format(ch.split(':')[0])).data()

            try:
                endIdx = tr.getNode('\electr_dens::top.camac:{}:endidx'.format(ch)).data()
            except MDS.TdiException, e:
                #print e
                pts = tr.getNode('\electr_dens::top.camac:{}:pts'.format(ch.split(':')[0])).data()
                endIdx = pts - startIdx - 1

            if digitiser_A14:
                Range = tr.getNode('\electr_dens::top.camac:{}:ext_clock_in'.format(ch.split(':')[0])).evaluate()
            else:
                Range = tr.getNode('\electr_dens::top.camac:{}:ext_clock'.format(ch.split(':')[0])).evaluate()

            win = MDS.Window(startIdx, endIdx, trigTime)
            dim = MDS.Dimension(win, Range)
            convExpr = MDS.Data.compile("0.35*$VALUE")
            convExpr.setUnits("1e18/m-3")
            pnode = tr.getNode('\electr_dens::top.ne_het:ne_'+str(n+1))
            rawMdsData = MDS.Float32Array(phs)
            rawMdsData.setUnits("rad")
            if digitiser_A14:
                signal = MDS.Signal(convExpr, rawMdsData, dim)
            else:
                signal = MDS.Signal(convExpr, rawMdsData, old_dim_of)
            pnode.putData(signal)
            if n==0:
                centrenode = tr.getNode('\electr_dens::top.ne_het:ne_centre')
                centrenode.putData(signal)
    except None:#Exception, reason:
        print('Exception on shot {s}, "{r}"'.format(s=shot,r=reason))




def error(test_vals, x, values, omega, mirror_depth):
    phi, a , b = test_vals
    art = mirror_depth * np.sin(omega*x + phi) + a*x + b
    return np.sqrt(np.sum((values - art)**2))/len(values)

def calc_art_big_fit(test_vals, t, values, x_locs, omega, mirror_depth):
    phi, k, a = test_vals
    term1 = t
    term2 = x_locs * float(k)
    art = mirror_depth*np.cos(omega*term1[np.newaxis,:] + term2[:,np.newaxis]+phi) + a * t
    return art

def error2(test_vals, t, values, x_locs, omega, mirror_depth):
    # phi, k = test_vals
    # term1 = t
    # term2 = x_locs * float(k)
    # art = np.sin(term1[np.newaxis,:] + term2[:,np.newaxis]+phi)
    art = calc_art_big_fit(test_vals, t, values, x_locs, omega, mirror_depth)
    return np.sqrt(np.sum((values - art)**2))#/np.sqrt(np.sum(values**2))


def plot_wobbly_mirror_results(big_fit_values, t, big_fit, art, mirror_depth, title, savefig_name):
    fig, ax = pt.subplots()
    cm_to_inch=0.393701
    fig.set_figwidth(8.48*cm_to_inch)
    fig.set_figheight(8.48*2*cm_to_inch)
    import matplotlib as mpl
    mpl.rcParams['font.size']=8.0
    mpl.rcParams['axes.titlesize']=8.0#'medium'                                                                                
    mpl.rcParams['xtick.labelsize']=8.0
    mpl.rcParams['ytick.labelsize']=8.0
    mpl.rcParams['lines.markersize']=5.0
    mpl.rcParams['savefig.dpi']=300


    for i, offset in enumerate(big_fit_values):
        for j in range(3):
            ax.plot(t + np.max(t)*j, big_fit[i,:]+offset*2*mirror_depth,'b')
            ax.plot(t + np.max(t)*j, art[i,:]+offset*2*mirror_depth,'g')
    ax.set_xlabel('time (3 periods have been pasted together)')
    ax.set_ylabel('Phase + i x 2.*np.pi where i is channel number')
    #ax.set_ylim([0,21*2.*np.pi])
    ax.hlines(np.arange(0,21)*2.*mirror_depth,ax.get_xlim()[0],ax.get_xlim()[1],linestyles='dashed')
    ax.set_xlim([0,np.max(t)*(j+1)])
    ax.set_title(title)
    fig.savefig(savefig_name)
    fig.canvas.draw(); fig.show()

def extract_useful_wobbly_mirror_data(inter_obj, start_loc, end_loc, window_length = 20, plot=False):
    if plot:
        fig, ax = pt.subplots(nrows = 2, sharey = True)
        fig2, ax2 = pt.subplots(nrows = 2, sharey = True)
    start_list = []
    output_new = inter_obj.output.copy()
    for i in range(inter_obj.output.shape[0]):
        #smooth the data
        tmp = np.convolve(inter_obj.output[i,:], np.ones(window_length)/float(window_length), mode= 'same')
        #remove dc offset by centering about the peaks
        #output_new[i,:] = tmp - (np.min(tmp) + np.max(tmp))/2
        output_new[i,:] = tmp - (np.mean(tmp))
        if (np.min(inter_obj.output[i,:]) + np.max(inter_obj.output[i,:]))<2.*np.pi:
            if plot:
                ax[0].plot(output_new[i,:])
                ax2[0].plot(output_new[i,start_loc:end_loc])
                ax2[1].plot(inter_obj.output[i,start_loc:end_loc])
        start_list.append(np.mean(output_new[i, 0:100]))
    if plot:
        tmp = ax[0].get_xlim()
        ax[0].hlines([np.pi, -np.pi], tmp[0], tmp[1])
        ax[0].vlines([start_loc, end_loc], -4, 4)
        ax[1].plot(start_list,'-o')
        ax[1].set_ylim([-np.pi*1.3, np.pi*1.3])
        ax2[0].set_ylim([-np.pi*1.3, np.pi*1.3])
        ax2[0].set_xlim([0, end_loc - start_loc])
        fig.canvas.draw(); fig.show()
        fig2.canvas.draw(); fig2.show()

    big_fit = []
    big_fit_xvals = []

    #remove measurements that are going to mess things up
    for i in range(inter_obj.output.shape[0]):
        print i
        values = inter_obj.output[i, start_loc:end_loc]
        #for j in range(3):
        if (-np.min(inter_obj.output[i,:]) + np.max(inter_obj.output[i,:]))<(2.*np.pi*1.5):
            #if j==0:
            big_fit.append(values - np.mean(values))
            big_fit_xvals.append(i)

    return output_new, np.array(big_fit), np.array(big_fit_xvals)



def wobbly_mirror_contour_plot(output_new, mirror_depth, n_contours=10):
    fig, ax = pt.subplots()
    #im = ax.imshow(output_new, cmap = 'RdBu', aspect='auto')#, interpolation = 'nearest')
    #im.set_clim([-np.pi, np.pi])
    cont = ax.contour(output_new,np.linspace(-mirror_depth,mirror_depth,10), cmap='RdBu')
    pt.colorbar(cont)
    ax.set_xlabel('time a.u')
    ax.set_ylabel('Interferometer channel')
    ax.set_title('Contour map of the phase shift due to moving sinusoidal mirror\nMoves up, stops, then moves down')
    fig.canvas.draw(); fig.show()



def fft_fit_wobbly_mirror(big_fit, big_fit_xvals, t, omega_mirror, mirror_depth, plot = False):
    #big_fit = np.array(big_fit); big_fit_values = np.array(big_fit_xvals)
    fft_freqs = np.fft.fftfreq(big_fit.shape[1], d=(t[1]-t[0]))
    fft_vals = np.fft.fft(big_fit)/big_fit.shape[1]
    valid_ones = fft_vals[:,np.argmin(np.abs(fft_freqs - omega_mirror/2./np.pi))]
    answer_phases = np.unwrap(np.angle(valid_ones))
    for i in range(len(answer_phases)-1):
        if (answer_phases[i+1] - answer_phases[i])<0:
            answer_phases[i+1:]=answer_phases[i+1:]+2.*np.pi
    [k_guess, phase_guess] = np.polyfit(big_fit_xvals, answer_phases, 1)
    if plot:
        fig, ax = pt.subplots(ncols = 2, sharex = True)
        ax[0].plot(big_fit_xvals, answer_phases ,'o-')
        ax[1].plot(big_fit_xvals, np.abs(valid_ones)/np.pi,'o-')
        ax[1].axhline(mirror_depth/2/np.pi)
        ax[1].set_ylim([0,mirror_depth/2/np.pi*1.2])
        ax[0].plot(big_fit_xvals, np.polyval([k_guess, phase_guess], big_fit_xvals))
        fig.canvas.draw(); fig.show()

        fig, ax = pt.subplots(nrows = 21, sharex = True)
        for i in range(big_fit.shape[0]):
            ax[i].plot(t, big_fit[i,:])
            #ax[i].plot(t, mirror_depth * np.cos(omega_mirror*t + np.angle(valid_ones[i])))
            ax[i].plot(t, mirror_depth * np.cos(omega_mirror*t + k_guess * big_fit_xvals[i] + phase_guess))
        ax[-1].set_xlim([np.min(t), np.max(t)])
        fig.canvas.draw(); fig.show()

    return answer_phases, k_guess, phase_guess, 
