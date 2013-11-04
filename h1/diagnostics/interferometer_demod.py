import scipy.special as spec
import MDSplus as MDS
import numpy as np
import matplotlib.pyplot as pt
import scipy.signal as sig
import scipy.optimize as opt


class interf_demod(object):
    def __init__(self,shot_number, f_m, dig_mult, start_samples = 10000):
        self.shot_number = shot_number
        self.f_m = f_m
        self.dig_mult = dig_mult
        self.start_samples = start_samples
        self.f_s = self.f_m * self.dig_mult

    def extract_mdsplus_data(self,):
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
        if full:
            stop_loc = len(signal)
        else:
            stop_loc = self.start_samples
        ac_rms = np.sqrt(1./self.start_samples*np.sum((signal[:stop_loc]-np.mean(signal[:stop_loc]))**2))
        dc = np.mean(signal[:stop_loc])
        return dc, ac_rms

    def plot_raw_signals(self, subtract_dc = 1):
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
        chan_num_dict = {'0':[22,5],
                         '1':[22,4],
                         '2':[22,3],
                         '3':[22,2],
                         '4':[22,1],
                         '5':[21,6],
                         '6':[21,5],
                         '7':[21,4],
                         '8':[21,3],
                         '9':[21,2],
                         '10':[21,1],
                         '11':[22,6],
                         '12':[23,1],
                         '13':[23,2],
                         '14':[23,3],
                         '15':[23,4],
                         '16':[23,5],
                         '17':[23,6],
                         '18':[24,1],
                         '19':[24,2],
                         '20':[24,3]}
        return chan_num_dict[str(chan_num)]

    def fft_signal(self,):
        self.signal_fft = np.fft.fft(self.raw_signals)
        self.fft_ax = np.fft.fftfreq(self.signal_fft.shape[1], 1./self.f_s)

        #Find the modulation frequency
        self.f_m_loc = np.argmin(np.abs(self.fft_ax - self.f_m))
        self.f_m2_loc = np.argmin(np.abs(self.fft_ax - 2*self.f_m))

    def calc_epsilon(self,):
        #calculate e
        self.e_recon = np.angle(self.signal_fft[:,self.f_m_loc]/1j)
        self.e_recon2 = np.angle(self.signal_fft[:,self.f_m_loc*2]/1j)
        self.e_recon3 = np.angle(self.signal_fft[:,self.f_m_loc*3]/1j)

    def calc_carriers_time(self, n_carriers = 6, force_symmetric = 0, individual_eps_calc = 0):
        self.carriers = np.zeros((n_carriers,self.signal_fft.shape[0], self.signal_fft.shape[1]),dtype = complex)
        #Width of the tophat
        self.width = int(self.f_m_loc)
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

    def calc_carriers_time_rfft(self, n_carriers = 6, force_symmetric = 0, individual_eps_calc = 0):
        self.carriers_rfft = np.zeros((n_carriers,self.signal_fft.shape[0], self.signal_fft.shape[1]),dtype = complex)
        self.signal_rfft = np.fft.rfft(self.raw_signals)
        #Width of the tophat
        self.width = int(self.f_m_loc)
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
        return np.abs(ratio - spec.jn(bessel_numerator,phi_1_test)/spec.jn(bessel_denomenator,phi_1_test))

    def phi1_min_brent(self,phi_1_test,ratio, bessel_numerator, bessel_denomenator):
        return ratio - spec.jn(bessel_numerator,phi_1_test)/spec.jn(bessel_denomenator,phi_1_test)

    def calc_phi1_means_std(self, harm1, harm2, start_point = 0, end_point = None):
        if end_point == None: end_point = self.start_samples
        self.ratio = np.abs(self.carriers[harm1,:,start_point:end_point]) / np.abs(self.carriers[harm2,:,start_point:end_point])
        self.ratio_means = np.mean(self.ratio, axis = -1)
        self.ratio_stds = np.std(self.ratio, axis = -1)
        fig, ax = pt.subplots()
        ax.errorbar(range(len(self.ratio_means)), self.ratio_means, yerr=self.ratio_stds)
        fig.canvas.draw(); fig.show()
        self.phi1_mean = []
        self.phi1_std_upper = []
        self.phi1_std_lower = []
        for i in range(self.ratio_means.shape[0]):
            print opt.fmin(self.phi1_min, np.deg2rad(135), args=(self.ratio_means[i],harm1, harm2),disp=0, xtol=0.00001)
            print opt.brentq(self.phi1_min_brent, 0, 2.*np.pi, args=(self.ratio_means[i],harm1, harm2),disp=0, xtol=0.00001)
            self.phi1_mean.append(opt.fmin(self.phi1_min, np.deg2rad(135), args=(self.ratio_means[i],harm1, harm2),disp=0, xtol=0.00001)[0])
            self.phi1_std_upper.append(opt.fmin(self.phi1_min, np.deg2rad(135), args=(self.ratio_means[i]+self.ratio_stds[i], harm1, harm2),disp=0, xtol=0.00001)[0])
            self.phi1_std_lower.append(opt.fmin(self.phi1_min, np.deg2rad(135), args=(self.ratio_means[i]-self.ratio_stds[i], harm1, harm2),disp=0, xtol=0.00001)[0])
        


    def extract_harmonic_amps(self,start_point = 0, end_point = None):
        if end_point == None: end_point = self.start_samples
        end_point += (start_point - end_point)%self.dig_mult
        tmp_fft = np.fft.fft(self.raw_signals[:,start_point:end_point])
        tmp_fft_ax = np.fft.fftfreq(tmp_fft.shape[1], 1./self.f_s)
        f_m_loc_tmp = np.argmin(np.abs(tmp_fft_ax-self.f_m))
        n_harmonics = int((self.f_s/2.)/self.f_m)
        if (self.f_m*n_harmonics)>=(self.f_s/2.): n_harmonics -= 1
        self.harmonic_list = np.arange(n_harmonics+1)*f_m_loc_tmp
        self.harmonics = tmp_fft[:,self.harmonic_list]
        fig, ax = pt.subplots(nrows=5, ncols = 5, sharex = True, sharey = True); ax = ax.flatten()
        for i in range(tmp_fft.shape[0]):
            ax[i].plot(tmp_fft_ax, np.abs(tmp_fft[i,:]))
            ax[i].plot(tmp_fft_ax[self.harmonic_list], np.abs(self.harmonics[i,:]), 'o')
        ax[-1].set_xlim([0,np.max(tmp_fft_ax)])
        ax[-1].set_ylim([0,np.max(np.abs(self.harmonics))])
        fig.subplots_adjust(hspace=0.015, wspace=0.015,left=0.10, bottom=0.10,top=0.95, right=0.95)
        fig.canvas.draw(); fig.show()

    def compare_harmonics(self,phi1 = None):
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

        fig, ax = pt.subplots(nrows=5, ncols = 5, sharex = True, sharey = True); ax = ax.flatten()
        for i in range(chans):
            ax[i].plot(range(1,max_harmonic),np.sign(np.real(self.harmonics[i,1:]))*np.abs(self.harmonics[i,1:])*np.abs(bessel_amps[i,1])/ np.abs(self.harmonics[i,1]),'o-')
            ax[i].plot(range(1,max_harmonic), bessel_amps[i,1:],'x-')
        ax[-1].set_xlim([0,max_harmonic])
        ax[-1].set_ylim([0,np.max(bessel_amps)])
        
        #ax[-1].set_ylim([0,np.max(np.abs(self.harmonics))])
        fig.subplots_adjust(hspace=0.015, wspace=0.015,left=0.10, bottom=0.10,top=0.95, right=0.95)
        fig.canvas.draw(); fig.show()




    def calculate_phase(self, phi1, odd, even, use_rfft = 0, clim = [0,5]):
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
        self.im_fig, self.im_ax = pt.subplots(ncols = 2, sharey=True)
        self.im =self.im_ax[0].imshow(np.abs(self.output),interpolation = 'nearest', aspect='auto',cmap = 'hot', extent=[self.t[0],self.t[-1],self.output.shape[0],0], origin = 'upper')
        times = [0.02,0.04,0.06]
        for t in times:
            t_loc = np.argmin(np.abs(self.t - t))
            plot_vals = np.abs(self.output[:,t_loc])<10
            self.im_ax[1].plot(np.abs(self.output[:,t_loc][plot_vals]),np.arange(self.output.shape[0])[plot_vals],'x-')
            self.im_ax[0].vlines(t, self.output.shape[0],0, colors = 'b')
        self.im_ax[1].grid()
        self.im_ax[1].set_xlim(clim)
        #self.im =self.im_ax.imshow(np.abs(np.flipud(self.output)), aspect='auto',cmap = 'hot', extent=[self.t[0],self.t[-1],self.output.shape[0],0])
        self.im.set_clim(clim)
        self.im_ax[0].set_xlim([-0.005,0.12])
        self.im_ax[1].set_ylim([0.,self.output.shape[0]])
        pt.colorbar(self.im,ax = self.im_ax[0])
        self.im_fig.canvas.draw(); self.im_fig.show()

            # for i in range(0,shot_start*phi1_calc_width,500):
            #     ratio = np.real(h_n[3])/np.real(h_n[1])
            #     ratio2 = np.real(h_n[4])/np.real(h_n[2])
            #     ratio3 = np.real(h_n[5])/np.real(h_n[3])

            #     ratio = np.abs(h_n[3])/np.abs(h_n[1])
            #     ratio2 = np.abs(h_n[4])/np.abs(h_n[2])
            #     ratio3 = np.abs(h_n[5])/np.abs(h_n[3])
            #     tmp1 = opt.fmin(phi_1_min,np.deg2rad(135), args=(ratio[i],3, 1),disp=0, xtol=0.00001)
            #     tmp2 = opt.fmin(phi_1_min,np.deg2rad(135), args=(ratio2[i],4, 2),disp=0, xtol=0.00001)
            #     tmp3 = opt.fmin(phi_1_min,np.deg2rad(135), args=(ratio3[i],5, 3),disp=0, xtol=0.00001)
            #     #phi1_vals1.append(tmp1)
            #     #phi1_vals2.append(tmp2)
            #     #i_vals.append(i)
            #     tmp1_error = np.abs(ratio[i] - spec.jn(3,tmp1)/spec.jn(1,tmp1))
            #     tmp2_error = np.abs(ratio2[i] - spec.jn(4,tmp2)/spec.jn(2,tmp2))
            #     tmp3_error = np.abs(ratio3[i] - spec.jn(5,tmp3)/spec.jn(3,tmp3))
            #     if tmp1_error<1.e-1: 
            #         phi1_vals1.append(tmp1)
            #         i_vals1.append(i)
            #     if tmp2_error<1.e4: 
            #         phi1_vals2.append(tmp2)
            #         i_vals2.append(i)
            #     if tmp3_error<1.e4: 
            #         phi1_vals3.append(tmp3)
            #         i_vals3.append(i)
            #     print tmp1, tmp1_error, tmp2, tmp2_error, tmp3, tmp3_error
