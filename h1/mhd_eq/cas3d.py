import numpy as np
import matplotlib.pyplot as pt
import h1.helper.generic_funcs as gen_funcs
import copy
import pylab as pl
from matplotlib.font_manager import FontProperties
from StringIO import StringIO
import re

class cas3d_results():
    def __init__(self, directory, single_eig = True, output_filename = 'out_cas3d'):
        self.directory = directory
        self.single_eig = single_eig
        self.output_filename = output_filename

    def get_ew_khz(self,):
        fname = '{}/{}'.format(self.directory, self.output_filename)
        print 'fname :', fname
        with file(fname, 'r') as file_handle:
            data = file_handle.readlines()
        found = False
        for i in range(len(data)):
            if data[i].find('hydrogen')>=0:
                print 'found!!', data[i]
                out_data = data[i:i+4]
                found = True
                break
        #print out_data
        if found:
            for i in out_data:
                if i.find('frequency')>=0:
                    dat = i[i.find('=')+1:i.find('[')]
                    break
            self.ew_freq_khz = float(dat)
        else:
            self.ew_freq_khz = 0
        print self.ew_freq_khz

    def plot_eigenvalues(self, ax = None, ylims = None, max_harms = None, sqrt_s = False, multiplier = 1, plot_kwargs = None,):

        self.get_ew_khz()
        with file(self.directory + '/'+ self.output_filename, 'r') as file_handle:
            file_lines = file_handle.readlines()
        if plot_kwargs == None: plot_kwargs = {}
        harmonic_strings = ['largest xis harmonics', 'largest eta harmonics', 'largest mu  harmonics']
        harmonic_strings = ['largest xis harmonics']
        coeff_strings = ['Fourier coefficients of     xis', 'Fourier coefficients of     eta','Fourier coefficients of      mu']
        titles = [r'$\xi^s$', r'$\xi^\eta$', r'$\xi^\mu$']
        own_plots = True if ax == None else False
        if own_plots: fig, ax = pt.subplots(nrows = len(harmonic_strings), sharex =True , sharey=False)
        self.identifiers = []; self.m_n_list = []; self.max_value = []
        if ylims == None: ylims = [None] * len(harmonic_strings)

        for i in range(0, len(harmonic_strings)):
            max_val_plot = 0.
            min_val_plot = 0.
            identifier, m_n_list, max_value = self.return_harmonics(file_lines, harmonic_strings[i])
            svals, coefficient_array = self.return_coeffs(file_lines, coeff_strings[i])
            dominant_ind = np.argmax(np.sum(np.abs(coefficient_array[3:-3,:]**2), axis = 0))
            max_mode = coefficient_array[3:-3,dominant_ind]
            max_val = np.min(max_mode) if (np.min(max_mode)**2>np.max(max_mode)**2) else np.max(max_mode)
            multiplier = multiplier/np.abs(multiplier) / max_val

            max_ind = np.argmin(np.abs(svals - 0.95))
            if max_harms == None: max_harms = len(m_n_list) - 1
            plotted_non_dom = False
            for j in range(0,max_harms+1):
                dominant = dominant_ind == j
                plot_kwargs['color'] = 'k' if dominant else 'r'
                text_string = 'CAS3D:({},{})'.format(-m_n_list[j][1],-m_n_list[j][0])
                plot_kwargs['label'] = text_string if dominant else None
                if not plotted_non_dom and not dominant:
                    plotted_non_dom = True; plot_kwargs['label'] = 'CAS3D:rest'
                xax = svals if not sqrt_s else np.sqrt(svals)
                #ax[i].plot(xax, coefficient_array[:,j]*multiplier, label = '{}'.format(m_n_list[j]), **plot_kwargs)
                ax[i].plot(xax, coefficient_array[:,j]*multiplier,  **plot_kwargs)
                max_loc = np.argmax(np.abs(coefficient_array[:max_ind,j]))
                max_val = np.max(coefficient_array[:max_ind,j])
                min_val = np.min(coefficient_array[:max_ind,j])
                if max_val>max_val_plot: max_val_plot = max_val
                if min_val<min_val_plot: min_val_plot = min_val
                text_string = m_n_list[j]
                #ax[i].text(xax[max_loc], coefficient_array[max_loc,j]*multiplier, '{}'.format('{},{}'.format(-text_string[1], -text_string[0])))


            ax[i].grid()
            #ax[i].legend(fontsize = 7, loc = 'best')
            if ylims[i]!=None: ax[i].set_ylim(ylims[i])
            print max_val_plot
            #ax[i].set_ylim([-max_val_plot*1.05, max_val_plot*1.05])
            extra = (max_val_plot - min_val_plot)*0.1
            ax[i].set_ylim([min_val_plot - extra, max_val_plot+extra])
            ax[i].set_ylabel(titles[i])
            self.identifiers.append(identifier); self.m_n_list.append(m_n_list); self.max_value.append(max_value)
        #ax[i].set_ylim([-0.0008,0.0008])
        ax[i].set_xlim([0,1])
        ax[i].set_xlabel('s')
        if own_plots: 
            fig.suptitle(self.directory.replace('_','\_'))
            fig.canvas.draw(); fig.show()

        #direct read in of an xis file
        # tmp_data = np.loadtxt(self.directory + '/xis_largest01-10.dat')
        # fig, ax = pt.subplots()
        # for j in range(1,tmp_data.shape[1]):
        #     xax = tmp_data[:,0] if not sqrt_s else np.sqrt(tmp_data[:,0])
        #     xlabel = 's' if not sqrt_s else '$\sqrt{s}$'
        #     #ax.plot(tmp_data[:,0], tmp_data[:,j],'.-')
        #     ax.plot(xax, tmp_data[:,j],'.-')
        # ax.set_xlabel(xlabel)
        # #ax.set_ylabel(r'$\xi^s$')
        # ax.set_ylabel(r'Normal displacement (a.u)')
        # fig.suptitle(self.directory.replace('_','\_'))
        # fig.canvas.draw();fig.show()
        
    def return_harmonics(self, file_lines, harmonic_string):
        #find the start point of the interesting stuff
        found = 0; i=0
        while found==0 and i<len(file_lines):
            if file_lines[i].find(harmonic_string)>=0:
                start_point = i; found = 1
            i+=1
        identifier = []; m_n_list = []; max_value = []
        for i in range(start_point+2, start_point+12):
            tmp = file_lines[i].split()
            identifier.append(int(tmp[0]))
            m_n_list.append([int(tmp[1]), int(tmp[2])])
            max_value.append(float(tmp[3]))
        return np.array(identifier), np.array(m_n_list), np.array(max_value)

    def return_coeffs(self, file_lines, harmonic_string):
        #find the start point of the interesting stuff
        found = 0; i=0
        while found==0 and i<=len(file_lines):
            if file_lines[i].find(harmonic_string)>=0:
                start_point = i; found = 1
            i+=1
        coefficient_array = []; svals = []

        current_loc = start_point + 2
        while len(file_lines[current_loc].split())>4:
            tmp = file_lines[current_loc].split()
            svals.append(float(tmp[0]))
            coefficient_array.append([float(x) for x in tmp[1:]])
            current_loc+=1
        return np.array(svals), np.array(coefficient_array)

    def plot_ne_Te_P(self,ne = True, Te = False, P = True, fname_ne = 'perturbed_density.dat', fname_Te = 'perturbed_temperature.dat', fname_P = 'perturbed_pressure.dat', max_harms = None, ax = None, ylims = None, divide = False, multiply = False, sqrt_s = False, multiplier = 1, ne_mult = None):
        #fname_ne = 'perturbed_ne_f.dat'
        #fname_Te = 'perturbed_pressure_f.dat'
        #fname_P = 'perturbed_ne_f.dat'
        print 'divide {}, multiply {}'.format(divide, multiply)
        fnames = []; yax = []
        divisors = []
        pert_quants = []
        multipliers = []
        if ne: 
            fnames.append(fname_ne)
            yax.append('Perturbed ne (a.u)')
            divisors.append(self.rho)
            if ne_mult == None:
                multipliers.append(self.rho)
            else:
                print 'using provided ne mult'
                multipliers.append(ne_mult)
            pert_quants.append('ne1')
        if Te: 
            fnames.append(fname_Te)
            yax.append('Perturbed Te (a.u)')
            divisors.append(self.Te)
            multipliers.append(self.Te)
            pert_quants.append('Te1')
        if P: 
            fnames.append(fname_P)
            yax.append('Perturbed P (a.u)')
            divisors.append(self.P)
            multipliers.append(self.P)
            pert_quants.append('P1')
        #yax = ['Perturbed Pressure', 'Perturbed Density']
        own_plots = True if ax == None else False
        if ylims == None: ylims = [None] * len(yax)
        if own_plots: fig, ax = pt.subplots(nrows = len(fnames), sharex = True)
        for fname, yname, ax_cur, ylim, div, mult, pert_name in zip(fnames, yax, ax, ylims, divisors, multipliers, pert_quants):

            max_val_plot = 0.
            min_val_plot = 0.
            print self.directory + '/'+fname
            with file(self.directory + '/'+fname,'r') as file_handle:
                lines = file_handle.readlines()
                for i_tmp,line in enumerate(lines):
                    if line.find('poloidal m')>=0:
                        start_line = i_tmp+1
                        break
                for i_tmp,line in enumerate(lines):
                    if line.find('radial grid')>=0:
                        end_line = i_tmp-1
                        break
                mode_nums = [i_tmp.lstrip('#').lstrip(' ').rstrip('\n').rstrip(' ') for i_tmp in lines[start_line:end_line]]
                mode_nums = [re.sub('\ +',',',i_tmp) for i_tmp in mode_nums]
            a = np.loadtxt(self.directory + '/'+fname, comments='#')
            tmp_amps = +a[:,1:] if not divide else a[:,1:]/div[:, np.newaxis]
            tmp_amps = tmp_amps if not multiply else tmp_amps*mult[:, np.newaxis]
            dominant_ind = np.argmax(np.sum(np.abs(tmp_amps[3:-3,:]), axis = 0))
            #mode_sums = np.sum(tmp_amps[5:-5,:]**2,axis = 1)
            max_mode = tmp_amps[:, dominant_ind]
            #dx = a[1,0] - a[0,0]
            #multiplier = multiplier/np.abs(multiplier) / np.max(np.sum(np.abs(tmp_amps), axis = 0)* dx)
            #max_loc = np.max(tmp_amps)
            #max_val = np.min(tmp_amps) if np.abs(np.min(tmp_amps))>np.max(tmp_amps) else np.max(tmp_amps)
            max_val = np.min(max_mode) if (np.min(max_mode)**2>np.max(max_mode)**2) else np.max(max_mode)
            #multiplier = multiplier/np.abs(multiplier) / tmp_amps[max_loc]
            multiplier = multiplier/np.abs(multiplier) / max_val
            #multiplier = multiplier/np.abs(multiplier) / np.max(tmp_amps)
            print np.sum(tmp_amps, axis = 0), np.sum(tmp_amps, axis = 0).shape
            print 'mult at end', multiplier
            self.s_cur = a[:,0]
            #max_ind = np.argmin(np.abs(self.s_cur - 1.0))
            max_ind = None
            if max_harms == None: max_harms = a.shape[1]-1
            setattr(self, pert_name, a[:,1:])
            plotted_non_dom = False
            for i in range(1, max_harms+1):
                try:
                    #self.identifiers[0].tolist().index(i)
                    x_ax = np.sqrt(self.s_cur) if sqrt_s else self.s_cur
                    plot_val = +a[:,i] if not divide else a[:,i]/div
                    plot_val = plot_val if not multiply else plot_val*mult
                    text_string = 'CAS3D:{},{}'.format(-int(mode_nums[i-1].split(',')[1]), -int(mode_nums[i-1].split(',')[0]))
                    dominant = dominant_ind == (i-1)
                    plot_style = '-k.' if dominant else '-r'
                    label = text_string if dominant else None
                    if not plotted_non_dom and not dominant:
                        plotted_non_dom = True; label = 'CAS3D:rest'
                    #if dominant and np.sum(plot_val)<0: plot_val *= -1
                    ax_cur.plot(x_ax, plot_val*multiplier, plot_style, label = label)
                    max_val = np.max(plot_val[:max_ind]*multiplier)
                    min_val = np.min(plot_val[:max_ind]*multiplier)
                    if max_val>max_val_plot: max_val_plot = max_val
                    if min_val<min_val_plot: min_val_plot = min_val
                except ValueError:
                    print 'errors....'
                    pass

            #ax_cur.set_xlabel('s')
            ax_cur.set_ylabel(yname)
            print 'hello'
            #if ylim!=None: ax_cur.set_ylim(ylim)
            print max_val_plot
            #ax_cur.set_ylim([-max_val_plot*1.05, max_val_plot*1.05])
            ax_cur.set_ylim([min_val_plot-0.2, max_val_plot+0.2])
            ax_cur.grid(True)
        if own_plots: fig.canvas.draw(); fig.show()

    def eq_quants_pnt(self, fname = 'equilibrium_profiles_pnt.dat', ax = None,plot = True):
        print self.directory + '/equilibrium_profiles_pnt.dat'
        data = np.loadtxt(self.directory + '/equilibrium_profiles_pnt.dat', comments='#')
        labels = ['P', 'dP/ds','rho','drho/ds','Te','dTe/ds']
        names = ['P', 'dPds','rho','drhods','Te','dTeds']
        own_ax = True if ax == None else False
        if own_ax and plot: 
            fig, ax = pt.subplots(ncols = 2, nrows = 3, sharex = True)
            gen_funcs.setup_publication_image(fig, height_prop = 1./1.618, single_col = False, replacement_kwargs = None)
        self.s = data[:,0]
        for i, name in enumerate(names):
            setattr(self,name,data[:,i+1])
        if plot:
            for i, (cur_ax, label, name) in enumerate(zip(ax.flatten(), labels, names)):
                cur_ax.plot(data[:,0], data[:,i+1])
                #setattr(self,name,data[:,i+1])
                cur_ax.set_ylabel(label)
            if own_ax:
                for j in ax[-1,:]: j.set_xlabel('s')
                fig.tight_layout(pad = 0.1)
                fig.savefig('/home/srh112/win_shared/tomo_referee_response/profiles.png')
                fig.canvas.draw(); fig.show()


class ModeTable:
    """ Recover info about mode table, a bit fussy about searches, for good 
        reason.
        Should be able to add methods to show the table, and check mode families
    """
    def __init__(self, target, buff, name=None):
        #target = 'User supplied perturbation spectrum'
        if name == None: name = target[0].strip().upper()
        self.name = name
        wst = np.where([target[0] in lin for lin in buff])[0]
        if len(wst) != 1: 
            raise IOError('Target {t} not found just once in {f}'.
                          format(t=target[0],f=outfn))
        wst = wst[0]
        wend = wst + np.where([target[1] in lin for lin in buff[wst:]])[0]
        if len(wend)<1: 
            raise IOError('Target{t1} not found in {f}'.
                          format(t1=target[1],f=outfn))
        else: 
            if target[1].strip() != '': 
                wend=wend[0]   # don't wish to read target
            else: wend = wend[2]
        wst += 1 # move past target
        while (buff[wst].strip() == ''): wst += 1
        # extra code for multi-line mode tables
        wblank = wst
        while buff[wblank].strip() != '': wblank += 1
        nlines = wblank - wst
        cbuff = (''.join(buff[wst:wblank])).replace('\n',' ')
        self.ms = np.loadtxt(StringIO(cbuff),dtype=int)
        #self.ms = np.loadtxt(StringIO(buff[wst]),dtype=int)  #used to be 1 line
        #extra code for multi-line
        first = wst + 2 + nlines-1
        cbuffs = [(''.join(buff[fr:fr+nlines])).replace('\n',' ') for 
                  fr in range(first,wend,nlines)]
        self.wtabl = np.loadtxt(StringIO('\n'.join(cbuffs)),dtype=int)
        #self.wtabl = np.loadtxt(StringIO(' '.join(buff[wst+2:wend])),dtype=int)
        self.ns = self.wtabl[:,0]
        self.inds = self.wtabl[:,1:]
        self.lut = {}  # make a look up table as a dict keyed on mode number
        nmin = min(self.ns) # each entry is a small dict with m and n values
        mmin = min(self.ms)
        for n in self.ns:
            for m in self.ms:
                self.lut.update({self.inds[n-nmin,m-mmin]: {'m':m, 'n':n}})
        self.rest = buff[wend:]
        active_Ns = [self.lut[nn]['n'] for nn in self.lut.keys()]
        self.Np = np.unique(np.abs(np.diff(np.unique(active_Ns))))
        self.info = str('{nm} mode table has {n} elements, '
                        '{n0}<=n<={n1} (step {Np}), {m0}<=m<={m1}'
                        .format(n = len(self.lut.keys()),
                                nm=self.name, Np = self.Np,
                                n0=min(self.ns),n1=max(self.ns),
                                m0=min(self.ms),m1=max(self.ms), 
                 ))



class cas3d_continua_data():
    def __init__(self, path = '/home/jbb105/short/030_n1_low_tinier/', suff = '', debug = 1,  minscan = 0, maxscan = 40, modes = 'polar', get_wkin = True, no_hash = True, exception = Exception, output_log = 'out_cas3d', load_eig_funcs = True, alf_sound_cut = 0.01, alf_sound_ratio = 'sum'):
        offs = 2 # allow for inserting both hash flags ahead of the rest. 
        self.m_ind = 2 + offs # the column containing the mode index
        self.A_sw=5 + offs     # column for the Aflven/sound flag (1/0)
        self.A_cpt=self.A_sw+1 
        self.s_ind = 0 + offs
        self.f_ind = self.s_ind+1
        hash_ind = 0

        if path[-1] != '/': path = path + '/'  # add a trailing / if missing
        if suff == '.bz2': 
            import bz2  # this import will only be attempted if suff='.bz2'
            OPEN=bz2.BZ2File
        elif suff == '.gz': 
            import gzip
            OPEN=gzip.GzipFile
        else: OPEN = open

        outfn = path+output_log+suff
        try:
            outcas3d = OPEN(outfn,'r')
            buff=outcas3d.readlines(100000)   # arg is # characters to read, approx.
        except IOError, details:
            print('**Warning: out_cas3d file {o} not found - cannot label modes'.format(o=outfn))
        self.XIS=ModeTable(buff=buff,target=[' xis\n','eta'])
        self.mode_table_info = self.XIS.info

        try:
            self.ETA=ModeTable(buff=self.XIS.rest,target=[' eta\n', 'mu'])
            self.MU=ModeTable(buff=self.ETA.rest,target=[' mu\n', ' \n'])
        except exception, details:
            print('Unable to read other mode tables {d}'.format(d=details))
            if verbose>0: print('{d}'.format(d=details))
        else: # here if successful
            for minfo in [self.ETA.info,self.MU.info]:
                if self.XIS.info.split(None,1)[1] == minfo.split(None,1)[1]:
                    minfo = '<as per XIS>'
                self.mode_table_info += ', '+minfo

        scan_arr = None
        all_lines=[]
        tmp1 = []
        tmp2 = []
        tmp3 = []
        import time
        def read_scan_old():
            overall = []
            for ind in range(1,4):
                tmp = []
                for n in range(minscan,maxscan):
                    xisn = path+'xis{}_{}.dat'.format(ind, n) + suff
                    with OPEN(xisn,"r") as xf: 
                        xfbuff = xf.readlines()
                        for i in xfbuff:
                            z = i.rstrip('\n').split(' ')
                            for j in z:
                                if j!='':
                                    if j.find('E')>=0:
                                        tmp.append(float(j))
                                    else:
                                        tmp.append(0)
                overall.append(tmp)
            return overall
        import fileinput

        def read_scan_new():
            overall = []
            for ind in range(1,4):
                data = np.loadtxt(fileinput.input([path+'xis{}_{}.dat'.format(ind, n) + suff for n in range(minscan,maxscan)]))
                print data
                print data.shape
                overall.append(data)
            return overall

        for n in range(minscan,maxscan):
            scfn = path+'scan_{0}.dat'.format(n) + suff
            wkfn = path+'wkin_{0}.dat'.format(n) + suff
            if get_wkin:
                try:
                    with OPEN(scfn,"r") as sf: 
                        scbuff = sf.readlines()
                    with OPEN(wkfn,"r") as wf: 
                        wkbuff = wf.readlines()
                    for (ln, s_line) in enumerate(scbuff[2:]):
                        hashes=(s_line[0]+wkbuff[ln+1][0])
                        hashes_cvt = hashes.replace(' ','0 ').replace('#','1 ')
                        all_lines.append(hashes_cvt+   # # replaced by 1, space by 0
                                         s_line[1:-1]  # omit first # and last \n
                                         +wkbuff[ln+1][1:]) # omit just first # 
                except IOError, details: 
                    print('File read error - perhaps reduce maxscan to {n}: {d}'
                          .format(n=n, d=details))
                    break

            else:
                a = np.loadtxt(scfn).T
                if len(np.shape(a)) == 1: a = np.array([a]) # make a 2D array if one line
                if scan_arr == None: scan_arr = a
                elif len(np.shape(a)) == 2 : scan_arr = np.append(scan_arr, a, 1)
                else: print('discarding {f}, shape is {s}'.format(f=scfn, s=np.shape(a)))

        #scan_arr - 0:scan #, 1:wkin#, 2:s, 3:w(kHz), 4:harm, 5:num_in_scan, 6:max xis, 7:type 0=sound
        #8:?, 9:?, 10:?, 11:s, 12: ew, 13:wkin, 14:wkin_norm, 15:wkin_binorm, 16:wkin_parallel, 17:sum
        #18:?, 19:?,20:?
        if get_wkin: self.scan_arr = np.loadtxt(StringIO(' '.join(all_lines))).T
        if load_eig_funcs:
            with file(path+'xis{}_{}.dat'.format(1, 0) + suff,'r') as file_handle:
                self.new_method = True if len(file_handle.readlines()[0])>100 else False
            if self.new_method == False:
                tmp1, tmp2, tmp3 = read_scan_old()
                self.mode_shapes_all = np.zeros((3,self.scan_arr.shape[1],len(tmp1)/self.scan_arr.shape[1]),dtype = float)
                for i, tmp_tmp in enumerate([tmp1,tmp2, tmp3]):
                    self.mode_shapes_all[i,:,:] = np.resize(tmp_tmp, (self.scan_arr.shape[1],len(tmp1)/self.scan_arr.shape[1]))
                self.mode_nums = False
            else:
                overall = read_scan_new()
                self.mode_shapes_all = np.zeros((len(overall),self.scan_arr.shape[1],overall[0].shape[1]-3),dtype = float)
                self.mode_nums = np.zeros((len(overall),self.scan_arr.shape[1]),dtype = float)
                print 'mode_shapes_all', self.mode_shapes_all
                print 'mode_nums', self.mode_nums
                for i, tmp_tmp in enumerate(overall):
                    print tmp_tmp[:,1:-2].shape
                    self.mode_shapes_all[i,:,:] = tmp_tmp[:,1:-2]
                    self.mode_nums[i,:] = tmp_tmp[:,0]
        if no_hash: 
            w_no_hash = np.where(self.scan_arr[hash_ind] == 0)[0]
            print('Ignoring {h} hashed of {t} instances'.
                  format(h=len(self.scan_arr[0])-len(w_no_hash), t=len(self.scan_arr[0])))
            self.scan_arr=self.scan_arr[:,w_no_hash]
            if load_eig_funcs:
                self.mode_shapes_all = self.mode_shapes_all[:,w_no_hash,:]
                if self.new_method: self.mode_nums = self.mode_nums[:,w_no_hash]
        print self.scan_arr.shape
        if alf_sound_ratio == 'sum':
            self.mratio = self.scan_arr[self.A_cpt]/np.sum(self.scan_arr[self.A_cpt:self.A_cpt+3].T,1)
        else:
            self.mratio = self.scan_arr[self.A_cpt]/np.std(self.scan_arr[self.A_cpt:self.A_cpt+3].T,1)

        self.allmodes = np.unique(self.scan_arr[self.m_ind]).astype(int)
        self.wA = np.where(self.scan_arr[self.A_sw] == 1)[0]
        self.Alfvenicmodes = np.unique(self.scan_arr[self.m_ind][self.wA]).astype(int)
        self.wMM = np.where(self.mratio>alf_sound_cut)[0]
        self.mmodes = np.unique(self.scan_arr[self.m_ind][self.wMM]).astype(int)

        def dither(dith=.001, arr=None):
            return dith*(np.random.random(len(arr))-0.5)

        if debug>0: print('Xis modes {m} found, {a} are classifed as '
                          'Alfvenic by "sound_test()"'.
                          format(m=self.allmodes,a=self.Alfvenicmodes))

        if pl.is_string_like(modes):
            if 'alf' in modes.lower(): self.modes = self.Alfvenicmodes
            elif 'pol' in modes.lower(): self.modes = self.mmodes
            else: 
                try:
                    self.modes = eval(modes)
                except:
                    raise ValueError,"modes = {m} is not understood".format(m=modes)

        # mode score will be used to allocate colors and position in legend
        self.mode_score = np.zeros(max([k for k in self.XIS.lut.keys()])+1, dtype=int)
        self.mode_score[self.mmodes] += 100  # give M-type modes the highest rank
        self.mode_score[self.Alfvenicmodes] += 50  # give Alf modes a high rank
        for k in self.XIS.lut.keys():    # bonus for small n,m
            self.mode_score[k] += -(abs(self.XIS.lut[k]['n']) + abs(self.XIS.lut[k]['m']))  
            # alternatively could rank by closeness to resonance?
        self.ct = np.arange(len(self.mode_score))
        self.ct[np.argsort(-self.mode_score)] = np.arange(len(self.mode_score))


        a = np.loadtxt(path+'physical_profiles.dat', comments = '#')
        self.s = a[:,0]; self.beta = a[:,1]; self.iota = a[:,2]; self.rho = a[:,3]; self.p = a[:,4]
        #list_of_lines = [] #SH added

    def plot_continua(self, ax, colorset = None, edgecolorset = None, dith = None, rest_same = True, size_scale = 50, sym_all = '.c', ms_all = 1, label_type = 'full', linewidths = None, ylims = None, plot_crosses = True,plot_styles = None):
        self.plot_styles = {} if plot_styles == None else plot_styles
        if dith==None: dith = [0.001, .1]
        if colorset == None: colorset = ('b,g,r,c,m,y,k,orange,purple,lightgreen,gray'.split(','))
        if edgecolorset == None: edgecolorset='k,b,darkgreen,r,purple,green,darkgray,darkblue,c,m'.split(',')
        edgecolors = []
        colors = 10*colorset
        for i in range(10): 
            for j in range(len(colorset)): 
                edgecolors.append(edgecolorset[i])

        def dither(dith=.001, arr=None):
            return dith*(np.random.random(len(arr))-0.5)
        tot_plot = 0
        if rest_same: 
            tmp, = ax.plot(self.scan_arr[self.s_ind]+dither(dith[0],self.scan_arr[0]), # all modes 
                       self.scan_arr[self.f_ind],sym_all,markersize=ms_all,
                           label='all ({n})'.format(n=len(self.scan_arr[0])), rasterized = True)
            
        else: 
            pass
            # for (m,mode)in enumerate(self.allmodes): # 
            #     wM = np.where((self.scan_arr[m_ind]==mode))[0]
            #     tmp, = ax.plot(self.scan_arr[s_ind][wM],self.scan_arr[f_ind][wM],'.',
            #                    markersize=msall,color =colors[self.ct[mode]],
            #                    label='allmode {m}'.format(m=mode))

        if plot_crosses:
            for (m,mode)in enumerate(self.Alfvenicmodes): # (tiny offset in s to separate)
                wM = np.where((self.scan_arr[self.m_ind]==mode) &(self.scan_arr[self.A_sw] == 1))[0]
                if label_type == 'full':
                    label = 'Amode {m} ({n})'.format(m=mode,n=len(wM))
                else:
                    label = None
                tmp,=ax.plot(self.scan_arr[self.s_ind][wM]+0.0001,self.scan_arr[self.f_ind][wM],'+',color =colors[self.ct[mode]],
                             label= label, rasterized = True)
                #list_of_lines.append(tmp)
        for (m,mode)in enumerate(self.modes):
            wMode = np.where((self.scan_arr[self.m_ind][self.wMM]==mode))[0]
            tot_plot += len(wMode)
            if len(wMode) != 0:
                if label_type == 'full':
                    label='M. mode {mode} {N}/{M} ({m}:{num})&{ord:03d}'.format(mode=mode,m=self.ct[mode],num=len(wMode),
                                                                            M=self.XIS.lut[mode]['m'],N=self.XIS.lut[mode]['n'],
                                                                            ord=self.ct[mode] )
                else: 
                    #label='{N}/{M}'.format(M=self.XIS.lut[mode]['m'],N=self.XIS.lut[mode]['n'])
                    label='({N},{M})'.format(M=-self.XIS.lut[mode]['m'],N=-self.XIS.lut[mode]['n'])
                tmp_x = self.scan_arr[self.s_ind][self.wMM[wMode]]+dither(dith[0],wMode)
                tmp_y = self.scan_arr[self.f_ind][self.wMM[wMode]]+dither(dith[1],wMode)
                tmp_z = size_scale*np.sqrt(self.mratio[self.wMM[wMode]])
                if ylims!=None: 
                    truth = (tmp_y > ylims[0]) * (tmp_y < ylims[1])
                else:
                    truth = np.ones(len(tmp_x),dtype = bool)
                pop_list = []
                for tmp_i in self.plot_styles.keys():
                    clr = self.plot_styles[tmp_i][0]
                    edge_clr = self.plot_styles[tmp_i][1]
                    for iii in range(len(colors)):
                        if (colors[iii]==clr) and (edgecolors[iii] == edge_clr):pop_list.append(iii)
                for iii in pop_list:
                    colors.pop(iii)
                    edgecolors.pop(iii)
                if np.sum(truth)>10:
                    plot_style_key = '{},{}'.format(self.XIS.lut[mode]['n'], self.XIS.lut[mode]['m'])
                    if plot_style_key in self.plot_styles.keys():
                        clr = self.plot_styles[plot_style_key][0]
                        edge_clr = self.plot_styles[plot_style_key][1]
                        #pop_list = []
                        #for iii in range(len(colors)):
                        #    if (colors[iii]==clr) and (edgecolors[iii] == edge_clr):pop_list.append(iii)
                        #for iii in pop_list:
                        #    colors.pop(iii)
                        #    edgecolors.pop(iii)
                    else:
                        clr = colors[self.ct[mode]]
                        edge_clr = edgecolors[self.ct[mode]]
                        self.plot_styles[plot_style_key] = [clr, edge_clr]
                    z=ax.scatter(tmp_x[truth], tmp_y[truth], s = tmp_z[truth],marker='o',
                                 c = clr,edgecolors=edge_clr,
                                 label=label, zorder = 100,rasterized = True, linewidths = linewidths)  
                    #z=ax.scatter(tmp_x[truth], tmp_y[truth], s = tmp_z[truth],marker='o',
                    #             c = colors[self.ct[mode]],edgecolors=edgecolors[self.ct[mode]],
                    #             label=label, zorder = 100,rasterized = True, linewidths = linewidths)  
        # the fine dot plots are dithered only in S

    def picker_widget(self, path = '/home/jbb105/short/030_n1_low_tinier/', suff = '', debug = 1, size_scale = 50, rest_same = True, sym_all = '.c', ms_all = 1, colorset = None, edgecolorset = None, dith = None, minscan = 0, maxscan = 40, modes = 'polar', get_wkin = True, no_hash = True, exception = Exception, output_log = 'out_cas3d', sharex = False):
        tot_plot = 0
        fig, ax = pt.subplots(nrows = 2, sharex = sharex)
        ax, ax2 = ax
        self.plot_continua(ax, colorset = colorset, edgecolorset = edgecolorset, dith = dith, rest_same = rest_same, size_scale = size_scale, sym_all = sym_all, ms_all = ms_all)
        leg_props = FontProperties(size='small')
        try:       # try to order legend by the codes after "&" in the label
            #ax=pl.gca() 
            handles, labels = ax.get_legend_handles_labels()
            lab_temp = [lab+'&999' for lab in labels]  # add a low priority code 
                                                       # so that all have some code.
            lab_text = np.array([ltemp.split('&')[0] for ltemp in lab_temp])
            lab_order = [ltemp.split('&')[1] for ltemp in lab_temp]
            lab_inds = np.argsort(lab_order)           # sort by increasing code
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
            leg_sh = ax.legend(np.array(handles)[lab_inds], lab_text[lab_inds],prop=leg_props, fancybox=True,loc='center left', bbox_to_anchor=(1, 0.5)) #SH
            box = ax2.get_position()
            ax2.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        except None, details:    
            print 'hello'
        ax.set_title(path.replace('_','\_'))
        ax.set_xlabel('s (flux)')
        ax.set_ylabel('Frequency (kHz)')
        am = ' '.join([ str(a) for a in self.allmodes])
        fig.suptitle('Total of {t} instances, {m} Xis modes used, {a} Alfvenic, {s}, selected indices: {am}\n{mti}'.
                    format(t=len(self.scan_arr[0]), a = len(np.where(self.scan_arr[self.A_sw]==1)[0]),
                           s=tot_plot, m=len(self.allmodes), am = am, 
                           mti=self.mode_table_info),fontsize='small')
        self.star = None
        norm = np.max(self.scan_arr[3,:])
        def onclick(event):
            print 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(
                event.button, event.x, event.y, event.xdata, event.ydata)
            min_loc = np.argmin((self.scan_arr[3,:]/norm - event.ydata/norm)**2 + (self.scan_arr[2,:] - event.xdata)**2)
            #min_loc = np.argmin(((scan_arr[3,:] - event.y)/np.max(scan_arr[3,:]))**2 + (scan_arr[2,:] - event.x)**2)
            print self.mode_shapes_all.shape, self.scan_arr.shape, self.scan_arr[2,min_loc], self.scan_arr[3,min_loc], min_loc
            tmp_lut = self.XIS.lut[self.scan_arr[4,min_loc]]

            if self.new_method:
                mode_str = ','.join(['({},{})'.format(self.XIS.lut[self.mode_nums[i,min_loc]]['n'], self.XIS.lut[self.mode_nums[i,min_loc]]['m']) for i in range(self.mode_nums.shape[0])])
            else:
                mode_str = ''
            print mode_str
            if event.inaxes==ax and not fig.canvas.manager.toolbar._active:
                mode_n, mode_m  = [tmp_lut['n'], tmp_lut['m']]
                en = np.sum(np.abs(self.scan_arr[self.A_cpt:self.A_cpt+3,min_loc]))
                title = 'freq:{}, ew:{:.4e}, n={} , m={}, {}, Shear:{:.1f},Sound:{:.1f},Fast:{:.1f},{:.3f}'.format(self.scan_arr[3,min_loc], self.scan_arr[12,min_loc], mode_n, mode_m, mode_str, self.scan_arr[self.A_cpt,min_loc]/en*100, self.scan_arr[self.A_cpt+1,min_loc]/en*100, self.scan_arr[self.A_cpt+2,min_loc]/en*100,self.mratio[min_loc])
                if self.new_method:
                    print min_loc,self.scan_arr[4,min_loc], self.mode_nums[0,min_loc], self.mode_nums[1,min_loc], self.mode_nums[2,min_loc]
                ax2.cla()
                #fig2.clf()
                if self.star!=None: self.star[0].remove()
                self.star = ax.plot(self.scan_arr[2,min_loc], self.scan_arr[3,min_loc],'*', markersize = 20)
                for i in range(self.mode_shapes_all.shape[0]):
                    ax2.plot(np.linspace(0,1,self.mode_shapes_all.shape[2]),self.mode_shapes_all[i,min_loc,:], zorder = 200)
                ax2.set_title(title)
                ax2.set_xlabel('xis')
                ax2.set_ylabel('s')
                #ax2.draw_artist(self.star[0])
                #ax2.draw_artist(ax2.patch)
                #ax2.draw_artist(ax2.xaxis)
                #ax2.draw_artist(ax2.yaxis)
                #fig.canvas.update()
                fig.canvas.draw(); fig.show()
        #self.scan_arr = scan_arr

        cid = fig.canvas.mpl_connect('button_press_event', onclick)
        #fig.canvas.mpl_connect('pick_event',onpick)
        ax.grid(); ax.set_xlim([0,1]);ax.set_ylim([np.min(self.scan_arr[3,:]),np.max(self.scan_arr[3,:])])
        fig.canvas.draw(); fig.show()


    def pub_plot(self, size_scale = 50, rest_same = True, sym_all = '.c', ms_all = 1, colorset = None, edgecolorset = None, dith = None, minscan = 0, maxscan = 40, modes = 'polar', get_wkin = True, no_hash = True, exception = Exception, output_log = 'out_cas3d', figname = None):
        tot_plot = 0
        fig, ax = pt.subplots(nrows = 2, sharex = True); ax,ax2 = ax
        gen_funcs.setup_publication_image(fig, height_prop = 1.5)
        self.plot_continua(ax, colorset = colorset, edgecolorset = edgecolorset, dith = dith, rest_same = rest_same, size_scale = size_scale, sym_all = sym_all, ms_all = ms_all, label_type = 'small', linewidths = 0.5)
        self.plot_continua(ax2, colorset = colorset, edgecolorset = edgecolorset, dith = dith, rest_same = rest_same, size_scale = size_scale, sym_all = sym_all, ms_all = ms_all, label_type = 'small', linewidths = 0.5)
        ax.grid(); ax.set_xlim([0,1]);ax.set_ylim([np.min(self.scan_arr[3,:]),np.max(self.scan_arr[3,:])])
        ax2.grid(); ax2.set_xlim([0,1]);ax2.set_ylim([0,50])
        ax.set_ylabel('Frequency (kHz)')
        ax2.set_ylabel('Frequency (kHz)')
        ax2.set_xlabel('s')
        leg_props = FontProperties(size='small')
        try:
            handles, labels = ax.get_legend_handles_labels()
            lab_temp = [lab+'&999' for lab in labels]  # add a low priority code 
                                                       # so that all have some code.
            lab_text = np.array([ltemp.split('&')[0] for ltemp in lab_temp])
            lab_order = [ltemp.split('&')[1] for ltemp in lab_temp]
            lab_inds = np.argsort(lab_order)           # sort by increasing code
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
            leg_sh = ax.legend(np.array(handles)[lab_inds], lab_text[lab_inds],prop=leg_props, fancybox=True,loc='center left', bbox_to_anchor=(1, 0.5)) #SH
            box = ax2.get_position()
            ax2.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        except None, details:    
            print 'hello'
        #fig.tight_layout()
        if figname!=None: 
            for i in ['.svg','.pdf','.eps']: fig.savefig(figname + i, bbox_extra_artists=(leg_sh,), bbox_inches='tight')
        fig.canvas.draw(); fig.show()





class plot_continua():
    def __init__(self, path = '/home/jbb105/short/030_n1_low_tinier/', suff = '', debug = 1, size_scale = 50, rest_same = True, sym_all = '.c', ms_all = 1, colorset = None, edgecolorset = None, dith = None, minscan = 0, maxscan = 40, modes = 'polar', get_wkin = True, no_hash = True, exception = Exception, output_log = 'out_cas3d'):

        offs = 2 # allow for inserting both hash flags ahead of the rest. 
        m_ind = 2 + offs # the column containing the mode index
        A_sw=5 + offs     # column for the Aflven/sound flag (1/0)
        A_cpt=A_sw+1 
        s_ind = 0 + offs
        f_ind = s_ind+1
        hash_ind = 0
        if colorset == None: colorset = ('b,g,r,c,m,y,k,orange,purple,lightgreen,gray'.split(','))
        if edgecolorset == None: edgecolorset='k,b,darkgreen,r,purple,green,darkgray,darkblue,c,m'.split(',')
        if dith==None: dith = [0.001, .1]
        edgecolors = []
        colors = 10*colorset
        for i in range(10): 
            for j in range(len(colorset)): 
                edgecolors.append(edgecolorset[i])
        fig, ax = pt.subplots(nrows = 2)
        ax, ax2 = ax
        #pl.clf()
        if path[-1] != '/': path = path + '/'  # add a trailing / if missing

        if suff == '.bz2': 
            import bz2  # this import will only be attempted if suff='.bz2'
            OPEN=bz2.BZ2File
        elif suff == '.gz': 
            import gzip
            OPEN=gzip.GzipFile
        else: OPEN = open

        outfn = path+output_log+suff
        try:
            outcas3d = OPEN(outfn,'r')
            buff=outcas3d.readlines(100000)   # arg is # characters to read, approx.
        except IOError, details:
            print('**Warning: out_cas3d file {o} not found - cannot label modes'.format(o=outfn))
        XIS=ModeTable(buff=buff,target=[' xis\n','eta'])
        mode_table_info = XIS.info

        try:
            ETA=ModeTable(buff=XIS.rest,target=[' eta\n', 'mu'])
            MU=ModeTable(buff=ETA.rest,target=[' mu\n', ' \n'])
        except exception, details:
            print('Unable to read other mode tables {d}'.format(d=details))
            if verbose>0: print('{d}'.format(d=details))
        else: # here if successful
            for minfo in [ETA.info,MU.info]:
                if XIS.info.split(None,1)[1] == minfo.split(None,1)[1]:
                    minfo = '<as per XIS>'
                mode_table_info += ', '+minfo

        scan_arr = None
        all_lines=[]
        tmp1 = []
        tmp2 = []
        tmp3 = []
        for n in range(minscan,maxscan):
            scfn = path+'scan_{0}.dat'.format(n) + suff
            wkfn = path+'wkin_{0}.dat'.format(n) + suff
            xisn1 = path+'xis1_{0}.dat'.format(n) + suff
            xisn2 = path+'xis2_{0}.dat'.format(n) + suff
            xisn3 = path+'xis3_{0}.dat'.format(n) + suff

            if get_wkin:
                try:
                    with OPEN(scfn,"r") as sf: 
                        scbuff = sf.readlines()
                    with OPEN(wkfn,"r") as wf: 
                        wkbuff = wf.readlines()
                    with OPEN(xisn1,"r") as xf: 
                        xfbuff = xf.readlines()
                        for i in xfbuff:
                            z = i.rstrip('\n').split(' ')
                            for j in z:
                                if j!='':
                                    if j.find('E')>=0:
                                        tmp1.append(float(j))
                                    else:
                                        tmp1.append(0)
                    with OPEN(xisn2,"r") as xf: 
                        xfbuff = xf.readlines()
                        for i in xfbuff:
                            z = i.rstrip('\n').split(' ')
                            for j in z:
                                if j!='':
                                    if j.find('E')>=0:
                                        tmp2.append(float(j))
                                    else:
                                        tmp2.append(0)
                    with OPEN(xisn3,"r") as xf: 
                        xfbuff = xf.readlines()
                        for i in xfbuff:
                            z = i.rstrip('\n').split(' ')
                            for j in z:
                                if j!='':
                                    if j.find('E')>=0:
                                        tmp3.append(float(j))
                                    else:
                                        tmp3.append(0)

                    for (ln, s_line) in enumerate(scbuff[2:]):
                        hashes=(s_line[0]+wkbuff[ln+1][0])
                        hashes_cvt = hashes.replace(' ','0 ').replace('#','1 ')
                        all_lines.append(hashes_cvt+   # # replaced by 1, space by 0
                                         s_line[1:-1]  # omit first # and last \n
                                         +wkbuff[ln+1][1:]) # omit just first # 
                except IOError, details: 
                    print('File read error - perhaps reduce maxscan to {n}: {d}'
                          .format(n=n, d=details))
                    break

            else:
                a = np.loadtxt(scfn).T
                #a = np.loadtxt(scfn, comments='!',usecols=range(1,10),skiprows=2).T
                #b = np.loadtxt(wkfn, comments='!',usecols=range(1,8),skiprows=1).T
                #zip_dim = len(np.shape(a))-1   # if 1D, zipper on the 0th dim, 2d, on the 1st
                #a = np.append(a,b,zip_dim)
                if len(np.shape(a)) == 1: a = np.array([a]) # make a 2D array if one line

                if scan_arr == None: scan_arr = a
                elif len(np.shape(a)) == 2 : scan_arr = np.append(scan_arr, a, 1)
                else: print('discarding {f}, shape is {s}'.format(f=scfn, s=np.shape(a)))

        #scan_arr - 0:scan #, 1:wkin#, 2:s, 3:w(kHz), 4:harm, 5:num_in_scan, 6:max xis, 7:type 0=sound
        #8:?, 9:?, 10:?, 11:s, 12: ew, 13:wkin, 14:wkin_norm, 15:wkin_binorm, 16:wkin_parallel, 17:sum
        #18:?, 19:?,20:?
        if get_wkin: scan_arr = np.loadtxt(StringIO(' '.join(all_lines))).T
        self.mode_shapes = np.abs(np.resize(tmp1, (scan_arr.shape[1],len(tmp1)/scan_arr.shape[1])))
        self.mode_shapes2 = np.abs(np.resize(tmp2, (scan_arr.shape[1],len(tmp1)/scan_arr.shape[1])))
        self.mode_shapes3 = np.abs(np.resize(tmp3, (scan_arr.shape[1],len(tmp1)/scan_arr.shape[1])))
        print len(tmp1), scan_arr.shape, len(tmp1)/scan_arr.shape[1], self.mode_shapes.shape

        if no_hash: 
            w_no_hash = np.where(scan_arr[hash_ind] == 0)[0]
            print('Ignoring {h} hashed of {t} instances'.
                  format(h=len(scan_arr[0])-len(w_no_hash), t=len(scan_arr[0])))
            scan_arr=scan_arr[:,w_no_hash]
            self.mode_shapes = self.mode_shapes[w_no_hash,:]
            self.mode_shapes2 = self.mode_shapes2[w_no_hash,:]
            self.mode_shapes3 = self.mode_shapes3[w_no_hash,:]
        mratio = scan_arr[A_cpt]/np.std(scan_arr[A_cpt:A_cpt+3].T,1)

        allmodes = np.unique(scan_arr[m_ind]).astype(int)
        wA = np.where(scan_arr[A_sw] == 1)[0]
        Alfvenicmodes = np.unique(scan_arr[m_ind][wA]).astype(int)
        wMM = np.where(mratio>0.01)[0]
        mmodes = np.unique(scan_arr[m_ind][wMM]).astype(int)

        def dither(dith=.001, arr=None):
            return dith*(np.random.random(len(arr))-0.5)

        if debug>0: print('Xis modes {m} found, {a} are classifed as '
                          'Alfvenic by "sound_test()"'.
                          format(m=allmodes,a=Alfvenicmodes))

        if pl.is_string_like(modes):
            if 'alf' in modes.lower(): modes = Alfvenicmodes
            elif 'pol' in modes.lower(): modes = mmodes
            else: 
                try:
                    modes = eval(modes)
                except:
                    raise ValueError,"modes = {m} is not understood".format(m=modes)

        # mode score will be used to allocate colors and position in legend
        mode_score = np.zeros(max([k for k in XIS.lut.keys()])+1, dtype=int)
        mode_score[mmodes] += 100  # give M-type modes the highest rank
        mode_score[Alfvenicmodes] += 50  # give Alf modes a high rank
        for k in XIS.lut.keys():    # bonus for small n,m
            mode_score[k] += -(abs(XIS.lut[k]['n']) + abs(XIS.lut[k]['m']))  
            # alternatively could rank by closeness to resonance?
        ct = np.arange(len(mode_score))
        ct[np.argsort(-mode_score)] = np.arange(len(mode_score))

        tot_plot = 0
        list_of_lines = [] #SH added
        for (m,mode)in enumerate(modes):
            wMode = np.where((scan_arr[m_ind][wMM]==mode))[0]
            tot_plot += len(wMode)
            if len(wMode) != 0:
                z=ax.scatter(scan_arr[s_ind][wMM[wMode]]+dither(dith[0],wMode),
                           scan_arr[f_ind][wMM[wMode]]+dither(dith[1],wMode),
                           size_scale*np.sqrt(mratio[wMM[wMode]]),marker='o',
                           c = colors[ct[mode]],edgecolors=edgecolors[ct[mode]],
                           label='M. mode {mode} {N}/{M} ({m}:{num})&{ord:03d}'.
                           format(mode=mode,m=ct[mode],num=len(wMode),
                                  M=XIS.lut[mode]['m'],N=XIS.lut[mode]['n'],
                                  ord=ct[mode] ))  
        # add an "order" code &123 at end.
                list_of_lines.append(z) #SH added

        # the fine dot plots are dithered only in S
        if rest_same: 
            tmp, = ax.plot(scan_arr[s_ind]+dither(dith[0],scan_arr[0]), # all modes 
                       scan_arr[f_ind],sym_all,markersize=ms_all,
                       label='all ({n})'.format(n=len(scan_arr[0])))
            list_of_lines.append(tmp)
        else: 
            for (m,mode)in enumerate(allmodes): # 
                wM = np.where((scan_arr[m_ind]==mode))[0]
                tmp, = ax.plot(scan_arr[s_ind][wM],scan_arr[f_ind][wM],'.',
                               markersize=msall,color =colors[ct[mode]],
                               label='allmode {m}'.format(m=mode))
                list_of_lines.append(tmp)

        for (m,mode)in enumerate(Alfvenicmodes): # (tiny offset in s to separate)
            wM = np.where((scan_arr[m_ind]==mode) &(scan_arr[A_sw] == 1))[0]
            tmp,=ax.plot(scan_arr[s_ind][wM]+0.0001,scan_arr[f_ind][wM],'+',color =colors[ct[mode]],
                    label='Amode {m} ({n})'.format(m=mode,n=len(wM)))
            list_of_lines.append(tmp)
        leg_props = FontProperties(size='small')

        try:       # try to order legend by the codes after "&" in the label
            #ax=pl.gca() 
            handles, labels = ax.get_legend_handles_labels()
            lab_temp = [lab+'&999' for lab in labels]  # add a low priority code 
                                                       # so that all have some code.
            lab_text = np.array([ltemp.split('&')[0] for ltemp in lab_temp])
            lab_order = [ltemp.split('&')[1] for ltemp in lab_temp]
            lab_inds = np.argsort(lab_order)           # sort by increasing code
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
            leg_sh = ax.legend(np.array(handles)[lab_inds], lab_text[lab_inds],prop=leg_props, fancybox=True,loc='center left', bbox_to_anchor=(1, 0.5)) #SH
            box = ax2.get_position()
            ax2.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        except None, details:    
            print 'hello'
            if verbose: print('label ordering failed - {d}'.format(details))
            leg_sh = pl.legend(prop=leg_prop,fancybox=True) #SH


        ########################SH mod###################
        def onpick(event):
            print 'hello in onpick'
            # on the pick event, find the orig line corresponding to the
            # legend proxy line, and toggle the visibility
            legline = event.artist
            print legline
            vis = not legline.get_visible()

            origline = lined[legline]
            vis = not origline.get_visible()
            origline.set_visible(vis)
            # Change the alpha on the line in the legend so we can see what lines
            # have been toggled
            if vis:
                legline.set_alpha(1.0)
            else:
                legline.set_alpha(0.2)
            fig.canvas.draw()

        lined = dict()
        handles_tmp, labels_tmp = pl.gca().get_legend_handles_labels()
        for legline, origline in zip(handles_tmp, list_of_lines):
            legline.set_picker(5)  # 5 pts tolerance
            lined[legline] = origline


        ax.set_title(path.replace('_','-'))
        ax.set_xlabel('s (flux)')
        ax.set_ylabel('frequency')
        am = ' '.join([ str(a) for a in allmodes])
        fig.suptitle('Total of {t} instances, {m} Xis modes used, {a} Alfvenic, {s}, selected indices: {am}\n{mti}'.
                    format(t=len(scan_arr[0]), a = len(np.where(scan_arr[A_sw]==1)[0]),
                           s=tot_plot, m=len(allmodes), am = am, 
                           mti=mode_table_info),fontsize='small')

        #fig = pl.gcf()
        #fig2, ax2 = pt.subplots()
        self.star = None
        norm = np.max(scan_arr[3,:])
        def onclick(event):
            print 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(
                event.button, event.x, event.y, event.xdata, event.ydata)
            min_loc = np.argmin((scan_arr[3,:]/norm - event.ydata/norm)**2 + (scan_arr[2,:] - event.xdata)**2)
            #min_loc = np.argmin(((scan_arr[3,:] - event.y)/np.max(scan_arr[3,:]))**2 + (scan_arr[2,:] - event.x)**2)
            print self.mode_shapes.shape, scan_arr.shape, scan_arr[2,min_loc], scan_arr[3,min_loc], min_loc
            tmp_lut = XIS.lut[scan_arr[4,min_loc]]
            if event.inaxes==ax and not fig.canvas.manager.toolbar._active:
                mode_n, mode_m  = [tmp_lut['n'], tmp_lut['m']]
                title = 'freq:{}, ew:{:.4e}, n={} , m={}'.format(scan_arr[3,min_loc], scan_arr[12,min_loc], mode_n, mode_m)
                ax2.cla()
                #fig2.clf()
                if self.star!=None: self.star[0].remove()
                self.star = ax.plot(scan_arr[2,min_loc], scan_arr[3,min_loc],'*', markersize = 20)
                ax2.plot(np.linspace(0,1,self.mode_shapes.shape[1]),self.mode_shapes[min_loc,:])
                ax2.plot(np.linspace(0,1,self.mode_shapes.shape[1]),self.mode_shapes2[min_loc,:])
                ax2.plot(np.linspace(0,1,self.mode_shapes.shape[1]),self.mode_shapes3[min_loc,:])
                ax2.set_title(title)
                #ax2.draw_artist(self.star[0])
                #ax2.draw_artist(ax2.patch)
                #ax2.draw_artist(ax2.xaxis)
                #ax2.draw_artist(ax2.yaxis)
                #fig.canvas.update()
                fig.canvas.draw(); fig.show()
        self.scan_arr = scan_arr
        cid = fig.canvas.mpl_connect('button_press_event', onclick)
        #fig.canvas.mpl_connect('pick_event',onpick)
        ax.grid(); ax.set_xlim([0,1]);ax.set_ylim([np.min(scan_arr[3,:]),np.max(scan_arr[3,:])])
        fig.canvas.draw(); fig.show()
        #pl.grid()
        #pl.show()


def plot_mode_shapes(path = '/home/srh112/code/python/cas3d_playground/input.kh0.850-kv1.000fixed_dir/cas3d_r_n1_free/', output_log = 'out_cas3d', std_dev_cut = 0.28, allowed_n = None, allowed_m = None):
    '''Meant to help find the interesting modes
    '''
    minscan = 1
    maxscan = 32
    #path = '/home/srh112/code/python/cas3d_playground/030_n1_full_full_flat_dens_100_surfaces/test_new_script/'
    suff = ''
    get_wkin = 1
    OPEN=open
    all_lines = []
    tmp = []
    limit_m = None

    #outfn = 'out_cas3d'
    with file(path + '/'+output_log,'r') as out_cas3d:
        buff=out_cas3d.readlines(100000)   # arg is # characters to read, approx.

    XIS=ModeTable(buff=buff,target=[' xis\n','eta'])
    for n in range(minscan,maxscan):
        scfn = path+'scan_{0}.dat'.format(n) + suff
        wkfn = path+'wkin_{0}.dat'.format(n) + suff
        xisn = path+'xis1_{0}.dat'.format(n) + suff
        if get_wkin:
            try:
                with OPEN(scfn,"r") as sf: 
                    scbuff = sf.readlines()
                with OPEN(wkfn,"r") as wf: 
                    wkbuff = wf.readlines()
                with OPEN(xisn,"r") as xf: 
                    xfbuff = xf.readlines()
                    for i in xfbuff:
                        z = i.rstrip('\n').split(' ')
                        for j in z:
                            if j!='':tmp.append(float(j))

                for (ln, s_line) in enumerate(scbuff[2:]):
                    hashes=(s_line[0]+wkbuff[ln+1][0])
                    hashes_cvt = hashes.replace(' ','0 ').replace('#','1 ')
                    all_lines.append(hashes_cvt+   # # replaced by 1, space by 0
                                     s_line[1:-1]  # omit first # and last \n
                                     +wkbuff[ln+1][1:]) # omit just first # 
            except IOError, details: 
                print('File read error - perhaps reduce maxscan to {n}: {d}'
                      .format(n=n, d=details))
                break

        else:
            a = np.loadtxt(scfn).T
            #a = np.loadtxt(scfn, comments='!',usecols=range(1,10),skiprows=2).T
            #b = np.loadtxt(wkfn, comments='!',usecols=range(1,8),skiprows=1).T
            #zip_dim = len(np.shape(a))-1   # if 1D, zipper on the 0th dim, 2d, on the 1st
            #a = np.append(a,b,zip_dim)
            if len(np.shape(a)) == 1: a = np.array([a]) # make a 2D array if one line

            if scan_arr == None: scan_arr = a
            elif len(np.shape(a)) == 2 : scan_arr = np.append(scan_arr, a, 1)
            else: print('discarding {f}, shape is {s}'.format(f=scfn, s=np.shape(a)))
    print all_lines[-1]
    if get_wkin: scan_arr = np.loadtxt(StringIO(' '.join(all_lines))).T
    scan_arr.shape[1]
    print len(tmp), scan_arr.shape, len(tmp)/scan_arr.shape[1]

    def unique(a):
        order = np.lexsort(a.T)
        a = a[order]
        diff = np.diff(a, axis=0)
        ui = np.ones(len(a), 'bool')
        ui[1:] = (np.abs(diff) > 0.01).any(axis=1) 
        return a[ui]

    mode_shapes = np.abs(np.resize(tmp, (scan_arr.shape[1],len(tmp)/scan_arr.shape[1])))
    print mode_shapes.shape, scan_arr.shape
    mode_shapes_copy = copy.deepcopy(mode_shapes)
    mode_shapes = unique(mode_shapes)
    std_dev_list = []
    fig_shapes, ax_shapes = pt.subplots(nrows=2,sharex = True)
    count = 0
    good = np.zeros(mode_shapes_copy.shape[0])
    #allowed_n = [-5]; allowed_m = [4]; 
    all_allowed_n = True if allowed_n == None else False
    all_allowed_m = True if allowed_m == None else False
    for i in range(0,mode_shapes.shape[0],2):
        if ((i%100)==0):print i
        tmp_data = np.abs(mode_shapes[i,:]).flatten()

        #find original data spot....
        max_val = np.max(tmp_data)
        mean_val = np.mean(tmp_data)
        std_dev = np.std(tmp_data/max_val)
        std_dev_list.append(std_dev)


        if std_dev>std_dev_cut:
            s_axis = np.linspace(0,1,len(tmp_data))
            orig_index = np.argmin(np.sum(np.abs(mode_shapes_copy - tmp_data),axis=1))

            tmp_lut = XIS.lut[scan_arr[4,orig_index]]
            mode_n, mode_m  = [tmp_lut['n'], tmp_lut['m']]
            if (mode_n in allowed_n or all_allowed_n) and (mode_m in allowed_m or all_allowed_m):
                mode_details = '{},{}'.format(mode_n, mode_m)
                if limit_m != None:
                    if mode_m!=limit_m:
                        plot_mode = False
                    else:
                        plot_mode = True
                else:
                    plot_mode = True
                if plot_mode:
                    ax_shapes[0].plot(s_axis,tmp_data)
                    max_loc = np.argmax(tmp_data)
                    #print scan_arr[:,orig_index]
                    #ax_shapes[0].text(s_axis[max_loc],np.max(tmp_data),'{:.2f}kHz,{},{}'.format(scan_arr[3,orig_index], int(scan_arr[4,orig_index]), mode_details))
                    print scan_arr[12,orig_index]
                    ax_shapes[0].text(s_axis[max_loc],np.max(tmp_data),'{:.2f}kHz,{},{}'.format(1000.*scan_arr[12,orig_index], int(scan_arr[4,orig_index]), mode_details))
                    count+=1
                    good[orig_index] = 1

    #fig, ax = pt.subplots()
    good = np.array(good)
    ax_shapes[1].plot(scan_arr[2,good==0], scan_arr[3,good==0],'.',markersize=4)
    ax_shapes[1].plot(scan_arr[2,good==1], scan_arr[3,good==1],'o',markersize=7)
    ax_shapes[1].set_xlabel('s')
    ax_shapes[1].set_ylabel('freq (kHz)')
    ax_shapes[0].set_ylabel('CAS3D eig func disp')
    #fig.canvas.draw(), fig.show()
    fig_shapes.suptitle('first attempt global mode finder')
    fig,ax = pt.subplots()
    ax.hist(std_dev_list)
    fig.canvas.draw(), fig.show()
    fig_shapes.canvas.draw(); fig_shapes.show()
    print count


def pub_plot(size_scale = 50, rest_same = True, sym_all = '.c', ms_all = 1, colorset = None, edgecolorset = None, dith = None, minscan = 0, maxscan = 40, modes = 'polar', get_wkin = True, no_hash = True, exception = Exception, output_log = 'out_cas3d', figname = None, kh_list = None, inc_iota = False, ylims = None, inc_legend = True, linewidths = 0.5, ellipses = None, ax = None, fig = None, direct_format = None, txt_details = None, plot_styles = None, adjust_for_leg = True, plot_a_b = False):
    #tot_plot = 0
    n_plots = len(kh_list)+inc_iota
    if direct_format == None: direct_format = r'/home/srh112/SSD_mount/raijin_short/whale_tail/input.kh0.{}0-kv1.000fixed_dir/cas3d_r_n1_free_lin/'
    if ax == None and fig == None:
        print 'creating a new figure'
        fig, ax = pt.subplots(nrows = n_plots, sharex = True); 
        gen_funcs.setup_publication_image(fig, height_prop = n_plots*0.65)
        sup_fig = False
    else:
        sup_fig = True
    if n_plots == 1: ax = [ax]
    if ylims == None: ylims = [[0,100] for i in kh_list]
    if ellipses == None: ellipses = [[] for i in kh_list]
    if inc_legend==True: inc_legend = [True]*len(kh_list)
    try:
        tmp = int(size_scale)
        size_scale = [size_scale]*len(kh_list)
    except (ValueError, TypeError) as e:
        pass
    try:
        tmp = float(linewidths)
        linewidths = [linewidths]*len(kh_list)
    except (ValueError, TypeError) as e:
        pass
    for i, (ax_tmp, kh, ylim, size_scale_tmp, linewidth, tmp_ellipse) in enumerate(zip(ax, kh_list, ylims, size_scale, linewidths, ellipses)):
        directory = '/home/srh112/SSD_mount/raijin_short/whale_tail/input.kh0.{}0-kv1.000fixed_dir/cas3d_r_n1_fixed/'.format(kh)
        directory = direct_format.format(kh)
        #directory = '/home/srh112/SSD_mount/raijin_short/whale_tail/input.kh0.{}0-kv1.000fixed_dir/cas3d_r_n1_tiniest_free_lin/'.format(kh)
        cont = cas3d_continua_data(path = directory, maxscan = 32,output_log='out_cas3d_r', load_eig_funcs = False)
        cont.plot_continua(ax_tmp, colorset = colorset, edgecolorset = edgecolorset, dith = dith, rest_same = rest_same, size_scale = size_scale_tmp, sym_all = sym_all, ms_all = ms_all, label_type = 'small', linewidths = linewidth, ylims = ylim, plot_crosses = False, plot_styles = plot_styles)
        plot_styles = cont.plot_styles
        ax_tmp.grid(); ax_tmp.set_xlim([0,1])
        for ii in tmp_ellipse:
            ii['ax']=ax_tmp
            ellipse(**ii)
        ax_tmp.set_ylim(ylim)#ax.set_ylim([np.min(self.scan_arr[3,:]),np.max(self.scan_arr[3,:])])
        if inc_iota:ax[-1].plot(cont.s, cont.iota,label = 'kh=0.{}'.format(kh))
        if i == 0:
            ax_tmp.set_title('$\kappa_H = '+ '0.{}$'.format(kh))
        else:
            if kh!=kh_list[i-1]:
                ax_tmp.set_title('$\kappa_H = '+ '0.{}$'.format(kh))
        if txt_details!=None:
            if txt_details['axes'][i] == True:
                for txt, loc in zip(txt_details['txts'],txt_details['locs']):ax_tmp.text(loc[0],loc[1],txt, **txt_details['txt_kwargs'])
    if inc_iota:
        for n,m in zip([5,4], [4,3]):
            ax[-1].axhline(float(n)/m)
            ax[-1].text(0.01,float(n)/m,'{}/{}'.format(n,m))
    for i in ax[:-1]: i.set_ylabel('Frequency (kHz)')
    if inc_iota:
        ax[-1].set_ylabel('rot transform')
    else:
        ax[-1].set_ylabel('Frequency (kHz)')
    ax[-1].set_xlabel('s')
    leg_props = FontProperties(size='small')
    for ax_tmp, leg_tmp in zip(ax, inc_legend):
        handles, labels = ax_tmp.get_legend_handles_labels()
        lab_temp = [lab+'&999' for lab in labels]  # add a low priority code 
                                                   # so that all have some code.
        lab_text = np.array([ltemp.split('&')[0] for ltemp in lab_temp])
        lab_order = [ltemp.split('&')[1] for ltemp in lab_temp]
        lab_inds = np.argsort(lab_order)           # sort by increasing code
        box = ax_tmp.get_position()
        if adjust_for_leg: ax_tmp.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        if leg_tmp: leg_sh = ax_tmp.legend(np.array(handles)[lab_inds], lab_text[lab_inds],prop=leg_props, fancybox=True,loc='center left', bbox_to_anchor=(1, 0.5))
    #box = ax[-1].get_position()
    #ax[-1].set_position([box.x0, box.y0, box.width * 0.8, box.height])
    if plot_a_b :gen_funcs.a_b_c(ax, x_locs = 0.3, y_locs = 0.85, kwargs = {'backgroundcolor':'white','zorder':1000})
    if figname!=None: 
        for i in ['.svg','.pdf']: fig.savefig(figname + i, bbox_extra_artists=(leg_sh,), bbox_inches='tight')
    if sup_fig==False:
        fig.canvas.draw(); fig.show()
    return plot_styles

def ellipse(w_radius=1, h_radius=1, c_x=0, c_y=0, n_pts=100, plot = False, ax = None, plot_kwargs = None):
    theta = np.linspace(-np.pi,np.pi,n_pts)
    if plot_kwargs == None: plot_kwargs = {}
    x = w_radius*np.cos(theta) + c_x
    y = h_radius*np.sin(theta) + c_y
    if plot:
        if ax != None:
            ax.plot(x,y, **plot_kwargs)
        else:
            fig, ax = pt.subplots(nrows = 1, ncols = 1, sharex = False, sharey = False)
            ax.plot(x,y, **plot_kwargs)
            #ax.set_xlim(-5,5)
            #ax.set_ylim(-5,5)
            fig.canvas.draw();fig.show()
