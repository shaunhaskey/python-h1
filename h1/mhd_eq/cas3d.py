import numpy as np
import matplotlib.pyplot as pt
import h1.helper.generic_funcs as gen_funcs
import copy
import pylab as pl
from matplotlib.font_manager import FontProperties
from StringIO import StringIO


class cas3d_results():
    def __init__(self, directory, single_eig = True, output_filename = 'out_cas3d'):
        self.directory = directory
        self.single_eig = single_eig
        self.output_filename = output_filename
    def plot_eigenvalues(self, ax = None, ylims = None, max_harms = None, sqrt_s = False, multiplier = 1):
        with file(self.directory + '/'+ self.output_filename, 'r') as file_handle:
            file_lines = file_handle.readlines()
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
            identifier, m_n_list, max_value = self.return_harmonics(file_lines, harmonic_strings[i])
            svals, coefficient_array = self.return_coeffs(file_lines, coeff_strings[i])
            max_ind = np.argmin(np.abs(svals - 0.95))
            if max_harms == None: max_harms = len(m_n_list) - 1
            for j in range(0,max_harms+1):
                xax = svals if not sqrt_s else np.sqrt(svals)
                ax[i].plot(xax, coefficient_array[:,j]*multiplier,'k-', label = '{}'.format(m_n_list[j]))
                max_loc = np.argmax(np.abs(coefficient_array[:max_ind,j]))
                max_val = np.abs(np.max(coefficient_array[:max_ind,j]))
                if max_val>max_val_plot: max_val_plot = max_val
                text_string = m_n_list[j]
                #ax[i].text(xax[max_loc], coefficient_array[max_loc,j]*multiplier, '{}'.format('{},{}'.format(-text_string[1], -text_string[0])))


            ax[i].grid()
            #ax[i].legend(fontsize = 7, loc = 'best')
            if ylims[i]!=None: ax[i].set_ylim(ylims[i])
            print max_val_plot
            ax[i].set_ylim([-max_val_plot*1.05, max_val_plot*1.05])
            ax[i].set_ylabel(titles[i])
            self.identifiers.append(identifier); self.m_n_list.append(m_n_list); self.max_value.append(max_value)
        #ax[i].set_ylim([-0.0008,0.0008])
        ax[i].set_xlim([0,1])
        ax[i].set_xlabel('s')
        if own_plots: 
            fig.suptitle(self.directory.replace('_','\_'))
            fig.canvas.draw(); fig.show()

        #direct read in of an xis file
        tmp_data = np.loadtxt(self.directory + '/xis_largest01-10.dat')
        fig, ax = pt.subplots()
        for j in range(1,tmp_data.shape[1]):
            xax = tmp_data[:,0] if not sqrt_s else np.sqrt(tmp_data[:,0])
            xlabel = 's' if not sqrt_s else '$\sqrt{s}$'
            #ax.plot(tmp_data[:,0], tmp_data[:,j],'.-')
            ax.plot(xax, tmp_data[:,j],'.-')
        ax.set_xlabel(xlabel)
        #ax.set_ylabel(r'$\xi^s$')
        ax.set_ylabel(r'Normal displacement (a.u)')
        fig.suptitle(self.directory.replace('_','\_'))
        fig.canvas.draw();fig.show()
        
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

    def plot_ne_Te_P(self,ne = True, Te = False, P = True, fname_ne = 'perturbed_density.dat', fname_Te = 'perturbed_temperature.dat', fname_P = 'perturbed_pressure.dat', max_harms = None, ax = None, ylims = None, divide = False, sqrt_s = False, multiplier = 1):
        #fname_ne = 'perturbed_ne_f.dat'
        #fname_Te = 'perturbed_pressure_f.dat'
        #fname_P = 'perturbed_ne_f.dat'
        fnames = []; yax = []
        divisors = []
        pert_quants = []
        if ne: 
            fnames.append(fname_ne)
            yax.append('Perturbed ne (a.u)')
            divisors.append(self.rho)
            pert_quants.append('ne1')
        if Te: 
            fnames.append(fname_Te)
            yax.append('Perturbed Te (a.u)')
            divisors.append(self.Te)
            pert_quants.append('Te1')
        if P: 
            fnames.append(fname_P)
            yax.append('Perturbed P (a.u)')
            divisors.append(self.P)
            pert_quants.append('P1')
        #yax = ['Perturbed Pressure', 'Perturbed Density']
        own_plots = True if ax == None else False
        if ylims == None: ylims = [None] * len(yax)
        if own_plots: fig, ax = pt.subplots(nrows = len(fnames), sharex = True)
        for fname, yname, ax_cur, ylim, div, pert_name in zip(fnames, yax, ax, ylims, divisors, pert_quants):
            max_val_plot = 0.
            a = np.loadtxt(self.directory + '/'+fname, comments='#')
            self.s_cur = a[:,0]
            max_ind = np.argmin(np.abs(self.s_cur - 0.95))
            if max_harms == None: max_harms = a.shape[1]-1
            setattr(self, pert_name, a[:,1:])
            for i in range(1, max_harms+1):
                try:
                    x_ax = np.sqrt(self.s_cur) if sqrt_s else self.s_cur
                    self.identifiers[0].tolist().index(i)
                    plot_val = a[:,i] if not divide else a[:,i]/div
                    ax_cur.plot(x_ax, plot_val*multiplier, 'b-')
                    max_val = np.max(np.abs(plot_val[:max_ind]))*multiplier
                    max_loc = np.argmax(np.abs(plot_val[:max_ind]))
                    if max_val>max_val_plot: max_val_plot = max_val
                    text_string = str(i)
                    text_string = self.m_n_list[0][self.identifiers[0].tolist().index(i)]
                    #ax_cur.text(s[max_loc], plot_val[max_loc]*multiplier, '{},{}'.format(-text_string[1], -text_string[0]))
                    #ax_cur.text(s[max_loc], a[max_loc,i], str(i))
                    #print 'errors....', 
                except ValueError:
                    print 'errors....'
                    pass

            #ax_cur.set_xlabel('s')
            ax_cur.set_ylabel(yname)
            print 'hello'
            #if ylim!=None: ax_cur.set_ylim(ylim)
            ax_cur.set_ylim([-max_val_plot*1.05, max_val_plot*1.05])
            ax_cur.grid(True)
        if own_plots: fig.canvas.draw(); fig.show()

    def eq_quants_pnt(self, fname = 'equilibrium_profiles_pnt.dat', ax = None,):
        print self.directory + '/equilibrium_profiles_pnt.dat'
        data = np.loadtxt(self.directory + '/equilibrium_profiles_pnt.dat', comments='#')
        labels = ['P', 'dP/ds','rho','drho/ds','Te','dTe/ds']
        names = ['P', 'dPds','rho','drhods','Te','dTeds']
        own_ax = True if ax == None else False
        if own_ax: 
            fig, ax = pt.subplots(ncols = 2, nrows = 3, sharex = True)
            gen_funcs.setup_publication_image(fig, height_prop = 1./1.618, single_col = False, replacement_kwargs = None)
        for i, (cur_ax, label, name) in enumerate(zip(ax.flatten(), labels, names)):
            cur_ax.plot(data[:,0], data[:,i+1])
            setattr(self,name,data[:,i+1])
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
        #for i in range(10): edgecolors.extend(len(colorset)*edgecolorset[i])

        pl.clf()
        # try:
        #     from warnings import warn  # this process_cmd is the latest and greatest (no pyfusion needed)
        #     #from bdb_utils import process_cmd_line_args
        #     #exec(process_cmd_line_args())
        # except None:
        #     warn('Importing process_cmd_line_args allows variables to be changed from the command line')

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

        #XIS=ModeTable(buff=buff,target=['User supplied perturbation spectrum','eta'])
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
                #a = np.loadtxt(scfn, comments='!',usecols=range(1,10),skiprows=2).T
                #b = np.loadtxt(wkfn, comments='!',usecols=range(1,8),skiprows=1).T
                #zip_dim = len(np.shape(a))-1   # if 1D, zipper on the 0th dim, 2d, on the 1st
                #a = np.append(a,b,zip_dim)
                if len(np.shape(a)) == 1: a = np.array([a]) # make a 2D array if one line

                if scan_arr == None: scan_arr = a
                elif len(np.shape(a)) == 2 : scan_arr = np.append(scan_arr, a, 1)
                else: print('discarding {f}, shape is {s}'.format(f=scfn, s=np.shape(a)))

        if get_wkin: scan_arr = np.loadtxt(StringIO(' '.join(all_lines))).T

        if no_hash: 
            w_no_hash = np.where(scan_arr[hash_ind] == 0)[0]
            print('Ignoring {h} hashed of {t} instances'.
                  format(h=len(scan_arr[0])-len(w_no_hash), t=len(scan_arr[0])))
            scan_arr=scan_arr[:,w_no_hash]

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
                z=pl.scatter(scan_arr[s_ind][wMM[wMode]]+dither(dith[0],wMode),
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
            tmp, = pl.plot(scan_arr[s_ind]+dither(dith[0],scan_arr[0]), # all modes 
                       scan_arr[f_ind],sym_all,markersize=ms_all,
                       label='all ({n})'.format(n=len(scan_arr[0])))
            list_of_lines.append(tmp)
        else: 
            for (m,mode)in enumerate(allmodes): # 
                wM = np.where((scan_arr[m_ind]==mode))[0]
                tmp, = pl.plot(scan_arr[s_ind][wM],scan_arr[f_ind][wM],'.',
                               markersize=msall,color =colors[ct[mode]],
                               label='allmode {m}'.format(m=mode))
                list_of_lines.append(tmp)

        for (m,mode)in enumerate(Alfvenicmodes): # (tiny offset in s to separate)
            wM = np.where((scan_arr[m_ind]==mode) &(scan_arr[A_sw] == 1))[0]
            tmp,=pl.plot(scan_arr[s_ind][wM]+0.0001,scan_arr[f_ind][wM],'+',color =colors[ct[mode]],
                    label='Amode {m} ({n})'.format(m=mode,n=len(wM)))
            list_of_lines.append(tmp)
        leg_props = FontProperties(size='small')

        try:       # try to order legend by the codes after "&" in the label
            ax=pl.gca() 
            handles, labels = ax.get_legend_handles_labels()
            lab_temp = [lab+'&999' for lab in labels]  # add a low priority code 
                                                       # so that all have some code.
            lab_text = np.array([ltemp.split('&')[0] for ltemp in lab_temp])
            lab_order = [ltemp.split('&')[1] for ltemp in lab_temp]
            lab_inds = np.argsort(lab_order)           # sort by increasing code
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
            leg_sh = ax.legend(np.array(handles)[lab_inds], lab_text[lab_inds],prop=leg_props, fancybox=True,loc='center left', bbox_to_anchor=(1, 0.5)) #SH
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


        pl.title(path.replace('_','-'))
        pl.xlabel('s (flux)')
        pl.ylabel('frequency')
        am = ' '.join([ str(a) for a in allmodes])
        pl.suptitle('Total of {t} instances, {m} Xis modes used, {a} Alfvenic, {s}, selected indices: {am}\n{mti}'.
                    format(t=len(scan_arr[0]), a = len(np.where(scan_arr[A_sw]==1)[0]),
                           s=tot_plot, m=len(allmodes), am = am, 
                           mti=mode_table_info),fontsize='small')

        #fig = pl.gcf()

        #fig.canvas.mpl_connect('pick_event',onpick)
        pl.grid()
        pl.show()



def plot_mode_shapes(path = '/home/srh112/code/python/cas3d_playground/input.kh0.850-kv1.000fixed_dir/cas3d_r_n1_free/', output_log = 'out_cas3d'):
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
    for i in range(0,mode_shapes.shape[0],2):
        if ((i%100)==0):print i
        tmp_data = np.abs(mode_shapes[i,:]).flatten()

        #find original data spot....
        max_val = np.max(tmp_data)
        mean_val = np.mean(tmp_data)
        std_dev = np.std(tmp_data/max_val)
        std_dev_list.append(std_dev)


        if std_dev>0.28:
            s_axis = np.linspace(0,1,len(tmp_data))
            orig_index = np.argmin(np.sum(np.abs(mode_shapes_copy - tmp_data),axis=1))

            tmp_lut = XIS.lut[scan_arr[4,orig_index]]
            mode_n, mode_m  = [tmp_lut['n'], tmp_lut['m']]
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
                ax_shapes[0].text(s_axis[max_loc],np.max(tmp_data),'{:.2}fkHz,{},{}'.format(scan_arr[3,orig_index], int(scan_arr[4,orig_index]), mode_details))
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
