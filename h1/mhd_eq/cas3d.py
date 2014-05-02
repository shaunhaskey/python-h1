import numpy as np
import matplotlib.pyplot as pt
import h1.helper.generic_funcs as gen_funcs

class cas3d_results():
    def __init__(self, directory, single_eig = True):
        self.directory = directory
        self.single_eig = single_eig

    def plot_eigenvalues(self, ax = None, ylims = None, max_harms = None):
        with file(self.directory + '/out_cas3d', 'r') as file_handle:
            file_lines = file_handle.readlines()
        harmonic_strings = ['largest xis harmonics', 'largest eta harmonics', 'largest mu  harmonics']
        coeff_strings = ['Fourier coefficients of     xis', 'Fourier coefficients of     eta','Fourier coefficients of      mu']
        titles = [r'$\xi^s$', r'$\xi^\eta$', r'$\xi^\mu$']
        own_plots = True if ax == None else False
        if own_plots: fig, ax = pt.subplots(nrows = len(harmonic_strings), sharex =True , sharey=False)
        self.identifiers = []; self.m_n_list = []; self.max_value = []
        if ylims == None: ylims = [None] * len(harmonic_strings)
        for i in range(0, len(harmonic_strings)):
            print i
            identifier, m_n_list, max_value = self.return_harmonics(file_lines, harmonic_strings[i])
            svals, coefficient_array = self.return_coeffs(file_lines, coeff_strings[i])
            if max_harms == None: max_harms = len(m_n_list) - 1
            for j in range(0,max_harms+1):
                ax[i].plot(svals, coefficient_array[:,j],'-', label = '{}'.format(m_n_list[j]))
                max_loc = np.argmax(np.abs(coefficient_array[:-100,j]))
                text_string = m_n_list[j]
                ax[i].text(svals[max_loc], coefficient_array[max_loc,j], '{}'.format('{},{}'.format(-text_string[1], -text_string[0])))


            ax[i].grid()
            #ax[i].legend(fontsize = 7, loc = 'best')
            if ylims[i]!=None: ax[i].set_ylim(ylims[i])
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
            ax.plot(tmp_data[:,0], tmp_data[:,j],'.-')
        fig.suptitle(self.directory.replace('_','\_'))
        fig.canvas.draw();fig.show()
        
    def return_harmonics(self, file_lines, harmonic_string):
        #find the start point of the interesting stuff
        found = 0; i=0
        while found==0 and i<=len(file_lines):
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

    def plot_ne_Te_P(self,ne = True, Te = False, P = True, fname_ne = 'perturbed_ne_f.dat', fname_Te = 'perturbed_te_f.dat', fname_P = 'perturbed_pressure_f.dat', max_harms = None, ax = None, ylims = None):
        #fname_ne = 'perturbed_ne_f.dat'
        #fname_Te = 'perturbed_pressure_f.dat'
        #fname_P = 'perturbed_ne_f.dat'
        fnames = []; yax = []
        if ne: fnames.append(fname_ne)
        if Te: fnames.append(fname_Te)
        if P: fnames.append(fname_P)
        if ne: yax.append('Perturbed ne')
        if Te: yax.append('Perturbed Te')
        if P: yax.append('Perturbed P')
        #yax = ['Perturbed Pressure', 'Perturbed Density']
        own_plots = True if ax == None else False
        if ylims == None: ylims = [None] * len(yax)
        if own_plots: fig, ax = pt.subplots(nrows = len(fnames), sharex = True)
        for fname, yname, ax_cur, ylim in zip(fnames, yax, ax, ylims):
            a = np.loadtxt(self.directory + '/'+fname, comments='#')
            s = a[:,0]
            if max_harms == None: max_harms = a.shape[1]-1
            for i in range(1, max_harms+1):
                try:
                    self.identifiers[0].tolist().index(i)
                    ax_cur.plot(s, a[:,i])
                    max_loc = np.argmax(np.abs(a[:-100,i]))
                    text_string = str(i)

                    text_string = self.m_n_list[0][self.identifiers[0].tolist().index(i)]
                    ax_cur.text(s[max_loc], a[max_loc,i], '{},{}'.format(-text_string[1], -text_string[0]))
                    #ax_cur.text(s[max_loc], a[max_loc,i], str(i))
                    print 'errors....', 
                except ValueError:
                    print 'errors....'
                    pass

            ax_cur.set_xlabel('s')
            ax_cur.set_ylabel(yname)
            print 'hello'
            if ylim!=None: ax_cur.set_ylim(ylim)
            ax_cur.grid(True)
        if own_plots: fig.canvas.draw(); fig.show()

directory = '/home/srh112/code/python/cas3d_playground/030_n1_full_full_free_wo_vac_dens_edge/eig_single_1265'
directory = '/home/srh112/code/python/cas3d_playground/replicate_jason_paper'


fig, ax = pt.subplots(nrows = 3, ncols = 2, sharex = True)
gen_funcs.setup_publication_image(fig, height_prop = 1./1.618, single_col = False, replacement_kwargs = None)
a = cas3d_results(directory, single_eig = True)
a.plot_eigenvalues(ax = ax[:,0], ylims = [[-1.4e-4,+1.4e-4], [-2.e-6,2.e-6], [-3e-7,3e-7]], max_harms = 5)
#a.plot_ne_Te_P(ne = True, Te = True, P = True, max_harms = None, ax = ax[:,1], ylims = [[-1.e15,1.e15], [-3.e-9,3.e-9], [-3.e-9,3.e-9]])
a.plot_ne_Te_P(ne = True, Te = False, P = True, max_harms = None, ax = ax[:,1], ylims = [[-1.e15,1.e15], [-3.e-9,3.e-9]])
fig.tight_layout()
fig.savefig('/home/srh112/win_shared/tomo_referee_response/test.png')
fig.canvas.draw(); fig.show()
