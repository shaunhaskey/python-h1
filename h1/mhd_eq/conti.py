from matplotlib.font_manager import FontProperties
import numpy as np
import matplotlib.pyplot as pt
import h1.helper.generic_funcs as gen_funcs
fname_list = []
kh_vals_list = range(25,40,1)
kh_vals_list.remove(36)

class conti_results():
    def __init__(self,directory, file_name,):
        self.file_name = file_name
        self.base_dir = directory

        with file(self.base_dir+self.file_name+'.dat','r') as dat_file_handle:
            lines = dat_file_handle.readlines()

        self.s_list = []; self.freq_list = []
        self.overal_s_list = []
        self.overal_freq_list = []
        for loc, tmp in enumerate(lines):
            line_split = tmp.split()
            if len(tmp.split())!=2:
                self.overal_s_list.append(self.s_list)
                self.overal_freq_list.append(self.freq_list)
                self.s_list = []
                self.freq_list = []
            else:
                self.s_list.append(float(line_split[0]))
                self.freq_list.append(float(line_split[1]))
        self.overal_s_list.pop(0)
        self.overal_freq_list.pop(0)

        #extract the mode numbers and the sound/alfven classification
        with file(self.base_dir+self.file_name+'.com','r') as com_file_handle:lines = com_file_handle.readlines()
        self.n_modes = []; self.m_modes=[];self.symbol = []
        for tmp in lines:
            if tmp.find('legend')>=0:
                start_point = tmp.find('"')
                tmp2 = tmp[start_point+1:-2].split()
                self.n_modes.append(int(tmp2[3]))
                self.m_modes.append(int(tmp2[1].rstrip(',')))
            if tmp.find('symbol'):
                if len(tmp.split())==3:
                    self.symbol.append(int(tmp.split()[-1]))

    def plot_sound_modes(self,allowed_n= None, allowed_m=None, ax = None, plot_props = None, sqrt_s = False):
        if allowed_n == None: allowed_n_list = self.n_modes
        if allowed_m == None: allowed_m_list = self.m_modes
        if ax == None: fig, ax = pt.subplots()
        #if plot_props == None: plot_props = {'color':(0.5,0.5,0.5), 'markersize':0.5,'marker':'.', 'linestyle':'None', 'rasterized':True, 'zorder':50}
        if plot_props == None: plot_props = {'markersize':0.5,'marker':'.', 'linestyle':'None', 'rasterized':True, 'zorder':50}
        print plot_props
        for i in range(1,len(self.overal_s_list)):
            if (self.symbol[i]==2) and (self.n_modes[i] in allowed_n_list) and (self.m_modes[i] in allowed_m_list):
                #ax.plot(self.overal_s_list[i], self.overal_freq_list[i],'.',color=(0.5,0.5,0.5),markersize=0.5)

                x_ax = np.sqrt(self.overal_s_list[i]) if sqrt_s else np.array(self.overal_s_list[i])
                y_ax = np.array(self.overal_freq_list[i])
                tmp = np.argsort(x_ax)
                plot_props['linestyle'] = 'None'
                ax.plot(x_ax, self.overal_freq_list[i],**plot_props)
                plot_props['linestyle'] = '-'
                ax.plot(x_ax[tmp], y_ax[tmp],**plot_props)

    def plot_alfven_modes(self,allowed_n= None, allowed_m=None, ax = None, plot_props = None, sqrt_s = False, ylim = None):
        if allowed_n == None: allowed_n_list = self.n_modes
        if allowed_m == None: allowed_m_list = self.m_modes
        if ax == None: fig, ax = pt.subplots()
        if plot_props == None: plot_props = {'markersize':2,'marker':'o', 'linestyle':'None', 'rasterized':True, 'zorder':100}
        print plot_props
        for i in range(1,len(self.overal_s_list)):
            if (self.symbol[i]==1) and (self.n_modes[i] in allowed_n_list) and (self.m_modes[i] in allowed_m_list):
                x_ax = np.sqrt(self.overal_s_list[i]) if sqrt_s else np.array(self.overal_s_list[i])
                y_ax = np.array(self.overal_freq_list[i])
                if ylim==None or (np.sum((y_ax>ylim[0]) * (y_ax<ylim[1]))>5):
                    #tmp, = ax.plot(x_ax, y_ax, label = 'm=%d,n=%d'%(self.m_modes[i],self.n_modes[i]),**plot_props)
                    tmp, = ax.plot(x_ax, y_ax, label = '({},{})'.format(-self.n_modes[i],-self.m_modes[i]),**plot_props)

def pub_plot(output_log = 'continuum_h1', figname = None, kh_list = None, inc_iota = False, ylims = None, inc_sound = True, alf_plot_props = None, inc_legend = None, ax = None, fig = None):
    n_plots = len(kh_list)+inc_iota
    if ax == None and fig == None:
        fig, ax = pt.subplots(nrows = n_plots, sharex = True); 
    if n_plots == 1: ax = [ax]
    if ylims == None: ylims = [[0,100] for i in kh_list]
    if inc_legend == None: inc_legend = [True for i in kh_list]
    gen_funcs.setup_publication_image(fig, height_prop = n_plots*0.65)
    if alf_plot_props == None: alf_plot_props = {'markersize':2,'marker':'o', 'linestyle':'None', 'rasterized':True, 'zorder':100}
    # try:
    #     tmp = int(size_scale)
    #     size_scale = [size_scale]*len(kh_list)
    # except (ValueError, TypeError) as e:
    #     pass
    # try:
    #     tmp = float(linewidths)
    #     linewidths = [linewidths]*len(kh_list)
    # except (ValueError, TypeError) as e:
    #     pass
    for i, (ax_tmp, kh, ylim) in enumerate(zip(ax, kh_list, ylims)):
        file_name = output_log
        base_dir = '/home/srh112/SSD_mount/raijin_short/whale_tail/input.kh0.{}0-kv1.000fixed_dir/CONTI_n1_fixed_expt_ne/'.format(kh)
        c = conti_results(base_dir, file_name)
        if inc_sound:c.plot_sound_modes(allowed_n= None, allowed_m=None, ax = ax_tmp, sqrt_s = False)
        c.plot_alfven_modes(allowed_n= None, allowed_m=None, ax = ax_tmp, sqrt_s = False, plot_props = alf_plot_props, ylim = ylim)
        ax_tmp.grid(); ax_tmp.set_xlim([0,1])
        ax_tmp.set_ylim(ylim)#ax.set_ylim([np.min(self.scan_arr[3,:]),np.max(self.scan_arr[3,:])])
        if inc_iota:ax[-1].plot(cont.s, cont.iota,label = 'kh=0.{}'.format(kh))
        if i == 0:
            ax_tmp.set_title('$\kappa_H = '+ '0.{}$'.format(kh))
        else:
            if kh!=kh_list[i-1]:
                ax_tmp.set_title('$\kappa_H = '+ '0.{}$'.format(kh))
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
        #handles, labels = ax_tmp.get_legend_handles_labels()
        #lab_temp = [lab+'&999' for lab in labels]  # add a low priority code 
                                                   # so that all have some code.
        #lab_text = np.array([ltemp.split('&')[0] for ltemp in lab_temp])
        #lab_order = [ltemp.split('&')[1] for ltemp in lab_temp]
        #lab_inds = np.argsort(lab_order)           # sort by increasing code
        box = ax_tmp.get_position()
        ax_tmp.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        if leg_tmp: leg_sh = ax_tmp.legend(prop=leg_props, fancybox=True,loc='center left', bbox_to_anchor=(1, 0.5))
    # #box = ax[-1].get_position()
    # #ax[-1].set_position([box.x0, box.y0, box.width * 0.8, box.height])
    if figname!=None: 
        for i in ['.svg','.pdf']: fig.savefig(figname + i, bbox_extra_artists=(leg_sh,),bbox_inches='tight')
    fig.canvas.draw(); fig.show()


#     def onpick(event):
#         print 'hello in onpick'
#         # on the pick event, find the orig line corresponding to the
#         # legend proxy line, and toggle the visibility
#         legline = event.artist
#         origline = lined[legline]
#         vis = not origline.get_visible()
#         origline.set_visible(vis)
#         # Change the alpha on the line in the legend so we can see what lines
#         # have been toggled
#         if vis:
#             legline.set_alpha(1.0)
#         else:
#             legline.set_alpha(0.2)
#         fig.canvas.draw()
#     #leg.get_frame().set_alpha(0.4)
#         ax.set_xlim([0,1.2])

#         lined = dict()
#         for legline, origline in zip(leg.get_lines(), lines):
#             legline.set_picker(5)  # 5 pts tolerance
#             lined[legline] = origline
#     a = np.loadtxt(base_dir + 'profiles_h1.dat', comments = '#')
#     for rat in [4./3, 5./4,6./5]:
#         ax.axvline(a[np.argmin(np.abs(a[:,2] - rat)),0])
#     fig.canvas.mpl_connect('pick_event',onpick)
#     fname_list.append('{}.png'.format(kh_tmp))
#     ax.set_ylim([0,60])
#     ax.set_title('{:.2f}'.format(kh_tmp/100.))
#     fig.savefig(fname_list[-1])
#     fig.canvas.draw(); fig.show()
# delay = 20
# import os
# os.system('convert -delay {} -loop 0 {} {}'.format(delay*32/float(len(fname_list)), ' {} {} '.format(fname_list[0],fname_list[0]) + ' '.join(fname_list), 'test.gif'))

