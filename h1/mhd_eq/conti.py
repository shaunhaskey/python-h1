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
        if plot_props == None: plot_props = {'color':(0.5,0.5,0.5), 'markersize':0.5,'marker':'.', 'linestyle':'None'}
        for i in range(1,len(self.overal_s_list)):
            if (self.symbol[i]==2) and (self.n_modes[i] in allowed_n_list) and (self.m_modes[i] in allowed_m_list):
                #ax.plot(self.overal_s_list[i], self.overal_freq_list[i],'.',color=(0.5,0.5,0.5),markersize=0.5)
                x_ax = np.sqrt(self.overal_s_list[i]) if sqrt_s else self.overal_s_list[i]
                ax.plot(x_ax, self.overal_freq_list[i],**plot_props)

    def plot_alfven_modes(self,allowed_n= None, allowed_m=None, ax = None, plot_props = None, sqrt_s = False):
        if allowed_n == None: allowed_n_list = self.n_modes
        if allowed_m == None: allowed_m_list = self.m_modes
        if ax == None: fig, ax = pt.subplots()
        if plot_props == None: plot_props = {'markersize':2,'marker':'o', 'linestyle':'None'}
        for i in range(1,len(self.overal_s_list)):
            if (self.symbol[i]==1) and (self.n_modes[i] in allowed_n_list) and (self.m_modes[i] in allowed_m_list):
                x_ax = np.sqrt(self.overal_s_list[i]) if sqrt_s else self.overal_s_list[i]
                tmp, = ax.plot(x_ax, self.overal_freq_list[i], label = 'm=%d,n=%d'%(self.m_modes[i],self.n_modes[i]),**plot_props)
                

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
