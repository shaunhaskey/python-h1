'''
This set of routines is for the tomographic inversion that is described in 

SH: 4May2014
'''
import matplotlib as mpl
import scipy.optimize as optimize
from scipy.interpolate import griddata as scipy_griddata
#import h1.mhd_eq.heliac_vmec_utils as hv_utils
import os,copy, time, scipy, pickle
import matplotlib.pyplot as pt
import numpy as np
from StringIO import StringIO
import h1.mhd_eq.heliac_worker_funcs as heliac
import scipy.interpolate as interp
import scipy.stats.distributions as dist
import scipy.optimize as opt
from scipy.stats import norm
import itertools
import multiprocessing
import h1.helper.generic_funcs as gen_funcs
import matplotlib.gridspec as gridspec

class Tomography():
    '''The main tomographic reconstruction class
    contains many plotting for the specific classes which are sub-classes of this one

    SRH : 4May2014
    '''
    def plot_amp_angle(self,):
        fig, ax = pt.subplots(nrows = 2)
        ax[0].plot(np.abs(self.T))
        ax[1].plot(np.angle(self.T))
        fig.suptitle(self.method)
        fig.canvas.draw(); fig.show()

    def plot_reprojection(self, ):
        fig, ax = pt.subplots(nrows = 2)
        ax[0].plot(np.abs(self.all_measurements[self.valid_channels]), label='data')
        ax[0].plot(np.abs(self.re_projection), label='reproj')
        ax[0].legend(loc='best')
        ax[1].plot(np.angle(self.all_measurements[self.valid_channels]))
        ax[1].plot(np.angle(self.re_projection))
        fig.suptitle(self.method)
        fig.canvas.draw(); fig.show()

    def plot_reprojection_comparison_extrap(self, n, m, LOS_object, valid_channels, all_measurements, geom_matrix, cut_values = None, multiplier = 1., pub_fig = 1, savefig_name = None):
        fig2, ax2 = pt.subplots(ncols = 8, sharey=True)
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
            fig2.set_figwidth(8.48*2.*cm_to_inch)
            fig2.set_figheight(8.48*1.25*cm_to_inch)

        if cut_values == None:
            cut_values = [33]
            cut_values = [130]
            cut_values=[LOS_object.valid_channels.shape[1]/2]
        if n.__class__ == int:
            n = [n]
            m = [m]
        measurements = all_measurements[valid_channels]
        valid_channels = valid_channels
        n_measurements, n_regions = geom_matrix.shape

        #If the geometry matrix or measurement matrix complex
        #Need to build a larger matrix that contains the linear system of equations
        #in real numbers only
        if geom_matrix.dtype == np.complex_ or measurements.dtype == np.complex:
            P = np.zeros((n_measurements*2, n_regions*2),dtype=float)
            #Build a real matrix 
            P[0::2,0::2] = +np.real(geom_matrix)
            P[0::2,1::2] = -np.imag(geom_matrix)
            P[1::2,0::2] = +np.imag(geom_matrix)
            P[1::2,1::2] = +np.real(geom_matrix)
            S = np.zeros(n_measurements*2, dtype=float)
            S[0::2] = +np.real(measurements)
            S[1::2] = +np.imag(measurements)
            T_tmp = np.zeros(n_regions*2, dtype=float)
            T_tmp[0::2] = +np.real(self.T)
            T_tmp[1::2] = +np.imag(self.T)
            input_type = 'complex'
            tmp = np.dot(P,T_tmp)- S
        else:
            print 'real input'
            P = geom_matrix
            S = measurements
            input_type = 'real'
            
        tmp = np.dot(P,T_tmp)
        re_projection = tmp[0::2]+1j*tmp[1::2]
        q = all_measurements*0
        q[valid_channels]= re_projection
        q[-valid_channels]= 0
        ax2 = np.array([[ax2[0],ax2[1], ax2[4], ax2[6]],[ax2[2], ax2[3], ax2[5],ax2[7]]])
        tmp = np.abs(all_measurements)*valid_channels
        tmp[-valid_channels] = np.nan
        im1_a = ax2[0,0].imshow(tmp,origin='upper',aspect='auto',interpolation='nearest', cmap='spectral')
        #im1_a = ax2[0,0].imshow(np.abs(all_measurements)*valid_channels,origin='upper',aspect='auto',interpolation='nearest', cmap='gist_stern')
        im1_p = ax2[1,0].imshow(np.angle(all_measurements)*valid_channels,origin='upper',aspect='auto',interpolation='nearest',cmap='RdBu')
        im1_p.set_clim([-np.pi,np.pi])

        tmp = np.abs(q)
        tmp[-valid_channels] = np.nan
        im2_a = ax2[0,1].imshow(tmp,origin='upper', aspect='auto',interpolation='nearest', cmap='spectral')
        
        im2_p = ax2[1,1].imshow(np.angle(q),origin='upper', aspect='auto',interpolation='nearest', cmap='RdBu')
        im2_a.set_clim(im1_a.get_clim())
        im2_p.set_clim(im1_p.get_clim())
        start = 0; increment = len(T_tmp)/len(n)

        tmp_y_axis = np.arange(len(all_measurements[:,0]))
        tmp_y_axis = np.max(tmp_y_axis)- tmp_y_axis
        tmp_reproj = all_measurements*0.
        tmp_reproj[valid_channels] = re_projection

        for i in cut_values:
            ax2[0,2].plot((np.abs(all_measurements[:,i]*valid_channels[:,i]))[::-1], tmp_y_axis, 'b-', label='Expt')
            ax2[0,2].plot((multiplier*np.abs(tmp_reproj[:,i]))[::-1], tmp_y_axis,'g-',label = 'Reproj')
            ax2[1,2].plot((np.angle(all_measurements[:,i])*valid_channels[:,i])[::-1], tmp_y_axis, 'b-')
            ax2[1,2].plot((np.angle(tmp_reproj[:,i])*valid_channels[:,i])[::-1], tmp_y_axis,'g-', label = '')

            ax2[0,3].plot((np.real(all_measurements[:,i]*valid_channels[:,i]))[::-1], tmp_y_axis, 'b-', label='Expt')
            ax2[0,3].plot((multiplier*np.real(tmp_reproj[:,i]))[::-1], tmp_y_axis,'g-',label = 'Reproj')
            ax2[1,3].plot((np.imag(all_measurements[:,i])*valid_channels[:,i])[::-1], tmp_y_axis, 'b-')
            ax2[1,3].plot((np.imag(tmp_reproj[:,i])*valid_channels[:,i])[::-1], tmp_y_axis,'g-', label = '')
            for j in [0,1]:
                for k in [0,1]:
                    tmp = ax2[j,k].get_ylim()
                    ax2[j,k].vlines(i,tmp[0], tmp[1])

        ax2[0,0].set_title('Amp Expt')
        ax2[0,1].set_title('Amp Reproj')
        ax2[1,0].set_title('Phase Expt')
        ax2[1,1].set_title('Phase Reproj')

        ax2[1,2].set_xlabel('Phase (rad)')
        #ax2[1,0].set_xlabel('Pixel')
        #ax2[1,1].set_xlabel('Pixel')
        #ax2[1,0].set_ylabel('Pixel')
        ax2[0,0].set_ylabel('Pixel')
        ax2[1,1].set_xlabel('Pixel')
        ax2[0,2].set_xlabel('Amp (a.u)')
        ax2[0,3].set_xlabel('Real (a.u)')
        ax2[1,3].set_xlabel('Imag (a.u)')
        ax2[0,2].legend(prop={'size':5}, loc = 'best')
        ax2[1,2].set_ylim([0,valid_channels.shape[0]])
        ax2[1,2].set_xlim([-np.pi,np.pi])
        center = valid_channels.shape[1]/2
        edge = center/4.
        for i in ax2:
            for j in i:
                j.tick_params(axis='both',which='both',labelbottom='off',labelleft='off')
                
        for i in [ax2[0,0], ax2[0,1], ax2[1,0], ax2[1,1]]:
            i.set_xlim([center+edge, center-edge])
            #i.tick_params(axis='both',which='both',labelbottom='off',labelleft='off')
            i.set_xlabel('Pixel')
        # ax2[0,0].set_xlim([center+edge, center-edge])
        # ax2[0,1].set_xlim([center+edge,center-edge])
        # ax2[1,0].set_xlim([center+edge, center-edge])
        # ax2[1,1].set_xlim([center+edge, center-edge])
        ax2[0,2].set_ylim([0,valid_channels.shape[0]])
        for i in ax2: 
            for j in i:
                j.grid()

        re_imag_limit = np.max(np.abs(np.array([ax2[1,3].get_xlim(), ax2[0,3].get_xlim()])))
        ax2[1,3].set_xlim([-re_imag_limit, re_imag_limit])
        ax2[0,3].set_xlim([-re_imag_limit, re_imag_limit])

        #fig2.subplots_adjust(hspace=0.0, wspace=0.03,left=0.05, bottom=0.05,top=0.95, right=0.95)
        fig2.subplots_adjust(hspace=0.0, wspace=0.03,left=0.05, bottom=0.05,top=0.95, right=0.95)
        #fig2.tight_layout()
        #fig2.savefig('reprojection.pdf')
        #fig2.savefig('reprojection.eps')
        if savefig_name!=None:
            fig2.savefig(savefig_name+'.pdf')
            fig2.savefig(savefig_name+'.eps')
        fig2.canvas.draw(); fig2.show()

    def plot_reprojection_comparison_extrap_diff(self, n, m, LOS_object, valid_channels, all_measurements, geom_matrix, cut_values = None, multiplier = 1., pub_fig = 1, savefig_name = None, include_lines = True, decimate_lines = 10, tomo_DC = None, dI_I = False):
        '''This generates the plots that are in the paper

        SRH: 4May2014
        '''
        #Setup axes and figure
        fig2 = pt.figure()
        if n.__class__ == int: n = [n];m = [m]
        if pub_fig: gen_funcs.setup_publication_image(fig2, height_prop = 1., fig_height = 8.48*1.25, single_col = True, fig_width = 8.48*2, replacement_kwargs = None)
        ax2 = []
        gs = gridspec.GridSpec(9,7)
        ax2.append(pt.subplot(gs[:8,2]))
        for i in range(3,7): ax2.append(pt.subplot(gs[:8,i], sharex = ax2[0],sharey = ax2[0]))
        for i,j in zip([0,2,4], [2,4,8]): ax2.append(pt.subplot(gs[i:j,0:2]))
        cbar_wave_ax, cbar_phase_ax, cbar_amp_ax = [pt.subplot(gs[8,i:j]) for i,j in zip([0,2,4], [2,4,7])]
        [phase_exp_ax, phase_best_ax, amp_exp_ax, amp_best_ax, error_ax, foo1, foo2, foo3] = ax2

        #wave in a poloidal cross section, and radial functions
        norm = False
        s_ax,s_x_label = (np.sqrt(LOS_object.segment_midpoints),'$\sqrt{s}$') if LOS_object.radial_s_spacing else (LOS_object.segment_midpoints, 's')
        extra_txt = '' if dI_I==False else 'dI/I '
        plot_radial_structure(self.T, s_ax, n, m, prov_ax = [ax2[5],ax2[6]], norm = norm, extra_txt = extra_txt, single_mode = None)#max_phase = 0)
        wave_field = self.calc_wave_field_phasors(n, m, LOS_object, LOS_object.segment_midpoints, tomo_DC = tomo_DC)
        wave_field_phasors = np.ma.array(wave_field, mask = LOS_object.grid_mask)
        image_extent = [LOS_object.r_grid[0,0], LOS_object.r_grid[0,-1], LOS_object.z_grid[0,0], LOS_object.z_grid[-1,0]]
        pol_clim = np.max(np.abs(wave_field_phasors))
        pol_clim = [0, +pol_clim] if (n==[0] and m==[0]) else [-pol_clim, +pol_clim]
        im = ax2[7].imshow(np.real(wave_field_phasors), interpolation = 'nearest', extent=image_extent, aspect='auto')
        im.set_clim(pol_clim)
        cbar_xsection = pt.colorbar(im, cax=cbar_wave_ax, orientation = 'horizontal')
        cbar_xsection.set_label('Real (a. u.)')
        ax2[7].set_xlabel('R (m)'); ax2[7].set_ylabel('Z (m)')
        ax2[5].set_title('Tomographic Reconstruction')
        ax2[7].set_xlim([1,1.4]); ax2[7].set_ylim([-0.25,0.25])
        cbar_xsection.set_ticks(np.round(np.linspace(im.get_clim()[0],im.get_clim()[1],7),1)[0:-1])
        if norm: ax2[5].set_ylim([0,0.12])
        for i in ax2[5:7]:i.set_xlabel(s_x_label); i.set_xlim([0,1])
        #ax2[5].legend(loc='center left')
        ax2[5].legend(loc='best')
        for i in ax2[5:8]: gen_funcs.setup_axis_publication(i, n_yticks = 5)
        gen_funcs.setup_axis_publication(ax2[5], n_xticks = 5)
        ax2[5].text(0.85,0.83*ax2[5].get_ylim()[1],'(a)')
        ax2[7].text(1.35,0.2,'(c)')
        ax2[6].text(0.85,(ax2[6].get_ylim()[1] - ax2[6].get_ylim()[0])*0.85 + ax2[6].get_ylim()[0],'(b)')
        measurements = all_measurements[valid_channels]
        valid_channels = valid_channels
        n_measurements, n_regions = geom_matrix.shape

        if include_lines:
            loc = LOS_object.valid_channels.shape[1]/2
            start_loc, end_loc, style = [LOS_object.start_indices[1], LOS_object.end_indices [1],'k-']
            LOS_r, LOS_z, meas_vals, reproj_vals, z_vals = self.LOS_start_end_cut(LOS_object, loc, start_loc, end_loc, None, r_val = np.array([0.5,1.5]))
            ax2[7].plot((np.array(LOS_r).T)[:,::8],(np.array(LOS_z).T)[:,::8], style, linewidth=0.3)

        #Phasors arrays for the image sets
        best_fit = self.all_measurements * 0
        best_fit[self.valid_channels] = self.re_projection
        max_val = np.max(np.abs(self.re_projection))
        im_clim = [-max_val, +max_val]
        meas_phasors = np.ma.array(self.all_measurements[:,:], mask = np.invert(self.valid_channels[:,:]))
        best_phasors = np.ma.array(best_fit[:,:], mask = np.invert(self.valid_channels[:,:]))
        error_phasors = meas_phasors - best_phasors
        amp_cmap = 'jet'

        im1_a = amp_exp_ax.imshow(np.abs(meas_phasors),origin='upper',aspect='auto',interpolation='nearest', cmap=amp_cmap)
        im2_a = amp_best_ax.imshow(np.abs(best_phasors),origin='upper', aspect='auto',interpolation='nearest', cmap=amp_cmap)
        im1_p = phase_exp_ax.imshow(np.arctan2(meas_phasors.imag,meas_phasors.real),origin='upper',aspect='auto',interpolation='nearest',cmap='RdBu')
        im2_p = phase_best_ax.imshow(np.arctan2(best_phasors.imag,best_phasors.real),origin='upper', aspect='auto',interpolation='nearest', cmap='RdBu')
        amp_diff = error_ax.imshow(np.abs(error_phasors),origin='upper', aspect='auto',interpolation='nearest', cmap=amp_cmap)

        for i in [im1_p, im2_p]: i.set_clim([-np.pi,np.pi])
        for i in [im1_a, im2_a, amp_diff]: i.set_clim([0,im1_a.get_clim()[1]])

        cbar2 = pt.colorbar(im1_a,cax = cbar_amp_ax, orientation = 'horizontal')
        cbar2.set_ticks(np.round(np.linspace(im1_a.get_clim()[0],im1_a.get_clim()[1],7),1)[0:-1])
        cbar2.set_label('Amp (a.u)')
        cbar3 = pt.colorbar(im1_p,cax = cbar_phase_ax, orientation = 'horizontal')
        cbar3.set_ticks(np.round(np.linspace(-np.pi,np.pi,7),1)[0:-1])
        cbar3.set_ticks([-np.pi,0,np.pi])
        cbar3.ax.set_xticklabels(['$-\pi$', '$0$','$\pi$'])
        cbar3.set_label('Phase (rad)')
        amp_exp_ax.set_title('(f)Experiment')
        amp_best_ax.set_title('(g)Best Fit')
        phase_exp_ax.set_title('(d)Experiment')
        phase_best_ax.set_title('(e)Best Fit')
        error_ax.set_title('(h)Error')
        center = valid_channels.shape[1]/2
        edge = center/4.
        for j in ax2[0:5]:j.tick_params(axis='both',which='both',labelbottom='off',labelleft='off')
        for j in ax2:j.grid()
        for i in ax2[0:5]: i.set_xlim([center+edge*0.85, center-edge*1.2])
        amp_exp_ax.set_ylim([0,valid_channels.shape[0]-40])
        if savefig_name!=None:
            gs.tight_layout(fig2, pad = 0.0001)
            for i in ['.pdf','.eps','.svg','.png']: fig2.savefig(savefig_name+i, bbox_inches='tight',pad_inches=0.05)
        fig2.canvas.draw(); fig2.show()

    def plot_reprojection_comparison_extrap_diff_hue_sat(self, n, m, LOS_object, valid_channels, all_measurements, geom_matrix, cut_values = None, multiplier = 1., pub_fig = 1, savefig_name = None, include_lines = True, decimate_lines = 10):
        '''Plots, but combining the phase and amplitude into hue and saturation

        SRH: 4May2014
        '''
        #fig2, ax2 = pt.subplots(ncols = 5, sharey=True, sharex = True)
        fig2 = pt.figure()
        if n.__class__ == int: n = [n];m = [m]
        if pub_fig: gen_funcs.setup_publication_image(fig2, height_prop = 1., fig_height = 8.48*1.25, single_col = True, fig_width = 8.48*2, replacement_kwargs = None)
        ax2 = []
        gs = gridspec.GridSpec(9,6)
        ax2.append(None)
        ax2.append(None)
        ax2.append(pt.subplot(gs[:8,3]))
        ax2.append(pt.subplot(gs[:8,4], sharex = ax2[2],sharey = ax2[2]))
        ax2.append(pt.subplot(gs[:8,5], sharex = ax2[2],sharey = ax2[2]))
        ax2.append(pt.subplot(gs[0:2,0:3]))
        ax2.append(pt.subplot(gs[2:4,0:3]))#, sharex = ax2[5]))
        ax2.append(pt.subplot(gs[4:8,0:3]))
        cbar_wave_ax = pt.subplot(gs[8,0:3])
        cbar_amp_ax = pt.subplot(gs[8,3:])
        norm = False
        if LOS_object.radial_s_spacing:
            plot_radial_structure(self.T, np.sqrt(LOS_object.segment_midpoints), n, m, prov_ax = [ax2[5],ax2[6]], norm = norm, extra_txt = '', single_mode = None)
        else:
            plot_radial_structure(self.T, LOS_object.segment_midpoints, n, m, prov_ax = [ax2[5],ax2[6]], norm = norm, extra_txt = '', single_mode = None)

        start = 0; increment = len(self.T)/len(n)
        s_vals = LOS_object.segment_midpoints
        #s_vals = (s_values[1:]+s_values[0:-1])/2
        wave_field = 0
        for n_cur, m_cur in zip(n,m):
            end = start + increment
            wave_amp = np.interp(LOS_object.s_cross_sect.flatten(), s_vals, np.real(self.T[start:end])).reshape(LOS_object.theta_cross_sect.shape) + 1j* np.interp(LOS_object.s_cross_sect.flatten(), s_vals, np.imag(self.T[start:end])).reshape(LOS_object.theta_cross_sect.shape)
            wave_field = wave_field + wave_amp * np.exp(1j*(n_cur*LOS_object.phi_cross_sect + m_cur * LOS_object.theta_cross_sect))
            start = +end

        wave_field_tmp = wave_field 
        image_extent = [LOS_object.r_grid[0,0], LOS_object.r_grid[0,-1], LOS_object.z_grid[0,0], LOS_object.z_grid[-1,0]]
        im = ax2[7].imshow(np.ma.array(np.real(wave_field_tmp), mask = LOS_object.grid_mask), interpolation = 'nearest', extent=image_extent, aspect='auto')
        tmp = im.get_clim()
        max_val = np.max(np.abs(tmp))
        if n==[0] and m==[0]:
            im.set_clim([0, max_val])
        else:
            im.set_clim([-max_val, max_val])
        #cbar_xsection = pt.colorbar(im, ax=[ax2[5],ax2[6],ax2[7]], orientation = 'horizontal', pad = 0.01, aspect = 20./3*2)
        cbar_xsection = pt.colorbar(im, cax=cbar_wave_ax, orientation = 'horizontal')
        cbar_xsection.set_label('Real (a. u.)')
        ax2[7].set_xlabel('R (m)')
        ax2[7].set_ylabel('Z (m)')
        ax2[5].set_title('Tomographic Reconstruction')
        ax2[7].set_xlim([1,1.4])
        ax2[7].set_ylim([-0.25,0.25])
        cbar_xsection.set_ticks(np.round(np.linspace(im.get_clim()[0],im.get_clim()[1],7),1)[0:-1])
        ax2[5].set_xlim([0,1])
        ax2[6].set_xlim([0,1])
        if norm:
            ax2[5].set_ylim([0,0.12])
        if LOS_object.radial_s_spacing:
            ax2[6].set_xlabel(r'$\sqrt{s}$')
            ax2[5].set_xlabel(r'$\sqrt{s}$')
        else:
            ax2[6].set_xlabel(r'$s$')
            ax2[5].set_xlabel(r'$s$')
        ax2[5].legend(loc='center left')
        ax2[7].set_xticks(ax2[7].get_xticks()[::2])
        ax2[6].set_yticks(ax2[6].get_yticks()[::2])
        ax2[5].set_yticks(ax2[5].get_yticks()[::2])
        ax2[5].text(0.85,0.83*ax2[5].get_ylim()[1],'(a)')
        ax2[7].text(1.35,0.2,'(c)')
        ax2[6].text(0.85,(ax2[6].get_ylim()[1] - ax2[6].get_ylim()[0])*0.85 + ax2[6].get_ylim()[0],'(b)')
        measurements = all_measurements[valid_channels]
        valid_channels = valid_channels
        n_measurements, n_regions = geom_matrix.shape
        if include_lines:
            loc = LOS_object.valid_channels.shape[1]/2
            #for start_loc, end_loc, style in zip(LOS_object.start_indices, LOS_object.end_indices,['b-','k-','r-']):
            start_loc, end_loc, style = [LOS_object.start_indices[1], LOS_object.end_indices [1],'k-']
            for i in range(start_loc, end_loc, decimate_lines):
                if LOS_object.valid_channels[i,loc]:
                    point1_z = LOS_object.intersection1[i,loc,2] 
                    point2_z = LOS_object.intersection2[i,loc,2]
                    point1_r = np.sqrt(LOS_object.intersection1[i,loc,0]**2+LOS_object.intersection1[i,loc,1]**2) 
                    point2_r = np.sqrt(LOS_object.intersection2[i,loc,0]**2+LOS_object.intersection2[i,loc,1]**2)
                    gradient = (point2_z - point1_z)/(point2_r - point1_r)
                    offset = point1_z - gradient*point1_r
                    r_val = np.array([0.5,1.5])
                    z_val = gradient*r_val + offset
                    #z_val = [point1_z, point2_z]
                    #r_val = [point1_r, point2_r]
                    ax2[7].plot(r_val,z_val,style, linewidth=0.3)

        #If the geometry matrix or measurement matrix complex
        #Need to build a larger matrix that contains the linear system of equations
        #in real numbers only
        if geom_matrix.dtype == np.complex_ or measurements.dtype == np.complex:
            P = np.zeros((n_measurements*2, n_regions*2),dtype=float)
            #Build a real matrix 
            P[0::2,0::2] = +np.real(geom_matrix)
            P[0::2,1::2] = -np.imag(geom_matrix)
            P[1::2,0::2] = +np.imag(geom_matrix)
            P[1::2,1::2] = +np.real(geom_matrix)
            S = np.zeros(n_measurements*2, dtype=float)
            S[0::2] = +np.real(measurements)
            S[1::2] = +np.imag(measurements)
            T_tmp = np.zeros(n_regions*2, dtype=float)
            T_tmp[0::2] = +np.real(self.T)
            T_tmp[1::2] = +np.imag(self.T)
            input_type = 'complex'
            tmp = np.dot(P,T_tmp)- S
        else:
            print 'real input'
            P = geom_matrix
            S = measurements
            input_type = 'real'
            
        tmp = np.dot(P,T_tmp)
        tmp = tmp[0::2]+1j*tmp[1::2]
        re_projection = all_measurements*0
        re_projection[valid_channels]= tmp
        re_projection[-valid_channels]= 0
        tmp_meas = np.abs(all_measurements)*valid_channels
        tmp_meas = all_measurements
        #tmp_meas[-valid_channels] = np.nan
        tmp_meas[-valid_channels] = 0
        [phase_exp_ax, phase_best_ax, amp_exp_ax, amp_best_ax, error_ax, foo1, foo2, foo3] = ax2
        amp_cmap = 'jet'
        im1_a = amp_exp_ax.imshow(complex_array_to_rgb(tmp_meas, rmax=1, theme = 'dark'),origin='upper',aspect='auto',interpolation='nearest', cmap=amp_cmap)
        tmp = np.abs(re_projection)
        tmp = re_projection
        tmp[-valid_channels] = 0
        #im2_a = amp_best_ax.imshow(tmp,origin='upper', aspect='auto',interpolation='nearest', cmap=amp_cmap)
        im2_a = amp_best_ax.imshow(complex_array_to_rgb(tmp, rmax=1, theme = 'dark'),origin='upper', aspect='auto',interpolation='nearest', cmap=amp_cmap)
        foo_tmp = np.angle(all_measurements)
        foo_tmp[-valid_channels] = np.nan
        average_magnitude = np.mean(np.abs(re_projection[valid_channels]))
        difference = re_projection*0
        difference[valid_channels] = re_projection[valid_channels]- all_measurements[valid_channels]
        difference[-valid_channels] = 0
        amp_diff = error_ax.imshow(complex_array_to_rgb(difference, rmax=1, theme = 'dark'),origin='upper', aspect='auto',interpolation='nearest', cmap=amp_cmap)
        amp_diff.set_clim(im1_a.get_clim())
        cbar2 = hue_sat_cbar(cbar_amp_ax, rmax = 1)
        start = 0; increment = len(T_tmp)/len(n)

        tmp_y_axis = np.arange(len(all_measurements[:,0]))
        tmp_y_axis = np.max(tmp_y_axis)- tmp_y_axis
        amp_exp_ax.set_title('(f)Experiment')
        amp_best_ax.set_title('(g)Best Fit')
        error_ax.set_title('(h)Error')
        center = valid_channels.shape[1]/2
        edge = center/4.
        for j in ax2[2:5]:j.tick_params(axis='both',which='both',labelbottom='off',labelleft='off')
        for j in ax2[2:]:j.grid()
        for i in ax2[2:5]: i.set_xlim([center+edge*0.85, center-edge*1.2])
        amp_exp_ax.set_ylim([0,valid_channels.shape[0]-40])
        
        #fig2.subplots_adjust(hspace=0.0, wspace=0.03,left=0.05, bottom=0.05,top=0.95, right=0.95)
        if savefig_name!=None:
            gs.tight_layout(fig2, pad = 0.0001)
            #fig2.tight_layout(pad = 0.05)
            fig2.savefig(savefig_name+'.pdf', bbox_inches='tight',pad_inches=0.05)
            fig2.savefig(savefig_name+'.eps', bbox_inches='tight',pad_inches=0.05)
        fig2.canvas.draw(); fig2.show()

    def plot_reprojection_comparison(self, n, m, LOS_object, cut_values = None, multiplier = 1., pub_fig = 1):
        '''This function is for making publication images comparing
        the reprojection

        SRH:29Nov2013 
        '''
        if cut_values == None:
            cut_values = [33]
            cut_values = [130]
            cut_values=[LOS_object.valid_channels.shape[1]/2]
        if n.__class__ == int:
            n = [n]
            m = [m]

        fig2, ax2 = pt.subplots(nrows = 2, ncols = 3, sharey=True)
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
            fig2.set_figwidth(8.48*2.*cm_to_inch)
            fig2.set_figheight(8.48*1.35*cm_to_inch)

        im1_a = ax2[0,0].imshow(np.abs(self.all_measurements)*self.valid_channels,origin='upper',aspect='auto',interpolation='nearest')
        im1_p = ax2[1,0].imshow(np.angle(self.all_measurements)*self.valid_channels,origin='upper',aspect='auto',interpolation='nearest')
        im1_p.set_clim([-np.pi,np.pi])
        q = self.all_measurements*0
        q[self.valid_channels]= self.re_projection
        q[-self.valid_channels]= 0
        im2_a = ax2[0,1].imshow(multiplier*np.abs(q),origin='upper', aspect='auto',interpolation='nearest')
        im2_p = ax2[1,1].imshow(np.angle(q),origin='upper', aspect='auto',interpolation='nearest')
        im2_a.set_clim(im1_a.get_clim())
        im2_p.set_clim(im1_p.get_clim())
        start = 0; increment = len(self.T)/len(n)

        #ax2[0,3].plot(np.abs(self.all_measurements[self.valid_channels]), label='|data|')
        #ax2[0,3].plot(np.abs(self.re_projection), label='|reproj|')
        #ax2[0,3].legend(loc='best')
        tmp_y_axis = np.arange(len(self.all_measurements[:,0]))
        tmp_y_axis = np.max(tmp_y_axis)- tmp_y_axis
        tmp_reproj = self.all_measurements*0.
        tmp_reproj[self.valid_channels] = self.re_projection

        for i in cut_values:
            ax2[0,2].plot((np.abs(self.all_measurements[:,i]*self.valid_channels[:,i]))[::-1], tmp_y_axis, '-', label='Measurement')
            ax2[0,2].plot((multiplier*np.abs(tmp_reproj[:,i]))[::-1], tmp_y_axis,'-',label = 'Reprojection')
            ax2[1,2].plot((np.angle(self.all_measurements[:,i])*self.valid_channels[:,i])[::-1], tmp_y_axis, '-')
            ax2[1,2].plot((np.angle(tmp_reproj[:,i])*self.valid_channels[:,i])[::-1], tmp_y_axis,'-', label = '')
            for j in [0,1]:
                for k in [0,1]:
                    tmp = ax2[j,k].get_ylim()
                    ax2[j,k].vlines(i,tmp[0], tmp[1])
        ax2[1,2].set_xlabel('Phase (rad)')
        ax2[1,0].set_xlabel('Pixel')
        ax2[1,1].set_xlabel('Pixel')
        ax2[1,0].set_ylabel('Pixel')
        ax2[0,0].set_ylabel('Pixel')
        ax2[0,2].set_xlabel('Amplitude (a.u)')
        ax2[0,2].legend(prop={'size':5}, loc = 'upper left')
        ax2[1,2].set_ylim([0,self.valid_channels.shape[0]])
        ax2[1,2].set_xlim([-np.pi,np.pi])
        center = self.valid_channels.shape[1]/2
        edge = center/2.
        ax2[0,0].set_xlim([center+edge, center-edge])
        ax2[0,1].set_xlim([center+edge,center-edge])
        ax2[1,0].set_xlim([center+edge, center-edge])
        ax2[1,1].set_xlim([center+edge, center-edge])
        ax2[0,2].set_ylim([0,self.valid_channels.shape[0]])

            
        #ax2[0,3].plot(np.abs(self.all_measurements[self.valid_channels]), label='|data|')
        #ax2[0,3].plot(np.abs(self.re_projection), label='|reproj|')


        #ax2[1,3].plot(np.angle(self.all_measurements[self.valid_channels]), label = 'Arg(data)')
        #ax2[1,3].plot(np.angle(self.re_projection), label = 'Arg(reproj)')
        #ax2[1,3].legend(loc='best')
        fig2.subplots_adjust(hspace=0.05, wspace=0.05,left=0.05, bottom=0.05,top=0.95, right=0.95)
        fig2.tight_layout()
        #fig2.suptitle('Tomo method: {} Top : (Meas amp, reproj amp, amp slice, amp vs norm flux), Bot:(Meas Phase, reproj Phase, phase slice, phase vs norm flux)'.format(self.method))
        fig2.savefig('reprojection.pdf')
        fig2.savefig('reprojection.eps')
        fig2.canvas.draw(); fig2.show()



    def plot_lots_of_things(self, n, m, LOS_object, cut_values = None, multiplier = 1.):
        if cut_values == None:
            cut_values=[LOS_object.valid_channels.shape[1]/2]
        if n.__class__ == int:
            n = [n]
            m = [m]
        fig2, ax2 = pt.subplots(nrows = 2, ncols = 4)
        im1_a = ax2[0,0].imshow(np.abs(self.all_measurements)*self.valid_channels,origin='upper', aspect = 'auto',interpolation='nearest')
        im1_p = ax2[1,0].imshow(np.angle(self.all_measurements)*self.valid_channels,origin='upper', aspect='auto',interpolation='nearest')
        im1_p.set_clim([-np.pi,np.pi])
        q = self.all_measurements*0
        q[self.valid_channels]= self.re_projection
        q[-self.valid_channels]= 0
        im2_a = ax2[0,1].imshow(multiplier*np.abs(q),origin='upper', aspect='auto',interpolation='nearest')
        im2_p = ax2[1,1].imshow(np.angle(q),origin='upper', aspect='auto',interpolation='nearest')
        im2_a.set_clim(im1_a.get_clim())
        im2_p.set_clim(im1_p.get_clim())
        start = 0; increment = len(self.T)/len(n)
        if len(self.T)%len(n)!=0: raise ValueError('Something wrong')
        for n_cur, m_cur in zip(n,m):
            end = start + increment
            cur_T = self.T[start:end]
            start = +end
            ax2[0,3].plot(LOS_object.segment_midpoints, np.abs(cur_T), label='|{},{}|'.format(n_cur, m_cur))
            ax2[1,3].plot(LOS_object.segment_midpoints, np.angle(cur_T), label='Arg({},{})'.format(n_cur, m_cur))
        ax2[0,3].legend(loc='best')
        ax2[0,3].set_xlim([0,1])
        ax2[1,3].set_xlim([0,1])
        ax2[1,3].set_ylim([-np.pi,np.pi])

        ax2[1,3].legend(loc='best')

        #ax2[0,3].plot(np.abs(self.all_measurements[self.valid_channels]), label='|data|')
        #ax2[0,3].plot(np.abs(self.re_projection), label='|reproj|')
        #ax2[0,3].legend(loc='best')
        tmp_y_axis = np.arange(len(self.all_measurements[:,0]))
        tmp_y_axis = np.max(tmp_y_axis)- tmp_y_axis
        tmp_reproj = self.all_measurements*0.
        tmp_reproj[self.valid_channels] = self.re_projection

        for i in cut_values:
            ax2[0,2].plot(np.abs(self.all_measurements[:,i]*self.valid_channels[:,i]), tmp_y_axis, '-')
            ax2[0,2].plot(multiplier*np.abs(tmp_reproj[:,i]), tmp_y_axis)
            ax2[1,2].plot(np.angle(self.all_measurements[:,i])*self.valid_channels[:,i], tmp_y_axis, '-')
            ax2[1,2].plot(np.angle(tmp_reproj[:,i])*self.valid_channels[:,i], tmp_y_axis,'-')
            for j in [0,1]:
                for k in [0,1]:
                    tmp = ax2[j,k].get_ylim()
                    ax2[j,k].vlines(i,tmp[0], tmp[1])
        ax2[1,2].set_ylim([0,self.valid_channels.shape[0]])
        ax2[1,2].set_xlim([-np.pi,np.pi])
        ax2[0,2].set_ylim([0,self.valid_channels.shape[0]])
        #ax2[0,3].plot(np.abs(self.all_measurements[self.valid_channels]), label='|data|')
        #ax2[0,3].plot(np.abs(self.re_projection), label='|reproj|')


        #ax2[1,3].plot(np.angle(self.all_measurements[self.valid_channels]), label = 'Arg(data)')
        #ax2[1,3].plot(np.angle(self.re_projection), label = 'Arg(reproj)')
        #ax2[1,3].legend(loc='best')
        fig2.subplots_adjust(hspace=0.0, wspace=0.0,left=0., bottom=0.05,top=0.95, right=0.95)
        fig2.suptitle('Tomo method: {} Top : (Meas amp, reproj amp, amp slice, amp vs norm flux), Bot:(Meas Phase, reproj Phase, phase slice, phase vs norm flux)'.format(self.method))
        fig2.canvas.draw(); fig2.show()

    def plot_wave_field(self, n, m, LOS_object,s_values):
        #Get the boozer locations to interpolate
        if n.__class__ == int:
            n = [n]
            m = [m]
        fig, ax = pt.subplots(nrows = 2, ncols = 3)
        im_s = ax[0,0].imshow(np.ma.array(LOS_object.s_cross_sect, mask=LOS_object.grid_mask))
        im_s.set_clim([0,1])
        im_theta = ax[0,1].imshow(np.ma.array(LOS_object.theta_cross_sect%(2.*np.pi), mask=LOS_object.grid_mask))
        im_theta.set_clim([0,2.*np.pi])
        im_phi = ax[0,2].imshow(np.ma.array(LOS_object.phi_cross_sect, mask=LOS_object.grid_mask))
        im_phi.set_clim([-np.pi,np.pi])
        s_vals = (s_values[1:]+s_values[0:-1])/2
        wave_field = 0
        start = 0; increment = len(self.T)/len(n)
        if len(self.T)%len(n)!=0: raise ValueError('Something wrong')
        for n_cur, m_cur in zip(n,m):
            end = start + increment
            wave_amp = np.interp(LOS_object.s_cross_sect.flatten(), s_vals, np.real(self.T[start:end])).reshape(LOS_object.theta_cross_sect.shape) + 1j* np.interp(LOS_object.s_cross_sect.flatten(), s_vals, np.imag(self.T[start:end])).reshape(LOS_object.theta_cross_sect.shape)
            wave_field = wave_field + wave_amp * np.exp(1j*(n_cur*LOS_object.phi_cross_sect + m_cur * LOS_object.theta_cross_sect))
            start = +end
        im = ax[1,0].imshow(np.ma.array(np.abs(wave_field), mask = LOS_object.grid_mask), interpolation = 'nearest')
        im.set_clim([0, np.mean(np.abs(wave_field))*3])
        ax[1,2].imshow(np.ma.array(np.real(wave_field), mask = LOS_object.grid_mask), interpolation = 'nearest')
        im = ax[1,1].imshow(np.ma.array(np.angle(wave_field), mask = LOS_object.grid_mask), interpolation = 'nearest')
        im.set_clim([-np.pi,np.pi])
        fig.subplots_adjust(hspace=0.0, wspace=0.0,left=0., bottom=0.05,top=0.95, right=0.95)
        fig.suptitle('{} top row : s, theta, phi; bottom row : abs, phase, real'.format(self.method))
        fig.canvas.draw(); fig.show()

    def calc_wave_field_phasors(self, n, m, LOS_object, s_vals, tomo_DC = None):
        '''This function sums the contribution from each mode for a
        poloidal cross section
        
        SRH : 4May2014 '''
        wave_field = 0
        start = 0; increment = len(self.T)/len(n)
        if len(self.T)%len(n)!=0: raise ValueError('Something wrong')
        DC_mult = 1. if tomo_DC == None else tomo_DC.T
        for n_cur, m_cur in zip(n,m):
            end = start + increment
            T_cur = self.T[start:end] * DC_mult
            wave_amp = np.interp(LOS_object.s_cross_sect.flatten(), s_vals, np.real(T_cur)).reshape(LOS_object.theta_cross_sect.shape) + 1j* np.interp(LOS_object.s_cross_sect.flatten(), s_vals, np.imag(T_cur)).reshape(LOS_object.theta_cross_sect.shape)
            wave_field = wave_field + wave_amp * np.exp(1j*(n_cur*LOS_object.phi_cross_sect + m_cur * LOS_object.theta_cross_sect))
            start = +end
        return wave_field

    def LOS_start_end_cut(self,LOS_object, loc, start_loc, end_loc, best_fit, r_val = None):
        if r_val == None: r_val = np.array([0.5,1.4])
        meas_vals = []; z_vals = []; reproj_vals = []
        LOS_r = []; LOS_z = []
        for j in range(start_loc, end_loc):
            if LOS_object.valid_channels[j,loc]:
                point1_z = LOS_object.intersection1[j,loc,2] 
                point2_z = LOS_object.intersection2[j,loc,2]
                point1_r = np.sqrt(LOS_object.intersection1[j,loc,0]**2+LOS_object.intersection1[j,loc,1]**2) 
                point2_r = np.sqrt(LOS_object.intersection2[j,loc,0]**2+LOS_object.intersection2[j,loc,1]**2)
                gradient = (point2_z - point1_z)/(point2_r - point1_r)
                offset = point1_z - gradient*point1_r
                #r_val = np.array([0.5,1.4])
                z_val = gradient*r_val + offset
                LOS_r.append(r_val)
                LOS_z.append(z_val)
                meas_vals.append(self.all_measurements[j , loc])
                if best_fit!=None: reproj_vals.append(best_fit[j , loc])
                z_vals.append(+z_val[1])
                #scale_fact = 0.05/np.max(np.abs(meas_vals))
        return LOS_r, LOS_z, meas_vals, reproj_vals, z_vals

    def plot_wave_field_animation(self, n, m, LOS_object,s_values, n_images = 10, save_name = None, delay = 10, inc_cbar=0, pub_fig = False, save_fig = None, inc_profile = False,kh = None, tomo_DC = None, dI_I = False):
        '''
        '''
        #Setup the figure
        if inc_profile:
            fig = pt.figure(); ax = []
            ax.append(pt.subplot2grid((8,6), (0,0), rowspan=4, colspan = 3))
            ax.append(pt.subplot2grid((8,6), (4,0), rowspan = 2, colspan = 3))
            ax[-1].grid(True)
            ax.append(pt.subplot2grid((8,6), (6,0), rowspan = 2, colspan = 3))
            ax[-1].grid(True)
            ax_im1 = pt.subplot2grid((8,6), (0,3), rowspan = 7)
            ax_im2 = pt.subplot2grid((8,6), (0,4), rowspan = 7)
            ax_im3 = pt.subplot2grid((8,6), (0,5), rowspan = 7)
            cbar_ax = pt.subplot2grid((8,6), (7,3), rowspan = 1, colspan = 3)
        else:
            fig, ax = pt.subplots()
            ax = [ax]
        if inc_cbar:
            from mpl_toolkits.axes_grid1 import make_axes_locatable
            divider = make_axes_locatable(ax[0])
            cax = divider.append_axes("right", size="10%", pad=0.05)
        gen_funcs.setup_publication_image(fig, height_prop = 1., single_col = True, fig_width = 8.48*1.6, replacement_kwargs = None)

        if n.__class__ == int: n = [n];m = [m]
        aspect = 'auto'
        tmp_cmap = 'RdBu'

        #calculate the wave field phasors for a poloidal cross section
        s_vals = (s_values[1:]+s_values[0:-1])/2
        wave_field = self.calc_wave_field_phasors(n, m, LOS_object, s_vals, tomo_DC=tomo_DC)
        wave_field_phasors = np.ma.array(wave_field, mask = LOS_object.grid_mask)
        pol_clim = np.max(np.abs(wave_field_phasors))
        pol_clim = [-pol_clim, +pol_clim]

        image_list = []
        image_extent = [LOS_object.r_grid[0,0], LOS_object.r_grid[0,-1], LOS_object.z_grid[0,0], LOS_object.z_grid[-1,0]]

        start_loc, end_loc = [LOS_object.start_indices[1], LOS_object.end_indices[1]]
        start_loc = 0; end_loc = self.all_measurements.shape[0]
        center = self.valid_channels.shape[1]/2; edge = center/4.

        #Phasors arrays for the image sets
        best_fit = self.all_measurements * 0
        best_fit[self.valid_channels] = self.re_projection
        max_val = np.max(np.abs(self.re_projection))
        im_clim = [-max_val, +max_val]
        meas_phasors = np.ma.array(self.all_measurements[start_loc:end_loc,:], mask = np.invert(self.valid_channels[start_loc:end_loc,:]))
        best_phasors = np.ma.array(best_fit[start_loc:end_loc,:], mask = np.invert(self.valid_channels[start_loc:end_loc,:]))
        error_phasors = meas_phasors - best_phasors

        #Plot the wave profile
        extra_txt = '' if dI_I==False else 'dI/I '
        if inc_profile:
            norm = False
            s_ax,x_label = (np.sqrt(LOS_object.segment_midpoints),'$\sqrt{s}$') if LOS_object.radial_s_spacing else (LOS_object.segment_midpoints, 's')
            plot_radial_structure(self.T, s_ax, n, m, prov_ax = [ax[1],ax[2]], norm = norm, extra_txt = extra_txt, single_mode = None, max_phase = 0)
            ax[1].legend(loc='best')
            ax[1].set_xlim([0,1]); ax[2].set_xlim([0,1]); ax[2].set_ylim([-7,0])
            ax[2].set_ylabel('Phase (rad)'); ax[1].set_ylabel('Amplitude (a.u)'); ax[2].set_xlabel('$\sqrt{s}$')

        loc = LOS_object.valid_channels.shape[1]/2
        start_loc, end_loc, style = [LOS_object.start_indices[1], LOS_object.end_indices [1],'k-']
        LOS_r, LOS_z, meas_vals, reproj_vals, z_vals = self.LOS_start_end_cut(LOS_object, loc, start_loc, end_loc, best_fit, r_val = None)
        scale_fact = 0.05/np.max(np.abs(meas_vals))

        #Make the animated parts
        for i, t_cur in enumerate(np.linspace(0,2.*np.pi,n_images,endpoint=False)):
            for tmp_ax in [ax[0], ax_im1, ax_im2, ax_im3]: tmp_ax.cla()
            #masked_array1 = np.ma.array(np.real(self.all_measurements* np.exp(1j*t_cur))[start_loc:end_loc,:], mask = np.invert(self.valid_channels[start_loc:end_loc,:]))
            for ax_im, tmp_data, tmp_title in zip([ax_im1, ax_im2, ax_im3],[meas_phasors, best_phasors, error_phasors], ['Expt', 'Fit', 'Error']):
                cam_im = ax_im.imshow(np.real(tmp_data * np.exp(1j*t_cur)), cmap = tmp_cmap, aspect = aspect, origin='upper',interpolation='nearest')
                ax_im.set_title(tmp_title)
                cam_im.set_clim(im_clim)
                ax_im.set_xlim([center+edge*0.90, center-edge*1.25])
                ax_im.set_ylim([0,self.valid_channels.shape[0]-40])
                ax_im.tick_params(axis='both',which='both',labelbottom='off',labelleft='off')
            xval = np.mean([center])
            L = float(self.valid_channels.shape[0]-40)
            for j in [ax_im1, ax_im2]: j.vlines(xval, L/2 - L/6, L/2 + L/6)
            #im = ax[0].imshow(np.ma.array(np.real(wave_field * np.exp(1j*t_cur)), mask = LOS_object.grid_mask), interpolation = 'nearest', extent=image_extent, aspect='auto',origin='lower')
            im = ax[0].imshow(np.real(wave_field_phasors * np.exp(1j*t_cur)), interpolation = 'nearest', extent=image_extent, aspect='auto',origin='lower')
            #tmp = im.get_clim()
            #print tmp
            #max_val = np.max(np.abs(tmp))
            #im.set_clim([-max_val, +max_val])
            im.set_clim(pol_clim)
            if inc_cbar and i==0:
                cbar = pt.colorbar(im, cax=cax)
                cbar.set_label('Wave Amplitude (a.u)')
                cbar2 = pt.colorbar(cam_im, cax = cbar_ax, orientation = 'horizontal')
                cbar2.set_label('LOS Amplitude (a.u)')
            ax[0].plot((np.array(LOS_r).T)[:,::8],(np.array(LOS_z).T)[:,::8], style, linewidth=0.3)
            tmp_val = np.real(np.array(meas_vals)*np.exp(1j*t_cur))*scale_fact
            tmp_val_reproj = np.real(np.array(reproj_vals)*np.exp(1j*t_cur))*scale_fact

            ax[0].plot(tmp_val+1.4, np.array(z_vals),'-')
            ax[0].plot(tmp_val_reproj+1.4, np.array(z_vals),'-')
            ax[0].plot([1.4]*len(z_vals),np.array(z_vals))
            if kh == None: kh = -1
            modifiers = {'set_xlabel':'R (m)', 'set_ylabel':'Z (m)','set_xlim':[1,1.45], 
                         'set_ylim':[-0.25,0.25], 'set_aspect':'equal', 'set_title': '$\kappa_h = '+'{:.2f}$ '.format(kh) + 'LOS projections,\nblue: experiment, green: resconstruction'}
            gen_funcs.modify_axis(ax[0], modifiers)
            if save_fig!=None:
                xticks = ax[0].xaxis.get_major_ticks()
                for tmp_tick in range(0,len(xticks),2):
                    xticks[tmp_tick].label1.set_visible(False)
                if i==0:
                    fig.tight_layout()
                    fig.savefig(save_fig+'.pdf')
                    fig.savefig(save_fig+'.eps')
            fig.savefig('{:02d}.png'.format(i),dpi=200)
            image_list.append('{:02d}.png'.format(i))

        if save_name!=None:
            image_string = ''
            for i in image_list: image_string+=i + ' '
            print image_string
            os.system('convert -delay {} -loop 0 {} {}'.format(delay*32/float(n_images), image_string, save_name))
            os.system('zip {} {}'.format(save_name.rstrip('gif')+'zip', image_string))
        fig.canvas.draw(); fig.show()

    def plot_convergence(self):
        fig, ax = pt.subplots(nrows = 2)
        for i, T in enumerate(self.T_list):
            real_part = T
            if self.input_type == 'complex':
                real_part = T[0::2]
                ax[1].plot(T[1::2])
            ax[0].plot(real_part)
            ax[0].text(np.argmax(real_part), np.max(real_part),str(i))
        fig.canvas.draw(); fig.show()

    def plot_convergence(self,scale_clim = 1):
        fig2, ax2 = pt.subplots(nrows = 2, ncols = 2)
        im1_a = ax2[0,0].imshow(np.abs(self.all_measurements)*self.valid_channels,origin='upper', aspect = 'auto',interpolation='nearest')
        im1_p = ax2[0,1].imshow(np.angle(self.all_measurements)*self.valid_channels,origin='upper', aspect='auto',interpolation='nearest')
        for i, T in enumerate(self.T_list):
            re_projection = np.dot(self.geom_matrix, T)
            im1_p.set_clim([-np.pi,np.pi])
            q = self.all_measurements*0
            q[self.valid_channels]= re_projection
            im2_a = ax2[1,0].imshow(np.abs(q),origin='upper', aspect='auto',interpolation='nearest')
            im2_p = ax2[1,1].imshow(np.angle(q),origin='upper', aspect='auto',interpolation='nearest')
            if i == 0:
                clim_lower = im1_a.get_clim()[0]
                clim_upper = im1_a.get_clim()[1] * scale_clim
                clim = [clim_lower, clim_upper]
            im2_a.set_clim(clim)
            im1_a.set_clim(clim)
            im2_p.set_clim(im1_p.get_clim())
            fig2.savefig('{:02d}.png'.format(i))


class DirectSolution(Tomography):
    def __init__(self, geom_matrix, measurements, valid_channels = None, tomo_DC = None):
        n_measurements, n_regions = geom_matrix.shape
        if tomo_DC == None:
            self.geom_matrix = geom_matrix
        else:
            print tomo_DC.T.shape, geom_matrix.shape
            if geom_matrix.shape[1]%tomo_DC.T.shape[0]!=0: raise(ValueError)
            n_modes = geom_matrix.shape[1]/tomo_DC.T.shape[0]
            print n_modes
            self.geom_matrix = np.dot(geom_matrix, np.diag(np.hstack((tomo_DC.T for i in range(n_modes)))))
        self.all_measurements = measurements
        self.method = 'DirectSolution'
        if valid_channels!=None:
            self.measurements = measurements[valid_channels]
            self.valid_channels = valid_channels
        else:
            self.measurements = measurements
            self.valid_channels = np.ones(measurements.shape, dtype=bool)
    def run(self):
        #This allows measurements that only go through one or two surfaces to be removed:
        #Currently disabled
        measures_to_remove = (np.sum(self.geom_matrix>0, axis = 1))>=0
        print self.geom_matrix.shape[0], np.sum(measures_to_remove)
        self.geom_matrix_pinv = np.linalg.pinv(self.geom_matrix[measures_to_remove,:])
        valid_meas = self.all_measurements[self.valid_channels]
        self.T = np.dot(self.geom_matrix_pinv, valid_meas[measures_to_remove])
        self.re_projection = np.dot(self.geom_matrix, self.T)

        #self.error = calc_RMSE(self, plot = False)
        #self.error = tomo_recon_error_calc(self.geom_matrix, self.T, valid_meas)
        self.error = tomo_recon_error_NRMSE(self.geom_matrix, self.T, valid_meas)

class DirectSolutionFixedPhase(Tomography):
    def __init__(self, geom_matrix, measurements, valid_channels = None, tomo_DC = None):
        n_measurements, n_regions = geom_matrix.shape
        if tomo_DC == None:
            self.geom_matrix = geom_matrix
        else:
            print tomo_DC.T.shape, geom_matrix.shape
            if geom_matrix.shape[1]%tomo_DC.T.shape[0]!=0: raise(ValueError)
            n_modes = geom_matrix.shape[1]/tomo_DC.T.shape[0]
            print n_modes
            self.geom_matrix = np.dot(geom_matrix, np.diag(np.hstack((tomo_DC.T for i in range(n_modes)))))
        self.all_measurements = measurements
        self.method = 'DirectSolutionFixedPhase'
        if valid_channels!=None:
            self.measurements = measurements[valid_channels]
            self.valid_channels = valid_channels
        else:
            self.measurements = measurements
            self.valid_channels = np.ones(measurements.shape, dtype=bool)
        self.P = np.zeros((n_measurements*2, n_regions*2),dtype=float)
        #Build a real matrix 
        self.P[0::2,0::2] = +np.real(self.geom_matrix)
        self.P[0::2,1::2] = -np.imag(self.geom_matrix)
        self.P[1::2,0::2] = +np.imag(self.geom_matrix)
        self.P[1::2,1::2] = +np.real(self.geom_matrix)
        self.S = np.zeros(n_measurements*2, dtype=float)
        self.S[0::2] = +np.real(self.measurements)
        self.S[1::2] = +np.imag(self.measurements)
        #self.measurements = measurements[self.valid_channels]

    def run(self,):
        res = opt.minimize(fixed_phase_error, np.pi, method='nelder-mead', args = (self.P, self.S), options={'xtol': 1e-1, 'disp': False})
        geom_matrix, T, re_projection = fixed_phase_calc(res['x'], self.P, self.S)
        self.phasing = res['x']
        self.T = (T*np.cos(self.phasing) + T*1j*np.sin(self.phasing))
        self.phasing = np.angle(np.sum(self.T))
        print 'PHASING : ', self.phasing
        #tmp = np.dot(self.geom_matrix, T)
        self.re_projection = np.dot(self.geom_matrix, self.T)
        #self.re_projection = tmp[0::2]+1j*tmp[1::2]
        #print res

def fixed_phase_calc(phasing, geom_matrix, all_measurements):
    geom_matrix = geom_matrix[:,::2]*np.cos(phasing) + geom_matrix[:,1::2]*np.sin(phasing)
    geom_matrix_pinv = np.linalg.pinv(geom_matrix)
    T = np.dot(geom_matrix_pinv, all_measurements)
    re_projection = np.dot(geom_matrix, T)
    return geom_matrix, T, re_projection

def fixed_phase_error(phasing, geom_matrix, all_measurements):
    geom_matrix, T, re_projection = fixed_phase_calc(phasing, geom_matrix, all_measurements)
    error = tomo_recon_error_calc(geom_matrix, T, all_measurements)
    print phasing, error, 'finished'
    return error[0]

def tomo_recon_error_calc(geom_matrix, tomo_recon, valid_meas):
    '''Function to calculate the reconstruction error given the geom_matrix, reconstruction and measurements

    SRH: 7Feb2014
    '''
    reproj = np.dot(geom_matrix, tomo_recon)- valid_meas
    return [np.sqrt(np.mean(np.real(reproj)**2 + np.imag(reproj)**2) / np.mean(np.real(valid_meas)**2 + np.imag(valid_meas)**2))]

def tomo_recon_error_NRMSE(geom_matrix, tomo_recon, valid_meas, plot = False):
    '''Function to calculate the reconstruction error given the geom_matrix, reconstruction and measurements

    SRH: 7Feb2014
    '''
    reproj = np.dot(geom_matrix, tomo_recon)- valid_meas
    M_reproj = np.dot(geom_matrix, tomo_recon)
    NRMSD_list = []; RMS_ratio_list = []
    if plot: fig, ax = pt.subplots()
    for i in np.linspace(0,2.*np.pi,16):
        M_reproj_curr = np.real(M_reproj * np.exp(1j*i))
        meas_cur = np.real(valid_meas * np.exp(1j*i))
        diff = M_reproj_curr - meas_cur
        RMSD = np.sqrt(np.mean(diff**2))
        RMS = np.sqrt(np.mean(meas_cur**2))
        NRMSD = RMSD/(np.max(meas_cur) - np.min(meas_cur))
        #print RMSD, RMS, RMSD/RMS, NRMSD, RMSD/np.mean(np.abs(meas_cur))
        RMS_ratio_list.append(RMSD/RMS)
        NRMSD_list.append(NRMSD)
    if plot:
        ax.plot(np.linspace(0,2.*np.pi,50), NRMSD_list, label = 'NRMSD')
        ax.plot(np.linspace(0,2.*np.pi,50), RMS_ratio_list, label = 'RMSE/RMS')
        ax.set_xlabel('Phase')
        ax.set_ylabel('Error Measure')
        ax.legend(loc='best')
        fig.canvas.draw(); fig.show()
    return [np.mean(NRMSD)]



class ART(Tomography):
    def __init__(self, geom_matrix, measurements, lamda, initial_guess = None, save_increment = 50, produce_plots = 0, random_ordering = 1, valid_channels = None):
        '''This is an implementation of the algebraic reconstruction
        technique It will accept real or complex inputs. If the inputs are
        complex, the arrays are expanded out to twice the width and
        length, to accomodate complex numbers, while leaving the arrays
        real.

        lamda is a mixing coefficient
        SRH: 21July2013
        '''
        self.method = 'ART'
        n_measurements, n_regions = geom_matrix.shape
        self.random_ordering = random_ordering
        self.lamda = lamda
        self.geom_matrix = geom_matrix
        self.all_measurements = measurements
        if valid_channels!=None:
            self.measurements = measurements[valid_channels]
            self.valid_channels = valid_channels
        else:
            self.measurements = measurements
            self.valid_channels = np.ones(measurements.shape, dtyp=bool)
        self.initial_guess = initial_guess
        self.save_increment = save_increment
        self.produce_plots = produce_plots
        #If the geometry matrix or measurement matrix complex
        #Need to build a larger matrix that contains the linear system of equations
        #in real numbers only
        if geom_matrix.dtype == np.complex_ or measurements.dtype == np.complex:
            print 'complex input'
            self.P = np.zeros((n_measurements*2, n_regions*2),dtype=float)
            #Build a real matrix 
            self.P[0::2,0::2] = +np.real(self.geom_matrix)
            self.P[0::2,1::2] = -np.imag(self.geom_matrix)
            self.P[1::2,0::2] = +np.imag(self.geom_matrix)
            self.P[1::2,1::2] = +np.real(self.geom_matrix)
            self.S = np.zeros(n_measurements*2, dtype=float)
            self.S[0::2] = +np.real(self.measurements)
            self.S[1::2] = +np.imag(self.measurements)
            self.input_type = 'complex'
        else:
            print 'real input'
            self.P = self.geom_matrix
            self.S = self.measurements
            self.input_type = 'real'

    def run(self, initial_guess = None, cycles = 1, skip_if_norm_below=0.0001):
        self.initial_guess = initial_guess
        self.T_list = []
        if self.initial_guess == None:
            self.T = np.zeros((self.P.shape[1]), dtype = float)
        else:
            print 'using the initial guess'
            if self.initial_guess.dtype == np.complex_:
                self.T = np.zeros((self.P.shape[1]), dtype = self.geom_matrix.dtype)
                self.T[::2] = np.real(self.initial_guess)
                self.T[1::2] = np.imag(self.initial_guess)
            else:
                self.T = self.initial_guess*1.
        if self.random_ordering:
            mixture = np.random.rand(self.P.shape[0]*cycles)
            ordering = np.argsort(mixture)
        else:
            ordering = np.range(self.P.shape[0]*cycles)
        for k in ordering:
            i = k%self.P.shape[0]
            a_i = self.P[i,:]
            norm = np.sum(a_i**2)
            if norm>skip_if_norm_below:
                self.T = self.T + self.lamda* (self.S[i] - np.dot(a_i, self.T)) * a_i / norm
                if (i%self.save_increment) == 0:
                    if self.input_type == 'complex':
                        self.T_list.append(self.T[0::2] + 1j*self.T[1::2])
                    else:
                        self.T_list.append(self.T)
        if self.input_type == 'complex':
            self.T = self.T[0::2] + self.T[1::2]*1j
        self.re_projection = np.dot(self.geom_matrix, self.T)




class SIRT(Tomography):
    def __init__(self, geom_matrix, measurements, lamda, initial_guess = None, save_increment = 50, produce_plots = 0, random_ordering = 1, valid_channels = None, tomo_DC = None):
        '''This is an implementation of the algebraic reconstruction
        technique It will accept real or complex inputs. If the inputs are
        complex, the arrays are expanded out to twice the width and
        length, to accomodate complex numbers, while leaving the arrays
        real.

        lamda is a mixing coefficient
        SRH: 21July2013
        '''
        self.method = 'SIRT'
        n_measurements, n_regions = geom_matrix.shape
        self.random_ordering = random_ordering
        self.lamda = lamda
        self.geom_matrix = geom_matrix
        self.all_measurements = measurements
        if valid_channels!=None:
            self.measurements = measurements[valid_channels]
            self.valid_channels = valid_channels
        else:
            self.measurements = measurements
            self.valid_channels = np.ones(measurements.shape, dtyp=bool)
        self.initial_guess = initial_guess
        self.save_increment = save_increment
        self.produce_plots = produce_plots
        #If the geometry matrix or measurement matrix complex
        #Need to build a larger matrix that contains the linear system of equations
        #in real numbers only
        if geom_matrix.dtype == np.complex_ or measurements.dtype == np.complex:
            print 'complex input'
            self.P = np.zeros((n_measurements*2, n_regions*2),dtype=float)
            #Build a real matrix 
            self.P[0::2,0::2] = +np.real(self.geom_matrix)
            self.P[0::2,1::2] = -np.imag(self.geom_matrix)
            self.P[1::2,0::2] = +np.imag(self.geom_matrix)
            self.P[1::2,1::2] = +np.real(self.geom_matrix)

            #self.P[0::2,0::2] = +np.real(self.geom_matrix)
            #self.P[0::2,1::2] = +np.imag(self.geom_matrix)
            #self.P[1::2,0::2] = +np.imag(self.geom_matrix)
            #self.P[1::2,1::2] = -np.real(self.geom_matrix)
            self.S = np.zeros(n_measurements*2, dtype=float)
            self.S[0::2] = +np.real(self.measurements)
            self.S[1::2] = +np.imag(self.measurements)
            self.input_type = 'complex'
        else:
            print 'real input'
            self.P = self.geom_matrix
            self.S = self.measurements
            self.input_type = 'real'

    def run(self, initial_guess = None, cycles = 1, skip_if_norm_below=0.0001):
        self.initial_guess = initial_guess
        self.T_list = []
        if self.initial_guess == None:
            self.T = np.zeros((self.P.shape[1]), dtype = float)
        else:
            print 'using the initial guess'
            if self.initial_guess.dtype == np.complex_:
                self.T = np.zeros((self.P.shape[1]), dtype = self.geom_matrix.dtype)
                self.T[::2] = np.real(self.initial_guess)
                self.T[1::2] = np.imag(self.initial_guess)
            else:
                self.T = self.initial_guess*1.
        #if self.random_ordering:
        #    mixture = np.random.rand(self.P.shape[0]*cycles)
        #    ordering = np.argsort(mixture)
        #else:
        #    ordering = np.range(self.P.shape[0]*cycles)
        #for k in ordering:
        self.error = []
        for k in range(cycles):
            if (k%20)==0:
                self.error.append(tomo_recon_error_calc(self.P, self.T, self.S)[0])
                #tmp = np.dot(self.P,self.T)- self.S
                print('cycle {} of {}, error {:.2f}'.format(k,cycles, self.error[-1]))
            current_modification = self.T*0
            curr_count = 0
            for i in range(self.P.shape[0]):
                a_i = self.P[i,:]
                norm = np.sum(a_i**2)
                if norm>skip_if_norm_below:
                    current_modification +=  self.lamda* (self.S[i] - np.dot(a_i, self.T)) * a_i / norm
                    curr_count += 1
            self.T = self.T + current_modification / curr_count
            if (k%self.save_increment) == 0:
                if self.input_type == 'complex':
                    self.T_list.append(self.T[0::2] + 1j*self.T[1::2])
                else:
                    self.T_list.append(self.T)
        if self.input_type == 'complex':
            tmp = np.dot(self.P,self.T)
            self.re_projection = tmp[0::2]+1j*tmp[1::2]
            self.T = (self.T[0::2] + self.T[1::2]*1j)
            print 'hello'
        else:
            self.re_projection = np.dot(self.geom_matrix, self.T)

def single_channel(boozer_pts, interp_pts, dl, segments, interp_fact, n, m, scaling_factor = None, scaling_factor_s = None):
    '''
    finds the weightings for each flux region for a single LOS.
    boozer_pts is expected to be (s, theta, phi)
    interp_fact is an optional increase in the number of points for each line of sight to ensure each region is well sampled
    SRH : 12Feb2014
    '''
    mode_amps = np.exp(1j*n*boozer_pts[:,2] + 1j*m*boozer_pts[:,1])
    dl_tmp = dl
    s_tmp = boozer_pts[:,0]
    if interp_fact != None:
        l = np.arange(0,len(interp_pts[:,0]))*dl_tmp
        dl_tmp = dl_tmp / interp_fact
        l_new = np.arange(0,len(l)*interp_fact)*dl_tmp
        mode_amps = np.interp(l_new, l, np.real(mode_amps)) + 1j* np.interp(l_new, l, np.imag(mode_amps))
        s_tmp = np.interp(l_new, l, boozer_pts[:,0])
    if scaling_factor != None:
        print('hello')
        scale_fact = np.interp(s_tmp, scaling_factor_s, scaling_factor)
        mode_amps = mode_amps / scale_fact
    bin_allocation = np.digitize(s_tmp, segments)
    output_list = []
    for i in range(1,len(segments)):
        #output_list.append(np.sum(bin_allocation==i)*dl)
        output_list.append(np.sum(mode_amps[bin_allocation==i])*dl)
    return np.array(output_list)

# def calculate_inverse_matrix2(LOS_obj, n_segments, s_min =0., s_max = 1.,n_list=None, m_list=None, interp_fact = None, mode_amps_input=None, channel_mask=None, scaling_factor = None, scaling_factor_s = None):
#     '''Calculate the inverse matrix for the calculation on density profile
#     based on a fixed number of segments, and some mode numbers
#     SH: 2013
#     '''
#     print 'calculating geometric matrix and its pseudo inverse'
#     if n_list.__class__ == int:
#         n_list = [n_list]
#         m_list = [m_list]
#     segments = np.linspace(s_min, s_max, n_segments)
#     LOS_obj.segments = segments
#     LOS_obj.segment_midpoints = (segments[1:]+segments[0:-1])/2.
#     n_measurements = np.sum(LOS_obj.valid_channels)
#     geom_mat_list = []
#     for i in n_list:
#         geom_mat_list.append(np.zeros((n_measurements,n_segments-1), dtype=complex))
#     #geom_mat2[channel_count, :] = output_array*1.
#     #s=0, theta = 1, phi = 2
#     pixel_list = []
#     for geom_mat, n, m in zip(geom_mat_list, n_list, m_list):
#         channel_count = 0
#         for pixel_y in range(LOS_obj.CCD_pixels_y):
#             for pixel_x in range(LOS_obj.CCD_pixels_x):
#                 if LOS_obj.valid_channels[pixel_y, pixel_x]:
#                     interp_pts = LOS_obj.interpolation_pts[pixel_y, pixel_x,:,:]
#                     boozer_pts = LOS_obj.interp_boozer[pixel_y, pixel_x,:,:]
#                     dl = LOS_obj.dl[pixel_y, pixel_x]
#                     output_array = single_channel(boozer_pts, interp_pts, dl, segments, interp_fact, n, m, scaling_factor = scaling_factor, scaling_factor_s = scaling_factor_s)
#                     geom_mat[channel_count, :] = output_array*1.
#                     channel_count += 1
#                     #print channel_count
#                     pixel_list.append([pixel_y, pixel_x])
#     geom_mat_comb = np.zeros((geom_mat_list[0].shape[0],geom_mat_list[0].shape[1]*len(geom_mat_list)) ,dtype=complex)
#     start = 0
#     for geom_mat in geom_mat_list:
#         end = start + n_segments - 1
#         geom_mat_comb[:, start:end] = +geom_mat
#         start = +end
#     return geom_mat_comb, pixel_list, segments


def calculate_inverse_matrix_sep_images(LOS_obj, n_segments, s_min =0., s_max = 1.,n=None, m=None, interp_fact = None, scaling_factor = None, scaling_factor_s = None, radial_s_spacing = False, dict_of_values= None):
    '''Calculate the geometry matrix for multiple views, for a given n and m

    SRH: 12Feb2014
    '''
    #Create the discrete flux surfaces
    segments = np.linspace(s_min, s_max, n_segments)
    if radial_s_spacing: segments = segments**2
    LOS_obj.segments = segments
    LOS_obj.radial_s_spacing = radial_s_spacing
    LOS_obj.segment_midpoints = (segments[1:]+segments[0:-1])/2.

    #Number of valid LOS for each view
    n_measurement_list = [np.sum(LOS_obj.valid_channels[start_ind:end_ind,:]) for start_ind, end_ind in zip(LOS_obj.start_indices, LOS_obj.end_indices)]

    geom_mat_list = []
    for i in n_measurement_list: geom_mat_list.append(np.zeros((i,n_segments-1), dtype=complex))

    #find the proportion in each region for each LOS that is valid
    for geom_mat, start_ind, end_ind, n_measurement, view in zip(geom_mat_list, LOS_obj.start_indices, LOS_obj.end_indices, n_measurement_list, LOS_obj.orientations):
        channel_count = 0
        for pixel_y in range(start_ind, end_ind):#LOS_obj.CCD_pixels_y):
            for pixel_x in range(LOS_obj.CCD_pixels_x):
                if LOS_obj.valid_channels[pixel_y, pixel_x]:
                    interp_pts = LOS_obj.interpolation_pts[pixel_y, pixel_x,:,:]
                    boozer_pts = LOS_obj.interp_boozer[pixel_y, pixel_x,:,:]
                    dl = LOS_obj.dl[pixel_y, pixel_x]
                    output_array = single_channel(boozer_pts, interp_pts, dl, segments, interp_fact, n, m, scaling_factor = scaling_factor, scaling_factor_s = scaling_factor_s)
                    geom_mat[channel_count, :] = output_array*1.
                    channel_count += 1
        if channel_count!=n_measurement: raise(Exception)
        if dict_of_values!=None:dict_of_values['{}_{}_{}'.format(n, m, view)]=geom_mat
    return geom_mat_list, segments


def intersection_points(R_values, Z_values, gradient, offset, plot = 0):
    '''Find the intersection of the line with gradient and offset
    and the surface defined by R_values, and Z_values
    '''
    a1 = gradient
    b1 = offset
    #Calculate the intersection of LOS with all the lines between
    #points (note R_values and Z_values must be in order around the
    #surface)
    a_vals = np.diff(Z_values)/np.diff(R_values)
    b_vals = Z_values[1:] - a_vals*R_values[1:]
    r_int = (b_vals -b1)/(a1 - a_vals)
    z_int = a1*r_int + b1
    #Check to see if the intersection happens where the line is valid (i.e on the surface)
    #Slight shift in points that have identical adjacent points (makes problems for truth1 and truth2)
    R_values[np.abs(R_values[1:] == R_values[0:-1])] +=0.00001 #<0.00001]+=0.00001
    Z_values[np.abs(Z_values[1:] == Z_values[0:-1])] +=0.00001 #<0.00001]+=0.00001
    truth1 = (r_int<R_values[1:]) * (r_int>R_values[0:-1]) + (r_int>R_values[1:]) * (r_int<R_values[0:-1])
    truth2 = (z_int<Z_values[1:]) * (z_int>Z_values[0:-1]) + (z_int>Z_values[1:]) * (z_int<Z_values[0:-1])
    truth3 = truth1 * truth2
    n_intersections = np.sum(truth3)
    if plot:
        fig, ax = pt.subplots()
        ax.plot(R_values, Z_values, '.-')
        ax.plot(rs, zs, '-')
        ax.plot(R_values[0], Z_values[0],'o')
        ax.plot(r_int[truth3], z_int[truth3],'x', markersize=10)
        fig.canvas.draw(); fig.show()
    if n_intersections>2:
        print 'possible problems - found 3 intersection points!!??'
    elif n_intersections == 1:
        print 'possible problems - found 1 intersection points!!??', np.sum(truth1), np.sum(truth2)
    elif n_intersections == 0:
        print 'no intersection with surface'
    return r_int[truth3], z_int[truth3]

class interferometer:
    def __init__(self, z_inter, heliac_object, gradient = None, method = 1, pts_along_line = 300):#, R_values, Z_values, label_values, psi_values):
        '''z_inter are the elevations of the channels
        heliac_object is a heliac ovject from heliac_worker_funcs file
        '''
        self.valid_LOS = np.zeros(len(z_inter), dtype = bool)
        print('hello')
        if gradient == None:
            print('gradient was none')
            self.gradient = np.zeros(z_inter.shape,dtype=float)
        else:
            self.gradient = gradient
        self.R_values = heliac_object.ordered_R
        self.Z_values = heliac_object.ordered_Z
        self.psi_values = heliac_object.valid_points_psi_norm#/ np.max(heliac_object.valid_points_psi_norm)
        self.label_values = heliac_object.ordered_labels
        self.z_inter = z_inter
        self.dl_list = [];self.s_list = []; self.z_list = []; self.r_list = []
        self.r_intersects = []; self.z_intersects = []
        for z, gradient, LOS_index in zip(self.z_inter, self.gradient, range(len(self.z_inter))):
            #extract the last surface and find intersection points
            r_cur, z_cur, label_cur, psi_values_cur  = zip(self.R_values, self.Z_values, self.label_values, self.psi_values)[-1]
            ans_x, ans_y, intersect, success = self.find_intersection(r_cur,z_cur,z)
            r_intersects, z_intersects = intersection_points(r_cur, z_cur, gradient, z, plot = 0)
            self.r_intersects.extend(r_intersects)
            self.z_intersects.extend(z_intersects)
            #make sure there is an intersection!
            if (len(intersect) == 2) and (method==1):
                #interpolation points, note slight offset 0.005 so that we are INTERpolating
                r_vals = np.linspace(np.min(intersect)+0.005,np.max(intersect)-0.005,pts_along_line)
                z_vals = r_vals*0+z
                self.dl_list.append(r_vals[1]-r_vals[0])
                self.s_list.append(scipy_griddata((np.array(self.R_values).flatten(), np.array(self.Z_values).flatten()), np.array(self.psi_values).flatten(),(r_vals,z_vals)))
                self.r_list.append(r_vals)
                self.z_list.append(z_vals)
                self.valid_LOS[LOS_index] = True
            elif len(r_intersects) == 2:
                r_vals = np.linspace(r_intersects[0], r_intersects[1], pts_along_line)
                z_vals = np.linspace(z_intersects[0], z_intersects[1], pts_along_line)
                self.dl_list.append(np.sqrt((r_vals[1]-r_vals[0])**2 + (z_vals[1]-z_vals[0])**2))
                #self.s_list.append(scipy_griddata((np.array(self.R_values).flatten(), np.array(self.Z_values).flatten()), np.array(self.psi_values).flatten(),(r_vals,z_vals)))
                self.r_list.append(r_vals)
                self.z_list.append(z_vals)
                self.valid_LOS[LOS_index] = True
        #self.s_list_poly = copy.deepcopy(self.s_list)
        self.dl_list_poly = copy.deepcopy(self.dl_list)


    def plot_eigenfunction(self,n,m,eig_func,im_ax=None, eig_ax=None):
        # booz_s_vals = np.array(self.s_list)
        # booz_th_vals = np.array(self.th_booz_list)
        # booz_phi_vals = np.array(self.phi_booz_list)
        # r_list = np.array(self.r_list_trans)
        # z_list = np.array(self.z_list_trans)

        # s_vals = self.cross_sect_booz[:,0,0]
        # s_amps = np.interp(booz_s_vals,eig_func[0],np.real(eig_func[1]))+ np.interp(booz_s_vals,eig_func[0],np.imag(eig_func[1]))*1j
        # s_amps = np.interp(s_vals,eig_func[0],np.real(eig_func[1]))+ np.interp(s_vals,eig_func[0],np.imag(eig_func[1]))*1j
        # wave_amp = np.zeros((self.cross_sect_booz.shape[0],self.cross_sect_booz.shape[2]),dtype=float)
        # for i in range(len(s_amps)):
        #     wave_amp[i,:] = np.real(s_amps[i] * np.exp(1j*n*self.cross_sect_booz[i,1,:] + 1j*m*self.cross_sect_booz[i,2,:]))
        # fig_tmp,ax_tmp = pt.subplots()
        # wave = np.real(s_amps*np.exp(1j*n*booz_phi_vals + 1j*m*booz_th_vals))

        #Get the boozer locations to interpolate
        phi_booz_list = self.cross_sect_booz[:,1,:].flatten()
        th_booz_list = self.cross_sect_booz[:,2,:].flatten()
        s_list = self.cross_sect_booz[:,0,:].flatten()
        r_list_trans = self.cross_sect_cyl[:,0,:].flatten()
        z_list_trans = self.cross_sect_cyl[:,2,:].flatten()

        #grid of points so we can use imshow
        r_vals = np.linspace(1.01,1.37,150)
        z_vals = np.linspace(-0.24,0.24,150)
        r_grid, z_grid = np.meshgrid(r_vals,z_vals)

        #find boozer coordinates of the grid
        s = interp.griddata((r_list_trans,z_list_trans),s_list,(r_grid.flatten(),z_grid.flatten()))

        phi_sin = interp.griddata((r_list_trans,z_list_trans),np.sin(phi_booz_list),(r_grid.flatten(),z_grid.flatten()))
        phi_cos = interp.griddata((r_list_trans,z_list_trans),np.cos(phi_booz_list),(r_grid.flatten(),z_grid.flatten()))
        phi = np.arctan2(phi_sin, phi_cos)

        theta_sin = interp.griddata((r_list_trans,z_list_trans),np.sin(th_booz_list),(r_grid.flatten(),z_grid.flatten()))
        theta_cos = interp.griddata((r_list_trans,z_list_trans),np.cos(th_booz_list),(r_grid.flatten(),z_grid.flatten()))
        theta = np.arctan2(theta_sin,theta_cos)

        #find the mode amplitudes at the grid points
        #calc the wave field
        s_amps = np.interp(s,eig_func[0],np.real(eig_func[1])) + np.interp(s,eig_func[0],np.imag(eig_func[1]))*1j
        wave = np.real(s_amps * np.exp(1j*n*phi + 1j*m*theta))
        wave = np.resize(wave,z_grid.shape)


        L,R = np.min(r_vals), np.max(r_vals)
        B,T = np.min(z_vals), np.max(z_vals)

        # r_vals = np.linspace(1.01,1.37,150)
        # z_vals = np.linspace(-0.24,0.24,150)
        # r_grid, z_grid = np.meshgrid(r_vals,z_vals)
        # tmp = interp.griddata((r_list,z_list),wave,(r_grid.flatten(),z_grid.flatten()))
        # tmp = np.resize(tmp,z_grid.shape)
        # L,R = np.min(r_vals), np.max(r_vals)
        # B,T = np.min(z_vals), np.max(z_vals)
        if im_ax!=None:
            clr =im_ax.imshow(wave,aspect='auto',origin='lower',extent=[L,R,B,T])
            #im_ax.plot(r_list,z_list,linewidth=0.1)
            im_ax.set_xlabel('r')
            im_ax.set_ylabel('z (simulated mode at t=0')
        if eig_ax!=None:
            eig_ax.plot(s, s_amps)
            eig_ax.set_xlabel('s')
            eig_ax.set_ylabel('disp (CAS3D eigenfunc)')


    #signal = np.sum(s_amps * np.exp(1j* n * phi_valid) *np.exp(1j* m* theta_valid)*dl)/(np.sum(valid_points)*dl)

    def simulate_channel_outputs(self, eig_func, s_vals, n = -4, m = 3):
        '''
        Given an eigenfunction as a function of s, simulate the LOS outputs

        SH:27June2013
        '''
        mode_output = []
        print 'hello'
        for i in range(len(self.LOS_coords)):
            #print 'new channel %d of %d'%(i,len(self.LOS_coords))
            curr_coords = self.LOS_coords[i]
            valid_points = np.isfinite(curr_coords[:,0]) * np.isfinite(curr_coords[:,1]) * np.isfinite(curr_coords[:,2])
            if np.sum(valid_points)>0:
                dl = self.dl_list[i]
                s_valid = curr_coords[valid_points,0]
                theta_valid = curr_coords[valid_points,1]
                phi_valid = curr_coords[valid_points,2]
                #print curr_coords, np.sum(valid_points)
                #s_amps = s_valid*0+1
                #eig_func = [s_valid,s_amps]
                s_amps = np.interp(s_valid,s_vals,np.real(eig_func)) + np.interp(s_valid,s_vals,np.imag(eig_func))*1j
                signal = np.sum(s_amps * np.exp(1j* n * phi_valid + 1j* m* theta_valid)*dl)#/(np.sum(valid_points)*dl)
            else:
                signal = 0.
            mode_output.append(signal)
        return mode_output

    def output_mode_signals(self,n = -4, m = 3, produce_plot=False, eig_func=None, pixels=512):
        '''
        Starting to try and get a simulated interferometer output for a mode
        based on the geometry of the line of sight
        SH:5May2013
        '''
        mode_output = []; measurements_complex = []
        for i in range(len(self.LOS_coords)):
            #print 'new channel %d of %d'%(i,len(self.LOS_coords))
            curr_coords = self.LOS_coords[i]
            valid_points = np.invert(np.isnan(curr_coords[:,0])) * np.invert(np.isnan(curr_coords[:,1])) * np.invert(np.isnan(curr_coords[:,2]))
            dl = self.dl_list[i]
            s_valid = curr_coords[valid_points,0]
            theta_valid = curr_coords[valid_points,1]
            phi_valid = curr_coords[valid_points,2]
            #print curr_coords, np.sum(valid_points)
            if eig_func == None:
                s_amps = s_valid*0+1
                eig_func = [s_valid,s_amps]
            s_amps = np.interp(s_valid,eig_func[0],np.real(eig_func[1])) + np.interp(s_valid,eig_func[0],np.imag(eig_func[1]))*1j
            signal = np.sum(s_amps * np.exp(1j* n * phi_valid + 1j* m* theta_valid)*dl)#/(np.sum(valid_points)*dl)
            mode_output.append(signal)
        #print mode_output
        #fit_mode(self,mode_output,n,m)
        self.mode_output = mode_output
        if produce_plot:
            fig, ax = pt.subplots(nrows = 2,ncols=2);ax = ax.flatten()
            if pixels!=0:
                x_axis = np.arange(pixels)[self.valid_LOS]
                y_axis = np.array(self.mode_output)
                xlim = [0,pixels]
            else:
                x_axis = [self.z_list[i][np.argmax(self.r_list[i])] for i in range(len(self.r_list))]
                y_axis = self.mode_output
                xlim=[-0.25,0.25]
            #ax[0].plot([i[0] for i in self.z_list],np.abs(self.mode_output),'o')
            #ax[1].plot([i[0] for i in self.z_list],np.angle(self.mode_output)/np.pi*180.,'o')
            ax[0].plot(x_axis,np.abs(y_axis),'o')
            ax[1].plot(x_axis,np.angle(y_axis)/np.pi*180.,'o')
            #ax[0].set_ylim(0,1)
            ax[1].set_ylim(-200,200)
            ax[0].set_xlim(xlim)
            ax[0].set_ylabel('Amplitude')
            ax[1].set_ylabel('Phase')
            ax[1].set_xlabel('z (m)')
            fig.suptitle('n:{n}, m:{m}'.format(n=n, m=m))
            self.plot_eigenfunction(n,m,eig_func,im_ax=ax[2], eig_ax=ax[3])
            fig.canvas.draw(); fig.show()

    def get_LOS_boozer_coords(self,boozer_filename=None, divide_by = 1, method='interpolate', surface_data_filename = None, read_surface_data = False, save_surface_data = False, save_surface_filename=None, no_theta = 100, create_plot = False, phi_val=0):
        '''
        First attempt at getting the probe coordinates in Boozer coordinates
        using Bernhards tools. Need to make this faster somehow
        
        if method is interpolate we try to interpolate the interferometer points onto the surface
        otherwise, we do the inversion for each individual point which is slow...
        SH:5May2013
        '''
        if method!='interpolate':
            import h1.mhd_eq.BOOZER as BOOZER
            self.boozer_object = BOOZER.BOOZER(boozer_filename,import_all=True,compute_spline_type=1,load_spline=False,save_spline=False,load_grid=False)
            output_data = []
            for i in range(len(self.r_list)):
                start_time = time.time()
                print i, len(self.r_list[i])
                RZ_input_data = np.zeros((len(self.r_list[i])/divide_by,2),dtype=float)
                RZ_input_data[:,0] = self.r_list[i][::divide_by]
                RZ_input_data[:,1] = self.z_list[i][::divide_by]
                a = self.boozer_object.real2Mag(0,RZ_input_data)
                print 'time to finish one chord : %.2f'%(time.time() - start_time)
                output_data.append(a)
            self.LOS_coords = output_data
        else:
            #either calculate the data or get it from a saved pickle file
            if read_surface_data:
                self.s_list, self.phi_booz_list, self.th_booz_list, self.r_list_trans, self.z_list_trans,self.cross_sect_booz, self.cross_sect_cyl = pickle.load(file(surface_data_filename,'r'))
            else:
                print 'getting the realspace poloidal cross-section in Boozer coords for interpolation'
                import h1.mhd_eq.BOOZER as BOOZER
                self.boozer_object = BOOZER.BOOZER(boozer_filename,import_all=True,compute_spline_type=1,load_spline=False,save_spline=False,load_grid=False)
                n_surfaces = len(self.boozer_object.es_b)
                start_time = time.time()
                self.cross_sect_cyl, self.cross_sect_booz = self.boozer_object.getCrossSectData(phi_val,s_b=None, s_ind=range(0,n_surfaces),no_theta=no_theta, phi_is_cyl=True,coordsys='cyl', return_booz_also=1)
                print 'time_taken {taken}'.format(taken=time.time() - start_time)
                #Flatten the output
                self.s_list = []; self.phi_booz_list = []; self.th_booz_list = []; self.r_list_trans = []; self.z_list_trans = []
                for i in range(self.cross_sect_booz.shape[0]):
                    self.s_list.extend(self.cross_sect_booz[i,0,:].tolist())
                    self.phi_booz_list.extend(self.cross_sect_booz[i,1,:].tolist())
                    self.th_booz_list.extend(self.cross_sect_booz[i,2,:].tolist())
                    self.r_list_trans.extend(self.cross_sect_cyl[i,0,:].tolist())
                    self.z_list_trans.extend(self.cross_sect_cyl[i,2,:].tolist())
                if save_surface_data:
                    pickle.dump((self.s_list, self.phi_booz_list, self.th_booz_list, self.r_list_trans, self.z_list_trans, self.cross_sect_booz, self.cross_sect_cyl), file(surface_data_filename,'w'))

            #now do the interpolation
            self.intersections_w_boozer(200)
            output_data =[]

            #build a single r_vals and z_vals
            z_vals_complete = []; r_vals_complete = []

            print 'trying new method'
            output_shape = np.array(self.r_list).shape
            z_flat = np.array(self.z_list).flatten()
            r_flat = np.array(self.r_list).flatten()
            interp_data_s = scipy_griddata((np.array(self.r_list_trans), np.array(self.z_list_trans)), np.array(self.s_list), (r_flat, z_flat))
            interp_data_phi = scipy_griddata((np.array(self.r_list_trans), np.array(self.z_list_trans)), np.array(self.phi_booz_list), (r_flat, z_flat))
            #note the sin and cosine are to get around a problem with interpolation at the wrapping points
            interp_data_theta_sin = scipy_griddata((np.array(self.r_list_trans), np.array(self.z_list_trans)), np.array(np.sin(self.th_booz_list)), (r_flat, z_flat))
            interp_data_theta_cos = scipy_griddata((np.array(self.r_list_trans), np.array(self.z_list_trans)), np.array(np.cos(self.th_booz_list)), (r_flat, z_flat))
            interp_data_theta = np.arctan2(interp_data_theta_sin, interp_data_theta_cos)
            interp_data_theta[interp_data_theta<-0.01]+=2.*np.pi
            interp_data_s_new = interp_data_s.reshape(output_shape)
            interp_data_theta_new = interp_data_theta.reshape(output_shape)
            interp_data_phi_new = interp_data_phi.reshape(output_shape)
            print 'finished'
            for i in range(len(self.r_list)):
                assembled_data = np.zeros((output_shape[1], 3),dtype=float)
                assembled_data[:,0] = interp_data_s_new[i,:]
                assembled_data[:,1] = interp_data_theta_new[i,:]
                assembled_data[:,2] = interp_data_phi_new[i,:]
                output_data.append(assembled_data)

            # for i in range(len(self.r_list)):
            #     print('interpolating Boozer {} of {}'.format(i,len(self.r_list)))
            #     r_vals = self.r_list[i]
            #     assembled_data = np.zeros((len(r_vals),3),dtype=float)
            #     z_vals = self.z_list[i]
            #     interp_data_s = scipy_griddata((np.array(self.r_list_trans), np.array(self.z_list_trans)), np.array(self.s_list), (r_vals, z_vals))
            #     interp_data_phi = scipy_griddata((np.array(self.r_list_trans), np.array(self.z_list_trans)), np.array(self.phi_booz_list), (r_vals, z_vals))
            #     #note the sin and cosine are to get around a problem with interpolation at the wrapping points
            #     interp_data_theta_sin = scipy_griddata((np.array(self.r_list_trans), np.array(self.z_list_trans)), np.array(np.sin(self.th_booz_list)), (r_vals, z_vals))
            #     interp_data_theta_cos = scipy_griddata((np.array(self.r_list_trans), np.array(self.z_list_trans)), np.array(np.cos(self.th_booz_list)), (r_vals, z_vals))
            #     interp_data_theta = np.arctan2(interp_data_theta_sin, interp_data_theta_cos)
            #     interp_data_theta[interp_data_theta<-0.01]+=2.*np.pi
            #     valid = np.isfinite(interp_data_s)
            #     print np.sum(valid)
            #     print np.allclose(interp_data_theta[valid], interp_data_theta_new[i,valid])
            #     print np.allclose(interp_data_s[valid], interp_data_s_new[i,valid])
            #     print np.allclose(interp_data_phi[valid], interp_data_phi_new[i,valid])
            #     assembled_data[:,0] = interp_data_s
            #     assembled_data[:,1] = interp_data_theta
            #     assembled_data[:,2] = interp_data_phi
            #     output_data.append(assembled_data)
            self.LOS_coords = output_data
        if create_plot:
            fig_tmp = pt.figure()
            ax_tmp = []
            ax_tmp.append(fig_tmp.add_subplot(1,2,1))
            ax_tmp.append(fig_tmp.add_subplot(3,2,2,sharex = ax_tmp[0]))
            ax_tmp.append(fig_tmp.add_subplot(3,2,4,sharex = ax_tmp[0]))
            ax_tmp.append(fig_tmp.add_subplot(3,2,6,sharex = ax_tmp[0]))
            ax_tmp[0].plot(self.r_list_trans, self.z_list_trans, '-', linewidth=0.2)
            for i in range(len(self.r_list)):
                r_vals = self.r_list[i]
                z_vals = self.z_list[i]
                ax_tmp[0].plot(r_vals, z_vals, 'k-',linewidth=0.2)
                ax_tmp[1].plot(r_vals, self.LOS_coords[i][:,0],linewidth=0.2)
                ax_tmp[2].plot(r_vals, self.LOS_coords[i][:,1], linewidth=0.2)
                ax_tmp[3].plot(r_vals, self.LOS_coords[i][:,2], linewidth=0.2)

            ax_tmp[0].set_ylabel('z(m)')
            ax_tmp[0].set_xlabel('r(m)')
            ax_tmp[1].set_ylabel('s')
            ax_tmp[2].set_ylabel('theta_B')
            ax_tmp[3].set_ylabel('phi_B')
            ax_tmp[3].set_xlabel('r (m)')

            fig_tmp.canvas.draw(); fig_tmp.show()
        #return output_data

    def intersections_w_boozer(self,pts_along_line):
        self.booz_r_LCF = self.cross_sect_cyl[-1,0,:]
        self.booz_z_LCF = self.cross_sect_cyl[-1,2,:]
        self.dl_list = [];self.z_list = []; self.r_list = []
        self.r_intersects = []; self.z_intersects = []
        fig, ax = pt.subplots()
        ax.plot(self.booz_r_LCF, self.booz_z_LCF,'.-')
        self.valid_LOS = np.zeros(len(self.z_inter), dtype = bool)
        for z, gradient, LOS_index in zip(self.z_inter, self.gradient, range(len(self.z_inter))):
            #extract the last surface and find intersection points
            r_intersects, z_intersects = intersection_points(self.booz_r_LCF, self.booz_z_LCF, gradient, z, plot = 0)
            self.r_intersects.extend(r_intersects)
            self.z_intersects.extend(z_intersects)
            ax.plot(r_intersects, z_intersects,'o-')
            #make sure there is an intersection!
            if len(r_intersects) == 2:
                r_vals = np.linspace(r_intersects[0], r_intersects[1], pts_along_line)
                z_vals = np.linspace(z_intersects[0], z_intersects[1], pts_along_line)
                self.dl_list.append(np.sqrt((r_vals[1]-r_vals[0])**2 + (z_vals[1]-z_vals[0])**2))
                self.r_list.append(r_vals)
                self.z_list.append(z_vals)
                self.valid_LOS[LOS_index] = True


    def plot_heliac_surface_channel_psi(self):
        fig, ax = pt.subplots(ncols = 2)
        for r_cur, z_cur, label_cur in zip(self.R_values, self.Z_values, self.label_values):
            ax[0].plot(r_cur, z_cur,'-')
        for z in self.z_inter:
            for r_cur, z_cur, label_cur in zip(self.R_values, self.Z_values, self.label_values):
                ans_x, ans_y, intersect, success = self.find_intersection(r_cur,z_cur,z)
                ax[0].plot(ans_x,ans_y,'x')
                ax[0].plot(intersect,intersect*0+z,'o')
        ax[0].hlines(self.z_inter,ax[0].get_xlim()[0], ax[0].get_xlim()[1])
        for i in range(len(self.s_list)):
            ax[1].plot(self.r_list[i],self.s_list[i],'o-')
        fig.canvas.draw(); fig.show()

    def plot_heliac_surface_channel_psi2(self):
        fig, ax = pt.subplots(ncols = 2)
        for r_cur, z_cur, label_cur in zip(self.R_values, self.Z_values, self.label_values):
            ax[0].plot(r_cur, z_cur,'-')
        for z, gradient, r_LOS, z_LOS in zip(self.z_inter, self.gradient, self.r_list, self.z_list):
            for r_cur, z_cur, label_cur in zip(self.R_values, self.Z_values, self.label_values):
                ans_x, ans_y, intersect, success = self.find_intersection(r_cur,z_cur,z)
            ax[0].plot(r_LOS, z_LOS,'-')
                #ax[0].plot(ans_x,ans_y,'x')
                #ax[0].plot(intersect,intersect*0+z,'o')
        #ax[0].hlines(self.z_inter,ax[0].get_xlim()[0], ax[0].get_xlim()[1])
        ax[0].plot(self.r_intersects,self.z_intersects,'o')
        for i in range(len(self.s_list)):
            ax[1].plot(self.r_list[i],self.s_list[i],'o-')
        fig.canvas.draw(); fig.show()


    def find_intersection(self, R, Z, z):
        ''' Finds the point of intersection between a bean surface defined by
        R and Z, and the horizontal line z (meant to represent an interferometer
        line of sight
        Returns the two closest points in flip_r and flip_z, as well as the intersection
        Also returns success to say if the surface is intersected by the line
        This works by forming a line between the two points closest to the intersection
        then working out the intersection point with z
        SH : 8Apr2013
        '''
        old_sign = np.copysign(1,Z[0]-z)
        flip_count = 0
        flip_r=[];flip_z=[]; intersect = []
        success=0
        for i in range(1,len(Z)):
            sign = np.copysign(1,Z[i]-z)
            if sign != old_sign:
                #print "FLIPPED"
                flip_count+=1
                #print flip_count, R[i],Z[i], R[i-1],Z[i-1]
                flip_r.append(R[i]); flip_r.append(R[i-1])
                flip_z.append(Z[i]); flip_z.append(Z[i-1])
                old_sign=sign
                a, b = np.polyfit([R[i],R[i-1]],[Z[i],Z[i-1]],1)
                intersect.append((z-b)/a)
        if len(intersect)==2:success=1
        return flip_r, flip_z, np.array(intersect), success

    def calculate_profile_bessel_fit_positive(self, amps):#, measurement):
        '''Calculate the interferometer output based on 
        a bessel representation of psi
        Need the svalues for the line of sight of the interferometer channels and the dl between them
        SH : 8Apr2013
        '''
        return np.min(self.calculate_profile_bessel(amps,np.linspace(0,1,50)))


    def calculate_profile_bessel(self, amps, s_vals):
        '''Calculate the interferometer output based on 
        a bessel representation of psi
        Need the svalues for the line of sight of the interferometer channels and the dl between them
        SH : 8Apr2013
        '''
        m_list, l_list = self.m_l
        count = 0
        for m in m_list:
            for l in l_list:
                if count ==0:
                    profile = amps[count]*scipy.special.jn(m, s_vals*self.bessel_zeros[m,l])
                else:
                    profile =profile +  amps[count]*scipy.special.jn(m, s_vals*self.bessel_zeros[m,l])
                count+=1
        return profile
        #return np.sum(np.polyval(poly, s_vals)*dl)/(dl*len(s_vals))

    def calculate_measurements_bessel(self, amps, channel):
        '''Calculate the interferometer output based on 
        a shifted zeroth order bessel representation of psi
        
        SH : 8Apr2013
        '''
        s_vals = self.s_list[channel]
        dl = self.dl_list[channel]
        m_list, l_list = self.m_l
        length = dl*len(s_vals)/dl
        total2 = np.sum(np.array(amps)*self.bessel_dict[channel]['sum']/length)
        return total2

    def setup_bessel(self, m_l):
        '''Calculate the common values for the bessel functions that are 
        repeatedly used - call this before doing any further bessel calculations

        SH: 9Apr2013
        '''
        self.m_l = m_l; self.s_list
        m_list, l_list = self.m_l
        #generate the list of bessel zeros for future calculations, include a few extra
        #m and l values
        bessel_zeros = []
        for i in range(np.max(m_list)+2):
            bessel_zeros.append(scipy.special.jn_zeros(i,np.max(l_list)+3))
        self.bessel_zeros = np.array(bessel_zeros)
        #generate the bessel basis functions for future use and store in a dictionary
        input_dict = {}; 
        for channel,s_vals in enumerate(self.s_list):
            bessel_array = np.zeros((len(m_list)*len(l_list),len(s_vals)),dtype=float)
            count=0            
            input_dict[channel]={}
            for m in m_list:
                input_dict[channel][m]={}
                for l in l_list:
                    tmp = scipy.special.jn(m, s_vals*self.bessel_zeros[m,l])
                    input_dict[channel][m][l] = tmp
                    bessel_array[count,:] = tmp
                    count+=1
            input_dict[channel]['mat']=bessel_array
            input_dict[channel]['sum']=np.sum(bessel_array,axis=1)
        self.bessel_dict = input_dict

    def calc_diff_bessel(self, amps, measurement):
        '''Calculate an error between the actual measurement and the
        estimated value from calculate_measurement
        measurement are the experimental values
        s_list is a list of the psi values for line of sight of each channel
        dl_list is a list of dl values for each channel
        amps are the current best fit bessel amplitudes to psi
        SH: 8Apr2013
        '''
        predicted = []
        for channel in range(0,len(self.s_list)):
            predicted.append(self.calculate_measurements_bessel(amps, channel))
        predicted = np.array(predicted)
        return np.sum((measurement - np.array(predicted))**2)

    def fit_bessel(self, index, increment, amps):
        '''Find the bessel amplitudes that fit best to the averaged data from index -> index+increment
        in the interferometer data
        SH: 9Apr2013
        '''
        start_time = time.time()
        measurement = np.average(self.dens_data[:,index:index+increment],axis = 1)
        q, fopt, func_calls, iterations, warnflag = optimize.fmin(self.calc_diff_bessel,amps,args=(measurement,),disp=0,full_output=1)
        if warnflag!=0:
            print 'Possible errors optimising, flag :', warnflag
        print 'time :', time.time() - start_time
        predictions = [self.calculate_measurements_bessel(q,channel) for channel in range(len(self.s_list))]
        profile = self.calculate_profile_bessel(q, np.linspace(0,1,50))
        return q, fopt, measurement, predictions, profile

    def fit_bessel2(self, index, increment, amps, method='Nelder-Mead', force_positive =  0):
        '''Find the bessel amplitudes that fit best to the averaged data from index -> index+increment
        in the interferometer data
        SH: 9Apr2013
        '''
        start_time = time.time()
        measurement = np.average(self.dens_data[:,index:index+increment],axis = 1)
        if force_positive:
            print 'force positive'
            constraints = ({'type': 'ineq',
                            'fun': self.calculate_profile_bessel_fit_positive})
            res = optimize.minimize(self.calc_diff_bessel, amps, args=(measurement,), method=method, options={'disp': False}, constraints = constraints)

        else:
            print 'not force positive'
            res = optimize.minimize(self.calc_diff_bessel, amps, args=(measurement,), method=method, options={'disp': False})



        #res = optimize.minimize(self.calc_diff_bessel, amps, args=(measurement,), method=method, options={'disp': False})#, constraints=cons2)
        #optimize.minimize(self.calc_diff_bessel,amps,args=(measurement,),disp=0,full_output=1)
        #q, fopt, func_calls, iterations, warnflag = optimize.fmin(self.calc_diff_bessel,amps,args=(measurement,),disp=0,full_output=1)
        if res.success!=True:
            print 'Possible errors optimising, flag :'
        print 'time :', time.time() - start_time
        predictions = [self.calculate_measurements_bessel(res.x,channel) for channel in range(len(self.s_list))]
        profile = self.calculate_profile_bessel(res.x, np.linspace(0,1,50))
        return res.x, res.fun, measurement, predictions, profile
    
    def fit_bessel_shot(self, dens_data, time_data, start_time, end_time, increment, plot_results = 0,method='Nelder-Mead', force_positive = 0):
        '''Perform the bessel fits over a range of values
        SH: 9Apr2013
        '''
        self.dens_data = dens_data
        self.time_data = time_data
        start_value = np.argmin(np.abs(time_data-start_time))
        end_value = np.argmin(np.abs(time_data-end_time))
        self.bessel_fit_errors = []; self.bessel_fit_predictions = []; self.bessel_fit_measurements = []
        self.bessel_fit_profiles = []; self.bessel_fit_times = []
        #create all zeros as a starting point 
        amps = [0 for i in range(0,len(self.m_l[0])*len(self.m_l[1]))]
        for i in range(start_value, end_value-increment,increment):
            q, fopt, measurement, predictions, profile = self.fit_bessel2(i,increment,amps, method=method, force_positive = force_positive)
            #q, fopt, measurement, predictions, profile = self.fit_bessel2(i,increment,amps)
            self.bessel_fit_errors.append(fopt)
            self.bessel_fit_measurements.append(measurement)
            self.bessel_fit_predictions.append(predictions)
            self.bessel_fit_profiles.append(profile)
            time_average = np.average(self.time_data[i:i+increment])
            self.bessel_fit_times.append([time_average for tmp_tmp in predictions])
            #use previous value to help covergence
            amps = q
        self.bessel_fit_errors=np.array(self.bessel_fit_errors)
        self.bessel_fit_measurements=np.array(self.bessel_fit_measurements)
        self.bessel_fit_predictions=np.array(self.bessel_fit_predictions)
        self.bessel_fit_profiles=np.array(self.bessel_fit_profiles)
        self.bessel_fit_times=np.array(self.bessel_fit_times)
        self.bessel_fit_profile_s_vals = np.linspace(0,1,50)
        if plot_results:
            colors = ['k','r','b','y','m','k','r','b','y','m']
            colors.extend(colors)
            colors.extend(colors)
            fig_new, ax_new = pt.subplots(nrows = 3)
            for i in range(0,self.bessel_fit_measurements.shape[1]):
                ax_new[0].plot(self.bessel_fit_times[:,i],self.bessel_fit_measurements[:,i],color=colors[i],marker='x',label=str(i))
                ax_new[0].plot(self.bessel_fit_times[:,i],self.bessel_fit_predictions[:,i],color=colors[i],marker='o')
            for i in range(self.bessel_fit_profiles.shape[0]):
                ax_new[1].plot(self.bessel_fit_profile_s_vals, self.bessel_fit_profiles[i,:],'-o')
            ax_new[0].legend(loc='best')
            im = ax_new[2].imshow(np.transpose(self.bessel_fit_profiles),extent=[np.min(self.bessel_fit_times),np.max(self.bessel_fit_times),0,1],aspect='auto',origin='lower',cmap='jet')
            clim = [0,4.5]
            im.set_clim(clim)
            fig_new.suptitle('bessel fit: m_values :%s, l_values:%s, force_positive:%d'%(self.m_l[0].__str__(), self.m_l[1].__str__(), force_positive))
            fig_new.canvas.draw(); fig_new.show()


    def calculate_measurements_poly(self, poly, channel):
        '''Calculate the interferometer output based on 
        a polynomial representation of psi
        Need the svalues for the line of sight of the interferometer and the dl between them
        SH : 8Apr2013
        '''
        s_vals = self.s_list_poly[channel]
        dl = self.dl_list_poly[channel]
        return np.sum(np.polyval(poly, s_vals)*dl)/(dl*len(s_vals))

    def calc_diff_poly(self, test_poly, measurement):
        '''Calculate an error between the actual measurement and the
        estimated value from calculate_measurement
        measurement are the experimental values
        s_list is a list of the psi values for line of sight of each channel
        dl_list is a list of dl values for each channel
        test_poly is the current best fit polynomial to psi
        SH: 8Apr2013
        '''
        predicted = []
        for channel in range(0,len(self.s_list_poly)):
            predicted.append(self.calculate_measurements_poly(test_poly, channel))
        predicted = np.array(predicted)
        return np.sum((measurement - np.array(predicted))**2)

    def calculate_profile_poly(self, amps, s_vals):
        '''Calculate the interferometer output based on 
        a bessel representation of psi
        Need the svalues for the line of sight of the interferometer channels and the dl between them
        SH : 8Apr2013
        '''
        return np.polyval(amps,np.linspace(0,1,50))

    def calculate_profile_poly_fit_positive(self, amps):#, measurement):
        '''Calculate the interferometer output based on 
        a bessel representation of psi
        Need the svalues for the line of sight of the interferometer channels and the dl between them
        SH : 8Apr2013
        '''
        return np.min(np.polyval(amps,np.linspace(0,1,50)))


    def fit_poly(self, index, increment,amps, method='Nelder-Mead', force_positive = 0):
        '''Find the bessel amplitudes that fit best to the averaged data from index -> index+increment
        in the interferometer data
        SH: 9Apr2013
        '''
        start_time = time.time()
        measurement = np.zeros(len(self.s_list_poly),dtype=float) #make this one larger than needed for the 0 measurement
        measurement[:self.dens_data.shape[0]] = np.average(self.dens_data[:,index:index+increment],axis = 1)
        if force_positive:
            print 'force positive'
            constraints = ({'type': 'ineq',
                            'fun': self.calculate_profile_poly_fit_positive})
            res = optimize.minimize(self.calc_diff_poly, amps, args=(measurement,), method=method, options={'disp': False}, constraints=constraints)
        else:
            print 'not force positive'
            res = optimize.minimize(self.calc_diff_poly, amps, args=(measurement,), method=method, options={'disp': False})
        #res = optimize.minimize(self.calc_diff_poly, amps, args=(measurement,), method=method, options={'disp': False}, constraints=constraints)
        if res.success!=True:
            print 'Possible errors optimising, flag :'
        print 'time :', time.time() - start_time
        predictions = [self.calculate_measurements_poly(res.x,channel) for channel in range(len(self.s_list_poly))]
        profile = self.calculate_profile_poly(res.x, np.linspace(0,1,50))
        return res.x, res.fun, measurement, predictions, profile

        # if warnflag!=0:
        #     print 'Possible errors optimising, flag :'
        # print 'time :', time.time() - start_time
        # predictions = [self.calculate_measurements_poly(q,channel) for channel in range(len(self.s_list_poly))]
        # profile = self.calculate_profile_poly(q, np.linspace(0,1,50))
        # return q, fopt, measurement, predictions, profile

    def fit_poly_shot(self, dens_data, time_data, start_time, end_time, increment, poly_order, include_dummy_channel=1, plot_results = 0, method='Nelder-Mead', force_positive = 0):
        '''Perform the bessel fits over a range of values
        SH: 9Apr2013
        '''
        self.dens_data = dens_data; self.time_data = time_data
        #include dummy channels to set the density at the edge to be zero
        if include_dummy_channel:
            self.s_list_poly = copy.deepcopy(self.s_list)
            self.dl_list_poly = copy.deepcopy(self.dl_list)
            self.dl_list_poly.append(self.dl_list[-1])
            self.s_list_poly.append(self.s_list[-1]*0+1)
        else:
            self.s_list_poly = self.s_list
            self.dl_list_poly = self.dl_list
        start_value = np.argmin(np.abs(time_data-start_time))
        end_value = np.argmin(np.abs(time_data-end_time))
        self.poly_fit_errors = []; self.poly_fit_predictions = []; self.poly_fit_measurements = []
        self.poly_fit_profiles = []; self.poly_fit_times = []
        #create all zeros as a starting point 
        amps = [0 for i in range(poly_order+1)]
        for i in range(start_value, end_value-increment,increment):
            q, fopt, measurement, predictions, profile = self.fit_poly(i,increment,amps,method=method,force_positive = force_positive)
            self.poly_fit_errors.append(fopt)
            self.poly_fit_measurements.append(measurement)
            self.poly_fit_predictions.append(predictions)
            self.poly_fit_profiles.append(profile)
            time_average = np.average(self.time_data[i:i+increment])
            self.poly_fit_times.append([time_average for tmp_tmp in predictions])
            #use previous value to help covergence
            amps = q
        self.poly_fit_errors=np.array(self.poly_fit_errors)
        self.poly_fit_measurements=np.array(self.poly_fit_measurements)
        self.poly_fit_predictions=np.array(self.poly_fit_predictions)
        self.poly_fit_profiles=np.array(self.poly_fit_profiles)
        self.poly_fit_times=np.array(self.poly_fit_times)
        self.poly_fit_profile_s_vals = np.linspace(0,1,50)
        if plot_results:
            colors = ['k','r','b','y','m','k','r','b','y','m']
            colors.extend(colors)
            colors.extend(colors)
            print 'start poly way'
            fig_new, ax_new = pt.subplots(nrows = 3)
            for i in range(0,self.poly_fit_measurements.shape[1]):
                ax_new[0].plot(self.poly_fit_times[:,i],self.poly_fit_measurements[:,i],color=colors[i],marker='x')
                ax_new[0].plot(self.poly_fit_times[:,i],self.poly_fit_predictions[:,i],color=colors[i],marker='o')
            for i in range(self.poly_fit_profiles.shape[0]):
                ax_new[1].plot(self.poly_fit_profile_s_vals, self.poly_fit_profiles[i,:],'-o')
            im = ax_new[2].imshow(np.transpose(self.poly_fit_profiles),extent=[np.min(self.poly_fit_times),np.max(self.poly_fit_times),0,1],aspect='auto',origin='lower',cmap='jet')
            clim = [0,4.5]
            im.set_clim(clim)
            fig_new.suptitle('poly fit order: %d, extra zero edge channel : %d, force positive : %d'%(poly_order, include_dummy_channel, force_positive))
            fig_new.canvas.draw(); fig_new.show()

    def calculate_inverse_matrix(self,n_segments):
        '''Calculate the inverse matrix for the calculation on density profile
        based on a fixed number of segments
        SH: 9Apr2013
        '''
        print 'calculating inverse matrix'
        self.segments = np.linspace(0,1,n_segments + 1)
        n_measurements = len(self.z_inter)
        geom_mat = np.zeros((n_measurements,n_segments))
        for i in range(0,len(self.z_inter)):
            s_tmp = self.s_list[i]
            dl_tmp = self.dl_list[i]
            cumul_points = 0
            for j in range(0,n_segments):
                points = np.sum((s_tmp>=self.segments[j]) * (s_tmp<self.segments[j+1]))
                cumul_points += points
                geom_mat[i,j]= points * dl_tmp/(dl_tmp*len(s_tmp))
        self.geom_mat = geom_mat
        self.geom_mat_pinv = np.linalg.pinv(geom_mat)


    def calculate_inverse_matrix2(self,n_segments):
        '''Calculate the inverse matrix for the calculation on density profile
        based on a fixed number of segments
        SH: 9Apr2013
        '''
        print 'calculating inverse matrix'
        self.segments = np.linspace(0,1,n_segments + 1)
        n_measurements = len(self.z_inter)
        geom_mat = np.zeros((n_measurements,n_segments))
        for i in range(0,len(self.z_inter)):
            s_tmp = self.LOS_coors[i][:,0]
            #s_tmp = self.s_list[i]
            dl_tmp = self.dl_list[i]
            cumul_points = 0
            for j in range(0,n_segments):
                points = np.sum((s_tmp>=self.segments[j]) * (s_tmp<self.segments[j+1]))
                cumul_points += points
                geom_mat[i,j]= points * dl_tmp/(dl_tmp*len(s_tmp))
        self.geom_mat = geom_mat
        self.geom_mat_pinv = np.linalg.pinv(geom_mat)


    def calculate_segment_fits(self,dens_data, time_data, start_time, end_time, increment, plot_results = 0):
        '''Run through a shot and find the density profile using the inverse
        matrix with a fixed number of segments method
        SH: 9Apr2013
        '''
        self.dens_data = dens_data; self.time_data = time_data
        self.segment_fit_measurements=[]
        self.segment_fit_predictions=[]
        self.segment_fit_profiles=[]
        self.segment_fit_times=[]
        start_value = np.argmin(np.abs(time_data-start_time))
        end_value = np.argmin(np.abs(time_data-end_time))
        for i in range(start_value, end_value-increment,increment):
            measurement = np.average(self.dens_data[:,i:i+increment],axis = 1)
            profile = np.dot(self.geom_mat_pinv, measurement)
            self.segment_fit_profiles.append(profile)
            self.segment_fit_predictions.append(np.dot(self.geom_mat,profile))
            self.segment_fit_measurements.append(measurement)
            time_average = np.average(self.time_data[i:i+increment])
            self.segment_fit_times.append([time_average for tmp_tmp in self.segment_fit_predictions[-1]])
        self.segment_fit_measurements=np.array(self.segment_fit_measurements)
        self.segment_fit_predictions=np.array(self.segment_fit_predictions)
        self.segment_fit_profiles=np.array(self.segment_fit_profiles)
        self.segment_fit_times=np.array(self.segment_fit_times)
        self.segment_fit_profile_s_vals = self.segments[:-1]
        if plot_results:
            colors = ['k','r','b','y','m','k','r','b','y','m']
            colors.extend(colors)
            colors.extend(colors)
            fig_new, ax_new = pt.subplots(nrows = 3)
            for i in range(0,self.segment_fit_measurements.shape[1]):
                ax_new[0].plot(self.segment_fit_times[:,i],self.segment_fit_measurements[:,i],color=colors[i],marker='o')
                ax_new[0].plot(self.segment_fit_times[:,i],self.segment_fit_predictions[:,i],color=colors[i],marker='x')
            for i in range(self.segment_fit_profiles.shape[0]):
                ax_new[1].plot(self.segment_fit_profile_s_vals, self.segment_fit_profiles[i,:],'-o')
            im = ax_new[2].imshow(np.transpose(self.segment_fit_profiles),extent=[np.min(self.segment_fit_times),np.max(self.segment_fit_times),0,1],aspect='auto',origin='lower',cmap='jet')
            clim = [0,4.5]
            im.set_clim(clim)
            fig_new.suptitle('n_segments: %d'%(len(self.segment_fit_profile_s_vals)))
            fig_new.canvas.draw(); fig_new.show()

    def generate_dummy_data(self, start_time, end_time, density_profile_poly, frequency=1000000., plot_data = 1):
        '''Generates dummy array output based on a polynomial of type density_profile_poly
        useful for checking everything is working okay
        SH: 9Apr2013
        '''
        #time data at 1MHz
        time_data = np.arange(start_time,end_time,1./frequency)
        predictions = [self.calculate_measurements_poly(density_profile_poly, channel) for channel in range(len(self.s_list))]
        dens_data = np.zeros((len(predictions),len(time_data)))
        for i in range(len(time_data)):
            dens_data[:,i] = predictions
        print dens_data.shape, time_data.shape
        if plot_data:
            fig, ax = pt.subplots()
            for i in range(dens_data.shape[0]):
                ax.plot(dens_data[i,:])
            fig.canvas.draw(); fig.show()
        return dens_data, time_data


def test_profile(args,s_vals):
    #[a1,mu1,var1,a2,mu2,var2,phase] = args
    [a1,mu1,var1,phase] = args
    #return a1*dist.norm.pdf(s_vals,loc=mu1,scale=var1)+a2*dist.norm.pdf(s_vals,loc=mu2,scale=var2)
    return a1*dist.norm.pdf(s_vals,loc=mu1,scale=var1)


def get_signals(args, coords, n, m):
    #[a1,mu1,var1,a2,mu2,var2,phase] = args
    [a1,mu1,var1,phase] = args
    predicted = []
    s_linspace = np.linspace(0,1,1000)
    predicted_amps = test_profile(args,s_linspace)

    signals = []
    for s_vals,theta_valid,phi_valid in coords:
        s_amps = np.interp(s_vals,s_linspace,predicted_amps)
        signal = np.sum(s_amps * np.exp(1j* n * phi_valid) *np.exp(1j* m* theta_valid)*np.exp(1j* phase))/(len(s_vals))
        signals.append(signal)
    return np.array(signals)
    #return np.sum(np.array(measurements) - np.array(predicted))**2

def fit_error(args,measurements,coords,n,m):
    return np.sum(np.abs((measurements - get_signals(args,coords,n,m))**2))

def fit_mode(self,measurements,n,m):
    print measurements
    print np.abs(measurements)
    print np.angle(measurements)
    print 'hello!!!!'
    coord_list = []
    for i in range(len(self.LOS_coords)):
        print 'new channel %d of %d'%(i,len(self.LOS_coords))
        curr_coords = self.LOS_coords[i]
        valid_points = np.invert(np.isnan(curr_coords[:,0])) * np.invert(np.isnan(curr_coords[:,1])) * np.invert(np.isnan(curr_coords[:,2]))
        s_valid = curr_coords[valid_points,0]
        theta_valid = curr_coords[valid_points,1]
        phi_valid = curr_coords[valid_points,2]
        coord_list.append([s_valid, theta_valid, phi_valid])

    fig,ax = pt.subplots(nrows = 3)
    # correct = [1,0.3,0.1,1.5,0.6,0.2,0.1]
    # correct = [1,0.3,0.1,0.1]
    # measurements = get_signals(correct, coord_list,n,m)
    # ax[0].plot(test_profile(correct,np.linspace(0,1,100)))

    initialise = [0.5,0.5,0.5,0.5,0.5,0.5,0.5]
    initialise = [15,0.7,0.3,0.]
    tmp = optimize.fmin(fit_error,initialise,args=(measurements, coord_list,-4,3),disp=0,full_output=1)
    print tmp
    ax[0].plot(test_profile(tmp[0],np.linspace(0,1,1000)))
    ax[1].plot(np.abs(measurements))
    ax[1].plot(np.abs(get_signals(tmp[0], coord_list,n,m)))
    ax[2].plot(np.angle(measurements))
    ax[2].plot(np.angle(get_signals(tmp[0], coord_list,n,m)))

    fig.canvas.draw();fig.show()

def intersection_gradient(focal_point, points_r, points_z):
    points1r = focal_point[0]
    points1z = focal_point[1]
    points2r = points_r
    points2z = points_z
    gradient_view2 = (points2z - points1z)/(points2r - points1r)
    intersect_view2 = points1z - gradient_view2*points1r
    return gradient_view2, intersect_view2

def imax_geometry(theta_values=None, make_single_plot = 0, make_multiple_plots = 0, plot_props = None, heliac_object=None):
    if theta_values == None:theta_values = [0]
    gradient = []; intersect = []
    #focal point is 75.5mm past tank weld
    #Lens 1 to pivot point 75.5mm + 802mm - 79mm = 798.5mm Rc
    #915mm from PFC centre to outside of tank weld
    #radius of pivot point = 1 + 0.915-(0.723) = 1.192 = a
    #centre line is 50mm above PFC center line z_offset
    a = 1.192
    Rc = 0.7985 #0.833
    z_offset = 0.05
    fc = 17./1000.
    CCD = 0.01575
    n_c = 512
    if make_multiple_plots: 
        fig2, ax2 = pt.subplots(ncols = len(theta_values), sharex=1, sharey=1)
        if len(theta_values)== 1: ax2 = [ax2]
    if make_single_plot: fig_single, ax_single = pt.subplots()
    if plot_props == None : plot_props = {'marker':'.','markersize':5,'linestyle':'None'}
    for i, theta in enumerate(theta_values):
        print theta
        theta = theta/180.*np.pi
        #focal_point is really the location of the lens1 (centre of focus from the CCD)
        #vertical motion
        focal_point = [a+(Rc-fc),(Rc-fc)*np.tan(theta)+z_offset]
        #This ordering goes from bottom of CCD to top
        #which means, top of plasma to bottom
        points_r = np.linspace(a+ Rc + CCD/2.*np.sin(theta), a + Rc - CCD/2.*np.sin(theta), n_c)
        points_z = np.linspace(Rc*np.tan(theta) - CCD/2.*np.cos(theta) + z_offset,  Rc*np.tan(theta) + CCD/2.*np.cos(theta) + z_offset, n_c)
        #radial motion?
        #points_r = np.linspace(a+ Rc*np.cos(theta) + CCD/2.*np.sin(theta), a + Rc*np.cos(theta) - CCD/2.*np.sin(theta), n_c)
        #points_z = np.linspace(Rc*np.sin(theta) - CCD/2.*np.cos(theta) + z_offset,  Rc*np.sin(theta) + CCD/2.*np.cos(theta) + z_offset, n_c)
        gradient_tmp, intersect_tmp = intersection_gradient(focal_point, points_r, points_z)
        gradient.extend(gradient_tmp)
        intersect.extend(intersect_tmp)
        if make_multiple_plots:
            ax2[i].plot(points_r, points_z,'.')
            ax2[i].plot([focal_point[0],a], [focal_point[1],0],'o-')
            if heliac_object!=None: heliac_object.puncture_plot(ax2[i],0, plot_dict = plot_props)

            ax2[i].set_xlim([0.95,Rc+a+0.15])
            ax2[i].set_ylim([-0.5,0.5])
        if make_single_plot:
            heliac_object.puncture_plot(ax_single,0, plot_dict = plot_props)
            ax_single.plot(points_r, points_z,'o')
            ax_single.plot([focal_point[0],a], [focal_point[1],0 + z_offset],'o-')
            ax_single.plot(0, 0,'o')
            if heliac_object!=None: heliac_object.puncture_plot(ax_single,0, plot_dict = plot_props)
        if make_single_plot or make_multiple_plots:
            r_interp = np.linspace(1.0,2.5)
            for a_lin, b_lin in zip(gradient_tmp, intersect_tmp):
                z_interp = a_lin*r_interp + b_lin
                if make_multiple_plots:ax2[i].plot(r_interp,z_interp, 'k',linewidth=0.1)
                if make_single_plot:ax_single.plot(r_interp,z_interp,'k',linewidth=0.1)
    if make_single_plot:
        ax_single.set_xlim([0.95,Rc+a+0.15])
        ax_single.set_ylim([-0.5,0.5])
        ax_single.vlines(0.915 + 1,ax_single.get_ylim()[0], ax_single.get_ylim()[1])
        fig_single.canvas.draw(); fig_single.show()
    if make_multiple_plots:
        fig2.canvas.draw(); fig2.show()
    return gradient, intersect


def return_correct_mode(mode_tuples, tomo_modes_n, tomo_modes_m, orientations, tomo_orient, overall_geom_list, valid_channels, start_indices, end_indices, fourier_data, harmonic):
    tomo_modes_n_indices = [mode_tuples.index([i,j]) for i,j in zip(tomo_modes_n,tomo_modes_m)]
    tomo_modes_m_indices = [mode_tuples.index([i,j]) for i,j in zip(tomo_modes_n,tomo_modes_m)]
    if tomo_modes_n_indices!=tomo_modes_m_indices: raise(Exception)
    tomo_view_indices = [orientations.index(i) for i in tomo_orient]
    rev_indices = tomo_view_indices
    z = np.hstack((np.vstack((overall_geom_list[j][i] for i in tomo_view_indices)) for j in tomo_modes_n_indices))
    tomo_valid_channels = np.vstack((valid_channels[start_indices[i]:end_indices[i]] for i in tomo_view_indices))
    tomo_measurements = np.vstack((fourier_data[harmonic, start_indices[i]:end_indices[i]] for i in tomo_view_indices))
    return z, tomo_valid_channels, tomo_measurements

#def single(fourier_data, valid_channels, overall_geom_list, mode_tuples, tomo_modes_n, tomo_modes_m, tomo_orient, lamda, orientations, start_indices, end_indices, method, cycles, count, total, return_tomo_inv, harmonic, tomo_DC, tomo_eval, LOS_obj):
def single(tomo_modes_n, tomo_modes_m, tomo_orient, lamda, method, cycles, count, total, return_tomo_inv, harmonic, tomo_DC, tomo_eval, LOS_obj):
    '''Run a single tomographic inversion using the combination of modes in tomo_modes_n, tomo_modes_m and orientatino sin tomo_orient

    SRH : 12Feb2014
    '''
    z, tomo_valid_channels, tomo_measurements = LOS_obj.return_combined_matrices(tomo_modes_n, tomo_modes_m, harmonic, tomo_orient)
    if method=='SIRT':
        tomo_sirt = SIRT(z, tomo_measurements, lamda, valid_channels = tomo_valid_channels, tomo_DC = tomo_DC)
        tomo_sirt.run(cycles = cycles)
        cur_tomo = tomo_sirt
    elif method=='Direct':
        tomo_direct = DirectSolution(z, tomo_measurements, valid_channels = tomo_valid_channels, tomo_DC = tomo_DC)
        tomo_direct.run()
        cur_tomo = tomo_direct
    elif method=='DirectFixedPhase':
        print 'using fixed phase'
        tomo_direct = DirectSolutionFixedPhase(z, tomo_measurements, valid_channels = tomo_valid_channels, tomo_DC = tomo_DC)
        tomo_direct.run()
        cur_tomo = tomo_direct
    if tomo_eval == tomo_orient:
        error = tomo_recon_error_calc(z, cur_tomo.T, tomo_measurements[tomo_valid_channels])
        error = tomo_recon_error_NRMSE(z, cur_tomo.T, tomo_measurements[tomo_valid_channels])
        #error = calc_RMSE(cur_tomo)
    else:
        z_tmp, tomo_valid_channels_tmp, tomo_measurements_tmp = LOS_obj.return_combined_matrices(tomo_modes_n, tomo_modes_m, harmonic, tomo_eval)
        error = tomo_recon_error_calc(z_tmp, cur_tomo.T, tomo_measurements_tmp[tomo_valid_channels_tmp])
        error = tomo_recon_error_NRMSE(z_tmp, cur_tomo.T, tomo_measurements_tmp[tomo_valid_channels_tmp])
        #error = calc_RMSE(cur_tomo)
    print 'n {}, m{}, error {:.4f} {} of {}'.format(tomo_modes_n, tomo_modes_m, error[0], count, total)
    if return_tomo_inv:
        return tomo_modes_n,tomo_modes_m,error[0],cur_tomo
    else:
        return tomo_modes_n,tomo_modes_m,error[0]


def _single_multiproc_wrapper(arguments):
    return single(*arguments)

class try_several_modes_multi():
    def __init__(self,n_list, m_list, tomo_orient,number_modes, answer, overall_geom_list, fourier_data,mode_tuples, lamda=0.5,cycles=300, pool_size = 1, method='SIRT', tomo_DC = None, fixed_mode = None, tomo_eval = None):
        #self.answer = answer
        if tomo_eval == None: tomo_eval = tomo_orient
        #self.overall_geom_list = overall_geom_list
        #self.fourier_data = fourier_data
        #self.mode_tuples = mode_tuples
        #self.valid_channels = answer.valid_channels

        n_m = [[n_cur,m_cur] for n_cur,m_cur in zip(n_list, m_list)]
        if fixed_mode == None:
            mode_comb_list = [i for i in itertools.combinations(n_m,number_modes)]
            if number_modes>1:
                for i in n_m: 
                    tmp = ()
                    for j in range(number_modes): tmp = tmp + (i,)
                    mode_comb_list.append(tmp)
        else:
            mode_comb_list = [i for i in itertools.combinations(n_m,1)]
            for i in range(len(mode_comb_list)):mode_comb_list[i]+=(fixed_mode,)
        print mode_comb_list

        tomo_modes_n_list = []
        tomo_modes_m_list = []

        tomo_orient_list = [tomo_orient for i in mode_comb_list]
        for i in mode_comb_list:
            tomo_modes_n_list.append([j[0] for j in i])
            tomo_modes_m_list.append([j[1] for j in i])
            #tomo_orient_list.append(tomo_orient)
        orientations = answer.orientations
        arglist = [[tomo_n, tomo_m, tomo_orient, lamda, method, cycles, counter, len(tomo_orient_list),False,1, tomo_DC, tomo_eval, answer] for tomo_n, tomo_m, tomo_orient,counter in zip(tomo_modes_n_list, tomo_modes_m_list, tomo_orient_list, range(len(tomo_orient_list)))]
        if pool_size > 1:
            print "Using multiproc"
            pool = multiprocessing.Pool(processes=pool_size)
            #errors = pool.map(self._single_multiproc_wrapper, itertools.izip(tomo_modes_n_list, tomo_modes_m_list, tomo_orient_list))
            errors = pool.map(_single_multiproc_wrapper, arglist)
            pool.close(); pool.join() # no more tasks
        else:
            #errors = map(self._single_multiproc_wrapper, itertools.izip(tomo_modes_n_list, tomo_modes_m_list, tomo_orient_list))
            errors = map(_single_multiproc_wrapper, arglist)
            #errors = self._single_multiproc_wrapper(tomo_modes_n_list[0], tomo_modes_m_list[0], tomo_orient_list[0])
        print errors
        print '  closing pool and waiting for pool to finish'
        print '  pool finished'
        self.errors = errors

def plot_error_multi_list2(input_errors, filename = None):
    fig, ax = pt.subplots()
    mode_list = []; errors = []
    if filename!=None:
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
        fig.set_figheight(8.48*0.8*cm_to_inch)
    for i in input_errors:
        mode_list.append([i[0][0],i[1][0]])
        mode_list.append([i[0][1],i[1][1]])
        errors.append(i[2])
    print mode_list
    mode_unique_list = []
    for i in mode_list:
        if i not in mode_unique_list:
            mode_unique_list.append(i)
    mode_unique_list.sort(key=lambda pair: pair[1])
    error_array = np.zeros((len(mode_unique_list), len(mode_unique_list)),dtype = float)
    errors = [i[2] for i in input_errors]
    print 'Best mode combination : {}'.format(input_errors[np.argmin(errors)])
    for i in input_errors:
        m1 = [i[0][0],i[1][0]]
        m2 = [i[0][1],i[1][1]]
        i1 = mode_unique_list.index(m1)
        i2 = mode_unique_list.index(m2)
        error_array[i1,i2] = i[2]
        error_array[i2,i1] = i[2]
    im = ax.imshow(error_array, interpolation = 'nearest',aspect='auto', cmap='hot')
    im.set_clim([np.min(error_array),np.min(error_array)*4])
    ax.set_xticks(range(len(mode_unique_list)))
    ax.set_yticks(range(len(mode_unique_list)))
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    def mjrFormatter(x, pos):
        x = int(x)
        return '{},{}'.format(mode_unique_list[x][0], mode_unique_list[x][1])
    ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(mjrFormatter))
    ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(mjrFormatter))
    cbar = pt.colorbar(im, ax = ax)
    cbar.set_label('Tomo Recon Error')
    for tick in ax.xaxis.iter_ticks():
        #tick[0].label2On = True
        tick[0].label1On = True
        tick[0].label1.set_rotation('vertical')
    #for tick in ax.yaxis.iter_ticks():
    #    tick[0].label2On = True
    #    tick[0].label1On = False
    ax.set_xlabel('Mode 1 (n,m)')
    ax.set_ylabel('Mode 2 (n,m)')
    #ax.set_title('Reconstruction error using 2 modes')
    if filename!=None:
        fig.tight_layout()
        fig.savefig(filename+'.pdf')
        fig.savefig(filename+'.eps')
        print filename
    fig.canvas.draw(); fig.show()
    return input_errors[np.argmin(errors)]



def plot_error_multi_list_single(input_errors, filename = None, single_m = None,clim = None,second_errors = None):
    print 'hello'
    if single_m!=None:
        fig, [ax,ax2] = pt.subplots(nrows = 2)
    else:
        fig, ax = pt.subplots()
    if filename!=None:
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
        fig.set_figheight(8.48*1.2*cm_to_inch)

    mode_list_n = []; mode_list_m = []; errors = []
    for i in input_errors:
        mode_list_n.append(i[0][0])
        mode_list_m.append(i[1][0])
        errors.append(i[2])
    mode_list_n = list(set(mode_list_n))
    mode_list_m = list(set(mode_list_m))
    mode_list_n.sort()
    mode_list_m.sort()
    error_array = np.zeros((len(mode_list_n), len(mode_list_m)),dtype = float)
    for i in input_errors:
        m1 = i[0][0]
        m2 = i[1][0]
        i1 = mode_list_n.index(m1)
        i2 = mode_list_m.index(m2)
        error_array[i1,i2] = i[2]
    if second_errors!=None: 
        second_error_array = np.zeros((len(mode_list_n), len(mode_list_m)),dtype = float)
        for i in second_errors:
            m1 = i[0][0]
            m2 = i[1][0]
            i1 = mode_list_n.index(m1)
            i2 = mode_list_m.index(m2)
            second_error_array[i1,i2] = i[2]
    im = ax.imshow(error_array, interpolation = 'nearest',aspect='auto',cmap='hot')
    if clim==None:clim = [np.min(error_array),np.max(error_array)*0.85]
    im.set_clim(clim)
    ax.set_xticks(range(len(mode_list_m)))
    ax.set_yticks(range(len(mode_list_n)))
    print mode_list_n
    print mode_list_m
    def mjrFormatter(x, pos):
        x = int(x)
        return '{}'.format(mode_list_n[x])
    def mjrFormatter2(x, pos):
        x = int(x)
        return '{}'.format(mode_list_m[x])
    ax.set_xlabel('m')
    x_min, x_max = ax.get_xlim(); y_min, y_max = ax.get_ylim()
    ax.text(x_min+ (x_max-x_min)*0.05, y_min+ (y_max-y_min)*0.9,'(a)')
    ax.set_ylabel('n')
    #ax.set_title('Error from Single mode fit')
    ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(mjrFormatter))
    ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(mjrFormatter2))
    cbar = pt.colorbar(im, ax = ax)
    cbar.set_label('Tomo Recon Error')
    if single_m!=None:
        if single_m.__class__=='int': 
            single_m = [single_m]
            print single_m.__class__
        colors = ['b','g','r','k']
        if second_errors != None:
            for tmp_i in range(len(single_m)):
            #fig2,ax2 = pt.subplots()
                tmp_single_m = single_m[tmp_i]
                color = colors[tmp_i]
                index = mode_list_m.index(tmp_single_m)
                ax2.plot(mode_list_n, second_error_array[:,index], 'x-', color = color, label='center view m={}'.format(tmp_single_m))
            #ax2.text(mode_list_n[0], error_array[0,index],'m={}'.format(tmp_single_m))
        colors = ['b','g','r','k']
        for tmp_i in range(len(single_m)):
        #fig2,ax2 = pt.subplots()
            tmp_single_m = single_m[tmp_i]
            color = colors[tmp_i]
            index = mode_list_m.index(tmp_single_m)
            ax2.plot(mode_list_n, error_array[:,index], 's-',color = color, label='all views m={}'.format(tmp_single_m))
            min_loc = np.argmin(error_array[:,index])

        ax2.set_xlabel('n')
        ax2.set_ylim([0.25,0.85])
        x_min, x_max = ax2.get_xlim(); y_min, y_max = ax2.get_ylim()
        ax2.text(x_min+ (x_max-x_min)*0.05, y_min+ (y_max-y_min)*0.9,'(b)')
        ax2.set_ylabel('Tomo Recon Error')
        ax2.legend(loc='best')
        ax2.grid()
        #ax2.set_title('Reconstruction error for m={} and various n values'.format(single_m))
        #fig2.canvas.draw(); fig2.show()
    # for tick in ax.xaxis.iter_ticks():
    #     tick[0].label2On = True
    #     tick[0].label1On = False
    #     tick[0].label2.set_rotation('vertical')
    # for tick in ax.yaxis.iter_ticks():
    #     tick[0].label2On = True
    #     tick[0].label1On = False
    if filename!=None:
        fig.tight_layout(pad=0.1)
        fig.savefig(filename + '.pdf')
        fig.savefig(filename + '.eps')
        print filename
    fig.canvas.draw(); fig.show()



def iota(r,kappa):
    a =  [[1.24403098, 0.29927867, -0.04178176, -0.0113835, 0.01371373], 
          [-0.06438457, 0.17743677, -0.00568132, 0.11426079, -0.0981305],
          [0.16757832, -0.41083898, 0.00136293, -0.23903926, 0.22891545], 
          [-0.21602304, 0.16208048, 0.05840499, 0.1875845, -0.21617175],
          [0.12705246, -0.00544844, -0.03210589,-0.05116255, 0.07173953]]
    temp=0
    for i in range(5):
        for j in range(5):
            temp = temp + a[i][j]*(r+0.25938664)**i*(kappa-0.34786773)**j
    return temp

def cas3d_mode_table(mummin, mummax, munmin, munmax, mutol, N, kh, plot_table = False):
    mummin = 0; mummax = 7
    munmin = -8; munmax = 0
    #etatol = 5
    mutol = 4
    n_table_list = []; m_table_list = []
    sample = np.arange(0,1.01,0.01)
    iotamin = min([iota(r,kh) for r in sample])
    iotamax = max([iota(r,kh) for r in sample])
    N = 2
    m_n_table = np.zeros((munmax - munmin+1, mummax - mummin+1), dtype = int)
    for i, n_tmp in enumerate(range(munmin,munmax+1)):
        for j, m_tmp in enumerate(range(mummin,mummax+1)):
            if m_tmp*(iotamin)-mutol<-n_tmp and -n_tmp<m_tmp*(iotamax)+mutol and (m_tmp!=0 or n_tmp<=0) and (N==2 or not(N^(n_tmp%3>0))):
            #if m_tmp*(iotamin)-mutol<n_tmp and n_tmp<m_tmp*(iotamax)+mutol and (m_tmp!=0) and (N==2 or not(N^(n_tmp%3>0))):
                n_table_list.append(n_tmp)
                m_table_list.append(m_tmp)
                n_table_list.append(-n_tmp)
                m_table_list.append(m_tmp)
                m_n_table[i, j] = 1
    if plot_table:
        fig, ax = pt.subplots()
        ax.imshow(m_n_table, aspect = 'auto',extent=[mummin, mummax, munmax, munmin], origin = 'upper', interpolation='nearest')
        ax.set_xlabel('m'); ax.set_ylabel('n')
        fig.canvas.draw(); fig.show()
    return n_table_list, m_table_list, m_n_table


def compare_reconstruction_method(tomo1, tomo2, LOS_object, n, m, filename = None):
    fig, ax = pt.subplots(nrows = 2, sharex = True)
    if filename!=None:
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
        fig.set_figheight(8.48*1.0*cm_to_inch)
    method_txt = tomo1.method
    if method_txt == 'DirectSolution': method_txt = 'Direct'
    plot_radial_structure(tomo1.T, LOS_object.segment_midpoints, n, m,  prov_ax = ax, norm = True, extra_txt = method_txt + ',', single_mode = None, marker = 'x')
    method_txt = tomo2.method
    if method_txt == 'DirectSolution': method_txt = 'Direct'
    plot_radial_structure(tomo2.T, LOS_object.segment_midpoints, n, m, prov_ax = ax, norm = True, extra_txt = method_txt + ',', single_mode = None, marker = 'o')
    ax[0].legend(loc='lower center',prop={'size':6.5})
    ax[1].set_xlim([0,1])
    ax[1].set_ylim([-2.*np.pi,0])
    ax[1].set_xlabel(r'$\Psi_N$')
    ax[0].set_yticks(ax[0].get_yticks()[::2])
    for i in ax: i.grid()
    if filename!=None:
        fig.tight_layout(pad=0.1)
        fig.savefig(filename + '.pdf')
        fig.savefig(filename + '.eps')
        fig.savefig(filename + '.svg')
        print filename
    fig.canvas.draw(); fig.show()

def compare_num_of_views(tomo_list, LOS_object, n, m,marker=None,labels = None, filename = None, colors = None):
    if marker == None: marker = ['x','d','s','o']
    if colors == None: colors = [None]*len(marker)
    if labels == None: labels = ['']*len(tomo_list)
    fig, ax = pt.subplots(nrows = 2, sharex = True)
    if filename!=None:
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
        fig.set_figheight(8.48*1.0*cm_to_inch)
    for i, tmp_marker, tmp_label, color in zip(tomo_list,marker, labels, colors):
        if LOS_object.radial_s_spacing:
            plot_radial_structure(i.T, np.sqrt(LOS_object.segment_midpoints), n, m,  prov_ax = ax, norm = True, extra_txt = tmp_label + ',', single_mode = None, marker = tmp_marker, color = color)
        else:
            plot_radial_structure(i.T, LOS_object.segment_midpoints, n, m,  prov_ax = ax, norm = True, extra_txt = tmp_label + ',', single_mode = None, marker = tmp_marker, color=color)
    ax[0].legend(loc='upper left',prop={'size':7})
    ax[1].set_xlim([0,1])
    ax[1].set_ylim([-2.*np.pi,0])
    if LOS_object.radial_s_spacing:
        ax[1].set_xlabel(r'$\sqrt{s}$')
    else:
        ax[1].set_xlabel(r's')
    ax[0].set_yticks(ax[0].get_yticks()[::2])
    for i in ax: i.grid()
    if filename!=None:
        fig.tight_layout(pad=0.1)
        fig.savefig(filename + '.pdf')
        fig.savefig(filename + '.eps')
        fig.savefig(filename + '.svg')
        print filename
    fig.canvas.draw(); fig.show()

def plot_radial_structure(T, segment_midpoints, n, m, prov_ax = None, norm = False, extra_txt = '', single_mode =None,marker ='x',color=None,linestyle='solid', max_phase = np.pi):

    plot_dict = {'linestyle':linestyle, 'marker':marker}
    if color!=None: plot_dict['color'] = color
    start = 0; increment = len(T)/len(n)
    if prov_ax==None:
        fig, ax = pt.subplots(nrows = 2)
    else:
        ax = prov_ax
    #Find biggest mode and norm factor    
    amp_list = []
    if norm:
        for n_cur, m_cur in zip(n,m):
            cor_run = True
            end = start + increment
            cur_T = T[start:end]
            start = +end
            amp_list.append(np.sum(np.abs(cur_T)))
        norm_fact = np.max(amp_list)
    start = 0
    for n_cur, m_cur in zip(n,m):
        cor_run = True
        if single_mode!=None:
            cor_run = single_mode == [n_cur, m_cur]
        end = start + increment
        cur_T = T[start:end]
        start = +end
        if cor_run:
            if norm:
                cur_T = cur_T/norm_fact
            ax[0].plot(segment_midpoints, np.abs(cur_T), label='{}({},{})'.format(extra_txt, n_cur, m_cur), **plot_dict)
            max_ind = np.argmax(np.abs(cur_T[:-5]))
            max_s = segment_midpoints[max_ind]
            max_amp = np.max(np.abs(cur_T[:-5]))
            #ax[0].text(max_s, max_amp, extra_txt)
            angs = np.unwrap(np.angle(cur_T))
            start_pt = angs.shape[0]/3
            while np.mean(angs[start_pt:])>max_phase: angs += -2.*np.pi
            while np.mean(angs[start_pt:])<(max_phase - 2.*np.pi): angs += 2.*np.pi
            print np.sum(cur_T), np.abs(np.sum(cur_T)), np.angle(np.sum(cur_T))
            ax[1].plot(segment_midpoints, angs, label='Arg({},{})'.format(n_cur, m_cur), **plot_dict)
            #ax[1].text(max_s, angs[max_ind], extra_txt)
            ax[0].legend(loc='best', prop={'size':6})
            ax[0].set_xlim([0,1])
            ax[1].set_xlim([0,1])
            ax[1].set_ylabel('Phase (rad)')
            ax[0].set_ylabel('Amplitude (a.u)')
            ax[1].set_ylim([np.mean(angs[start_pt:])-np.pi,np.mean(angs[start_pt:])+np.pi])
    if prov_ax==None:
        fig.canvas.draw();fig.show()
    return max_s, max_amp, max_ind


def run_inversion(tomo_modes_n, tomo_modes_m, tomo_orient, tomo_orient_extrap, filename, answer, harmonic = 1, method = 'Direct', cycles=300, lamda=0.5, cut_values = None, plot_wave_fields = False, plot_old_reproj_comparison = False, plot_old_combo = False, plot_reprojection_comp1 = False, plot_reprojection_comp2 = False, plot_wave_animation = False, tomo_DC = None, dI_I = False):
    '''Run the tomographic inversion for the combination of modes in tomo_modes_n, tomo_modes_m

    SRH: 12Feb2014
    '''

    tomo_inv = single(tomo_modes_n, tomo_modes_m, tomo_orient, lamda, method, cycles, 0, 1,True, harmonic, tomo_DC, tomo_orient_extrap, answer)[3]
    if plot_old_combo:
        tomo_inv.plot_lots_of_things(tomo_modes_n, tomo_modes_m, answer, cut_values=cut_values)
    if plot_old_reproj_comparison:
        tomo_inv.plot_reprojection_comparison(tomo_modes_n, tomo_modes_m, answer, cut_values=None, multiplier=1.0, pub_fig=1)
    if plot_wave_fields:
        tomo_inv.plot_wave_field(tomo_modes_n, tomo_modes_m, answer,answer.s_values)
    z_extrap, tomo_valid_channels_extrap, tomo_measurements_extrap = answer.return_combined_matrices(tomo_modes_n, tomo_modes_m, harmonic, tomo_orient_extrap)
    if tomo_DC != None:
        if z_extrap.shape[1]%tomo_DC.T.shape[0]!=0: raise(ValueError)
        n_modes = z_extrap.shape[1]/tomo_DC.T.shape[0]
        z_extrap = np.dot(z_extrap, np.diag(np.hstack((tomo_DC.T for i in range(n_modes)))))
    if plot_reprojection_comp1:
        tomo_inv.plot_reprojection_comparison_extrap(tomo_modes_n, tomo_modes_m, answer, tomo_valid_channels_extrap, tomo_measurements_extrap, z_extrap, cut_values = None, multiplier = 1., pub_fig = 1, savefig_name='reproj_'+filename)
    if plot_wave_animation:
        tomo_inv.plot_wave_field_animation(tomo_modes_n, tomo_modes_m, answer,answer.s_values, n_images = 1, inc_cbar = 1, pub_fig=True, save_fig='wave_bean_'+filename, inc_profile = 1, tomo_DC = tomo_DC)

    if plot_reprojection_comp2:
        tomo_inv.plot_reprojection_comparison_extrap_diff(tomo_modes_n, tomo_modes_m, answer, tomo_valid_channels_extrap, tomo_measurements_extrap, z_extrap, cut_values = None, multiplier = 1., pub_fig = 1, savefig_name='reproj_diff_'+filename, tomo_DC = tomo_DC, dI_I = dI_I)
    return tomo_inv

def run_inv_and_save(kh, tomo_modes_n_best, tomo_modes_m_best, tomo_orient, answer, cut_values, tomo_DC, make_animation = False, run_dI_I = True, method = 'Direct', plot_reprojection_comp1 = True, plot_reprojection_comp2 = True, dI_anim = False ):
    kh_string = '{:.2f}'.format(kh).replace('.','_')
    filename = 'kh_{}_'.format(kh_string)
    for i, j in zip(tomo_modes_n_best, tomo_modes_m_best):
        filename+='{}_{}__'.format(i,j)

    #This makes a nice plot of m vs n
    #tomo.plot_error_multi_list_single(error_multi_list.errors,filename='kh_{}_helicity_check'.format(kh_string),single_m=-3)
    print tomo_orient
    if run_dI_I:
        tomo_dI = run_inversion(tomo_modes_n_best, tomo_modes_m_best, tomo_orient, tomo_orient, filename+'_dI-I', answer, harmonic = 1, cut_values = cut_values, method = method, cycles=300, lamda=0.5, plot_wave_fields = False, plot_old_reproj_comparison = False, plot_old_combo = False, plot_reprojection_comp1 = plot_reprojection_comp1, plot_reprojection_comp2 = plot_reprojection_comp2, plot_wave_animation = False, tomo_DC = tomo_DC, dI_I = True)
    else:
        tomo_dI = None
    tomo_norm = run_inversion(tomo_modes_n_best, tomo_modes_m_best, tomo_orient, tomo_orient, filename+'_norm',  answer, harmonic = 1, cut_values = cut_values, method = method, cycles=300, lamda=0.5, plot_wave_fields = False, plot_old_reproj_comparison = False, plot_old_combo = False, plot_reprojection_comp1 = plot_reprojection_comp1, plot_reprojection_comp2 = plot_reprojection_comp2, plot_wave_animation = False, tomo_DC = None, dI_I = False)
    if make_animation:
        print 'start animation'
        if dI_anim:
            tomo_dI.plot_wave_field_animation(tomo_modes_n_best, tomo_modes_m_best, answer,answer.s_values, n_images = 24, inc_cbar = 1, pub_fig=True, save_fig='wave_bean_'+filename, inc_profile = 1, save_name = filename+'_animation_dI-I.gif',delay=10, kh = kh, tomo_DC = tomo_DC, dI_I = True)
        else:
            tomo_norm.plot_wave_field_animation(tomo_modes_n_best, tomo_modes_m_best, answer,answer.s_values, n_images = 24, inc_cbar = 1, pub_fig=True, save_fig='wave_bean_'+filename, inc_profile = 1, save_name = filename+'_animation.gif',delay=10, kh = kh, tomo_DC = None)
    mode_output_data = {'T':tomo_norm.T, 'segment_midpoints':answer.segment_midpoints,'n':tomo_modes_n_best,'m':tomo_modes_m_best}
    inv_prof_fname = filename + 'mode_norm.pickle'
    pickle.dump(mode_output_data, file(inv_prof_fname,'w'))
    if run_dI_I:
        mode_output_data = {'T':tomo_dI.T, 'segment_midpoints':answer.segment_midpoints,'n':tomo_modes_n_best,'m':tomo_modes_m_best}
        inv_prof_fname = filename + 'mode_dII.pickle'
        pickle.dump(mode_output_data, file(inv_prof_fname,'w'))
    return tomo_dI, tomo_norm


def many_measurements(LOS_object, tomo_modes_n, tomo_modes_m, tomo_orient, harmonic, peak_loc = 0.5, peak_width = 0.1, noise_strength = 0.00001, filename = ''):
    filename_list = []
    for i, peak_loc in enumerate(np.linspace(0.1,0.9,5)):
        filename = str(i)
        create_measurements(LOS_object, tomo_modes_n, tomo_modes_m, tomo_orient, harmonic, peak_loc = peak_loc, peak_width = 0.1, noise_strength = noise_strength, filename = filename)
        filename_list.append('reproj_diff_{}'.format(filename))
    print filename_list
    for i in filename_list:
        os.system('convert {}.pdf {}.png'.format(i, i))

    os.system('convert -delay {} -loop 0 {} {}'.format(100, '.png '.join(filename_list)+'.png', 'moving_gaussian_{}_noise.gif'.format('{:.2f}'.format(noise_strength).replace('.','_'))))
    #os.system('zip {} {}'.format(save_name.rstrip('gif')+'zip', image_string))



def many_measurements_noise(LOS_object, tomo_modes_n, tomo_modes_m, tomo_orient, harmonic, peak_loc = 0.5, peak_width = 0.1, noise_strength = 0.00001, filename = ''):
    filename_list = []
    fig, ax2 = pt.subplots(nrows = 2, ncols = 2, sharex = True, sharey = True)
    ax = ax2.flatten()
    cm_to_inch=0.393701
    fig.set_figwidth(8.48*cm_to_inch)
    fig.set_figheight(8.48*0.75*cm_to_inch)
    count = 0
    markers = ['.','x']
    for j, peak_loc in enumerate([0.2, 0.4, 0.6, 0.8]):
        for i, (noise_strength, cur_marker) in enumerate(zip([0.5,1.5], markers)):
            filename = str(i)
            T_inv, T = create_measurements(LOS_object, tomo_modes_n, tomo_modes_m, tomo_orient, harmonic, peak_loc = peak_loc, peak_width = 0.1, noise_strength = noise_strength, filename = filename, plot_diff_images = False)
            filename_list.append('reproj_diff_{}'.format(filename))
            ax[count].plot(np.sqrt(LOS_object.segment_midpoints), np.abs(T_inv), label='{:.2f} noise'.format(noise_strength), marker = cur_marker)
            #ax[1,j].plot(np.sqrt(LOS_object.segment_midpoints), np.angle(T_inv))
        ax[count].plot(np.sqrt(LOS_object.segment_midpoints), np.abs(T),label='Correct', marker = '')
        count += 1
        #ax[1,j].plot(np.sqrt(LOS_object.segment_midpoints), np.angle(T))
    
    ax[0].legend(loc='best')
    ax[2].set_xlabel('$\sqrt{s}$')
    ax[3].set_xlabel('$\sqrt{s}$')
    ax[0].set_ylabel('Amp (a.u)')
    ax[2].set_ylabel('Amp (a.u)')
    for i in ax:i.grid()
    for i in [ax2[1,0]]:i.set_xticks(i.get_xticks()[::2])
    for i in [ax2[0,0]]:i.set_yticks(i.get_yticks()[::2])
    fig.tight_layout(pad = 0.01)
    fig.savefig('noise_tollerance.pdf')
    fig.canvas.draw(); fig.show()


def add_noise(LOS_object, tomo_modes_n, tomo_modes_m, tomo_orient, harmonic, noise_strength = None, filename = None):
    if noise_strength == None: noise_strength = 0.0001
    filename_list = []

    fig = pt.figure()
    cm_to_inch=0.393701
    fig.set_figwidth(8.48*cm_to_inch)
    fig.set_figheight(8.48*0.9*cm_to_inch)
    gs = gridspec.GridSpec(20,4)
    ax2 = []
    ax2.append(pt.subplot(gs[:9,3]))
    ax2.append(pt.subplot(gs[9:17,3], sharex = ax2[0],sharey = ax2[0]))
    ax2.append(pt.subplot(gs[:10,:3]))
    ax2.append(pt.subplot(gs[10:,:3], sharex = ax2[-1]))
    #ax2.append(pt.subplot(gs[:8,6], sharex = ax2[-1],sharey = ax2[-1]))

    #ax2.append(pt.subplot(gs[0:2,0:2]))
    #ax2.append(pt.subplot(gs[2:4,0:2]))#, sharex = ax2[5]))
    #ax2.append(pt.subplot(gs[4:8,0:2]))

    #cbar_wave_ax = pt.subplot(gs[8,0:2])
    #cbar_phase_ax = pt.subplot(gs[8,2:4])
    cbar_amp_ax = pt.subplot(gs[17:,3])

    #fig, ax2 = pt.subplots(nrows = 2, ncols = 2)
    ax = ax2#.flatten()


    tmp_meas = np.abs(LOS_object.fourier_data[harmonic,:,:])*LOS_object.valid_channels
    tmp_meas[-LOS_object.valid_channels] = np.nan
    amp_cmap = 'jet'
    [orig_im_ax, noise_im_ax, amp_ax, phase_ax] = ax
    orig_im = orig_im_ax.imshow(np.abs(tmp_meas[LOS_object.start_indices[1]:LOS_object.end_indices[1]]),origin='upper',aspect='auto',interpolation='nearest', cmap=amp_cmap)

    orig_im_ax.set_xlim([160,82])

    orig_im_ax.set_ylim([0,  (LOS_object.end_indices[1] - LOS_object.start_indices[1])*0.9])
    foo1, foo2, foo3, orig_inv = single(tomo_modes_n, tomo_modes_m, tomo_orient, 0.5, 'Direct', 300, 0, 1, True, 1, None, tomo_orient, LOS_object)
    plot_radial_structure(orig_inv.T, np.sqrt(LOS_object.segment_midpoints), tomo_modes_n, tomo_modes_m, prov_ax = [amp_ax, phase_ax], norm = False, extra_txt = 'Orig ', single_mode = None, marker = '.')

    #sig_strength = np.sqrt(np.sum(np.abs(re_projection)**2))
    fourier_data_tmp = LOS_object.fourier_data.copy()
    valid_data = LOS_object.fourier_data[harmonic,LOS_object.valid_channels]

    print 'max amp :{}, mean amp:{}'.format(np.max(np.abs(valid_data)), np.mean(np.abs(valid_data)))
    print 'RMS real :{}, RMS imag:{}'.format(np.sqrt(np.mean((np.real(valid_data))**2)),  np.sqrt(np.mean((np.imag(valid_data))**2)))

    noise1 = np.random.normal(0,noise_strength,LOS_object.fourier_data[harmonic,:,:].shape)
    noise2 = np.random.normal(0,noise_strength,LOS_object.fourier_data[harmonic,:,:].shape)
    print 'noise1 rms :{}, std_dev : {}'.format(np.sqrt(np.mean(noise1**2)), noise_strength)

    LOS_object.fourier_data[harmonic,:,:] += (noise1 + 1j*noise2)
    tmp_meas = np.abs(LOS_object.fourier_data[harmonic,:,:])*LOS_object.valid_channels
    tmp_meas[-LOS_object.valid_channels] = np.nan
    amp_cmap = 'jet'
    noise_im = noise_im_ax.imshow(np.abs(tmp_meas[LOS_object.start_indices[1]:LOS_object.end_indices[1]]),origin='upper',aspect='auto',interpolation='nearest', cmap=amp_cmap)
    orig_im.set_clim(noise_im.get_clim())
    cbar = pt.colorbar(orig_im, cax = cbar_amp_ax, orientation = 'horizontal')
    cbar.set_ticks([1,2,3,4])
    cbar.set_label('Amp (a.u.)')
    phase_ax.set_xlabel('$\sqrt{s}$')
    for i in [phase_ax, amp_ax]:i.set_yticks(i.get_yticks()[::2])
    for i in [phase_ax, amp_ax]:i.grid()
    phase_ax.set_xticks(phase_ax.get_xticks()[::2])
    orig_im_ax.axes.get_yaxis().set_visible(False)
    noise_im_ax.axes.get_yaxis().set_visible(False)
    orig_im_ax.axes.get_xaxis().set_visible(False)
    noise_im_ax.axes.get_xaxis().set_visible(False)
    orig_im_ax.set_title('Original')
    noise_im_ax.set_title('Noisy')
    
    tomo_orient_extrap = tomo_orient

    foo1, foo2, foo3, noise_inv = single(tomo_modes_n, tomo_modes_m, tomo_orient, 0.5, 'Direct', 300, 0, 1, True, 1, None, tomo_orient, LOS_object)
    plot_radial_structure(noise_inv.T, np.sqrt(LOS_object.segment_midpoints), tomo_modes_n, tomo_modes_m, prov_ax = [amp_ax, phase_ax], norm = False, extra_txt = 'Noisy ', single_mode = None)

    #restore the data
    LOS_object.fourier_data = +fourier_data_tmp
    gs.tight_layout(fig, pad = 0.2)
    if filename!=None:
        fig.savefig(filename + '.pdf')
        fig.savefig(filename + '.eps')
    fig.canvas.draw(); fig.show()



def many_measurements_shift(LOS_object, tomo_modes_n, tomo_modes_m, tomo_orient, harmonic, peak_loc = 0.5, peak_width = 0.1, noise_strength = 0.00001, filename = '', shift_left=False, shift_down = False):
    filename_list = []
    fig, ax2 = pt.subplots(nrows = 2, ncols = 2, sharex = True, sharey = True)
    ax = ax2.flatten()
    cm_to_inch=0.393701
    fig.set_figwidth(8.48*cm_to_inch)
    fig.set_figheight(8.48*0.75*cm_to_inch)
    count = 0
    markers = ['.','x','+']
    left_amount = None; down_amount = None
    for j, peak_loc in enumerate([0.2, 0.4, 0.6, 0.8]):
        for i, (shift, cur_marker) in enumerate(zip([1,5,10], markers)):
            filename = str(i)
            if shift_left: left_amount = shift
            if shift_down: down_amount = shift
            T_inv, T = create_measurements(LOS_object, tomo_modes_n, tomo_modes_m, tomo_orient, harmonic, peak_loc = peak_loc, peak_width = 0.1, noise_strength = 0.001,shift_left = left_amount, shift_down = down_amount, filename = filename, plot_diff_images = False)

            filename_list.append('reproj_diff_{}'.format(filename))
            ax[count].plot(np.sqrt(LOS_object.segment_midpoints), np.abs(T_inv), label='{:.2f} shift'.format(shift), marker = cur_marker)
            #ax[1,j].plot(np.sqrt(LOS_object.segment_midpoints), np.angle(T_inv))
        ax[count].plot(np.sqrt(LOS_object.segment_midpoints), np.abs(T),label='Correct', marker = '')
        count += 1
        #ax[1,j].plot(np.sqrt(LOS_object.segment_midpoints), np.angle(T))
    
    ax[0].legend(loc='best')
    ax[2].set_xlabel('$\sqrt{s}$')
    ax[3].set_xlabel('$\sqrt{s}$')
    ax[0].set_ylabel('Amp (a.u)')
    ax[2].set_ylabel('Amp (a.u)')
    for i in ax:i.grid()
    for i in [ax2[1,0]]:i.set_xticks(i.get_xticks()[::2])
    for i in [ax2[0,0]]:i.set_yticks(i.get_yticks()[::2])
    fig.tight_layout(pad = 0.01)
    if shift_left:
        fig.savefig('shift_left_tollerance.pdf')
    else:
        fig.savefig('shift_down_tollerance.pdf')
    fig.canvas.draw(); fig.show()


def many_measurements_vary_m_n(LOS_object, tomo_modes_n, tomo_modes_m, tomo_orient, harmonic, peak_loc = 0.5, peak_width = 0.1, noise_strength = 0.00001, filename = ''):
    filename_list = []
    for m_new in tomo_modes_m:
        for n_new in tomo_modes_n:
            filename = 'm_{}_n_{}'.format(m_new, n_new, )
            create_measurements(LOS_object, [n_new], [m_new], tomo_orient, harmonic, peak_loc = peak_loc, peak_width = 0.1, noise_strength = noise_strength, filename = filename)
            filename_list.append('reproj_diff_m_{}_n_{}'.format(m_new, n_new))
    print filename_list
    for i in filename_list:
        os.system('convert {}.pdf {}.png'.format(i, i))

    os.system('convert -delay {} -loop 0 {} {}'.format(20, '.png '.join(filename_list)+'.png', 'scan_n_m.gif'.format('{:.2f}'.format(noise_strength).replace('.','_'))))
    #os.system('zip {} {}'.format(save_name.rstrip('gif')+'zip', image_string))


def create_measurements(LOS_object,tomo_modes_n, tomo_modes_m, tomo_orient, harmonic, peak_loc = 0.5, peak_width = 0.1, noise_strength = 0.00001, shift_left = 0, shift_down = 0, filename = '', plot_image = False, plot_diff_images = True):
    '''Generates fake measurement data to see how the tomographic reconstruction behaves.
    Data is a gaussian with peak_loc = mean and peak_width=std 
    The phase is a constant value across the peak

    SRH : 12Feb2014
    '''
    
    geom_matrix, tomo_valid_channels, tomo_measurements = LOS_object.return_combined_matrices(tomo_modes_n, tomo_modes_m, harmonic, tomo_orient)
    s_vals = LOS_object.segment_midpoints
    rv = norm(loc = peak_loc, scale = peak_width)
    amps = rv.pdf(np.sqrt(LOS_object.segment_midpoints))
    phases = amps*0+2.13
    T = amps*np.cos(phases) + 1j*amps*np.sin(phases)
    re_projection = np.dot(geom_matrix, T)

    new_meas = np.zeros(tomo_valid_channels.shape, dtype=complex)
    new_valid_channels = np.zeros(tomo_valid_channels.shape, dtype=bool)
    sig_strength = np.sqrt(np.sum(np.abs(re_projection)**2))
    
    print 'max amp :{}, mean amp:{}'.format(np.max(np.abs(re_projection)), np.mean(np.abs(re_projection)))
    noise1 = np.random.normal(0,noise_strength,len(re_projection))
    noise2 = np.random.normal(0,noise_strength,len(re_projection))

    #noise1=noise_strength * noise1*np.sqrt(np.sum(np.abs(re_projection)**2))/np.sqrt(np.sum(np.abs(noise1)**2))
    #noise2=noise_strength * noise1*np.sqrt(np.sum(np.abs(re_projection)**2))/np.sqrt(np.sum(np.abs(noise2)**2))
    #noise2 = np.random.normal(0,0,len(re_projection))

    valid_channels_bak = LOS_object.valid_channels.copy()
    fourier_data_tmp = LOS_object.fourier_data.copy()

    new_meas[tomo_valid_channels] = re_projection + noise1 + 1j*noise2
    if shift_left!=None:
        new_meas[:,:-shift_left] = new_meas[:, shift_left:]
        new_valid_channels[:, :-shift_left] = +LOS_object.valid_channels[:, shift_left:]
    if shift_down!=None:
        new_meas[:-shift_down,:] = new_meas[shift_down:, :]
        new_valid_channels[:-shift_down, :] = +LOS_object.valid_channels[shift_down:, :]
    #new_meas[new_meas<0]=0

    if plot_image:
        fig, ax = pt.subplots(ncols = 2)
        ax[0].imshow(np.abs(new_meas), aspect = 'auto')
        ax[1].imshow(np.angle(new_meas), aspect = 'auto')

        fig.canvas.draw(); fig.show()
    tomo_orient_extrap = tomo_orient

    #fourier_data[harmonic,:,:] = +new_meas
    LOS_object.fourier_data[harmonic,:,:] = +new_meas    
    LOS_object.valid_channel = +new_valid_channels
    tomo_inv = run_inversion(tomo_modes_n, tomo_modes_m, tomo_orient, tomo_orient_extrap,  filename, LOS_object, harmonic = harmonic, method = 'Direct', cycles=300, lamda=0.5, cut_values = None, plot_wave_fields = False, plot_old_reproj_comparison = False, plot_old_combo = False, plot_reprojection_comp1 = False, plot_reprojection_comp2 = plot_diff_images, plot_wave_animation = False, tomo_DC = None)


    #restore the data
    LOS_object.fourier_data = +fourier_data_tmp
    LOS_object.valid_channels = +valid_channels_bak
    return tomo_inv.T, T
    
def complex_array_to_rgb(X, theme='dark', rmax=None):
    '''Takes an array of complex number and converts it to an array of [r, g, b],
    where phase gives hue and saturaton/value are given by the absolute value.
    Especially for use with imshow for complex plots.

    Taken from here:
    http://stackoverflow.com/questions/15207255/is-there-any-way-to-use-bivariate-colormaps-in-matplotlib

    SRH: 12Feb2014
    '''
    absmax = rmax or np.abs(X).max()
    Y = np.zeros(X.shape + (3,), dtype='float')
    Y[..., 0] = np.angle(X) / (2 * np.pi) % 1
    if theme == 'light':
        Y[..., 1] = np.clip(np.abs(X) / absmax, 0, 1)
        Y[..., 2] = 1
    elif theme == 'dark':
        Y[..., 1] = 1
        Y[..., 2] = np.clip(np.abs(X) / absmax, 0, 1)
    Y = mpl.colors.hsv_to_rgb(Y)
    return Y

def hue_sat_cbar(cbar_wave_ax, rmax = 1):
    '''Create a colorbar for the hue/sat images

    SRH: 12Feb2014
    '''
    phase = np.linspace(-np.pi,np.pi,50)
    amp = np.linspace(0,1,50)
    phase_mesh, amp_mesh  = np.meshgrid(phase, amp)
    new_array = amp_mesh * np.cos(phase_mesh) + 1j*amp_mesh *np.sin(phase_mesh)
    im = cbar_wave_ax.imshow(complex_array_to_rgb(new_array, rmax=rmax, theme = 'dark'), aspect='auto', origin = 'lower', extent=[-np.pi, np.pi, 0, 1])
    cbar_wave_ax.set_xlabel('phase')
    cbar_wave_ax.set_ylabel('amp')
    cbar_wave_ax.set_xticks([-3,0,3])
    cbar_wave_ax.set_yticks([0,1])


def calc_RMSE(b, plot = False):
    M_reproj = np.dot(b.geom_matrix, b.T)
    NRMSD_list = []; RMS_ratio_list = []
    if plot: fig, ax = pt.subplots()
    for i in np.linspace(0,2.*np.pi,16):
        M_reproj_curr = np.real(M_reproj * np.exp(1j*i))
        meas_cur = np.real(b.measurements * np.exp(1j*i))
        diff = M_reproj_curr - meas_cur
        RMSD = np.sqrt(np.mean(diff**2))
        RMS = np.sqrt(np.mean(meas_cur**2))
        NRMSD = RMSD/(np.max(meas_cur) - np.min(meas_cur))
        #print RMSD, RMS, RMSD/RMS, NRMSD, RMSD/np.mean(np.abs(meas_cur))
        RMS_ratio_list.append(RMSD/RMS)
        NRMSD_list.append(NRMSD)
    if plot:
        ax.plot(np.linspace(0,2.*np.pi,50), NRMSD_list, label = 'NRMSD')
        ax.plot(np.linspace(0,2.*np.pi,50), RMS_ratio_list, label = 'RMSE/RMS')
        ax.set_xlabel('Phase')
        ax.set_ylabel('Error Measure')
        ax.legend(loc='best')
        fig.canvas.draw(); fig.show()
    return [np.mean(NRMSD)]


def new_error(b):
    M_reproj = np.dot(b.geom_matrix, b.T)
    print np.sqrt(np.mean(np.abs(M_reproj - b.measurements)))/(np.max(np.abs(M_reproj)))


