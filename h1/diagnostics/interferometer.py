'''
This interferometer class is used to calculate density profiles using one of three methods at the moment:
Polyfitting, Bessel function fitting, and a linear algebra inverse matrix method.

The polynomial fit is not the best because it can be negative and has a tendency to behave poorly

Zernig basis, legendre polynomials
Search for Bessel function tomography Nagayama plasma paper done 
tomography of m-1 mode structure in tokamak plasma using leastsquare fitting method and Fourier-Bessel expansions
Yoshio Nagayama, J. Appl. Phys 64, 2702 (1987)

It is also supposed to produce the predicted interferometer output due to a particular mode, and in future
calculate the best fit to interferometer data.

SH: 5May2013
'''

import scipy.optimize as optimize
from scipy.interpolate import griddata as scipy_griddata
#import heliac_vmec_utils as hv_utils
import os,copy, time, scipy, pickle
import matplotlib.pyplot as pt
import numpy as np
from StringIO import StringIO
#import heliac_worker_funcs as heliac
import scipy.interpolate as interp
import scipy.stats.distributions as dist

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
        s_amps = np.interp(s,eig_func[0],np.real(eig_func[1])) + np.interp(s,eig_func[0],np.imag(eig_func[1]))*1j
        
        #calc the wave field
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
                print a
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


def calculate_inverse_matrix2(obj,n_segments, eig_func = None, n=None, m=None, interpolation = None, mode_amps_input=None, channel_mask=None):
    '''Calculate the inverse matrix for the calculation on density profile
    based on a fixed number of segments, and some mode numbers
    SH: 2013
    '''
    print 'calculating geometric matrix and its pseudo inverse'
    obj.segments = np.linspace(0,1,n_segments + 1)
    #n_measurements = len(obj.LOS_coords)
    n_measurements = np.sum(channel_mask)
    geom_mat = np.zeros((n_measurements,n_segments), dtype=complex)
    s_values = (obj.segments[1:]+obj.segments[0:-1])/2.
    valid_channel_count = 0
    #for i in range(0,len(obj.LOS_coords)):
    start_pt = np.argmax(obj.valid_LOS)
    for i in range(0,512):
        if channel_mask[i]==True:
            count_list = []
            s_tmp = obj.LOS_coords[-start_pt+i][:,0]
            theta_tmp = obj.LOS_coords[-start_pt+i][:,1]
            phi_tmp = obj.LOS_coords[-start_pt+i][:,2]
            dl_tmp = obj.dl_list[-start_pt+i]
            cumul_points = 0
            mode_amps = np.exp(1j*n*phi_tmp + 1j*m*theta_tmp)
            #Increase the dl spacing by the factor interpolation
            #Requires interpolating the mode output at these places
            if interpolation != None:
                l = np.arange(0,len(obj.LOS_coords[-start_pt+i][:,0]))*dl_tmp
                dl_tmp = dl_tmp / interpolation
                l_new = np.arange(0,len(l)*interpolation)*dl_tmp
                mode_amps = np.interp(l_new, l, np.real(mode_amps)) + 1j* np.interp(l_new, l, np.imag(mode_amps))
                s_tmp = np.interp(l_new, l, s_tmp)
            #For each solution segment, add up the dl's in that segment
            for j in range(0,n_segments):
                rel_pts = (s_tmp>=obj.segments[j]) * (s_tmp<obj.segments[j+1])
                points = np.sum(rel_pts)
                cumul_points += points
                #geom_mat[i,j]= points * dl_tmp/(dl_tmp*len(s_tmp))
                geom_mat[valid_channel_count, j] = np.sum(mode_amps[rel_pts]  * dl_tmp)# / (dl_tmp*len(s_tmp))
                #geom_mat[i,j] = np.sum(np.exp(1j*n*phi_tmp[rel_pts] + 1j*m*theta_tmp[rel_pts]) * dl_tmp)# / (dl_tmp*len(s_tmp))
                count_list.append(points)
            #print count_list
            valid_channel_count += 1
    obj.geom_mat = geom_mat
    obj.geom_mat_pinv = np.linalg.pinv(geom_mat)
    
    if (eig_func != None) or (mode_amps!=None): 
        if mode_amps==None:
            profile = np.dot(obj.geom_mat_pinv, obj.mode_output)
        else:
            #profile = np.dot(obj.geom_mat_pinv, mode_amps_input[obj.valid_LOS])
            profile = np.dot(obj.geom_mat_pinv, mode_amps_input[channel_mask])
        #print 'hello', profile.shape, obj.segments.shape
        xaxis = (obj.segments[0:-1] + obj.segments[1:])/2.
        fig, ax = pt.subplots(nrows = 2)
        ax[0].plot(xaxis, np.abs(profile), 'o', label='tomo inversion')
        ax[1].plot(xaxis, np.angle(profile, deg = True), '-o')
        if eig_func!=None:
            ax[0].plot(eig_func[0], np.abs(eig_func[1]), '-x', label='Eig func')
            ax[1].plot(eig_func[0], np.angle(eig_func[1],deg=True), '-x', label='Eig func')
        ax[1].set_xlabel('s')
        ax[0].set_ylabel('Amplitude')
        ax[0].legend(loc='best')
        ax[1].set_ylabel('Phase')
        fig.canvas.draw(); fig.show()
        return profile, s_values
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
