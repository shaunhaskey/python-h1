import h1.h1model.plot_functions as h1_plot
from scipy.interpolate import griddata as scipy_griddata
import numpy as np
import matplotlib.pyplot as pt
import BOOZER
#Note mayavi.mlab is imported for certain functions

'''
This set of classes provides a way to determine the geometry of LOS
diagnostics in Boozer co-ordinates

SurfacePatch : generates a patch in realspace that is based on a
regular grid in Boozer space. This is to provide points to interpolate
the LOS onto to find its Boozer coordinates

LOS : Class for a line of sight diagnostic. Determines the lines of
sight, intersection with LCFS, and can also provide the interpolation
onto a SurfacePatch grid

imax_camera : contains the geometry for the imax_camera

interferometer : contains the geometry for the interferometer

SRH: 12July2013
'''


class SurfacePatch():
    def __init__(self, phi_min, phi_max, n_phi = 10, no_theta = 50, boozer_filename = None, boozer_object = None):
        '''
        This creates a last closed surface mesh between phi_min and
        phi_max. It performs some calculations to make the LOS
        intersections with the surface easier.
        phi_min : boozer phi min for surfaces for interpolating grid
        phi_max : boozer phi max for surfaces for interpolating grid
        n_phi : number of phi slices in the grid
        no_theta : number of theta_boozer angles
        boozer_filename : location of boozer netCDF file to use
        boozer_object : BOOZER class object to use

        SRH : 11July2013
        '''
        self.phi_min = phi_min
        self.phi_max = phi_max
        self.n_phi = n_phi
        self.no_theta = no_theta
        if boozer_object==None and boozer_filename!=None:
            boozer_object = BOOZER.BOOZER(boozer_filename,import_all=True,compute_spline_type=0,load_spline=False,save_spline=False,load_grid=False)
        else:
            raise ValueError('Need to supply boozer_object or boozer_filename')
        self.boozer_object = boozer_object
        self.get_points_for_mesh()
        self.describe_surfaces()
        self.obtained_grid = 0

    def get_points_for_mesh(self,):
        '''
        Creates a regular grid in BOOZER coords for the LCFS, and
        transforms it into realspace to provide a grid to calculate
        the intersection with lines of sight

        SRH: 12July2013
        '''
        self.phi_vals = np.linspace(self.phi_min, self.phi_max, self.n_phi)
        self.cross_sect_x = []; self.cross_sect_y = []; self.cross_sect_z = []
        for phi_val in self.phi_vals:
            phi_val = phi_val / 180.*np.pi
            cross_sect_cyl2, cross_sect_booz2 = self.boozer_object.getCrossSectData(phi_val, s_b=None, s_ind=[-1],no_theta=self.no_theta, phi_is_cyl=False, coordsys='cart', return_booz_also=1)
            self.cross_sect_x.extend(cross_sect_cyl2[:,0,:].flatten())
            self.cross_sect_y.extend(cross_sect_cyl2[:,1,:].flatten())
            self.cross_sect_z.extend(cross_sect_cyl2[:,2,:].flatten())
        self.create_triangular_mesh()
        self.vertices = np.array([self.cross_sect_x, self.cross_sect_y, self.cross_sect_z]).T

    def get_grid_for_interpolation(self, no_theta = 50, s_increment = 1):
        '''
        Get a grid in Boozer as defined by s_increment, no_theta,
        phi_min(deg), phi_max(deg), n_phi. transform grid to realspace so
        there is a grid to interpolate onto

        SRH : 10July2013
        '''
        self.grid_no_theta = no_theta
        self.grid_s_increment = s_increment
        n_surfaces = self.boozer_object.es_b.shape[0]
        s_ind = range(0,n_surfaces, self.grid_s_increment)
        #cartesian
        self.grid_x = []; self.grid_y = []; self.grid_z = []
        #Boozer
        self.grid_s = []; self.grid_phi = []; self.grid_theta = []
        #self.cart_points = np.zeros((len(s_ind),self.grid_no_theta,3),dtype=float)
        #self.booz_grid = self.cart_points * 0
        for phi_val in self.phi_vals:
            print('{:.2f}deg'.format(phi_val))
            phi_val = np.deg2rad(phi_val)
            cross_sect_cyl2, cross_sect_booz2 = self.boozer_object.getCrossSectData(phi_val,s_b=None, s_ind=s_ind,no_theta=self.grid_no_theta, phi_is_cyl=False,coordsys='cart', return_booz_also=1)
            self.grid_x.extend(cross_sect_cyl2[:,0,:].flatten())
            self.grid_y.extend(cross_sect_cyl2[:,1,:].flatten())
            self.grid_z.extend(cross_sect_cyl2[:,2,:].flatten())
            self.grid_s.extend(cross_sect_booz2[:,0,:].flatten())
            self.grid_phi.extend(cross_sect_booz2[:,1,:].flatten())
            self.grid_theta.extend(cross_sect_booz2[:,2,:].flatten())
        self.obtained_grid = 1

    def create_triangular_mesh(self,):
        '''
        Assembles a triangular mesh that can be used to find the
        intersection locations

        SRH: 12July2013
        '''
        cut_start_points = np.arange(self.n_phi)*self.no_theta
        self.faces = np.zeros((2*(self.no_theta-1)*(self.n_phi-1),3),dtype=int)
        for i in range(self.n_phi-1):
            ind1 = i*(self.no_theta-1)
            ind2 = ind1 + self.no_theta - 1
            self.faces[ind1:ind2,:] = np.array([range(cut_start_points[i],cut_start_points[i]+self.no_theta-1), 
                                           range(cut_start_points[i+1],cut_start_points[i+1]+self.no_theta-1), 
                                           range(cut_start_points[i+1]+1,cut_start_points[i+1]+self.no_theta)]).T
            #faces[ind1:ind2,:] = np.array([range(cut_start_points[i],cut_start_points[i]+no_theta-1), range(cut_start_points[i+1]+1,cut_start_points[i+1]+no_theta), range(cut_start_points[i+1],cut_start_points[i+1]+no_theta-1)]).T
            ind1 = (i+self.n_phi-1)*(self.no_theta-1)
            ind2 = ind1 + self.no_theta-1
            self.faces[ind1:ind2,:] = np.array([range(cut_start_points[i],cut_start_points[i]+self.no_theta-1), 
                                                range(cut_start_points[i+1]+1,cut_start_points[i+1]+self.no_theta),  
                                                range(cut_start_points[i]+1,cut_start_points[i]+self.no_theta)]).T
            #faces[ind1:ind2,:] = np.array([range(cut_start_points[i],cut_start_points[i]+no_theta-1), range(cut_start_points[i+1]+1,cut_start_points[i+1]+no_theta),  range(cut_start_points[i]+1,cut_start_points[i]+no_theta)]).T

    def describe_surfaces(self,):
        '''Calculate the normals, d values, and relevant triangles for a
        mesh of a plasma Vertices are a list of points for the mesh, faces
        contains the mesh points for each triangle

        SRH : 10July2013
        '''
        #find the normal to the faces and normalise
        #https://sites.google.com/site/dlampetest/python/calculating-normals-of-a-triangle-mesh-using-numpy
        #find d = (\vec{n} \dot \vec{x}) by substituting in one of the points on the triangle
        self.tris = self.vertices[self.faces]
        self.n = np.cross(self.tris[:,1,:] - self.tris[:,0,:], self.tris[:,2,:] - self.tris[:,0,:] )
        lens = np.sqrt( self.n[:,0]**2 + self.n[:,1]**2 + self.n[:,2]**2 )
        self.n[:,0] /= lens; self.n[:,1] /= lens; self.n[:,2] /= lens                
        self.d = np.sum(self.n*self.tris[:,0,:],axis=1)

        #pt1, pt2, pt3 are precalculated for the in-out test later on
        self.pt1 = self.tris[:,1,:]-self.tris[:,0,:]
        self.pt2 = self.tris[:,2,:]-self.tris[:,1,:]
        self.pt3 = self.tris[:,0,:]-self.tris[:,2,:]

    def plot_mesh(self,):
        '''
        Plot the mesh that has been created

        SRH: 12 July2013
        '''
        from mayavi import mlab
        mlab.triangular_mesh(self.vertices[:,0], self.vertices[:,1], self.vertices[:,2],self.faces)
        mlab.triangular_mesh(self.vertices[:,0], self.vertices[:,1], self.vertices[:,2],self.faces,representation='wireframe',color=(0,0,0))


class LOS():
    def __init__(self, CCD_position, CCD_x, CCD_y, CCD_focal_distance, CCD_pixels_x, CCD_pixels_y, u_hat, patch, v_hat = None, w_hat = None, CCD_focal_point = None):
        '''
        CCD_position : location of the centre of the CCD in cart coords (m)
        CCD_x, CCD_y : length and width of the CCD in (m)
        CCD_pixels_x, CCD_pixels_y : number of pixels on the CCD
        CCD_focal_distance : focal distance infront of CCD in (m)
        u_hat : unit vector in vertical direction of CCD side
        v_hat : unit vector along CCD
        w_hat : unit vector normal to CCD
        u_hat x v_hat = w_hat
        must provide either v_hat or w_hat

        SRH : 11July2013
        '''
        #CCD_x = CCD_R*np.cos(np.deg2rad(CCD_phi))
        #CCD_y = CCD_R*np.sin(np.deg2rad(CCD_phi))
        #vector for location of CCD
        self.CCD_position = CCD_position
        self.CCD_length = CCD_y
        self.CCD_width = CCD_x
        self.f = CCD_focal_distance #17.3/1000
        self.u_hat = u_hat
        self.CCD_pixels_x = CCD_pixels_x
        self.CCD_pixels_y = CCD_pixels_y
        #u hat - verical direction on CCD
        if w_hat == None:
            self.v_hat = v_hat
            self.w_hat = np.cross(self.u_hat, self.v_hat)
        else:
            self.v_hat = np.cross(self.w_hat, self.u_hat)
        self.focal_point = self.f*self.w_hat + self.CCD_position
        self.calc_point1_gradient()
        self.patch = patch
        self.find_intersections()

    def calc_point1_gradient(self,):
        self.LOS_point = np.zeros((self.CCD_pixels_x, self.CCD_pixels_y, 3),dtype=float)
        self.LOS_gradient = self.LOS_point * 0
        for pixel_x in range(self.CCD_pixels_x):
            for pixel_y in range(self.CCD_pixels_y):
                print pixel_x,pixel_y
                self.LOS_point[pixel_x, pixel_y, :] = (float(pixel_x)/self.CCD_pixels_x - 0.5) * self.CCD_width * self.v_hat + (float(pixel_y)/self.CCD_pixels_y - 0.5) * self.CCD_length * self.u_hat + self.CCD_position
                self.LOS_gradient[pixel_x, pixel_y, :] = self.focal_point - self.LOS_point[pixel_x, pixel_y, :]

    def find_intersections(self,):
        '''Find the intersection for each triangle plane
        SRH: 10July2013
        '''
        print('finding intersections')
        self.valid_channels = np.zeros((self.CCD_pixels_x,self.CCD_pixels_y),dtype=bool)
        self.intersection1 = self.LOS_point * 0
        self.intersection2 = self.LOS_point * 0
        for pixel_x in range(self.CCD_pixels_x):
            for pixel_y in range(self.CCD_pixels_y):
                P0 = self.LOS_point[pixel_x, pixel_y, :]
                dP = self.LOS_gradient[pixel_x, pixel_y, :]
                t = (self.patch.d - np.sum(self.patch.n * P0,axis=1))/np.sum(self.patch.n*dP,axis=1)
                R = P0 + t[:,np.newaxis]*dP[np.newaxis,:]
                tmp1 = np.sum(self.patch.n*np.cross(self.patch.pt1,R - self.patch.tris[:,0,:]),axis=1)>=0
                tmp1[tmp1] = np.sum(self.patch.n[tmp1]*np.cross(self.patch.pt2[tmp1,:],R[tmp1,:] - self.patch.tris[tmp1,1,:]),axis=1)>=0
                tmp1[tmp1] = np.sum(self.patch.n[tmp1]*np.cross(self.patch.pt3[tmp1,:],R[tmp1,:] - self.patch.tris[tmp1,2,:]),axis=1)>=0
                #a, intersect_pts = find_intersections2(P0,dP,d,n,pt1,pt2,pt3,tris)
                if (np.sum(tmp1)) >= 2 and ((np.sum(tmp1)%2)==0):
                    #mlab.points3d(intersect_pts[:,0], intersect_pts[:,1], intersect_pts[:,2],scale_factor=0.02)
                    intersection_order = np.argsort(t[tmp1])
                    self.valid_channels[pixel_x,pixel_y] = True
                    self.intersection1[pixel_x, pixel_y] = R[tmp1][intersection_order[0]]
                    self.intersection2[pixel_x, pixel_y] = R[tmp1][intersection_order[1]]
                elif np.sum(tmp1)%2 == 1:
                    print('Warning odd number of intersections!!!, {}'.format(np.sum(tmp1)))
        return tmp1, R[tmp1]

    def create_interpolation_points(self, n_interp_pts):
        self.interpolation_pts = np.zeros((self.valid_channels.shape[0], self.valid_channels.shape[1],n_interp_pts,3),dtype=float)
        self.dl = np.zeros((self.valid_channels.shape[0], self.valid_channels.shape[1]),dtype=float)
        for pixel_x in range(self.CCD_pixels_x):
            for pixel_y in range(self.CCD_pixels_y):
                if self.valid_channels[pixel_x, pixel_y]:
                    self.interpolation_pts[pixel_x, pixel_y,:,0] = np.linspace(self.intersection1[pixel_x,pixel_y,0],self.intersection2[pixel_x, pixel_y,0],n_interp_pts)
                    self.interpolation_pts[pixel_x, pixel_y,:,1] = np.linspace(self.intersection1[pixel_x,pixel_y,1],self.intersection2[pixel_x, pixel_y,1],n_interp_pts)
                    self.interpolation_pts[pixel_x, pixel_y,:,2] = np.linspace(self.intersection1[pixel_x,pixel_y,2],self.intersection2[pixel_x, pixel_y,2],n_interp_pts)
                    self.dl[pixel_x, pixel_y] = np.sqrt(np.sum((self.interpolation_pts[pixel_x, pixel_y,1,:] - self.interpolation_pts[pixel_x, pixel_y,0,:])**2))

    def perform_interpolation(self,no_theta = 50, s_increment = 10):
        ''' 
        valid_channels : bool array listing the channels that intersect the plasma
        cross_sect_s : 1D array of s values (corresponding to the interpolation_pts
        cross_sect_theta : 1D array of theta_b values (corresponding to the interpolation_pts
        cross_sect_phi : 1D array of phi_b values (corresponding to the interpolation_pts

        interpolation_pts : grid of points in realspace [channel_num, interp_point, x/y/z]
        SRH: 10July2013
        '''
        valid_channels = self.valid_channels
        interp_pts_x = self.interpolation_pts[valid_channels,:,0].flatten()
        interp_pts_y = self.interpolation_pts[valid_channels,:,1].flatten()
        interp_pts_z = self.interpolation_pts[valid_channels,:,2].flatten()
        required_shape = self.interpolation_pts[valid_channels,:,0].shape
        self.interp_boozer = self.interpolation_pts*0
        if self.patch.obtained_grid != 1:
            print 'getting patch grid'
            self.patch.get_grid_for_interpolation(no_theta = no_theta, s_increment = s_increment)
        points_tuple = (self.patch.grid_x, self.patch.grid_y, self.patch.grid_z)
        cross_sect_s = self.patch.grid_s
        cross_sect_phi = self.patch.grid_phi
        cross_sect_theta = self.patch.grid_theta
        print 's'
        self.interp_boozer[valid_channels,:,0] = scipy_griddata(points_tuple, np.array(cross_sect_s), (interp_pts_x, interp_pts_y, interp_pts_z)).reshape(required_shape)
        print 'theta'
        interp_data_theta_sin = scipy_griddata(points_tuple, np.sin(cross_sect_theta), (interp_pts_x, interp_pts_y, interp_pts_z))
        interp_data_theta_cos = scipy_griddata(points_tuple, np.cos(cross_sect_theta), (interp_pts_x, interp_pts_y, interp_pts_z))
        self.interp_boozer[valid_channels,:,1] = np.arctan2(interp_data_theta_sin, interp_data_theta_cos).reshape(required_shape)
        print 'phi'
        interp_data_phi_sin = scipy_griddata(points_tuple, np.sin(cross_sect_phi), (interp_pts_x, interp_pts_y, interp_pts_z))
        interp_data_phi_cos = scipy_griddata(points_tuple, np.cos(cross_sect_phi), (interp_pts_x, interp_pts_y, interp_pts_z))
        self.interp_boozer[valid_channels,:,2] = np.arctan2(interp_data_phi_sin, interp_data_phi_cos).reshape(required_shape)

    def plot_intersections(self,):
        '''Plot the intersection points on the plasma surface
        SRH: 12July2013
        '''
        from mayavi import mlab
        mlab.points3d(self.intersection1[self.valid_channels,0], self.intersection1[self.valid_channels,1], self.intersection1[self.valid_channels,2],scale_factor=0.02, color=(0,0,1))
        mlab.points3d(self.intersection2[self.valid_channels,0], self.intersection2[self.valid_channels,1], self.intersection2[self.valid_channels,2],scale_factor=0.02,color=(1,0,0))


    def plot_LOS(self, to_focal_point = 0, grad_mult = 50, plot_arguments = None):
        '''
        Plot the line of sights on an mlab figure
        Can either plot from CCD - to focal point or give grad_mult value
        which extends the line grad_mult times past the focal point
        SRH: 12July2013
        '''
        from mayavi import mlab
        x = []; y = []; z = []; connections = []
        index = 0
        for x_tmp, y_tmp, z_tmp, x_grad, y_grad, z_grad in zip(self.LOS_point[:,:,0].flatten(), 
                                                               self.LOS_point[:,:,1].flatten(), 
                                                               self.LOS_point[:,:,2].flatten(),
                                                               self.LOS_gradient[:,:,0].flatten(),
                                                               self.LOS_gradient[:,:,1].flatten(),
                                                               self.LOS_gradient[:,:,2].flatten()):
            if to_focal_point:
                x.append(x_tmp); x.append(self.focal_point[0])
                y.append(y_tmp); y.append(self.focal_point[1])
                z.append(z_tmp); z.append(self.focal_point[2])
            else:
                x.append(x_tmp); x.append(x_grad*grad_mult+x_tmp)
                y.append(y_tmp); y.append(y_grad*grad_mult+y_tmp)
                z.append(z_tmp); z.append(z_grad*grad_mult+z_tmp)
                connections.append([index,index+1])
            index+=2

        src = mlab.pipeline.scalar_scatter(np.array(x), np.array(y), np.array(z))
        src.mlab_source.dataset.lines = np.array(connections)
        lines = mlab.pipeline.stripper(src)
        line_args = {'line_width':1.5,'opacity':0.5}
        if plot_arguments!=None:
            line_args.update(plot_arguments)
        mlab.pipeline.surface(lines,**line_args)#line_width=1.5, opacity=.1)

    def orient_camera(self,mayavi_fig,render=1):
        '''
        Orient the mayavi camera to the location of the actual camera
        Pass the mayavi figure

        SRH: 12July2013
        '''
        cam = mayavi_fig.scene.camera
        ren = mayavi_fig.scene.renderer
        cam.view_angle=np.rad2deg(np.arctan2(self.CCD_width/2,self.f))*2
        cam.position = self.CCD_position
        cam.focal_point = self.focal_point
        cam.compute_view_plane_normal()
        ren.reset_camera_clipping_range()
        if render:
            mayavi_fig.scene.render()

def imax_camera(boozer_filename = '/home/srh112/code/python/h1_eq_generation/results7/kh0.350-kv1.000fixed/boozmn_wout_kh0.350-kv1.000fixed.nc', plot_LOS = 0, plot_patch = 0, plot_intersections = 0, plot_pfc = 0, plot_tfc = 0):
    '''Convenience function containing the required geometry for the
    imax camera

    SRH: 12July2013
    '''
    # LOS(CCD, CCD_width, CCD_length)
    # #Nandi tomo17
    # CCD = np.array([-1.00239,1.673,0.0150])
    # w_hat = np.array([0.564,-0.82468,0.037002])
    # v_hat = np.cross(w_hat,u_hat)
    # f = 17./1000
    #phi_min =60; phi_max = 200;n_phi = 30;no_theta = 50
    #boozer_filename = '/home/srh112/code/python/h1_eq_generation/results7/kh0.350-kv1.000fixed/boozmn_wout_kh0.350-kv1.000fixed.nc'
    CCD_phi = 120; CCD_R = 1.946
    CCD_z = 5./100; CCD_x = CCD_R*np.cos(np.deg2rad(CCD_phi)); CCD_y = CCD_R*np.sin(np.deg2rad(CCD_phi))
    CCD_position = np.array([CCD_x, CCD_y, CCD_z])
    v_hat = np.array([-np.cos(np.deg2rad(90-CCD_phi)), np.sin(np.deg2rad(90-CCD_phi)),0])
    u_hat = np.array([0,0,1])
    CCD_L = 0.01575
    CCD_x = CCD_L; CCD_y = CCD_L
    n_pixels = 512
    CCD_pixels_x = 512/16; CCD_pixels_y = 512/16
    CCD_focal_distance = 17.3/1000.
    phi_min =CCD_phi - 50; phi_max = CCD_phi + 50;n_phi = 50;no_theta = 50;n_interp_pts = 100
    phi_min =CCD_phi - 180; phi_max = CCD_phi + 180;n_phi = 60;no_theta = 50;n_interp_pts = 100
    boozer_filename = '/home/srh112/code/python/h1_eq_generation/results7/kh0.350-kv1.000fixed/boozmn_wout_kh0.350-kv1.000fixed.nc'
    plot_length = 50
    patch =  SurfacePatch(phi_min, phi_max, n_phi = n_phi, no_theta = no_theta, boozer_filename = boozer_filename)
    answer = LOS(CCD_position, CCD_x, CCD_y, CCD_focal_distance, CCD_pixels_x, CCD_pixels_y, u_hat, patch, v_hat = v_hat, w_hat = None, CCD_focal_point = None)
    if plot_LOS or plot_patch or plot_intersections or plot_pfc or plot_tfc:
        import mayavi.mlab as mlab
        mayavi_fig = mlab.figure(1, fgcolor=(0, 0, 0), bgcolor=(1, 1, 1), size=(512, 512))
        if plot_patch : answer.patch.plot_mesh()
        if plot_pfc : h1_plot.plot_pfc()
        if plot_tfc : h1_plot.plot_tfc([11,12,13])
        if plot_intersections : answer.plot_intersections()
        if plot_LOS : answer.plot_LOS(grad_mult = plot_length,plot_arguments={'color':(0,1,1), 'line_width':2.5})
    return answer

def interferometer(boozer_filename = '/home/srh112/code/python/h1_eq_generation/results7/kh0.350-kv1.000fixed/boozmn_wout_kh0.350-kv1.000fixed.nc', plot_LOS = 0, plot_patch = 0, plot_intersections = 0, plot_pfc = 0, plot_tfc = 0):
    '''Convenience function containing the required geometry for the
    interferometer

    SRH: 12July2013
    '''
    CCD_phi = 240; CCD_R = 2.4
    CCD_z = 0./100; CCD_x = CCD_R*np.cos(np.deg2rad(CCD_phi)); CCD_y = CCD_R*np.sin(np.deg2rad(CCD_phi))
    CCD_position = np.array([CCD_x, CCD_y, CCD_z])
    v_hat = np.array([-np.cos(np.deg2rad(90-CCD_phi)), np.sin(np.deg2rad(90-CCD_phi)),0])
    u_hat = np.array([0,0,1])
    CCD_L = 40./100
    CCD_x = 0.001; CCD_y = CCD_L
    n_pixels = 512
    CCD_pixels_x = 1; CCD_pixels_y = 50
    CCD_focal_distance = 10000
    plot_length = 2./CCD_focal_distance
    phi_min =CCD_phi - 30; phi_max = CCD_phi+30;n_phi = 30;no_theta = 50;n_interp_pts = 100
    patch =  SurfacePatch(phi_min, phi_max, n_phi = n_phi, no_theta = no_theta, boozer_filename = boozer_filename)
    answer = LOS(CCD_position, CCD_x, CCD_y, CCD_focal_distance, CCD_pixels_x, CCD_pixels_y, u_hat, patch, v_hat = v_hat, w_hat = None, CCD_focal_point = None)
    if plot_LOS or plot_patch or plot_intersections or plot_pfc or plot_tfc:
        import mayavi.mlab as mlab
        mayavi_fig = mlab.figure(1, fgcolor=(0, 0, 0), bgcolor=(1, 1, 1), size=(512, 512))
        if plot_patch : answer.patch.plot_mesh()
        if plot_pfc : h1_plot.plot_pfc()
        if plot_tfc : h1_plot.plot_tfc([11,12,13])
        if plot_intersections : answer.plot_intersections()
        if plot_LOS : answer.plot_LOS(grad_mult = plot_length,plot_arguments={'color':(0,1,1), 'line_width':2.5})
    return answer
