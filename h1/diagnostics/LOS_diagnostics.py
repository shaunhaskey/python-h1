import MDSplus as MDS
import h1.h1model.plot_functions as h1_plot
import h1.diagnostics.imax as imax

from scipy.interpolate import griddata as scipy_griddata
import numpy as np
import matplotlib.pyplot as pt
import h1.mhd_eq.BOOZER as BOOZER
import cPickle as pickle
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

class GeneralSurface():
    def __init__(self, x_array, y_array, z_array):
        pass
        # x = []; y = []; z = []
        # self.edges, self.pts_per_edge = x_array.shape
        # for i in range(self.edges):
        #     x.extend(x_array[i,:])
        #     y.extend(y_array[i,:])
        #     z.extend(z_array[i,:])
        # self.vertices = np.array([x, y, z]).T
        # self.create_triangular_mesh()
        # self.describe_surfaces()

    def create_triangular_mesh(self,):
        '''
        Assembles a triangular mesh that can be used to find the
        intersection locations

        SRH: 12July2013
        '''
        x = []; y = []; z = []
        for i in range(self.edges):
            x.extend(self.x_array[i,:])
            y.extend(self.y_array[i,:])
            z.extend(self.z_array[i,:])
        self.vertices = np.array([x, y, z]).T
        cut_start_points = np.arange(self.edges)*self.pts_per_edge
        self.faces = np.zeros((2*(self.pts_per_edge-1)*(self.edges-1),3),dtype=int)
        for i in range(self.edges-1):
            ind1 = i*(self.pts_per_edge-1)
            ind2 = ind1 + self.pts_per_edge - 1
            self.faces[ind1:ind2,:] = np.array([range(cut_start_points[i],cut_start_points[i]+self.pts_per_edge-1), 
                                                range(cut_start_points[i+1],cut_start_points[i+1]+self.pts_per_edge-1), 
                                                range(cut_start_points[i+1]+1,cut_start_points[i+1]+self.pts_per_edge)]).T
            ind1 = (i+self.edges-1)*(self.pts_per_edge-1)
            ind2 = ind1 + self.pts_per_edge-1
            self.faces[ind1:ind2,:] = np.array([range(cut_start_points[i],cut_start_points[i]+self.pts_per_edge-1), 
                                                range(cut_start_points[i+1]+1,cut_start_points[i+1]+self.pts_per_edge),  
                                                range(cut_start_points[i]+1,cut_start_points[i]+self.pts_per_edge)]).T

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

    def find_intersections(self, P0, dP):
        '''Find the intersection for each triangle plane
        SRH: 10July2013
        '''
        t = (self.d - np.sum(self.n * P0,axis=1))/np.sum(self.n*dP,axis=1)
        R = P0 + t[:,np.newaxis]*dP[np.newaxis,:]
        tmp1 = np.sum(self.n*np.cross(self.pt1,R - self.tris[:,0,:]),axis=1)>=0
        tmp1[tmp1] = np.sum(self.n[tmp1]*np.cross(self.pt2[tmp1,:],R[tmp1,:] - self.tris[tmp1,1,:]),axis=1)>=0
        tmp1[tmp1] = np.sum(self.n[tmp1]*np.cross(self.pt3[tmp1,:],R[tmp1,:] - self.tris[tmp1,2,:]),axis=1)>=0
        #a, intersect_pts = find_intersections2(P0,dP,d,n,pt1,pt2,pt3,tris)
        return t, R, tmp1


class TFC(GeneralSurface):
    def __init__(self, **kwargs):
        self.x_array, self.y_array, self.z_array = h1_plot.tfc_points(**kwargs)
        self.edges, self.pts_per_edge = self.x_array.shape
        self.create_triangular_mesh()
        self.describe_surfaces()
    def plot(self,):
        import mayavi.mlab as mlab
        mlab.triangular_mesh(self.vertices[:,0], self.vertices[:,1], self.vertices[:,2],self.faces, color=(0.5,0.5,0.5))
        mlab.triangular_mesh(self.vertices[:,0], self.vertices[:,1], self.vertices[:,2],self.faces,representation='wireframe',color=(0,0,0))


class BoozerSurfacePatch(GeneralSurface):
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
            print(boozer_filename)
            boozer_object = BOOZER.BOOZER(boozer_filename,import_all=True,compute_spline_type=0,load_spline=False,save_spline=False,load_grid=False)
        else:
            raise ValueError('Need to supply boozer_object or boozer_filename')
        self.boozer_object = boozer_object
        self.get_points_for_mesh()
        self.create_triangular_mesh()
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
        self.x_array = np.zeros((len(self.phi_vals), self.no_theta),dtype=float)
        self.edges, self.pts_per_edge = self.x_array.shape
        self.y_array= +self.x_array
        self.z_array= +self.x_array
        for i, phi_val in enumerate(self.phi_vals):
            phi_val = phi_val / 180.*np.pi
            cross_sect_cyl2, cross_sect_booz2 = self.boozer_object.getCrossSectData(phi_val, s_b=None, s_ind=[-1],no_theta=self.no_theta, phi_is_cyl=False, coordsys='cart', return_booz_also=1)
            self.x_array[i,:] = cross_sect_cyl2[:,0,:].flatten()
            self.y_array[i,:] = cross_sect_cyl2[:,1,:].flatten()
            self.z_array[i,:] = cross_sect_cyl2[:,2,:].flatten()
            #self.cross_sect_x.extend(cross_sect_cyl2[:,0,:].flatten())
            #self.cross_sect_y.extend(cross_sect_cyl2[:,1,:].flatten())
            #self.cross_sect_z.extend(cross_sect_cyl2[:,2,:].flatten())
        #self.create_triangular_mesh()
        #self.vertices = np.array([self.cross_sect_x, self.cross_sect_y, self.cross_sect_z]).T

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
    def __init__(self,):
        pass

    def input_data(self, CCD_position, CCD_x, CCD_y, CCD_focal_distance, CCD_pixels_x, CCD_pixels_y, u_hat, patch, v_hat = None, w_hat = None, CCD_focal_point = None):
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
        self.LOS_point = np.zeros((self.CCD_pixels_y, self.CCD_pixels_x, 3),dtype=float)
        self.LOS_gradient = self.LOS_point * 0
        for pixel_y in range(self.CCD_pixels_y):
            for pixel_x in range(self.CCD_pixels_x):
                #print pixel_y,pixel_x
                self.LOS_point[pixel_y, pixel_x, :] = (float(pixel_x)/self.CCD_pixels_x - 0.5) * self.CCD_width * self.v_hat + (-float(pixel_y)/self.CCD_pixels_y + 0.5) * self.CCD_length * self.u_hat + self.CCD_position
                self.LOS_gradient[pixel_y, pixel_x, :] = self.focal_point - self.LOS_point[pixel_y, pixel_x, :]

    def find_intersections(self,):
        '''Find the intersection for each triangle plane
        SRH: 10July2013
        '''
        print('finding intersections - V2')
        self.valid_channels = np.zeros((self.CCD_pixels_y,self.CCD_pixels_x),dtype=bool)
        self.intersection1 = self.LOS_point * 0
        self.intersection2 = self.LOS_point * 0
        for pixel_y in range(self.CCD_pixels_y):
            for pixel_x in range(self.CCD_pixels_x):
                P0 = self.LOS_point[pixel_y, pixel_x, :]
                dP = self.LOS_gradient[pixel_y, pixel_x, :]
                t, R, tmp1 = self.patch.find_intersections(P0, dP)
                #t = (self.patch.d - np.sum(self.patch.n * P0,axis=1))/np.sum(self.patch.n*dP,axis=1)
                #R = P0 + t[:,np.newaxis]*dP[np.newaxis,:]
                #tmp1 = np.sum(self.patch.n*np.cross(self.patch.pt1,R - self.patch.tris[:,0,:]),axis=1)>=0
                #tmp1[tmp1] = np.sum(self.patch.n[tmp1]*np.cross(self.patch.pt2[tmp1,:],R[tmp1,:] - self.patch.tris[tmp1,1,:]),axis=1)>=0
                #tmp1[tmp1] = np.sum(self.patch.n[tmp1]*np.cross(self.patch.pt3[tmp1,:],R[tmp1,:] - self.patch.tris[tmp1,2,:]),axis=1)>=0
                #a, intersect_pts = find_intersections2(P0,dP,d,n,pt1,pt2,pt3,tris)
                if (np.sum(tmp1)) >= 2:# and ((np.sum(tmp1)%2)==0):
                    #mlab.points3d(intersect_pts[:,0], intersect_pts[:,1], intersect_pts[:,2],scale_factor=0.02)
                    intersection_order = np.argsort(t[tmp1])
                    self.valid_channels[pixel_y,pixel_x] = True
                    self.intersection1[pixel_y, pixel_x] = R[tmp1][intersection_order[0]]
                    self.intersection2[pixel_y, pixel_x] = R[tmp1][intersection_order[1]]
                if np.sum(tmp1)%2 == 1:
                    print('Warning odd number of intersections!!!, {}'.format(np.sum(tmp1)))
        return tmp1, R[tmp1]


    def find_intersections_TFC(self,tfc_kwargs, plot_first_intersection = 0, plot_tfc = 0):
        '''Find the intersection of the LOS with TFC's and remove the
        channels where this occurs tfc_coil_number : which TFC...

        SRH: 10July2013
        '''
        if tfc_kwargs==None:tfc_kwargs = {'coil_num':tfc_coil_number}
        print('finding TFC intersections - V2')
        TFC_intersection = np.zeros((self.CCD_pixels_y,self.CCD_pixels_x),dtype=bool)
        TFC_patch = TFC(**tfc_kwargs)
        TFC_intersection1 = self.LOS_point * 0
        TFC_intersection2 = self.LOS_point * 0
        for pixel_y in range(self.CCD_pixels_y):
            for pixel_x in range(self.CCD_pixels_x):
                P0 = self.LOS_point[pixel_y, pixel_x, :]
                dP = self.LOS_gradient[pixel_y, pixel_x, :]
                t, R, tmp1 = TFC_patch.find_intersections(P0, dP)
                if (np.sum(tmp1)) >= 2:# and ((np.sum(tmp1)%2)==0):
                    #mlab.points3d(intersect_pts[:,0], intersect_pts[:,1], intersect_pts[:,2],scale_factor=0.02)
                    intersection_order = np.argsort(t[tmp1])
                    dist_tfc = np.sum((P0 - R[tmp1][intersection_order[0]])**2)
                    dist_plas = np.sum((P0 - self.intersection1[pixel_y,pixel_x,:])**2)
                    if dist_tfc< dist_plas:
                        TFC_intersection[pixel_y,pixel_x] = True
                        TFC_intersection1[pixel_y, pixel_x] = R[tmp1][intersection_order[0]]
                        TFC_intersection2[pixel_y, pixel_x] = R[tmp1][intersection_order[1]]
                if np.sum(tmp1)%2 == 1:
                    print('Warning odd number of intersections!!!, {}'.format(np.sum(tmp1)))
        self.valid_channels *= np.invert(TFC_intersection)
        if plot_first_intersection:
            import mayavi.mlab as mlab
            mlab.points3d(TFC_intersection1[TFC_intersection,0], TFC_intersection1[TFC_intersection,1], TFC_intersection1[TFC_intersection,2],scale_factor=0.02, color=(0,0,1))
        if plot_tfc:
            TFC_patch.plot()
        return tmp1, R[tmp1]

    def create_interpolation_points(self, n_interp_pts):
        self.n_interp_pts = n_interp_pts
        self.interpolation_pts = np.zeros((self.valid_channels.shape[0], self.valid_channels.shape[1],n_interp_pts,3),dtype=float)
        self.dl = np.zeros((self.valid_channels.shape[0], self.valid_channels.shape[1]),dtype=float)
        for pixel_y in range(self.CCD_pixels_y):
            for pixel_x in range(self.CCD_pixels_x):
                if self.valid_channels[pixel_y, pixel_x]:
                    self.interpolation_pts[pixel_y, pixel_x,:,0] = np.linspace(self.intersection1[pixel_y,pixel_x,0],self.intersection2[pixel_y, pixel_x,0],n_interp_pts)
                    self.interpolation_pts[pixel_y, pixel_x,:,1] = np.linspace(self.intersection1[pixel_y,pixel_x,1],self.intersection2[pixel_y, pixel_x,1],n_interp_pts)
                    self.interpolation_pts[pixel_y, pixel_x,:,2] = np.linspace(self.intersection1[pixel_y,pixel_x,2],self.intersection2[pixel_y, pixel_x,2],n_interp_pts)
                    self.dl[pixel_y, pixel_x] = np.sqrt(np.sum((self.interpolation_pts[pixel_y, pixel_x,1,:] - self.interpolation_pts[pixel_y, pixel_x,0,:])**2))


    def perform_interpolation(self,no_theta = 50, s_increment = 2):
        ''' 
        valid_channels : bool array listing the channels that intersect the plasma
        cross_sect_s : 1D array of s values (corresponding to the interpolation_pts
        cross_sect_theta : 1D array of theta_b values (corresponding to the interpolation_pts
        cross_sect_phi : 1D array of phi_b values (corresponding to the interpolation_pts

        interpolation_pts : grid of points in realspace [channel_num, interp_point, x/y/z]
        SRH: 10July2013
        '''
        if self.patch.obtained_grid != 1:
            print 'getting patch grid'
            self.patch.get_grid_for_interpolation(no_theta = no_theta, s_increment = s_increment)
        points_tuple = (self.patch.grid_x, self.patch.grid_y, self.patch.grid_z)
        cross_sect_s = self.patch.grid_s
        cross_sect_phi = self.patch.grid_phi
        cross_sect_theta = self.patch.grid_theta

        valid_channels = self.valid_channels
        interp_pts_x = self.interpolation_pts[self.valid_channels,:,0].flatten()
        interp_pts_y = self.interpolation_pts[self.valid_channels,:,1].flatten()
        interp_pts_z = self.interpolation_pts[self.valid_channels,:,2].flatten()
        required_shape = self.interpolation_pts[self.valid_channels,:,0].shape
        self.interp_boozer = self.interpolation_pts*0
        print 's_pt1 - '
        #self.interp_boozer[self.valid_channels,:,0] = scipy_griddata(points_tuple, np.array(cross_sect_s), (interp_pts_x, interp_pts_y, interp_pts_z)).reshape(required_shape)
        self.s_pt1 = scipy_griddata(points_tuple, np.array(cross_sect_s), (interp_pts_x, interp_pts_y, interp_pts_z)).reshape(required_shape)

        print('before2 {}'.format(np.sum(np.isnan(self.s_pt1))))
        #tmp = 

        self.start_pt = np.argmax(np.isfinite(self.s_pt1),axis = 1)
        self.end_pt = np.argmax(np.isfinite(self.s_pt1)[:,::-1],axis = 1)
        #tmp1 = self.interpolation_pts[self.valid_channels,:,:]
        #tmp2 = self.intersection1[self.valid_channels,:]
        #tmp3 = self.intersection2[self.valid_channels,:]
        self.intersection1_old = +self.intersection1
        self.intersection2_old = +self.intersection2
        self.interpolation_pts_old = +self.interpolation_pts
        a = +self.intersection1[self.valid_channels,:]
        b = +self.intersection2[self.valid_channels,:]
        c = +self.interpolation_pts[self.valid_channels,:,:]
        d = +self.valid_channels[self.valid_channels]
        self.start_list = []
        self.end_list = []
        for i in range(len(self.start_pt)):
            if float(np.sum(np.isnan(self.s_pt1[i,:])))/self.s_pt1.shape[1] > 0.4:
                d[i]= False
                print('bad channel - nan content above 0.4')
            else:
                itemindex,=np.where(np.isfinite(self.s_pt1[i,:]))
                start_pt = np.min(itemindex)
                end_pt = np.max(itemindex)
                self.start_list.append(start_pt)
                self.end_list.append(end_pt)
                a[i,:] = (c[i, start_pt,:])*1.
                b[i,:] = (c[i, end_pt,:])*1.
            #self.intersection1[self.valid_channels,:][i,:] = self.interpolation_pts[self.valid_channels,:,:][i, self.start_pt[i],:]*1.
            #self.intersection2[self.valid_channels,:][i,:] = self.interpolation_pts[self.valid_channels,:,:][i, -self.end_pt[i],:]*1.
        self.intersection1[self.valid_channels,:] = +a
        self.intersection2[self.valid_channels,:] = +b
        self.valid_channels[self.valid_channels] = +d
        #self.intersection1[self.valid_channels,:] = self.interpolation_pts[self.valid_channels,:,:][:,self.start_pt,:]
        #self.intersection2[self.valid_channels,:] = self.interpolation_pts[self.valid_channels,:,:][:,self.end_pt,:]

        
        self.create_interpolation_points(self.n_interp_pts)

        interp_pts_x = self.interpolation_pts[self.valid_channels,:,0].flatten()
        interp_pts_y = self.interpolation_pts[self.valid_channels,:,1].flatten()
        interp_pts_z = self.interpolation_pts[self.valid_channels,:,2].flatten()
        required_shape = self.interpolation_pts[self.valid_channels,:,0].shape

        print 's_pt2'
        self.interp_boozer[self.valid_channels,:,0] = scipy_griddata(points_tuple, np.array(cross_sect_s), (interp_pts_x, interp_pts_y, interp_pts_z)).reshape(required_shape)
        print('after {}'.format(np.sum(np.isnan(self.interp_boozer[self.valid_channels,:,0]))))
        print 'theta'
        interp_data_theta_sin = scipy_griddata(points_tuple, np.sin(cross_sect_theta), (interp_pts_x, interp_pts_y, interp_pts_z))
        interp_data_theta_cos = scipy_griddata(points_tuple, np.cos(cross_sect_theta), (interp_pts_x, interp_pts_y, interp_pts_z))
        self.interp_boozer[self.valid_channels,:,1] = np.arctan2(interp_data_theta_sin, interp_data_theta_cos).reshape(required_shape)
        print('after {}'.format(np.sum(np.isnan(self.interp_boozer[self.valid_channels,:,1]))))
        print 'phi'
        interp_data_phi_sin = scipy_griddata(points_tuple, np.sin(cross_sect_phi), (interp_pts_x, interp_pts_y, interp_pts_z))
        interp_data_phi_cos = scipy_griddata(points_tuple, np.cos(cross_sect_phi), (interp_pts_x, interp_pts_y, interp_pts_z))
        self.interp_boozer[self.valid_channels,:,2] = np.arctan2(interp_data_phi_sin, interp_data_phi_cos).reshape(required_shape)
        print('after {}'.format(np.sum(np.isnan(self.interp_boozer[self.valid_channels,:,2]))))



    def grid_for_wave(self, phi_value=120., z_min = -0.4, z_max = 0.4, r_min = 0.9, r_max = 1.45, n_pts = 100, plot_mask = 0):
        ''' 

        SRH: 22July2013
        '''
        points_tuple = (self.patch.grid_x, self.patch.grid_y, self.patch.grid_z)
        cross_sect_s = self.patch.grid_s
        cross_sect_phi = self.patch.grid_phi
        cross_sect_theta = self.patch.grid_theta
        
        z_values = np.linspace(z_min,z_max, n_pts)
        r_values = np.linspace(r_min, r_max, n_pts)
        self.r_grid, self.z_grid = np.meshgrid(r_values, z_values)
        self.x_grid = self.r_grid * np.cos(np.deg2rad(phi_value))
        self.y_grid = self.r_grid * np.sin(np.deg2rad(phi_value))
        self.grid_mask = np.ones(self.y_grid.shape,dtype=bool)
        self.edge_points_rz = []; self.edge_points_xyz = []
        for i in range(self.r_grid.shape[0]):
            dP = np.array([self.x_grid[i,1]-self.x_grid[i,0], self.y_grid[i,1]-self.y_grid[i,0], self.z_grid[i,1]-self.z_grid[i,0]])
            dP = -dP/np.sqrt(np.sum(dP**2))
            P0 = np.array([self.x_grid[i, 0], self.y_grid[i, 0], self.z_grid[i,0]])
            t, R, tmp1 = self.patch.find_intersections(P0, dP)
            if (np.sum(tmp1)) >= 2:# and ((np.sum(tmp1)%2)==0):
                #mlab.points3d(intersect_pts[:,0], intersect_pts[:,1], intersect_pts[:,2],scale_factor=0.02)
                intersection_order = np.argsort(t[tmp1])
                int1_point = R[tmp1][intersection_order[0]]
                int2_point = R[tmp1][intersection_order[1]]
                tmp_points = np.array([self.x_grid[i,:],self.y_grid[i,:], self.z_grid[i,:]]).T
                int1 = np.argmin(np.sum((tmp_points - int1_point)**2,axis=1))
                int2 = np.argmin(np.sum((tmp_points - int2_point)**2,axis=1))
                self.edge_points_rz.append([self.r_grid[i,int1], self.z_grid[i,int1]])
                self.edge_points_rz.append([self.r_grid[i,int2], self.z_grid[i,int2]])
                self.edge_points_xyz.append([self.x_grid[i,int1], self.y_grid[i,int1], self.z_grid[i,int1]])
                self.edge_points_xyz.append([self.x_grid[i,int2], self.y_grid[i,int2], self.z_grid[i,int2]])
                if int1>int2:
                    self.grid_mask[i, int2:int1] = False
                else:
                    self.grid_mask[i, int1:int2] = False
            if np.sum(tmp1)%2 == 1:
                print('Warning odd number of intersections!!!, {}'.format(np.sum(tmp1)))
        self.edge_points_rz = np.array(self.edge_points_rz)
        self.edge_points_xyz = np.array(self.edge_points_xyz)
        
        if plot_mask:
            fig,ax = pt.subplots()
            ax.imshow(np.ma.array(self.grid_mask, mask=self.grid_mask), extent=[r_min,r_max,z_min,z_max], origin='lower')
            ax.plot(self.edge_points_rz[:,0], self.edge_points_rz[:,1], 'k.')
            fig.canvas.draw(); fig.show()

        self.valid_pts = self.grid_mask==False
        interp_pts_x = self.x_grid[self.valid_pts].flatten()
        interp_pts_y = self.y_grid[self.valid_pts].flatten()
        interp_pts_z = self.z_grid[self.valid_pts].flatten()
        print interp_pts_x.shape, n_pts**2
        required_shape = self.z_grid[self.valid_pts].shape
        self.s_cross_sect = self.x_grid *0
        self.theta_cross_sect = self.x_grid *0
        self.phi_cross_sect = self.x_grid *0
        print 's'
        #self.interp_boozer[self.valid_channels,:,0] = scipy_griddata(points_tuple, np.array(cross_sect_s), (interp_pts_x, interp_pts_y, interp_pts_z)).reshape(required_shape)
        self.s_cross_sect[self.valid_pts] = scipy_griddata(points_tuple, np.array(cross_sect_s), (interp_pts_x, interp_pts_y, interp_pts_z)).reshape(required_shape)
        print('theta')
        interp_data_theta_sin = scipy_griddata(points_tuple, np.sin(cross_sect_theta), (interp_pts_x, interp_pts_y, interp_pts_z))
        interp_data_theta_cos = scipy_griddata(points_tuple, np.cos(cross_sect_theta), (interp_pts_x, interp_pts_y, interp_pts_z))
        self.theta_cross_sect[self.valid_pts] = np.arctan2(interp_data_theta_sin, interp_data_theta_cos).reshape(required_shape)
        print('phi')
        interp_data_phi_sin = scipy_griddata(points_tuple, np.sin(cross_sect_phi), (interp_pts_x, interp_pts_y, interp_pts_z))
        interp_data_phi_cos = scipy_griddata(points_tuple, np.cos(cross_sect_phi), (interp_pts_x, interp_pts_y, interp_pts_z))
        self.phi_cross_sect[self.valid_pts] = np.arctan2(interp_data_phi_sin, interp_data_phi_cos).reshape(required_shape)

    def plot_boozer_LOS(self,y_vals = None, x_vals = None):
        mask = np.zeros(self.valid_channels.shape,dtype=bool)
        if y_vals == None and x_vals == None:
            mask = self.valid_channels
        elif y_vals != None and x_vals != None:
            mask[y_vals,x_vals] = True
        elif y_vals != None: 
            mask[y_vals,:] = True
        elif x_vals != None:
            mask[:,x_vals] = True
        else:
            raise ValueError('some problem! y_vals and x_vals must an integer or None')
        mask *= self.valid_channels
        fig, ax = pt.subplots(ncols = 3)
        ax[0].plot(self.interp_boozer[mask,:,0].T)
        ax[1].plot(self.interp_boozer[mask,:,1].T)
        ax[2].plot(self.interp_boozer[mask,:,2].T)
        fig.canvas.draw(); fig.show()


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
        count = 0
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
            if count in [0,1, np.prod(self.LOS_point.shape[0:2])-1]:#,np.prod(self.LOS_point.shape[0:2])):
                mlab.plot3d([x[-1],x[-2]],[y[-1],y[-2]],[z[-1],z[-2]],color=(0,0,0),line_width=1.5, tube_radius = None)
            count+=1
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




def imax_camera(boozer_filename = '/home/srh112/code/python/h1_eq_generation/results7/kh0.350-kv1.000fixed/boozmn_wout_kh0.350-kv1.000fixed.nc', plot_LOS = 0, plot_patch = 0, plot_intersections = 0, plot_pfc = 0, plot_tfc = 0,phi_range = 30, n_phi = 30, decimate_pixel=16, measurements = None, no_theta = 50, patch_pickle = None, elevation_angle = 0., elevation = 0., patch_object = None, n_pixels_x = 512, n_pixels_y=512, CCD_L=0.01575):
    '''Convenience function containing the required geometry for the
    imax camera

    SRH: 12July2013
    '''
    # LOS(CCD, CCD_width, CCD_length)
    u_hat = np.array([0.,0.,1.])
    CCD_focal_distance = 17.0/1000.
    if measurements == 'nandi_measurements':
        print 'nandi measurements'
        CCD_position = np.array([-1.01355,1.66046,0.04998])
        CCD_phi = 121.40
        w_hat = np.array([0.52084,-0.85327,0.025684])
        v_hat = np.cross(w_hat,u_hat)
    elif measurements == 'nandi_tomo17':
        print 'nandi tomo17'
        CCD_position = np.array([-1.00239,1.67301,0.01501])
        CCD_phi = 120.93
        w_hat = np.array([0.56439,-0.82468,0.037002])
        v_hat = np.cross(w_hat,u_hat)
    elif measurements == 'nandi_tomo19':
        print 'nandi tomo19'
        CCD_position = np.array([-0.982017,1.66487,0.01506])
        CCD_phi = 120.534
        w_hat = np.array([0.54081,-0.84041,0.035282])
        v_hat = np.cross(w_hat,u_hat)
    elif measurements == 'nandi_tomo20':
        print 'nandi tomo20'
        CCD_phi = 120.842
        CCD_position = np.array([-0.98496,1.64951,0.01500])
        w_hat = np.array([0.56164,-0.82642,0.03988])
        v_hat = np.cross(w_hat,u_hat)
    elif measurements == 'nandi_combined_1':
        print 'nandi combined 1'
        CCD_position = np.array([-0.98331,1.66879,0.01500])
        CCD_phi = 120.508
        w_hat = np.array([0.54205,-0.83946,0.03846])
        v_hat = np.cross(w_hat,u_hat)

    elif measurements=='pimax4':
        print 'pimax4'
        CCD_phi = 120.; CCD_R = 1.946
        CCD_R = 1.1475+0.963+0.0175 
        CCD_z = 7.0/100; 
        CCD_z = 8.0/100;
        CCD_x = CCD_R*np.cos(np.deg2rad(CCD_phi)); CCD_y = CCD_R*np.sin(np.deg2rad(CCD_phi))
        CCD_position = np.array([CCD_x, CCD_y, CCD_z])

        l1 = 1./np.cos(np.deg2rad(elevation_angle))
        l2 = np.tan(np.deg2rad(elevation_angle))
        print np.sqrt(l1**2 - l2**2)
        CCD_position = CCD_position + np.array([0,0,elevation])
        z_component = -np.sin(np.deg2rad(elevation_angle))
        w_hat = -np.array([CCD_x, CCD_y, 0])
        w_hat = w_hat / np.sqrt(np.sum(w_hat**2))
        w_hat = w_hat * np.cos(np.deg2rad(elevation_angle)) + np.array([0,0,-np.sin(np.deg2rad(elevation_angle))])
        print 'w_hat', np.sum(w_hat**2)
        #w_hat = -CCD_position/(np.sqrt(np.sum(CCD_position**2)))
        #print np.sum(w_hat**2)
        #print np.sqrt(np.sum((w_hat * l2)**2)), l2
        u_hat = w_hat * l2 + np.array([0.,0.,l1])
        print 'u hat', u_hat, np.sqrt(np.sum(u_hat**2))
        #1/0
        v_hat2 = np.cross(w_hat, u_hat)
        #v_hat2 = v_hat2 / (np.sqrt(np.sum(v_hat**2)))
        v_hat = np.array([-np.cos(np.deg2rad(90.-CCD_phi)), np.sin(np.deg2rad(90.-CCD_phi)),0])
        print 'v comparisons', v_hat2, v_hat, np.sqrt(np.sum(v_hat2**2)), np.sqrt(np.sum(v_hat**2))

        #CCD_focal_distance = 13/1000.
        CCD_focal_distance = 17.3/1000.

    else:
        print 'default'
        CCD_phi = 119.5; CCD_R = 1.946
        CCD_z = 5./100; CCD_x = CCD_R*np.cos(np.deg2rad(CCD_phi)); CCD_y = CCD_R*np.sin(np.deg2rad(CCD_phi))
        CCD_position = np.array([CCD_x, CCD_y, CCD_z])

        l1 = 1./np.cos(np.deg2rad(elevation_angle))
        l2 = np.tan(np.deg2rad(elevation_angle))
        print np.sqrt(l1**2 - l2**2)
        CCD_position = CCD_position + np.array([0,0,elevation])
        z_component = -np.sin(np.deg2rad(elevation_angle))
        w_hat = -np.array([CCD_x, CCD_y, 0])
        w_hat = w_hat / np.sqrt(np.sum(w_hat**2))
        w_hat = w_hat * np.cos(np.deg2rad(elevation_angle)) + np.array([0,0,-np.sin(np.deg2rad(elevation_angle))])
        print 'w_hat', np.sum(w_hat**2)
        #w_hat = -CCD_position/(np.sqrt(np.sum(CCD_position**2)))
        #print np.sum(w_hat**2)
        #print np.sqrt(np.sum((w_hat * l2)**2)), l2
        u_hat = w_hat * l2 + np.array([0.,0.,l1])
        print 'u hat', u_hat, np.sqrt(np.sum(u_hat**2))
        #1/0
        v_hat2 = np.cross(w_hat, u_hat)
        #v_hat2 = v_hat2 / (np.sqrt(np.sum(v_hat**2)))
        v_hat = np.array([-np.cos(np.deg2rad(90.-CCD_phi)), np.sin(np.deg2rad(90.-CCD_phi)),0])
        print 'v comparisons', v_hat2, v_hat, np.sqrt(np.sum(v_hat2**2)), np.sqrt(np.sum(v_hat**2))

        CCD_focal_distance = 17.3/1000.

    # #Nandi tomo17
    # CCD = np.array([-1.00239,1.673,0.0150])
    # w_hat = np.array([0.564,-0.82468,0.037002])
    # v_hat = np.cross(w_hat,u_hat)
    # f = 17./1000
    #phi_min =60; phi_max = 200;n_phi = 30;no_theta = 50

    CCD_x = CCD_L; CCD_y = CCD_L
    CCD_pixels_x = n_pixels_x/decimate_pixel; CCD_pixels_y = n_pixels_y/decimate_pixel
    phi_min =CCD_phi - phi_range; phi_max = CCD_phi + phi_range;no_theta = no_theta
    plot_length = 50
    if patch_object != None:
        print('using the patch object...')
        patch = patch_object
    elif patch_pickle == None:
        print('Generating patch surface...')
        patch =  BoozerSurfacePatch(phi_min, phi_max, n_phi = n_phi, no_theta = no_theta, boozer_filename = boozer_filename)
    else:
        print('Loading patch surface from pickle file...')
        patch = pickle.load(file(patch_pickle,'r'))

    print CCD_x, CCD_y, CCD_focal_distance, CCD_pixels_x, CCD_pixels_y


    answer = LOS()
    answer.input_data(CCD_position, CCD_x, CCD_y, CCD_focal_distance, CCD_pixels_x, CCD_pixels_y, u_hat, patch, v_hat = v_hat, w_hat = None, CCD_focal_point = None)

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
    phi_min =CCD_phi - 30; phi_max = CCD_phi+30;n_phi = 15;no_theta = 35;n_interp_pts = 100

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


def plot_imax_mask(LOS_object, shot_number = 71235, cal_file=None, image_array=None, cut_values = None):
    if cut_values==None:
        cut_values=[32]
    fig, ax = pt.subplots(ncols = 4)
    #fig, ax = pt.subplots(ncols = 4, sharex = True , sharey = True)
    if image_array != None:
        pass
    elif cal_file == None:
        if shot_number.__class__!=list: shot_number = [shot_number]
        for i, shot in enumerate(shot_number):
            image_array_tmp = MDS.Tree('imax',shot_number).getNode('.images').data()[0,:,:]
            if i == 0:
                len1, len2 = image_array_tmp.shape
                image_array = np.zeros(len1*len(shot_number),len2,dtype=float)
            image_array[len1*i:len1*(i+1),:] = +image_array_tmp
    else:
        if cal_file.__class__!=list: cal_file = [cal_file]
        for i, cal in enumerate(cal_file):
            tmp2 = imax.ImaxData(SPE_file= cal, plot_data = 0, plot_mirnov_triggers = 0,get_mirnov_triggers = 0)
            tmp2.calibrate_data(dark_shot=1103, white_shot=1107, plot_calibration_im = 0, clip = 0.2, plot_cal_images = 0, cal_sum = 1, med_filt = 3, mode_amp_filt = None)
            image_array_tmp = tmp2.image_array_cal[0,:,:]
            if i == 0:
                print 'hello'
                len1, len2 = image_array_tmp.shape
                image_array = np.zeros((len1*len(cal_file),len2),dtype=float)
            image_array[len1*i:len1*(i+1),:] = +image_array_tmp
    extent = image_array.shape
    #extent = [0, extent[0], 0, extent[1]]
    extent = [0, extent[1], 0, extent[0]]

    img3 = ax[2].imshow(image_array, interpolation='nearest', origin='upper', alpha=1.0,extent=extent,aspect='auto')
    #img3.set_clim([0,1.5])
    img3 = ax[0].imshow(image_array, interpolation='nearest', origin='upper', alpha=1.0,extent=extent,aspect='auto')
    #img3.set_clim([0,1.5])
    #ax.imshow(np.fliplr(np.flipud(valid)), origin='lower',alpha=0.6,extent=[0,512,0,512])
    #ax[1].imshow(np.flipud(np.fliplr(LOS_object.valid_channels)),origin='lower',alpha=0.6,extent=[0,512,0,512])
    #ax[2].imshow(np.flipud(np.fliplr(LOS_object.valid_channels)),origin='lower',alpha=0.6,extent=[0,512,0,512])
    ax[1].imshow(LOS_object.valid_channels,origin='upper',alpha=0.3,extent=extent,aspect='auto', cmap = pt.cm.bone)
    ax[2].imshow(LOS_object.valid_channels,origin='upper',alpha=0.3,extent=extent,aspect='auto', cmap = pt.cm.bone)
    tmp_y_axis = np.arange(len(LOS_object.valid_channels[:,0]))
    tmp_y_axis = np.max(tmp_y_axis) - tmp_y_axis
    for i in cut_values:
        ax[3].plot(LOS_object.valid_channels[:,i], tmp_y_axis)
        ax[3].plot(image_array[:,i]/np.max(image_array[:,i]), tmp_y_axis)
        
    #ax[1].imshow(LOS_object.valid_channels, origin='lower',alpha=0.6,extent=[0,512,0,512])
    #ax[2].imshow(LOS_object.valid_channels, origin='lower',alpha=0.6,extent=[0,512,0,512])
    #img3 = plt.imshow(zvals2, interpolation='nearest', cmap=cmap2, origin='lower', alpha=0.6)
    for i in range(3):
        ax[i].set_xlim([extent[0], extent[1]])
        ax[i].set_ylim([extent[2], extent[3]])
        ax[i].vlines(cut_values,extent[2], extent[3])
    ax[3].set_xlim([0,2])
    ax[3].set_ylim([extent[2], extent[3]])
    fig.canvas.draw(); fig.show()
