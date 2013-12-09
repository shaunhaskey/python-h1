import numpy as np
import mayavi.mlab as mlab
import mayavi

def _plot_coil_horizontal(centre, thickness, width, r, phi_samples=50,**kw):
    '''Plot a horizontal coil like hte PFC or OVC
    '''
    rmin = r - width/2
    rmax = r + width/2
    phi = np.linspace(0,2.*np.pi,phi_samples)
    heights = [-thickness/2,thickness/2,thickness/2,-thickness/2,-thickness/2]
    r_values = [rmax,rmax,rmin,rmin,rmax]
    circle_x = np.zeros((len(heights),len(phi)),dtype = float)
    circle_y = np.zeros((len(heights),len(phi)),dtype = float)
    circle_z = np.zeros((len(heights),len(phi)),dtype = float)
    for i in range(0,len(heights)):
        circle_x[i,:] = r_values[i]*np.sin(phi)+centre[0]
        circle_y[i,:] = r_values[i]*np.cos(phi)+centre[1]
        circle_z[i,:] = np.cos(phi)*0+heights[i]+centre[2]
    mlab.mesh(circle_x, circle_y, circle_z,**kw)#color=color,opacity=opacity)

def _plot_coil_vertical(center, radius, thickness, width, phi,  phi_samples = 50, ID_OD_offset=0, **kw):
    '''Plot a horizontal coil like hte PFC or OVC
    '''
    print('ID_OD_offset plot_tfc: {}'.format(ID_OD_offset))
    circle_x1, circle_y1, circle_z1 = _coil_vertical_points(center, radius, thickness, width, phi, phi_samples = phi_samples, ID_OD_offset = ID_OD_offset)
    mlab.mesh(circle_x1, circle_y1, circle_z1,**kw)#color=(1,0,0),opacity=opacity)

def _coil_vertical_points(center, radius, thickness, width, phi, phi_samples = 50, ID_OD_offset = 0):
    '''Plot a horizontal coil like hte PFC or OVC
    '''
    print('ID_OD_offset plot_tfc: {}'.format(ID_OD_offset))
    print('center {}, radius {}, thickness {}, width {}, phi {}, phi_samples {}'.format(center, radius, thickness, width, phi, phi_samples))
    rmin = radius - width/2
    rmax = radius + width/2
    print('rmin:{}, rmax:{}'.format(rmin,rmax))
    heights = [-thickness/2,thickness/2,thickness/2,-thickness/2,-thickness/2]
    r_values = [rmax,rmax,rmin,rmin,rmax]
    phi = phi/180.*np.pi
    centre_x = center[0]
    centre_y = center[1]
    centre_z = center[2]
    centre_z_values = np.array([-ID_OD_offset/2, -ID_OD_offset/2, +ID_OD_offset/2, +ID_OD_offset/2, -ID_OD_offset/2])+center[2]
    phi_tmp = np.linspace(0,2.*np.pi,phi_samples)
    circle_x1 = np.zeros((len(heights),len(phi_tmp)),dtype = float)
    circle_y1 = np.zeros((len(heights),len(phi_tmp)),dtype = float)
    circle_z1 = np.zeros((len(heights),len(phi_tmp)),dtype = float)
    for i in range(0,len(heights)):
        circle_r = r_values[i]*np.cos(phi_tmp)
        circle_z = r_values[i]*np.sin(phi_tmp)#+centre_z
        circle_x = circle_r*np.cos(phi)#+centre_x
        circle_y = circle_r*np.sin(phi)#+centre_y
        circle_x1[i,:] = circle_x + heights[i]*np.cos(phi+np.pi/2.) #r_values[i]*np.sin(phi)+centre[0]
        circle_y1[i,:] = circle_y + heights[i]*np.sin(phi+np.pi/2.)#r_values[i]*np.cos(phi)+centre[1]
        print('Using r:{}, Using z:{}, ID_OD_offset:{}'.format(r_values[i], centre_z_values[i], ID_OD_offset))
        circle_z1[i,:] = circle_z + centre_z_values[i] #np.cos(phi)*0+heights[i]+centre[2]
    circle_x1 = circle_x1 + centre_x
    circle_y1 = circle_y1 + centre_y
    #circle_z1 = circle_z1 + centre_z
    return circle_x1, circle_y1, circle_z1


def plot_pfc(pfc_thickness = 0.11, pfc_width = 0.11, pfc_radius = 1.0, pfc_mesh_props=None):
    '''plot the pfc

    SRH: 12July2013
    '''

    pfc_centre = [0,0,0]
    if pfc_mesh_props == None:
        pfc_mesh_props = {'opacity':1.,'color':(0.5,0.5,0.5)}
    _plot_coil_horizontal(pfc_centre, pfc_thickness, pfc_width, pfc_radius, **pfc_mesh_props)

def plot_ovc(ovc_mesh_props=None):
    '''plot the ovc

    SRH: 12July2013
    '''
    ovc_centre1 = [0,0,0.73]; ovc_centre2 = [0,0,-0.73]; ovc_thickness=0.1; ovc_width=0.2; ovc_radius = 2.13
    if ovc_mesh_props == None:
        ovc_mesh_props = {'opacity':1,'color':(255/255., 255/255., 51/255.)}
    _plot_coil_horizontal(ovc_centre1, ovc_thickness, ovc_width, ovc_radius, **ovc_mesh_props)
    _plot_coil_horizontal(ovc_centre2, ovc_thickness, ovc_width, ovc_radius, **ovc_mesh_props)

def plot_ivc(ivc_mesh_props=None):
    '''plot the ivc

    SRH: 12July2013
    '''
    ivc_centre1=[0,0,-1.07]; ivc_centre2=[0,0,1.07]; ivc_thickness = 0.1; ivc_width=0.2; ivc_radius = 0.72
    if ivc_mesh_props == None:
        ivc_mesh_props = {'opacity':1,'color':(0.5,0.5,0.5)}
    _plot_coil_horizontal(ivc_centre1, ivc_thickness, ivc_width, ivc_radius, **ivc_mesh_props)
    _plot_coil_horizontal(ivc_centre2, ivc_thickness, ivc_width, ivc_radius, **ivc_mesh_props)



def plot_tfc(include_coils, tfc_thickness=0.075, tfc_width=0.15, tfc_radius = 0.383, ID_OD_offset = 0, tfc_mesh_props = None, use_real = False):
    '''plot the tfc coils in the include_coils list

    SRH: 12July2013
    '''
    #TFC details
    print('ID_OD_offset plot_tfc: {}'.format(ID_OD_offset))
    if use_real:
        tfc = tfc_details_real()
    else:
        tfc = tfc_details()
    #plot the TFC's 
    #tfc_mesh_props = {'opacity':0.3,'color':(1,0.,0.)}
    if tfc_mesh_props==None:
        tfc_mesh_props = {'opacity':1.0,'color':(0.5,0.5,0.5)}
    #include_coils = [11,12,13]#include_coils.extend(range(1,15,4))
    #include_coils.extend(range(26,36,4))
    #for i in range(1,36,4):
    for i in include_coils:
         _plot_coil_vertical([tfc.x[i],tfc.y[i],tfc.z[i]], tfc_radius, tfc_thickness, tfc_width, np.arctan2(tfc.y[i],tfc.x[i])*180./np.pi, ID_OD_offset=ID_OD_offset, **tfc_mesh_props)

class tfc_details():
    def __init__(self,):
        self.e2 = [-86.482994, -79.050995, -70.533997, -60.533997, -48.550999,
        -36.482998, -23.517000, -10.948996, 0.533993, 10.533997,
        19.050997, 25.742990, 33.516987, 40.449001, 49.466000, 59.465992,
        71.448997, 83.516991, 83.516991, 70.949013, 59.466000, 49.466000,
        40.949005, 33.804981, 26.482990, 19.051008, 10.534008, 0.534011,
        -11.449009, -23.576992, -36.482990, -49.051003, -60.533997,
        -70.533997, -79.050995, -86.483002];

        self.x = [1.210616, 1.159895, 1.050229, 0.877091, 0.666271, 0.475905,
        0.319576, 0.167694, -0.009389, -0.203642, -0.386273, -0.525719,
        -0.669468, -0.769832, -0.845447, -0.869596, -0.841766, -0.795282,
        -0.798461, -0.835960, -0.865548, -0.846207, -0.773947, -0.674262,
        -0.540737, -0.385621, -0.203642, -0.009391, 0.176344, 0.320145,
        0.475311, 0.669492, 0.876394, 1.054943, 1.160091, 1.21115];

        self.y = [0.074405, 0.224389, 0.371204, 0.495547, 0.588410, 0.643548,
        0.734379, 0.866828, 1.007356, 1.095127, 1.118583, 1.090269,
        1.010805, 0.902983, 0.722948, 0.512927, 0.282484, 0.090372,
        -0.090733, -0.288677, -0.510539, -0.723598, -0.891926, -1.007011,
        -1.085358, -1.116693, -1.095127, -1.007556, -0.870722, -0.733585,
        -0.642744, -0.580935, -0.495154, -0.372871, -0.224427, -0.074435];

        self.z = [0.039700, 0.111500, 0.174000, 0.204500, 0.166500, 0.067200,
        -0.074000, -0.176000, -0.206500, -0.177000, -0.115000, -0.047000,
        0.038500, 0.109500, 0.174500, 0.203500, 0.164500, 0.069000,
        -0.071500, -0.175000, -0.208000, -0.177000, -0.115000, -0.043000,
        0.039000, 0.115500, 0.175500, 0.204500, 0.163500, 0.067500,
        -0.074000, -0.173500, -0.207500, -0.177000, -0.115500, -0.043000];


class tfc_details_real():
    def __init__(self,):

        self.phi = np.array([3.517, 10.949, 19.466, 29.466, 41.449, 53.517,
        66.483, 79.051, 90.534,100.534,109.051,115.743, 123.517,
        130.449, 139.466, 149.466, 161.449, 173.517, 186.483, 199.051,
        210.534, 220.534, 229.051, 236.195, 243.517, 250.949, 259.466,
        269.466, 281.449, 293.577, 306.483, 319.051, 330.534, 340.534,
        349.051, 356.483])

        #self.phi_correction = np.array([0, 0, 0, 0, 0.5, 0, 0, 0, 0,
        #0, 0, -0.7, 0, -0.5, 0, 0, 0.5, 0, 0, 0, 0, 0, -0.288, 0, 0,
        #0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0])
        self.phi_correction = self.phi*0
        self.phi = self.phi+self.phi_correction

        self.r = np.array([1212.9, 1181.4, 1113.9, 1007.4, 888.9, 800.4, 800.9,
        882.9, 1007.4, 1113.9, 1183.4, 1210.4, 1212.4, 1186.6, 1112.4,
        1009.6, 887.9, 800.4, 803.6, 884.4, 1004.9, 1113.4, 1180.9,
        1211.9, 1212.6, 1181.4, 1113.9, 1007.6, 888.4, 800.4, 799.4,
        886.4, 1006.6, 1118.9, 1181.6, 1213.4])/1000.

        self.z = np.array([39.7, 111.5, 174.0, 204.5, 166.5, 67.2,
        -74.0, -176.0, -206.5, -177.0, -115.0, -47.0, 38.5, 109.5,
        174.5, 203.5, 164.5, 69.0, -71.5, -175.0, -208.0, -177.0,
        -115.0, -43.0, 39.0, 115.5, 175.5, 204.5, 163.5, 67.5, -74.0,
        -173.5, -207.5, -177.0, -115.5, -43.0])/1000.

        self.x = self.r*np.cos(np.deg2rad(self.phi))
        self.y = self.r*np.sin(np.deg2rad(self.phi))
        self.e2 = self.phi-90
        # self.e2 = [-86.482994, -79.050995, -70.533997, -60.533997, -48.550999,
        # -36.482998, -23.517000, -10.948996, 0.533993, 10.533997,
        # 19.050997, 25.742990, 33.516987, 40.449001, 49.466000, 59.465992,
        # 71.448997, 83.516991, 83.516991, 70.949013, 59.466000, 49.466000,
        # 40.949005, 33.804981, 26.482990, 19.051008, 10.534008, 0.534011,
        # -11.449009, -23.576992, -36.482990, -49.051003, -60.533997,
        # -70.533997, -79.050995, -86.483002];

        # self.x = [1.210616, 1.159895, 1.050229, 0.877091, 0.666271, 0.475905,
        # 0.319576, 0.167694, -0.009389, -0.203642, -0.386273, -0.525719,
        # -0.669468, -0.769832, -0.845447, -0.869596, -0.841766, -0.795282,
        # -0.798461, -0.835960, -0.865548, -0.846207, -0.773947, -0.674262,
        # -0.540737, -0.385621, -0.203642, -0.009391, 0.176344, 0.320145,
        # 0.475311, 0.669492, 0.876394, 1.054943, 1.160091, 1.21115];

        # self.y = [0.074405, 0.224389, 0.371204, 0.495547, 0.588410, 0.643548,
        # 0.734379, 0.866828, 1.007356, 1.095127, 1.118583, 1.090269,
        # 1.010805, 0.902983, 0.722948, 0.512927, 0.282484, 0.090372,
        # -0.090733, -0.288677, -0.510539, -0.723598, -0.891926, -1.007011,
        # -1.085358, -1.116693, -1.095127, -1.007556, -0.870722, -0.733585,
        # -0.642744, -0.580935, -0.495154, -0.372871, -0.224427, -0.074435];

        # self.z = [0.039700, 0.111500, 0.174000, 0.204500, 0.166500, 0.067200,
        # -0.074000, -0.176000, -0.206500, -0.177000, -0.115000, -0.047000,
        # 0.038500, 0.109500, 0.174500, 0.203500, 0.164500, 0.069000,
        # -0.071500, -0.175000, -0.208000, -0.177000, -0.115000, -0.043000,
        # 0.039000, 0.115500, 0.175500, 0.204500, 0.163500, 0.067500,
        # -0.074000, -0.173500, -0.207500, -0.177000, -0.115500, -0.043000];

def tfc_points(coil_num = None, tfc_thickness=0.075, tfc_width=0.15, tfc_radius = 0.383, ID_OD_offset = 0,phi_samples = 50,coords=None, use_real=False):
    '''return the points describing a tfc

    SRH: 12July2013
    '''
    #TFC details
    if use_real:
        tfc = tfc_details_real()
    else:
        tfc = tfc_details()
    if coords==None:
        coords = [tfc.x[coil_num],tfc.y[coil_num],tfc.z[coil_num]]
    return _coil_vertical_points(coords, tfc_radius, tfc_thickness, tfc_width, np.arctan2(coords[1], coords[0])*180./np.pi, phi_samples = phi_samples, ID_OD_offset = ID_OD_offset)

# def tfc_points_triangular_mesh(coil, tfc_thickness=0.075, tfc_width=0.15, tfc_radius = 0.383, phi_samples = 50, plot=0):
#     '''return the points describing a tfc

#     SRH: 12July2013
#     '''
#     #TFC details
#     tfc = tfc_details()
#     coil_x, coil_y, coil_z =  _coil_vertical_points([tfc.x[coil],tfc.y[coil],tfc.z[coil]], tfc_radius, tfc_thickness, tfc_width, np.arctan2(tfc.y[coil],tfc.x[coil])*180./np.pi, phi_samples = phi_samples)
#     edges, pts_per_edge = coil_x.shape
#     x = []; y = []; z = []
#     for i in range(edges):
#         x.extend(coil_x[i,:])
#         y.extend(coil_y[i,:])
#         z.extend(coil_z[i,:])
#     vertices = np.array([x, y, z]).T
#     cut_start_points = np.arange(edges)*pts_per_edge
#     faces = np.zeros((2*(pts_per_edge-1)*(edges-1),3),dtype=int)
#     for i in range(edges-1):
#         ind1 = i*(pts_per_edge-1)
#         ind2 = ind1 + pts_per_edge - 1
#         faces[ind1:ind2,:] = np.array([range(cut_start_points[i],cut_start_points[i]+pts_per_edge-1), 
#                                        range(cut_start_points[i+1],cut_start_points[i+1]+pts_per_edge-1), 
#                                        range(cut_start_points[i+1]+1,cut_start_points[i+1]+pts_per_edge)]).T
#         ind1 = (i+edges-1)*(pts_per_edge-1)
#         ind2 = ind1 + pts_per_edge-1
#         faces[ind1:ind2,:] = np.array([range(cut_start_points[i],cut_start_points[i]+pts_per_edge-1), 
#                                        range(cut_start_points[i+1]+1,cut_start_points[i+1]+pts_per_edge),  
#                                        range(cut_start_points[i]+1,cut_start_points[i]+pts_per_edge)]).T
#     if plot:
#         import mayavi.mlab as mlab
#         mlab.triangular_mesh(vertices[:,0], vertices[:,1], vertices[:,2],faces)
#         mlab.triangular_mesh(vertices[:,0], vertices[:,1], vertices[:,2],faces,representation='wireframe',color=(0,0,0))
#     return vertices, faces

def extract_VMEC_surface_data(filename, s_ind=-1, phi_min = 0, phi_max = 2.*np.pi):
    '''extracts x, y, z and B for the surface index s_ind
    using Berhnard's VMEC utilities
    '''
    import h1.mhd_eq.VMEC as VMEC
    theta_vmec = np.linspace(0,2.*np.pi,100)
    phi_real = np.linspace(phi_min, phi_max,100)
    theta_grid, phi_grid = np.meshgrid(theta_vmec, phi_real)
    vmec_tmp = VMEC.VMEC(filename, import_all=True, load_spline=False, save_spline=False, compute_spline_type=0, compute_grid=False, load_grid=False, save_grid=False)
    Rmn_vmec = vmec_tmp.rmnc[s_ind,:]
    Zmn_vmec = vmec_tmp.zmns[s_ind,:]
    Bmn_vmec = vmec_tmp.bmnc[s_ind,:]
    m_vmec = vmec_tmp.xm
    n_vmec = vmec_tmp.xn
    R_vmec = theta_grid *0; Z_vmec = theta_grid * 0
    B_vmec = theta_grid *0;
    for i in range(0,len(n_vmec)):
        argv = (m_vmec[i]*theta_grid - n_vmec[i]*phi_grid) 
        sinv = np.sin(argv)
        cosv = np.cos(argv)
        R_vmec+= Rmn_vmec[i] * cosv
        Z_vmec+= Zmn_vmec[i] * sinv
        B_vmec+= Bmn_vmec[i] * cosv
    phi_vmec_cyl = np.array(phi_grid)
    x_vmec = R_vmec * np.cos(phi_vmec_cyl)
    y_vmec = R_vmec * np.sin(phi_vmec_cyl)
    z_vmec = np.array(Z_vmec)
    return x_vmec, y_vmec, z_vmec, B_vmec




def extract_BOOZER_surface_data(filename, s_ind=-1, phi_min = 0, phi_max = 2.*np.pi):
    '''extracts x, y, z and B for the surface index s_ind
    using Berhnard's VMEC utilities
    '''
    import h1.mhd_eq.BOOZER as BOOZER
    theta_boozer = np.linspace(0,2.*np.pi,100)
    phi_b = np.linspace(phi_min, phi_max,100)
    theta_grid, phi_b_grid = np.meshgrid(theta_boozer, phi_b)
    boozer_tmp = BOOZER.BOOZER(filename, import_all=True, load_spline=False, save_spline=False, compute_spline_type=0, compute_grid=False, load_grid=False, save_grid=False)
    print s_ind
    Rmn_boozer = boozer_tmp.rmnc_b[s_ind,:]
    Zmn_boozer = boozer_tmp.zmns_b[s_ind,:]
    Bmn_boozer = boozer_tmp.bmnc_b[s_ind,:]
    Pmn_boozer = boozer_tmp.pmns_b[s_ind,:]
    m_boozer = boozer_tmp.ixm_b
    n_boozer = boozer_tmp.ixn_b
    R = theta_grid *0; Z = theta_grid * 0
    Phi = theta_grid *0+phi_b_grid
    B = theta_grid *0; 
    for i in range(0,len(n_boozer)):
        argv = (m_boozer[i]*theta_grid - n_boozer[i]*phi_b_grid) 
        sinv = np.sin(argv)
        cosv = np.cos(argv)
        R+= Rmn_boozer[i] * cosv
        Z+= Zmn_boozer[i] * sinv
        Phi += Pmn_boozer[i]*sinv
        B += Bmn_boozer[i] * cosv
    x = R * np.cos(Phi)
    y = R * np.sin(Phi)
    z = +Z
    return x, y, z, B


def extract_BOOZER_fieldline(filename, s_ind=-1, theta_b_start = 0, no_cycles=1, pts_per_cycle=1000, phi_b_start = 0, s_is_ind=True):
    '''extracts x, y, z and B for the surface index s_ind
    using Berhnard's VMEC utilities
    '''
    import h1.mhd_eq.BOOZER as BOOZER
    boozer_tmp = BOOZER.BOOZER(filename, import_all=True, load_spline=False, save_spline=False, compute_spline_type=0, compute_grid=False, load_grid=False, save_grid=False)
    if not s_is_ind:
        print s_ind
        s_ind = np.argmin(np.abs(boozer_tmp.es-s_ind))
        print s_ind, boozer_tmp.es[s_ind]
        print('iota at this surface: {:.2f}'.format(boozer_tmp.iota_b[s_ind]))
    print s_ind, boozer_tmp.es[s_ind], boozer_tmp.rmnc_b.shape

    phi_b = np.linspace(phi_b_start, 2.*np.pi*no_cycles+phi_b_start,pts_per_cycle*np.abs(no_cycles))
    print len(phi_b)
    theta_boozer = phi_b * boozer_tmp.iota_b[s_ind]+theta_b_start
    theta_boozer = np.mod(theta_boozer, 2.*np.pi)
    phi_boozer = np.mod(phi_b, 2.*np.pi)
    theta_grid = theta_boozer
    phi_b_grid = phi_boozer
    Rmn_boozer = boozer_tmp.rmnc_b[s_ind,:]
    Zmn_boozer = boozer_tmp.zmns_b[s_ind,:]
    Bmn_boozer = boozer_tmp.bmnc_b[s_ind,:]
    Pmn_boozer = boozer_tmp.pmns_b[s_ind,:]
    m_boozer = boozer_tmp.ixm_b
    n_boozer = boozer_tmp.ixn_b
    R = theta_grid *0; Z = theta_grid * 0
    Phi = theta_grid *0+phi_b_grid
    B = theta_grid *0; 
    for i in range(0,len(n_boozer)):
        argv = (m_boozer[i]*theta_grid - n_boozer[i]*phi_b_grid) 
        sinv = np.sin(argv)
        cosv = np.cos(argv)
        R+= Rmn_boozer[i] * cosv
        Z+= Zmn_boozer[i] * sinv
        Phi += Pmn_boozer[i]*sinv
        B += Bmn_boozer[i] * cosv
    x = R * np.cos(Phi)
    y = R * np.sin(Phi)
    z = +Z
    return x, y, z, B, Phi


def mirnov_locations():
    #HMA coil locations R, Z, phi(deg)
    loc = np.array([[0.9118, 0.0287, 40.7600],
                    [0.9359, 0.0688, 31.5200],
                    [0.9746, 0.0926, 22.5200],
                    [1.0182, 0.0966, 13.8800],
                    [1.0575, 0.0821, 5.6000],
                    [1.0865, 0.0535, -2.4400],
                    [1.1011, 0.0167, -10.2400],
                    [1.0997, -0.0234, -18.1600],
                    [1.0822, -0.0595, -26.0800],
                    [1.0511, -0.0859, -34.1200],
                    [1.0108, -0.0973, -42.4000],
                    [0.9671, -0.0898, -51.1600],
                    [0.9304, -0.0628, -60.1600],
                    [0.9098, -0.0209, -69.4000],
                    [0.9111, 0.0264, -78.7600],
                    [0.9342, 0.0670, -88.0000]])

    HMA_x = loc[:,0]*np.cos(loc[:,2]/180.*np.pi)
    HMA_y = loc[:,0]*np.sin(loc[:,2]/180.*np.pi)
    HMA_z = loc[:,1]

    #Poloidal array1 R, phi(rad), Z (NOTE DIFFERENT TO HMA!!)
    pol_array1 = np.array([[1.114, 0.7732, 0.355],
                           [1.185, 0.7732, 0.289],
                           [1.216, 0.7732, 0.227],
                           [1.198, 0.7732, 0.137],
                           [1.129, 0.7732, 0.123],
                           [1.044, 0.7732, 0.128],
                           [0.963, 0.7732, 0.112],
                           [0.924, 0.7732, 0.087],
                           [0.902, 0.7732, 0.052],
                           [0.900, 0.7732, -0.008],
                           [0.925, 0.7732, -0.073],
                           [0.964, 0.7732, -0.169],
                           [0.897, 0.7732, -0.238],
                           [0.821, 0.7732, -0.221],
                           [0.696, 0.7732, -0.106],
                           [0.652, 0.7732, 0.036],
                           [0.676, 0.7732, 0.193],
                           #[0.790, 0.7732, 0.326],
                           [0.806, 0.7732, 0.336],
                           [0.934, 0.7732, 0.383],
                           [1.114, 0.7732, 0.355]])
    #Poloidal array1 R, phi(rad), Z (NOTE DIFFERENT TO HMA!!)
    pol_array2 = np.array([[1.114, 4.962, 0.355],
                           [1.185, 4.962, 0.289],
                           [1.216, 4.962, 0.227],
                           [1.198, 4.962, 0.137],
                           [1.129, 4.962, 0.123],
                           [1.044, 4.962, 0.128],
                           [0.963, 4.962, 0.112],
                           [0.924, 4.962, 0.087],
                           [0.902, 4.962, 0.052],
                           [0.900, 4.962, -0.008],
                           [0.925, 4.962, -0.073],
                           [0.964, 4.962, -0.169],
                           [0.897, 4.962, -0.238],
                           [0.821, 4.962, -0.221],
                           [0.696, 4.962, -0.106],
                           [0.652, 4.962, 0.036],
                           [0.676, 4.962, 0.193],
                           #[0.790, 4.962, 0.326],
                           [0.806, 4.962, 0.336],
                           [0.934, 4.962, 0.383],
                           [1.114, 4.962, 0.355]])
    #Convert to Cartesian
    pol_array1_x = pol_array1[:,0]*np.cos(pol_array1[:,1])
    pol_array1_y = pol_array1[:,0]*np.sin(pol_array1[:,1])
    pol_array1_z = pol_array1[:,2]
    pol_array2_x = pol_array2[:,0]*np.cos(pol_array2[:,1])
    pol_array2_y = pol_array2[:,0]*np.sin(pol_array2[:,1])
    pol_array2_z = pol_array2[:,2]
    return HMA_x, HMA_y, HMA_z, pol_array1_x, pol_array1_y, pol_array1_z, pol_array2_x, pol_array2_y, pol_array2_z


def plot_vmec(vmec_filename = '/home/srh112/code/python/h1_eq_generation/results7/kh0.100-kv1.000fixed/wout_kh0.100-kv1.000fixed.nc', phi_min=0, phi_max=2.*np.pi, s_ind = -1, vmec_args=None):
    #vmec_filename = '/home/srh112/code/python/h1_eq_generation/results7/kh0.100-kv1.000fixed/wout_kh0.100-kv1.000fixed.nc'
    x, y, z, B = extract_VMEC_surface_data(vmec_filename, s_ind=s_ind, phi_min = phi_min, phi_max = phi_max)
    vmec_kwargs={'opacity':1.0,'scalars':B,'colormap':'hot','representation':'surface' }
    if vmec_args!=None: 
        for i in vmec_args.keys(): vmec_kwargs[i] = vmec_args[i]

    #pts = mlab.mesh(x[:,:], y[:,:], z[:,:], opacity = 1.0, scalars = B, colormap = 'hot', representation='surface')
    pts = mlab.mesh(x[:,:], y[:,:], z[:,:], **vmec_kwargs)
    return x, y, z, B

def plot_boozer(boozer_filename = '/home/srh112/code/python/h1_eq_generation/results7/kh0.100-kv1.000fixed/boozmn_wout_kh0.100-kv1.000fixed.nc', phi_min=0, phi_max=2.*np.pi, s_ind = -1, boozer_args=None):
    #vmec_filename = '/home/srh112/code/python/h1_eq_generation/results7/kh0.100-kv1.000fixed/wout_kh0.100-kv1.000fixed.nc'
    x, y, z, B = extract_BOOZER_surface_data(boozer_filename, s_ind=s_ind, phi_min = phi_min, phi_max = phi_max)
    boozer_kwargs={'opacity':1.0,'scalars':B,'colormap':'hot','representation':'surface' }
    if boozer_args!=None: 
        for i in boozer_args.keys(): boozer_kwargs[i] = boozer_args[i]

    #pts = mlab.mesh(x[:,:], y[:,:], z[:,:], opacity = 1.0, scalars = B, colormap = 'hot', representation='surface')
    pts = mlab.mesh(x[:,:], y[:,:], z[:,:], **boozer_kwargs)
    return x, y, z, B

def make_plot(phi_min = 0, phi_max = 2.*np.pi):
    vmec_filename = '/home/srh112/code/python/h1_eq_generation/results7/kh0.100-kv1.000fixed/wout_kh0.100-kv1.000fixed.nc'
    f = mlab.figure(1, fgcolor=(0, 0, 0), bgcolor=(1, 1, 1))
    from mayavi.sources.poly_data_reader import PolyDataReader 
    a = PolyDataReader()
    a.initialize("/home/srh112/code/python/python-h1/h1/h1model/helical.stl")
    e = mlab.get_engine()
    e.add_source(a)
    from mayavi.modules.surface import Surface
    s = Surface()
    e.add_module(s)

    x, y, z, B = plot_vmec(vmec_filename=vmec_filename, phi_min = phi_min, phi_max = phi_max)
    #x, y, z, B = extract_VMEC_surface_data(vmec_filename, s_ind=-1, phi_min = phi_min, phi_max = phi_max)
    #pts = mlab.mesh(x[:,:], y[:,:], z[:,:], opacity = 1.0, scalars = B, colormap = 'hot', representation='surface')

    #plot the TFC's 
    include_coils = range(5,27,1)
    tfc_thickness=0.075;tfc_width=0.15; tfc_radius = 0.383
    #tfc_mesh_props = {'opacity':0.3,'color':(1,0.,0.)}
    tfc_mesh_props = {'opacity':1.0,'color':(0.5,0.5,0.5)}
    plot_tfc(include_coils, tfc_thickness=0.075, tfc_width=0.15, tfc_radius = 0.383, tfc_mesh_props = None)

    pfc_mesh_props = {'opacity':1.,'color':(0.5,0.5,0.5)}
    plot_pfc(pfc_thickness = 0.11, pfc_width = 0.11, pfc_radius = 1.0, pfc_mesh_props=pfc_mesh_props)

    ovc_mesh_props = {'opacity':1,'color':(255/255., 255/255., 51/255.)}
    #plot_ovc(ovc_mesh_props=ovc_mesh_props)

    ivc_mesh_props = {'opacity':1,'color':(0.5,0.5,0.5)}
    #plot_ivc(ivc_mesh_props=ivc_mesh_props)

    #Plot the HMA and poloidal Mirnov arrays as cubes joined by a line
    show_hma = 1; show_pol_array1 = 1; show_pol_array2 = 1
    HMA_x, HMA_y, HMA_z, pol_array1_x, pol_array1_y, pol_array1_z, pol_array2_x, pol_array2_y, pol_array2_z = mirnov_locations()
    if show_hma:
        mlab.plot3d(HMA_x, HMA_y, HMA_z,line_width=1,tube_radius=0.01)
        #mlab.points3d(HMA_x, HMA_y, HMA_z, scale_mode='none', scale_factor = 0.04, color=(0.5,0.5,0.5),mode='cube')
        mlab.points3d(HMA_x, HMA_y, HMA_z, scale_mode='none', scale_factor = 0.04, color=(0,0.,1),mode='cube')
    if show_pol_array1:
        mlab.plot3d(pol_array1_x, pol_array1_y, pol_array1_z,line_width=1,tube_radius=0.02)
        mlab.points3d(pol_array1_x, pol_array1_y, pol_array1_z, scale_mode='none', scale_factor = 0.04, color=(0,1,0),mode='cube')
    if show_pol_array2:
        mlab.plot3d(pol_array2_x, pol_array2_y, pol_array2_z,line_width=1,tube_radius=0.02)
        mlab.points3d(pol_array2_x, pol_array2_y, pol_array2_z, scale_mode='none', scale_factor = 0.04, color=(0.,1.,0.),mode='cube')

    return x, y, z
    #mlab.show()


def make_plot2():
    vmec_filename = '/home/srh112/code/python/h1_eq_generation/results7/kh0.100-kv1.000fixed/wout_kh0.100-kv1.000fixed.nc'
    mlab.figure(1, fgcolor=(0, 0, 0), bgcolor=(1, 1, 1))
    x, y, z, B = plot_vmec(vmec_filename=vmec_filename, phi_min = phi_min, phi_max = phi_max)
    # vmec_filename = '/home/srh112/code/python/h1_eq_generation/results7/kh0.100-kv1.000fixed/wout_kh0.100-kv1.000fixed.nc'
    # x, y, z, B = extract_VMEC_surface_data(vmec_filename, s_ind=-1)
    # pts = mlab.mesh(x, y, z, opacity = 1.0, scalars = B, colormap = 'hot', representation='surface')

    #plot the TFC's 
    include_coils = range(5,27)
    tfc_thickness=0.075;tfc_width=0.15; tfc_radius = 0.383
    #tfc_mesh_props = {'opacity':0.3,'color':(1,0.,0.)}
    tfc_mesh_props = {'opacity':1.0,'color':(0.5,0.5,0.5)}
    plot_tfc(include_coils, tfc_thickness=0.075, tfc_width=0.15, tfc_radius = 0.383, tfc_mesh_props = None)

    pfc_mesh_props = {'opacity':1.,'color':(0.5,0.5,0.5)}
    plot_pfc(pfc_thickness = 0.11, pfc_width = 0.11, pfc_radius = 1.0, pfc_mesh_props=pfc_mesh_props)

    ovc_mesh_props = {'opacity':1,'color':(255/255., 255/255., 51/255.)}
    plot_ovc(ovc_mesh_props=ovc_mesh_props)

    ivc_mesh_props = {'opacity':1,'color':(0.5,0.5,0.5)}
    plot_ivc(ivc_mesh_props=ivc_mesh_props)

    #Plot the HMA and poloidal Mirnov arrays as cubes joined by a line
    show_hma = 1; show_pol_array1 = 1; show_pol_array2 = 1
    HMA_x, HMA_y, HMA_z, pol_array1_x, pol_array1_y, pol_array1_z, pol_array2_x, pol_array2_y, pol_array2_z = mirnov_locations()
    if show_hma:
        mlab.plot3d(HMA_x, HMA_y, HMA_z,line_width=1,tube_radius=0.01)
        #mlab.points3d(HMA_x, HMA_y, HMA_z, scale_mode='none', scale_factor = 0.04, color=(0.5,0.5,0.5),mode='cube')
        mlab.points3d(HMA_x, HMA_y, HMA_z, scale_mode='none', scale_factor = 0.04, color=(0,0.,1),mode='cube')
    if show_pol_array1:
        mlab.plot3d(pol_array1_x, pol_array1_y, pol_array1_z,line_width=1,tube_radius=0.02)
        mlab.points3d(pol_array1_x, pol_array1_y, pol_array1_z, scale_mode='none', scale_factor = 0.04, color=(0,1,0),mode='cube')
    if show_pol_array2:
        mlab.plot3d(pol_array2_x, pol_array2_y, pol_array2_z,line_width=1,tube_radius=0.02)
        mlab.points3d(pol_array2_x, pol_array2_y, pol_array2_z, scale_mode='none', scale_factor = 0.04, color=(0.,1.,0.),mode='cube')

    mlab.show()
