import numpy as np
import mayavi.mlab as mlab

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

def _plot_coil_vertical(center, radius, thickness, width, phi, phi_samples = 50,**kw):
    '''Plot a horizontal coil like hte PFC or OVC
    '''
    rmin = radius - width/2
    rmax = radius + width/2
    heights = [-thickness/2,thickness/2,thickness/2,-thickness/2,-thickness/2]
    r_values = [rmax,rmax,rmin,rmin,rmax]
    phi = phi/180.*np.pi
    centre_x = center[0]
    centre_y = center[1]
    centre_z = center[2]
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
        circle_z1[i,:] = circle_z #np.cos(phi)*0+heights[i]+centre[2]
    circle_x1 = circle_x1 + centre_x
    circle_y1 = circle_y1 + centre_y
    circle_z1 = circle_z1 + centre_z
    mlab.mesh(circle_x1, circle_y1, circle_z1,**kw)#color=(1,0,0),opacity=opacity)


def plot_pfc(pfc_thickness = 0.11, pfc_width = 0.11, pfc_radius = 1.0, pfc_mesh_props=None):
    '''plot the pfc

    SRH: 12July2013
    '''

    pfc_centre = [0,0,0]
    if pfc_mesh_props == None:
        pfc_mesh_props = {'opacity':1.,'color':(0.5,0.5,0.5)}
    _plot_coil_horizontal(pfc_centre, pfc_thickness, pfc_width, pfc_radius, **pfc_mesh_props)

def plot_tfc(include_coils, tfc_thickness=0.075, tfc_width=0.15, tfc_radius = 0.383, tfc_mesh_props = None):
    '''plot the tfc coils in the include_coils list

    SRH: 12July2013
    '''
    #TFC details
    e2 = [-86.482994, -79.050995, -70.533997, -60.533997, -48.550999,
    -36.482998, -23.517000, -10.948996, 0.533993, 10.533997,
    19.050997, 25.742990, 33.516987, 40.449001, 49.466000, 59.465992,
    71.448997, 83.516991, 83.516991, 70.949013, 59.466000, 49.466000,
    40.949005, 33.804981, 26.482990, 19.051008, 10.534008, 0.534011,
    -11.449009, -23.576992, -36.482990, -49.051003, -60.533997,
    -70.533997, -79.050995, -86.483002];

    x = [1.210616, 1.159895, 1.050229, 0.877091, 0.666271, 0.475905,
    0.319576, 0.167694, -0.009389, -0.203642, -0.386273, -0.525719,
    -0.669468, -0.769832, -0.845447, -0.869596, -0.841766, -0.795282,
    -0.798461, -0.835960, -0.865548, -0.846207, -0.773947, -0.674262,
    -0.540737, -0.385621, -0.203642, -0.009391, 0.176344, 0.320145,
    0.475311, 0.669492, 0.876394, 1.054943, 1.160091, 1.21115];

    y = [0.074405, 0.224389, 0.371204, 0.495547, 0.588410, 0.643548,
    0.734379, 0.866828, 1.007356, 1.095127, 1.118583, 1.090269,
    1.010805, 0.902983, 0.722948, 0.512927, 0.282484, 0.090372,
    -0.090733, -0.288677, -0.510539, -0.723598, -0.891926, -1.007011,
    -1.085358, -1.116693, -1.095127, -1.007556, -0.870722, -0.733585,
    -0.642744, -0.580935, -0.495154, -0.372871, -0.224427, -0.074435];

    z = [0.039700, 0.111500, 0.174000, 0.204500, 0.166500, 0.067200,
    -0.074000, -0.176000, -0.206500, -0.177000, -0.115000, -0.047000,
    0.038500, 0.109500, 0.174500, 0.203500, 0.164500, 0.069000,
    -0.071500, -0.175000, -0.208000, -0.177000, -0.115000, -0.043000,
    0.039000, 0.115500, 0.175500, 0.204500, 0.163500, 0.067500,
    -0.074000, -0.173500, -0.207500, -0.177000, -0.115500, -0.043000];

    #plot the TFC's 
    #tfc_mesh_props = {'opacity':0.3,'color':(1,0.,0.)}
    if tfc_mesh_props==None:
        tfc_mesh_props = {'opacity':1.0,'color':(0.5,0.5,0.5)}
    #include_coils = [11,12,13]#include_coils.extend(range(1,15,4))
    #include_coils.extend(range(26,36,4))
    #for i in range(1,36,4):
    for i in include_coils:
         _plot_coil_vertical([x[i],y[i],z[i]], tfc_radius, tfc_thickness, tfc_width, np.arctan2(y[i],x[i])*180./np.pi, **tfc_mesh_props)

