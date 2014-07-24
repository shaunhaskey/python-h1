import h1.mhd_eq.BOOZER as BOOZER
import time
import scipy.optimize as optimize
import numpy as np
import matplotlib.pyplot as pt

class probe_array():
    def __init__(self,):
        pass
    def loc_boozer_coords(self, filename=None, boozer_phi = None, boozer_theta = None, distance = None):
        if filename == None:
            self.boozer_phi = boozer_phi
            self.boozer_theta = boozer_theta
            self.distance = distance
        else:
            booz_obj2 = BOOZER.BOOZER(filename, import_all=True, load_spline=False, save_spline=False, compute_spline_type=0, compute_grid=False, load_grid=False, save_grid=False)
            booz_obj = boozer_object(booz_obj2)
            kernel = -1.; a = 1.#2.*np.pi/3
            self.boozer_phi, self.boozer_theta, self.distance = find_coil_locs(self.cart_x, self.cart_y, self.cart_z, booz_obj, kernel, a,func)
        
    def plot_fig(self, ax = None, ax2 = None, mask = None, ax2_xaxis=None, ax3 = None):
        no_ax = True if  ax == None else False
        if no_ax: fig, ax = pt.subplots(nrows = 1, ncols = 1, sharex = False, sharey = False)
        if mask == None: mask = np.ones(len(self.boozer_phi), dtype = bool)
        if ax2_xaxis==None: ax2_xaxis = np.unwrap(np.deg2rad(self.boozer_phi[mask]))
        #ax.imshow(np.abs(vals), interpolation = 'nearest',cmap = 'binary', origin = 'upper',extent = [m_list[0],m_list[-1],n_list[-1],n_list[0]])
        im = ax.pcolormesh(self.m_rec, self.n_rec, np.abs(self.vals), cmap = 'binary')
        im.set_clim([0,0.75])
        #ax.plot(m_real, n_real,'bo')
        #pt.colorbar(im, ax = ax)
        max_ind = np.argsort(np.abs(self.vals.flatten()))
        txt = ''
        for i in range(-1,-5,-1):txt+='({},{})'.format(self.n_rec.flatten()[max_ind[i]],self.m_rec.flatten()[max_ind[i]])
        print txt
        ax.text(1,-8,txt, color = 'r')
        if ax2!=None:
            n = self.n_rec.flatten()[max_ind[-1]]
            m = self.m_rec.flatten()[max_ind[-1]]
            new_booz_th = np.linspace(np.min(self.boozer_theta[mask]),np.max(self.boozer_theta[mask]),300)
            new_booz_phi = np.interp(new_booz_th, self.boozer_theta[mask], self.boozer_phi[mask])
            kernel = 1j*m*np.deg2rad(self.boozer_theta[mask]) + self.inc_phi*1j*n*np.deg2rad(self.boozer_phi[mask])
            ax2.plot(ax2_xaxis, np.real(self.vals.flatten()[max_ind[-1]]*np.exp(kernel)), 'b--')
            kernel_interp = 1j*m*np.deg2rad(new_booz_th) + self.inc_phi * 1j*n*np.deg2rad(new_booz_phi)
            ax2.plot(np.deg2rad(new_booz_th), np.real(self.vals.flatten()[max_ind[-1]]*np.exp(kernel_interp)), 'b-.')
            #if ax3!=None: ax3.plot(self.m_rec[0,:], np.abs(self.vals[0,:]),'-bo')
            if ax3!=None: 
                foo1 = self.vals.flatten()[max_ind[-1]]* np.exp(kernel)
                foo1 = foo1/np.abs(foo1)
                ax3.plot(np.real(foo1), np.imag(foo1),'bo')
                for j in range(foo1.shape[0]): 
                    ax3.text(np.real(foo1[j]), np.imag(foo1[j]),str(j+1))

        if no_ax: fig.canvas.draw();fig.show()

    def perform_fit(self, data, mask = None, inc_phi = True):
        n_list = np.arange(-10,11)
        m_list = np.arange(-10,11)
        self.inc_phi = inc_phi
        self.vals = np.zeros((n_list.shape[0], m_list.shape[0]),dtype=complex)
        self.n_rec = np.zeros((n_list.shape[0], m_list.shape[0]),dtype=int)
        self.m_rec = np.zeros((n_list.shape[0], m_list.shape[0]),dtype=int)
        if mask == None: mask = np.ones(len(self.boozer_phi), dtype = bool)
        for i, n in enumerate(n_list):
            for j, m in enumerate(m_list):
                self.vals[i,j] = np.sum(data * np.exp(-1j*m*np.deg2rad(self.boozer_theta[mask]) - self.inc_phi*1j*n*np.deg2rad(self.boozer_phi[mask])))/np.sum(mask)
                self.n_rec[i,j] = +n
                self.m_rec[i,j] = +m

    def dummy_data(self,):
        probe_locations_theta = np.deg2rad(self.boozer_theta)
        probe_locations_phi = np.deg2rad(self.boozer_phi)
        n_real = 8
        m_real = 4
        data = np.exp(1j*m_real*probe_locations_theta + 1j*n_real*probe_locations_phi)
        self.perform_fit(data)
        self.plot_fig()

    def test_fit(self,):
        probe_locations_theta1 = [np.linspace(0,2.*np.pi,10, endpoint = False), np.linspace(0,2.*np.pi,10, endpoint = False), np.linspace(0,2.*np.pi,10, endpoint = False)]
        probe_locations_phi1 = [np.linspace(0,1./3*2.*np.pi,10, endpoint = False), np.linspace(0,0,10, endpoint = False), np.linspace(0.5,0.5,10, endpoint = False)]
        probe_locations_theta = np.hstack(probe_locations_theta1)
        probe_locations_phi = np.hstack(probe_locations_phi1)
        probe_locations_theta = np.deg2rad(self.boozer_theta)
        probe_locations_phi = np.deg2rad(self.boozer_phi)
        n_real = 8
        m_real = 4
        data = np.exp(1j*m_real*probe_locations_theta + 1j*n_real*probe_locations_phi)
        n_list = np.arange(-10,10)
        m_list = np.arange(-10,10)
        vals = np.zeros((n_list.shape[0], m_list.shape[0]),dtype=complex)
        n_rec = np.zeros((n_list.shape[0], m_list.shape[0]),dtype=int)
        m_rec = np.zeros((n_list.shape[0], m_list.shape[0]),dtype=int)

        for i, n in enumerate(n_list):
            for j, m in enumerate(m_list):
                vals[i,j] = np.sum(data * np.exp(-1j*m*probe_locations_theta - 1j*n*probe_locations_phi))
                n_rec[i,j] = n
                m_rec[i,j] = m

        fig, ax = pt.subplots(nrows = 1, ncols = 1, sharex = False, sharey = False)
        #ax.imshow(np.abs(vals), interpolation = 'nearest',cmap = 'binary', origin = 'upper',extent = [m_list[0],m_list[-1],n_list[-1],n_list[0]])
        im = ax.pcolormesh(m_rec, n_rec, np.abs(vals), cmap = 'binary')
        ax.plot(m_real, n_real,'bo')
        pt.colorbar(im, ax = ax)
        max_ind = np.argmax(np.abs(vals.flatten())) 
        max_ind = np.argsort(np.abs(vals.flatten()))

        ax.set_title('max:n,m = {}, {} and {},{}'.format(n_rec.flatten()[max_ind[-1]],m_rec.flatten()[max_ind[-1]],n_rec.flatten()[max_ind[-2]],m_rec.flatten()[max_ind[-2]]))
        fig.canvas.draw();fig.show()



def func(tmp, booz_obj, coil_loc, kernel, a):
    tmp2 = booz_obj.return_values_single2(tmp[0],tmp[1],3,coil_loc=coil_loc,kernel=kernel,a=a)
    #print tmp, tmp2
    return tmp2

def find_coil_locs(x_coil, y_coil, z_coil, booz_obj,kernel,a,opt_func):
    answer_list = []
    distance_list = []
    min_locations_phi = []
    min_locations_theta = []
    start_time= time.time()
    for i in range(0,len(x_coil)):
        coil_loc=[x_coil[i],y_coil[i],z_coil[i]]
        if i==0:
            q = optimize.fmin(opt_func,[0,0],args=(booz_obj,coil_loc, kernel, a))
        else:
            q = optimize.fmin(opt_func,q,args=(booz_obj,coil_loc, kernel, a))
        distance_list.append(booz_obj.distance)
        min_locations_phi.append(q[1]*180/np.pi)
        min_locations_theta.append(q[0]*180/np.pi)
    print 'finished all coils in : %.4fs'%(time.time() - start_time)
    #print min_locations_phi
    #print min_locations_theta
    min_locations_phi = np.array(min_locations_phi)
    min_locations_theta = np.array(min_locations_theta)
    for j in range(1,len(min_locations_phi)):
        if min_locations_phi[j]>min_locations_phi[j-1]:
            min_locations_phi[j]-=360
        if min_locations_theta[j]>min_locations_theta[j-1]:
            min_locations_theta[j]-=360
    if np.min(min_locations_theta)<-360:
        min_locations_theta = min_locations_theta+360
    return min_locations_phi, min_locations_theta, distance_list


#Mirnov coil locations
class HMA(probe_array):
    '''Returns the locations of the HMA 
    [R(m), phi(rad), z(m)], [x(m), y(m), z(m)]
    '''
    def __init__(self,):
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
        self.cart_x = loc[:,0]*np.cos(loc[:,2]/180.*np.pi)
        self.cart_y = loc[:,0]*np.sin(loc[:,2]/180.*np.pi)
        self.cart_z = loc[:,1]

        self.cyl = np.array([loc[:,0], np.deg2rad(loc[:,2]), loc[:,1]]).T
        self.cart = np.array([self.cart_x,self.cart_y,self.cart_z]).T

        #These are the outputs from each coil per AMP in the field coil
        #rows are former # starting at 1 and columns are x (Blue), y (Black axial), z (Grey)
        #The x, y, and z correspond with the markings on the copper amplifier box
        self.hfc_sens = np.array([[ 0.24680044,  0.38671494,  1.99371736],
                                  [ 1.8162571 ,  0.3139969 ,  1.11030899],
                                  [ 0.77506133,  0.34716982, -2.05733878],
                                  [-1.98664767,  0.27446152, -0.63596707],
                                  [-0.98941079,  0.24829954,  1.57720848],
                                  [ 0.07520408,  0.43666601,  1.64201736],
                                  [ 1.48546853,  0.28634358,  0.20937708],
                                  [ 0.96485716,  0.22613518, -1.01469717],
                                  [-0.84246615,  0.17836883, -1.27135465],
                                  [-1.6200651 ,  0.28164495, -0.15445505],
                                  [-0.64235141,  0.3044398 ,  1.42710179],
                                  [ 0.14169116,  0.27369179,  1.42189146],
                                  [ 1.40695573,  0.3618185 ,  0.9818182 ],
                                  [ 1.67037656,  0.41032717, -0.99520859],
                                  [-0.09590273,  0.38489415, -1.58732869],
                                  [-1.60452575,  0.3077779 ,  0.06969922]])

        #rows are former # starting at 1 and columns are x (Blue), y (Black axial), z (Grey)
        #The x, y, and z correspond with the markings on the copper amplifier box
        self.ovc_sens = np.array([[ 0.32148735,  0.14624306,  0.33290093],
                                  [ 0.47703375,  0.10144618,  0.12750015],
                                  [ 0.1966587 ,  0.0413601 , -0.45012207],
                                  [-0.37909552, -0.02401509, -0.31848575],
                                  [-0.47986755, -0.05278911,  0.11920395],
                                  [-0.46650299, -0.04302656,  0.18505433],
                                  [-0.10538589, -0.1053388 ,  0.47133019],
                                  [ 0.20085195, -0.13279013,  0.43445564],
                                  [ 0.47393702, -0.0971154 ,  0.12427098],
                                  [ 0.48078232, -0.07264135, -0.10010139],
                                  [ 0.33590748, -0.0291113 , -0.37420902],
                                  [ 0.28142761,  0.05629237, -0.40164629],
                                  [ 0.08206854,  0.13642537, -0.46844489],
                                  [-0.17726428,  0.11452079, -0.45364346],
                                  [-0.37137112,  0.14474828, -0.307499  ],
                                  [-0.41907554,  0.11838274,  0.24172061]])

        #rows are former # starting at 1 and columns are x (Blue), y (Black axial), z (Grey)
        #The x, y, and z correspond with the markings on the copper amplifier box
        self.pfc_sens = np.array([[ -7.99749218,  -3.05858866,  -4.5654661 ],
                                  [ -8.65408697,  -3.33189089,   3.81550182],
                                  [  5.98135272,  -2.95477332,   6.63519082],
                                  [  6.1890899 ,  -1.94152672,  -6.18416383],
                                  [ -4.85713685,  -1.59959204,  -7.15566348],
                                  [ -7.21754124,  -1.69387596,  -2.48902452],
                                  [ -3.33778916,  -1.62587849,   5.92913341],
                                  [  3.60105246,  -1.9112014 ,   5.92601954],
                                  [  7.28507464,  -2.08076584,  -2.91896605],
                                  [  1.03232588,  -2.26124645,  -8.1866978 ],
                                  [ -6.73277036,  -2.1534624 ,  -4.88527605],
                                  [ -9.00798967,  -2.89434804,  -1.71041815],
                                  [ -7.62557428,  -4.064121  ,   6.87788172],
                                  [  1.79385294,  -2.74557257,  10.05864949],
                                  [  8.09985463,  -3.08717648,   5.33294451],
                                  [  5.49885075,  -3.42119878,  -8.3755888 ]])

        #rows are former # starting at 1 and columns are x (Blue), y (Black axial), z (Grey)
        #The x, y, and z correspond with the markings on the copper amplifier box
        self.tfc_sens = np.array([[ 0.65129772, -7.49218624,  1.62549694],
                                  [ 1.63950254, -7.99839531, -0.24498995],
                                  [-0.42725304, -8.41199265, -1.50359276],
                                  [-1.09843213, -9.00588666,  0.1373979 ],
                                  [-0.42053863, -9.23841439,  0.80302916],
                                  [ 0.12297388, -8.88902202,  2.16134651],
                                  [ 1.48875708, -8.56707863, -0.30746247],
                                  [-0.11581799, -8.00095416, -1.48414972],
                                  [-1.2693579 , -7.47963575,  0.15921665],
                                  [-1.10589928, -6.99862441,  0.72853838],
                                  [ 0.09033151, -6.78518809,  1.32462539],
                                  [ 1.10077608, -6.79465548,  0.99386562],
                                  [ 1.78605707, -6.84787579, -0.42116655],
                                  [ 1.01115844, -7.21598888, -0.97219593],
                                  [-0.55493877, -7.52592757, -1.70943072],
                                  [-1.15746014, -8.02568961,  1.29325897]])
        
        #rows are former # starting at 1 and columns are x (Blue), y (Black axial), z (Grey)
        #The x, y, and z correspond with the markings on the copper amplifier box
        self.resistance_comp = np.array([[ 1.12688612,  1.06576173,  1.00143256],
                                      [ 1.14329642,  1.0242431 ,  1.05996262],
                                      [ 1.12909609,  1.01320953,  1.05431591],
                                      [ 1.12499777,  1.01664142,  1.0512358 ],
                                      [ 1.11117414,  1.00448864,  1.03784448],
                                      [ 1.12136575,  1.00807551,  1.04840498],
                                      [ 1.12423223,  1.01082688,  1.05108459],
                                      [ 1.12519099,  1.01134343,  1.04474993],
                                      [ 1.12768874,  1.01376638,  1.04724921],
                                      [ 1.12221355,  1.06012822,  1.00314682],
                                      [ 1.12394626,  1.06213138,  1.01811234],
                                      [ 1.12510437,  1.06304063,  0.99999999],
                                      [ 1.12976077,  1.07367974,  1.00379342],
                                      [ 1.12298962,  1.06854099,  1.01712037],
                                      [ 1.12435199,  1.0632905 ,  1.01356087],
                                      [ 1.1201823 ,  1.06012822,  1.00936083]])




class PMA1(probe_array):
    '''Returns the locations of the PMA1 
    [R(m), phi(rad), z(m)], [x(m), y(m), z(m)]
    '''
    def __init__(self,):
        self.cyl = np.array([[1.114, 0.7732, 0.355],
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
                             [0.790, 0.7732, 0.326],
                             [0.806, 0.7732, 0.336],
                             [0.934, 0.7732, 0.383]])
        self.cart_x = self.cyl[:,0]*np.cos(self.cyl[:,1])
        self.cart_y = self.cyl[:,0]*np.sin(self.cyl[:,1])
        self.cart_z = self.cyl[:,2]
        self.cart = np.array([self.cart_x, self.cart_y, self.cart_z]).T

class PMA1_reduced(probe_array):
    '''Returns the locations of the PMA1 
    [R(m), phi(rad), z(m)], [x(m), y(m), z(m)]
    '''
    def __init__(self,):
        self.cyl = np.array([[1.114, 0.7732, 0.355],
                             [1.185, 0.7732, 0.289],
                             [1.216, 0.7732, 0.227],
                             [1.198, 0.7732, 0.137],
                             #[1.129, 0.7732, 0.123],
                             #[1.044, 0.7732, 0.128],
                             [0.963, 0.7732, 0.112],
                             [0.924, 0.7732, 0.087],
                             [0.902, 0.7732, 0.052],
                             [0.900, 0.7732, -0.008],
                             #[0.925, 0.7732, -0.073],
                             #[0.964, 0.7732, -0.169],
                             #[0.897, 0.7732, -0.238],
                             #[0.821, 0.7732, -0.221],
                             [0.696, 0.7732, -0.106],
                             #[0.652, 0.7732, 0.036],
                             #[0.676, 0.7732, 0.193],
                             #[0.790, 0.7732, 0.326],
                             #[0.806, 0.7732, 0.336],
                             #[0.934, 0.7732, 0.383]
                             ])
        self.cart_x = self.cyl[:,0]*np.cos(self.cyl[:,1])
        self.cart_y = self.cyl[:,0]*np.sin(self.cyl[:,1])
        self.cart_z = self.cyl[:,2]
        self.cart = np.array([self.cart_x, self.cart_y, self.cart_z]).T

class PMA2(probe_array):
    '''Returns the locations of the PMA2 
    [R(m), phi(rad), z(m)], [x(m), y(m), z(m)]
    '''
    #Poloidal array1 R, phi(rad), Z (NOTE DIFFERENT TO HMA!!)
    def __init__(self,):
        self.cyl = np.array([[1.114, 4.962, 0.355],
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
                             [0.934, 4.962, 0.383]])

        self.cyl = np.array([[1.114, 4.962, 0.355],
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
                             [0.790, 4.962, 0.326],
                             [0.806, 4.962, 0.336],
                             [0.934, 4.962, 0.383]])
        #Convert to Cartesian
        self.cart_x = self.cyl[:,0]*np.cos(self.cyl[:,1])
        self.cart_y = self.cyl[:,0]*np.sin(self.cyl[:,1])
        self.cart_z = self.cyl[:,2]
        self.cart = np.array([self.cart_x, self.cart_y, self.cart_z]).T


class boozer_object:
    def __init__(self, input_obj):
        self.m = input_obj.ixm_b
        self.n = input_obj.ixn_b
        self.Rmn = input_obj.rmnc_b[-1,:]
        self.Zmn = input_obj.zmns_b[-1,:]
        self.deltaphimn = input_obj.pmns_b[-1,:]
        self.bmod_mn = input_obj.bmnc_b[-1,:]
        self.kernelSign = input_obj.kernelSign
        self.phi_b_fact = input_obj.phi_b_fact
        self.Nfp = 3

    def return_values(self, theta_b, phi_b, Nfp):
        self.theta_b_grid, self.phi_b_grid = np.meshgrid(theta_b,phi_b)
        self.R = self.theta_b_grid*0; 
        self.Z = self.theta_b_grid*0; 
        self.B = self.theta_b_grid*0; 
        self.phi = self.theta_b_grid*0;
        for i in range(0,len(self.n)):
            argv = 2.*np.pi*(self.m[i]*self.theta_b_grid + self.n[i]*self.phi_b_grid)
            sinv = np.sin(argv)
            cosv = np.cos(argv)
            self.R+=(self.Rmn[i] * cosv)
            self.Z+=(self.Zmn[i] * sinv)
            self.B+=(self.bmod_mn[i] * cosv)
            self.phi += self.deltaphimn[i]*sinv
        self.phi = (self.phi + self.phi_b_grid)*2.*np.pi/Nfp
        self.x = self.R * np.cos(self.phi)
        self.y = self.R * np.sin(self.phi)
        self.z = self.Z

    def return_values_booz_xform(self, theta_b, phi_b, Nfp, kernel, a):
        '''theta_b and phi_b are 0->2pi
        '''
        print 'hello booz_xform'
        self.theta_b_grid, self.phi_b_grid = np.meshgrid(theta_b,phi_b)
        self.R = self.theta_b_grid*0; 
        self.Z = self.theta_b_grid*0; 
        self.B = self.theta_b_grid*0; 
        self.phi = self.theta_b_grid*0;
        for i in range(len(self.n)):
            #argv = self.m[i]*theta_b + kernel*self.n[i]*phi_b
            #argv = 2.*np.pi*(self.m[i]*self.theta_b_grid + self.n[i]*self.phi_b_grid)
            argv = self.m[i]*self.theta_b_grid + kernel*self.n[i]*self.phi_b_grid
            sinv = np.sin(argv)
            cosv = np.cos(argv)
            self.R+=(self.Rmn[i] * cosv)
            self.Z+=(self.Zmn[i] * sinv)
            self.B+=(self.bmod_mn[i] * cosv)
            self.phi += self.deltaphimn[i]*sinv
        #self.phi = (self.phi + self.phi_b_grid)*2.*np.pi/Nfp
        self.phi = self.phi*a + self.phi_b_grid
        print 'min, max phi : ', np.max(self.phi), np.min(self.phi)
        self.x = self.R * np.cos(self.phi)
        self.y = self.R * np.sin(self.phi)
        self.z = self.Z

    def return_values_single(self, theta_b, phi_b, Nfp, coil_loc = None):
        self.R=0;self.Z=0;self.B=0;self.phi=0
        for i in range(len(self.n)):
            argv = 2.*np.pi*(self.m[i]*theta_b + self.n[i]*phi_b)
            sinv = np.sin(argv)
            cosv = np.cos(argv)
            self.R+=(self.Rmn[i] * cosv)
            self.Z+=(self.Zmn[i] * sinv)
            self.B+=(self.bmod_mn[i] * cosv)
            self.phi += self.deltaphimn[i]*sinv
        self.phi = (self.phi + phi_b)*2.*np.pi/Nfp
        self.x = self.R * np.cos(self.phi)
        self.y = self.R * np.sin(self.phi)
        self.z = self.Z
        if coil_loc!=None:
            self.distance = np.sqrt((coil_loc[0]-self.x)**2+(coil_loc[1]-self.y)**2+(coil_loc[2]-self.z)**2)
            
            return self.distance

    def return_values_single2(self, theta_b, phi_b, Nfp, coil_loc = None, kernel=-1, a=2*np.pi/3.):
        '''This is for BOOZ_XFORM
        '''
        argv = self.m*theta_b + kernel*self.n*phi_b
        #argv = self.m*theta_b + self.n*phi_b
        cosv = np.cos(argv)
        sinv = np.sin(argv)
        self.R = np.sum(self.Rmn * cosv)
        #a = self.phi_b_fact
        #a = np.pi/3.
        #a = 2.*np.pi/3.
        self.phi = phi_b + a*np.sum(self.deltaphimn * sinv)
        self.Z = np.sum(self.Zmn * sinv)
        self.x = self.R * np.cos(self.phi)
        self.y = self.R * np.sin(self.phi)
        self.z = self.Z
        if coil_loc!=None:
            self.distance = np.sqrt((coil_loc[0]-self.x)**2+(coil_loc[1]-self.y)**2+(coil_loc[2]-self.z)**2)
            return self.distance







def get_coil_orientation_data():
    '''
    Extracts data out of the cro files for the coil orientation
    measurements using each of the field coil sets individually to
    give information to figure out the 3 coil orientations on each of
    the formers.

    SRH: 24 July 2014
    '''

    def extract_data(coil_set):
        list_items = ['Amplitude(1)','Amplitude(M)']
        results = {}
        for jjj in range(0,48):
            file_name = '/home/srh112/Desktop/Orientation/{}_{}.txt'.format(coil_set,jjj)
            with file(file_name,'r') as file_opened:contents = file_opened.read().split('\n')
            #file_opened = open(file_name, 'r')
            #contents = file_opened.read().split('\n')
            #file_opened.close()
            answer_list = []
            for iii in range(0,len(list_items)):
                for i in range(0,len(contents)):
                    #print contents[i][0:len(list_items[iii])]
                    if contents[i][0:len(list_items[iii])]==list_items[iii]:
                        temp_line = contents[i]
                        start = temp_line.find('Cur ') + 4
                        end = temp_line[start:].find(',')+start
                        answer = temp_line[start:end]
                        if answer[-2:]=='mV':
                            number_answer = float(answer[:-2])/1000.
                        elif answer[-1]=='V':
                            number_answer = float(answer[:-1])
                        else:
                            print 'error'
                        answer_list.append(number_answer)
            a = np.loadtxt(file_name[:-3]+'csv',skiprows=4,delimiter=',')
            coil_ptp = np.ptp(a[:,1])
            field_current = np.ptp(a[:,2])
            correlation = np.correlate(a[:,1],a[:,2])

            #multiply by 100 to make a resonable answer (I think......)
            coil_ptp = coil_ptp * correlation / np.abs(correlation) *100/field_current
            answer_list.append(float(coil_ptp))
            answer_list.append(float(field_current))
            answer_list.append(float((np.abs(coil_ptp)-np.abs(answer_list[0]))/answer_list[0]*100))
            answer_list.append(float((field_current-answer_list[1])/answer_list[1]*100))
            results[jjj] = answer_list
        return results

    def rearange_results(hfc_results,ovf_results,pfc_results,tfc_results):
        coil_number = 1
        colour = 'blue'
        reordered_results = {}
        reordered_results[coil_number]={}
        reordered_results[coil_number][colour]={}

        for i in range(0,48):
            reordered_results[coil_number][colour]['pfc']=pfc_results[i]
            reordered_results[coil_number][colour]['hfc']=hfc_results[i]
            reordered_results[coil_number][colour]['ovf']=ovf_results[i]
            reordered_results[coil_number][colour]['tfc']=tfc_results[i]
            if colour == 'blue':
                colour = 'black'
                reordered_results[coil_number][colour]={}
            elif colour == 'black':
                colour = 'grey'
                reordered_results[coil_number][colour]={}
            elif colour == 'grey':
                colour = 'blue'
                coil_number+=1
                reordered_results[coil_number]={}
                reordered_results[coil_number][colour]={}
        return reordered_results

    def print_results(coil_set,results):
        print '\n\n',coil_set
        data = np.ones((16,3),dtype=float)

        #Organise data in Blue (x), Black(y), Grey (z)
        for i in range(1,17):
            #blue = results[i]['blue'][coil_set][2]
            #black = results[i]['black'][coil_set][2]
            #grey = results[i]['grey'][coil_set][2]
            #data[i-1,0]=black
            #data[i-1,1]=grey
            #data[i-1,2]=blue
            data[i-1,0]=results[i]['blue'][coil_set][2]
            data[i-1,1]=results[i]['black'][coil_set][2]
            data[i-1,2]=results[i]['grey'][coil_set][2]
        hma_tmp = HMA()
        #This is in Blue (x), Black(y), Grey (z)
        ResistanceComp = hma_tmp.resistance_comp
        #ResistanceComp blue, grey, black
        # ResistanceComp=[[168.6, 189.7212121, 178.2696774],
        #                 [166.180,179.245, 185.496],
        #                 [168.270,180.205,187.516],
        #                 [168.883,180.733,186.883],
        #                 [170.984,183.065,189.144],
        #                 [169.430,181.221,188.471],
        #                 [168.998,180.759,187.958],
        #                 [168.854,181.855,187.862],
        #                 [168.480,181.421,187.413],
        #                 [169.302,189.397,179.217],
        #                 [169.041,186.613,178.879],
        #                 [168.867,189.993,178.726],
        #                 [168.171,189.275,176.955],
        #                 [169.185,186.795,177.806],
        #                 [168.980,187.451,178.684],
        #                 [169.609,188.231,179.217]]
        #ResistanceComp = np.array(ResistanceComp)
        #reorder to make black, grey, blue
        #ResistanceComp_temp = ResistanceComp*1.0
        #ResistanceComp[:,0] = ResistanceComp_temp[:,2]
        #ResistanceComp[:,2] = ResistanceComp_temp[:,0]
        #ResistanceComp = 200./ResistanceComp
        data = data*ResistanceComp
        #for i in range(0,16):
        #    data[i,:]=data[i,:]/(np.sum(data[i,:]**2)**0.5)
        print data
        return data

    pfc_results = extract_data('pfc')
    hfc_results = extract_data('hfc')
    ovf_results = extract_data('ovf')
    tfc_results = extract_data('tfc')
    results = rearange_results(hfc_results,ovf_results,pfc_results,tfc_results)

    data_tfc = print_results('tfc',results)
    data_pfc = print_results('pfc',results)
    data_hfc = print_results('hfc',results)
    data_ovf = print_results('ovf',results)
    return data_tfc, data_pfc, data_hfc, data_ovf

# hma = HMA()
# kh = 0.35
# filename  = '/home/srh112/code/python/h1_eq_generation/results7/kh%.3f-kv1.000fixed/boozmn_wout_kh%.3f-kv1.000fixed.nc'%(kh, kh)
# hma.loc_boozer_coords(filename)
# hma.test_fit()

# pma1 = PMA1()
# kh = 0.6
# filename  = '/home/srh112/code/python/h1_eq_generation/results7/kh%.3f-kv1.000fixed/boozmn_wout_kh%.3f-kv1.000fixed.nc'%(kh, kh)
# pma1.loc_boozer_coords(filename)
# pma1.test_fit()

# pma2 = PMA2()
# kh = 0.6
# filename  = '/home/srh112/code/python/h1_eq_generation/results7/kh%.3f-kv1.000fixed/boozmn_wout_kh%.3f-kv1.000fixed.nc'%(kh, kh)
# pma2.loc_boozer_coords(filename)
# pma2.test_fit()

# fig, ax = pt.subplots(nrows = 2, ncols = 1, sharex = False, sharey = False)
# for i, style in zip([hma, pma1, pma2], ['x','o','d']):
#     ax[0].plot(i.boozer_phi%(360.), i.boozer_theta%(360.), style)
#     ax[1].plot(i.distance, style)
# fig.canvas.draw();fig.show()
