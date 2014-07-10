from StringIO import StringIO
import matplotlib.patches as mpatches
import matplotlib.pyplot as pt
import sys, re,os, copy,datetime, fpformat
import pylab as pl
import numpy as np
import matplotlib as mpl
import shutil
import datetime, os, tempfile, shutil, gzip, subprocess
from multiprocessing import Queue, Process, Pool

class heliac:
    def __init__(self, filename,Nfp):
        self.Nfp = Nfp
        self.filename = filename
        self.extract_heliac_plt_data()

    def extract_heliac_plt_data(self, ):
        # Guess from the filename whether the file has gzip of bzip2 compression

        filename = self.filename
        if filename.endswith('.gz'):
            import gzip
            file = gzip.open(filename)
        elif filename.endswith('.bz2'):
            import bz2
            file = bz2.BZ2File(filename)
        else:
            file = open(filename)

        print 'heliac opened :',filename
        pattern_string  = r"'CPLOT1'.*?\n(?P<row_of_numbers>.*?)\n"

        # The next three lines all have the same form. We want to extract the
        # text from these lines, which are the x and y labels, and the figure
        # title. We're following on from the previous string, so the starting
        # point here is a new line. \s*?' matches any whitespace up to the
        # first single quote ('), we then extract the text, e.g.
        # (?P<x_label>.*?) up to \$?'\n. \$? means zero or one instances of $
        # This is required as some of the labels have a $ before the closing
        # quote
        pattern_string += r"\s*?'(?P<x_label>.*?)\$?'\n"
        pattern_string += r"\s*?'(?P<y_label>.*?)\$?'\n"
        pattern_string += r"\s*?'(?P<title>.*?)\$?'\n"

        # We then grab the next row of numbers, which we probably won't use...
        pattern_string += r"(?P<data_info>.*?)\n"

        # The next lines have the data points, which are terminated by 'END POINTS'
        pattern_string += r"(?P<data>.+?)'END POINTS'"

        # Compiling the pattern makes it more efficient if it will be used multiple
        # times, but there's likely to be no real difference in performace in our
        # case.
        plot_pattern = re.compile(pattern_string, re.DOTALL)

        # This returns a list of all the above patterns in the file
        plt_contents = file.read()
        plot_data_strings = plot_pattern.findall(plt_contents)
        start_point = plt_contents.find('0OUTPUT SUMMARY:'); end_point = plt_contents.find('0MISCELLANEOUS SUMMARIES')
        interesting_text = plt_contents[start_point:end_point]
        interesting_text_list = interesting_text.split('\n')[3:-1]
        out_string = ''
        for i in interesting_text_list[1:]:
            out_string += i[0:16+6]+' '+i[85:85+6]+'\n'
        tmp=np.loadtxt(StringIO(out_string))
        if len(tmp.shape)==1: tmp = tmp[np.newaxis,:]
        self.trace_list_order = map(int, tmp[:,0])
        self.psi_list = tmp[:,1]
        self.R_list = tmp[:,2]
        self.Ra_list = tmp[:,3]
        self.iota_list = -tmp[:,4]

        # And now we're done with the datafile.
        file.close()

        # The datapoints are still in a single string. We need to convert them
        # into numbers. Let's make a dictionary for each plot and put them in a
        # list (plot_data)
        plot_data = []

        for heliac_plot in plot_data_strings:
            plot_dict = {'x_label':heliac_plot[1],
                         'y_label':heliac_plot[2],
                         'title':heliac_plot[3]}
            # split the long datapoint string into a list of strings, splitting where the
            # newlines \n occur. Because the long string ends with a newline, we ignore the
            # last (empty) element of the list (using [:-1])
            datapoint_strings = heliac_plot[5].split('\n')[:-1]
            # now, for each datapoint line, we strip off the whitespace, then split the numbers
            # and convert the strings to float
            plot_dict['data'] = [map(float,i.strip().split()) for i in datapoint_strings]

            # finally, add it to the list of plot dictionaries
            plot_data.append(plot_dict)

        self.plot_data = plot_data
        success = 0
        for i in range(0, len(self.plot_data)):
            if self.plot_data[i]['y_label']=='IOTA-BAR/N':
                tmp = np.array(self.plot_data[i]['data'])
                self.iota_bar_rad = tmp[:,0]
                self.iota_bar = -self.Nfp*tmp[:,1]
                success += 1
            if self.plot_data[i]['y_label']=='PSI (WB)':
                tmp = np.array(self.plot_data[i]['data'])
                self.psi_rad = tmp[1:,0]
                self.psi = tmp[1:,1]
                self.psi_max = np.max(self.psi)
                self.psi_norm = self.psi/np.max(self.psi)
                success += 1
        if success==2:self.iota_poly =np.polyfit(self.psi_norm,self.iota_bar,4)

        

    def read_trace_output2(self, infile):
        '''Read in a trace file and return the ordered surfaces
        Also returns psi as a function of the surfaces
        SH: 8Apr2013
        '''
        with open(infile,'r') as fin:
            file_data = fin.read()
        valid_traces = self.trace_list_order
        self.trace_data = {}
        for curr_trace,curr_psi, curr_R in zip(self.trace_list_order, self.psi_list,self.R_list):
            search_string = 'Start trace %2d\n'%(curr_trace)
            start_point = file_data.find(search_string)+ len(search_string)
            end_point = file_data[start_point:].find('Start trace')
            if (end_point==-1) and (start_point==-1):
                print 'trace %d does not exist!!'%(curr_trace)
                success = 0
            elif (end_point==-1):
                valid_data = file_data[start_point:]
                success = 1
            else:
                valid_data = file_data[start_point:start_point+end_point]#.lstrip(search_string)
                success = 1
            if success:
                self.trace_data[curr_trace] = np.loadtxt(StringIO(valid_data))
                #print curr_trace, self.trace_data[curr_trace][0,0], curr_R, self.trace_data[curr_trace][-1,0], self.trace_data[curr_trace].shape
                if np.abs(self.trace_data[curr_trace][0,0] - curr_R)>0.01:
                    print 'Possibly got the wrong trace!!!!'
        #At this point we have all of the traces in self.trace_data dictionary
        #It is ready to use as an input to descur at this stage
        #However, we want to do more manipulation for the interferometer stuff
        largest_label = np.max(self.trace_data.keys())
        Rtemp = self.trace_data[largest_label][:-1,0]; Ztemp= self.trace_data[largest_label][:-1,1]; PHItemp = self.trace_data[largest_label][:-1,2]
        truth_values = (PHItemp%(2*np.pi/3)<0.01) + ((2*np.pi/3-PHItemp%(2*np.pi/3))<0.01)
        self.trace_npuncs = np.sum(truth_values)
        self.trace_nsurf=Rtemp.shape[0]/self.trace_npuncs

    def plot_heliac_iota_bar(self,ax, Nfp, x_axis='r',label='', plot_style='o-', plot_hlines = True):
        if x_axis=='psi_norm':
            ax.plot(self.psi_norm,self.iota_bar*Nfp,plot_style,label=label)
            ax.plot(self.psi_norm, np.polyval(self.iota_poly*Nfp,self.psi_norm),'s')
            ax.set_xlabel(r'$\Psi_{norm}$')
        elif x_axis=='psi':
            ax.plot(self.psi,self.iota_bar*Nfp,plot_style,label=label)
            ax.set_xlabel(r'$\Psi$')
        elif x_axis=='r':
            ax.plot(self.iota_bar_rad,self.iota_bar * Nfp,plot_style,label=label)
            ax.set_xlabel('<r>')
        if plot_hlines:
            ax.hlines([4./3, 5./4],ax.get_xlim()[0],ax.get_xlim()[1])
    def plot_heliac_psi(self,ax):
        ax.plot(self.psi_rad, self.psi, 'o-')
                

    def puncture_plot(self,plot_ax, desired_phi=100251, desired_num = -10000, include_PFC=1, include_HFC=1, decimation = 1, all_surfaces = 1, plot_dict = None, select_surfaces = None):
        if plot_dict == None:
            plot_dict = {'marker':',','linestyle':'None'}

        print ' heliac puncture plot'
        is_puncture_plot=False
        phi_loc=10000
        num_punc_plots = 0
        count = 0
        while num_punc_plots!=desired_num and phi_loc!=desired_phi and count<=len(self.plot_data):
            heliac_plot = self.plot_data[count]
            if heliac_plot['x_label'] == 'R' and heliac_plot['y_label'] == 'Z':
                is_puncture_plot = True
                phi_loc = float(heliac_plot['title'][heliac_plot['title'].find('=')+1:heliac_plot['title'].find(',')])
                self.phi_loc = phi_loc
                print '##################', desired_phi, phi_loc
                num_punc_plots+=1
                print num_punc_plots
            else:
                is_puncture_plot = False
            count+=1
        print num_punc_plots!=desired_num, phi_loc!=desired_phi, count<=len(self.plot_data)
        print('!!############{},{}'.format(count, len(heliac_plot['data'])))
        if num_punc_plots==desired_num or phi_loc==desired_phi:
            count += 1
            x_data = np.array([heliac_plot['data'][i][0] for i in range(0,len(heliac_plot['data']),all_surfaces)])
            y_data = np.array([heliac_plot['data'][i][1] for i in range(0,len(heliac_plot['data']),all_surfaces)])
            number_good_surfaces = len(self.iota_bar)
            points_per_surface = x_data.shape[0]/number_good_surfaces
            x_data = x_data.reshape((number_good_surfaces,points_per_surface))
            y_data = y_data.reshape((number_good_surfaces,points_per_surface))
            #for i in range(0,len(heliac_plot['data']))
            #y_data = [i[1] for i in heliac_plot['data']]
            current_axis=plot_ax
            if select_surfaces!=None:
                x_data_plot = x_data[select_surfaces,::decimation].T
                y_data_plot = y_data[select_surfaces,::decimation].T
            else:
                x_data_plot = x_data[::all_surfaces,::decimation].T
                y_data_plot = y_data[::all_surfaces,::decimation].T
                
            current_axis.plot(x_data_plot, y_data_plot, **plot_dict)
            ring_conductor=mpatches.Rectangle((0.95,-0.05),0.1,0.1,axes=current_axis,facecolor='red',edgecolor='black')
            if include_PFC:
                current_axis.add_patch(ring_conductor)
            if include_HFC:
                hfc_coords = []
                for i in range(1,5):
                    hfc_coords.append(np.loadtxt('/home/srh112/code/python/heliac/HFC_%d.txt'%(i)))
                for hfc_cur in hfc_coords:
                    hfc_z_values = hfc_cur[:,2]
                    hfc_r_values = (hfc_cur[:,0]**2+hfc_cur[:,1]**2)**0.5
                    hfc_phi_values = np.arctan2(hfc_cur[:,1], hfc_cur[:,0])

                    tmp = np.argmin(np.abs(hfc_phi_values %(2.*np.pi) - np.deg2rad(desired_phi)%(2.*np.pi)))
                    current_axis.plot(hfc_r_values[tmp],hfc_z_values[tmp],'o',markerfacecolor='k', markeredgecolor='k', markersize = mpl.rcParams['lines.markersize']+1)



    def read_trace_output(self, infile, z_values):
        '''Read in a trace file and return the ordered surfaces
        Also returns psi as a function of the surfaces
        SH: 8Apr2013
        '''
        fin = open(infile,'r')
        lines = fin.readlines()
        i=0
        while lines[i].find('POINTS')<0:
            i+=1
        i+=1
        #get the flux as a function of surface
        psi_norm = self.psi_norm
        output_data = {}
        #go through the trace file finding the start of the traces and their numbers
        line_list = []; 
        for i in range(i,len(lines)):
            if lines[i].find('trace')>0: line_list.append([i,int(lines[i].split(' ')[-1])])
        if len(line_list)==1:
            line_list.append([len(lines),-1])

        #R, Z, PHI
        for i in range(0,len(line_list)-1):
            label = line_list[i][1]
            #make sure not to include lines after they have already been included - but in HELIAC?
            if label in output_data.keys():
                print 'label %d aleady counted'%(label)
            else:
                output_data[label] = np.loadtxt(StringIO(''.join(lines[line_list[i][0]+1:line_list[i+1][0]])))

        #Get the valid trace list, get the largest key and make sure it is long enough
        current_key = np.max(output_data.keys())
        while output_data[current_key].shape[0]<1000: current_key-=1
        print os.getpid(), 'descur using largest key:', current_key
        valid_list = [];max_radius_list= []; psi_norm_output_list = []

        #create the list of valid traces, and the corresponding psi norm values
        for j in output_data.keys():
            if j<=current_key:
                valid_list.append(j)
                max_radius_list.append(np.max(output_data[j][:,0]))
                #Not sure why this happens sometimes, but there aren't enough psi_norm values
                if (j)>len(psi_norm):
                    print '!!!past psi_norm boundary.. duplicating last value'
                    psi_norm_output_list.append(psi_norm[-1])
                else:
                    psi_norm_output_list.append(psi_norm[j-1])
        radius_list, valid_list, psi_norm_output_list = zip(*sorted(zip(max_radius_list,valid_list, psi_norm_output_list)))
        overall_result = []; valid_points_R = []; valid_points_Z = []; valid_points_label=[]; valid_points_psi_norm = []
        #Obtain the surface information
        for surf_loc in range(0,len(valid_list)):
            j = valid_list[surf_loc]
            data = output_data[j]
            Rtemp = data[:,0]; Ztemp= data[:,1]; PHItemp=data[:,2]
            truth_values = (PHItemp%(2*np.pi/3)<0.01) + ((2*np.pi/3-PHItemp%(2*np.pi/3))<0.01)
            R2=Rtemp[truth_values];Z2=Ztemp[truth_values];PHI2=PHItemp[truth_values]
            R_order, Z_order, distances = self.order_points_by_closeness(R2.tolist(), Z2.tolist())
            #only keep surfaces that have distances <0.05 between points 
            #this removes close to rational surfaces that behave strangely (not necessary for the new way??)
            if np.max(distances)<0.08:
                valid_points_R.append(R_order)
                valid_points_Z.append(Z_order)
                valid_points_label.append([surf_loc for i in range(len(R_order))])
                valid_points_psi_norm.append([psi_norm_output_list[surf_loc] for i in range(len(R_order))])
                #ax.plot(R_order[-1],Z_order[-1],'-')

        self.ordered_R = valid_points_R
        self.ordered_Z = valid_points_Z
        self.ordered_labels = valid_points_label
        self.valid_points_psi_norm = valid_points_psi_norm
        #return valid_points_R, valid_points_Z, valid_points_label, valid_points_psi_norm

    def order_points_by_closeness(self,R, Z):
        '''This function takes R and Z and orders them
        Note it makes a copy of them so that it doesn't mess up R and Z
        Returns R_order, Z_order which are ordered by closest point
        Also returns distances between neighbouring points
        SH : 8Apr2013
        '''
        R_copy = copy.deepcopy(R)
        Z_copy = copy.deepcopy(Z)
        R_order=[R_copy.pop(0)]
        Z_order=[Z_copy.pop(0)]
        #print len(R_copy), len(Z_copy)
        for i in range(len(R_copy)):
            next_pt = np.argmin((np.array(R_copy)-R_order[i])**2 + (np.array(Z_copy)-Z_order[i])**2)
            R_order.append(R_copy.pop(next_pt))
            Z_order.append(Z_copy.pop(next_pt))
        R_order.append(R_order[0])
        Z_order.append(Z_order[0])

        R_order = np.array(R_order); Z_order = np.array(Z_order)
        distances = np.sqrt((R_order[1:] - R_order[:-1])**2+(Z_order[1:] - Z_order[:-1])**2)
        return R_order, Z_order, distances

    def order_trace_pts_on_surface(self,):
        '''Order all of the trace points on a particular phi surface (currently forced to be 0)
        This is based on closest distance, returns a list of lists for each item that are the points for each surface

        SH: 11Apr2013
        '''
        valid_points_R=[]; valid_points_Z=[]; valid_points_label=[]; valid_points_psi_norm = []
        for curr_trace,curr_iota,curr_psi in zip(self.trace_list_order, self.iota_list, self.psi):#range(0,len(valid_list)):
            data = self.trace_data[curr_trace]
            truth_values = (data[:,2]%(2*np.pi/3)<0.01) + ((2*np.pi/3-data[:,2]%(2*np.pi/3))<0.01)
            R_order, Z_order, distances = self.order_points_by_closeness(data[truth_values,0].tolist(), data[truth_values,1].tolist())
            #print curr_trace, len(R_order), np.sum(truth_values)
            #only keep surfaces that have distances <0.05 between points 
            #this removes close to rational surfaces that behave strangely (not necessary for the new way??)
            if np.max(distances)<0.08:
                valid_points_R.append(R_order)
                valid_points_Z.append(Z_order)
                valid_points_label.append([curr_trace for i in R_order])
                valid_points_psi_norm.append([curr_psi for i in R_order])
                #ax.plot(R_order[-1],Z_order[-1],'-')

        self.ordered_R = valid_points_R
        self.ordered_Z = valid_points_Z
        self.ordered_labels = valid_points_label
        self.valid_points_psi_norm = (np.array(valid_points_psi_norm)/np.max(valid_points_psi_norm)).tolist()

        return valid_points_R, valid_points_Z, valid_points_label, valid_points_psi_norm

    def format_for_descur_input(self,force_trace=None,filename=None):
        '''Output a DESCUR input file based on the heliac trace file
        that has already been read in
        returns a string that should be put into the descur input file

        Use force_trace to pick a trace other than the largest one
        SH: 11Apr2013
        '''
        if force_trace==None:
            current_key = np.max(self.trace_list_order)
            print os.getpid(), 'descur using largest key:', current_key
        else:
            current_key = force_trace
            print os.getpid(), 'descur using a forced key!:', current_key
        data = self.trace_data[current_key][:-1,:]
        Rtemp = data[:,0]; Ztemp= data[:,1]; PHItemp=data[:,2]
        R=[];Z=[];PHI=[]
        count=0
        while count < self.trace_nsurf:
            if count==0:
                truth_values = (PHItemp%(2*np.pi/3)<0.01) + ((2*np.pi/3-PHItemp%(2*np.pi/3))<0.01)
            else:
                truth_values = (np.abs(PHItemp%(2*np.pi/3) - PHItemp[count])<0.01)
            R.extend(Rtemp[truth_values].tolist())
            Z.extend(Ztemp[truth_values].tolist())
            PHI.extend(PHItemp[truth_values].tolist())
            count +=1
        descur_input_string = ''
        descur_input_string += '%1d   %1d   3\n' % (self.trace_npuncs,self.trace_nsurf)
        for i in range(len(R)):
            descur_input_string += '%17.11f %20.11f %17.11f \n' % (R[i],PHI[i],Z[i])
        if filename!=None:
            with file(filename,'w') as f: f.write(descur_input_string)
        return descur_input_string



def iota_dave_fit(r,kappa):
    #Dave paramaterisation
    #Where r is a normalised radius r = r_a / 0.2
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


def generate_heliac_input_files(template_filename, base_output_directory,r0_offset = 0, r1_offset = 0.25, phi_values = [0,30,60,90], kh_values = [0.1,0.2,0.3], n_surfaces = 10, descur_surfaces=24, descur_points_per_surface=100,trace_increment=0.125, template_type = None):
    '''Generate the heliac input file
    SH: 2Apr2013
    '''
    # linear k_h parameterisation of axis location in the phi=0 plane, for k_v=1
    get_axis_r = lambda kappa_h: 1.22646966 + kappa_h*0.024363175
    get_axis_z = lambda kappa_h: -.01046828 + kappa_h*0.007136942

    # how to use these depends on input file (# of coils described) - make better.


    PFC_windings=36
    OVC_windings=8
    HC_windings=4
    current=500000.0/36


    #descur_surfaces = 24
    #descur_points_per_surface = 100
    #trace_increment = 0.125
    trace_save_increment = 360./(descur_surfaces*3)/trace_increment
    trace_iterations = descur_surfaces*descur_points_per_surface*trace_save_increment
    if np.abs(int(trace_save_increment)-trace_save_increment)>0.0001:
        print 'Trace_save_increment is not an integer rethink trace_increment and and'
        raise ValueError
    print 'descur_surfaces: %d, points_per_surface: %d, trace_increment:%.3f, trace iterations:%.2f, trace_save_increment:%.2f'%(descur_surfaces, descur_points_per_surface, trace_increment, trace_save_increment, trace_iterations)
    
    phi_placeholder = "<<phi values>>"
    helical_current_placeholder = "<<helical current>>"
    outer_vertical_current_placeholder = "<<outer vertical current>>"
    kappa_h_placeholder = "<<kappa helical>>"
    kappa_v_placeholder = "<<kappa vertical>>"
    inner_launch_r_placeholder = "<<r0>>"
    outer_launch_r_placeholder = "<<r1>>"
    launch_z_placeholder = "<<z>>"
    number_of_surfaces_placeholder = "<<N>>"
    timestamp_placeholder = "<<timestamp>>"
    ibline_placeholder="<<ibline>>"
    kappa_values = {'h':{},'v':{}}

    ## fix k_v = 1 so we can use the pre-calculated axis r,z
    kappa_values['v']['values'] = [1.0]
    kappa_values['h']['values']=kh_values 
    phi_values_strings=[]
    for i in phi_values: phi_values_strings.append('%.2f'%i)
    phi_values_str = ','.join(phi_values_strings)
    #print phi_values_str
    ibline_value_str = '1'
    output_filename_list = []
    for kappa_h in kappa_values['h']['values']:
        for kappa_v in kappa_values['v']['values']:
            axis_r = get_axis_r(kappa_h); axis_z = get_axis_z(kappa_h)

            launch_r0_str = "%.4f" %(axis_r + r0_offset); launch_r1_str = "%.4f" %(axis_r + r1_offset)
            if axis_z <0:
                sign_str = "-"
            else:
                sign_str = " "
            launch_z_str = sign_str+("%.3f" %(axis_z))[-4:]
            output_file_timestamp = datetime.datetime.now().strftime("%y%m%d")

            kappa_h_string=fpformat.fix(kappa_h,3)
            kappa_v_string=fpformat.fix(kappa_v,3)

            output_directory=base_output_directory + "kh"+kappa_h_string+"-kv"+kappa_v_string+'fixed/'
            os.system('mkdir '+output_directory)
            output_filename="kh"+kappa_h_string+"-kv"+kappa_v_string+'-desc%2d_%d'%(descur_surfaces,descur_points_per_surface)+'-fixed'
            output_filename_list.append([output_directory,output_filename, kappa_h, kappa_v])
            kappa_v_input=fpformat.sci(OVC_windings*current*kappa_v,4) #will -ve stuff up column allignment?

            kappa_h_input=fpformat.sci(-1*current*kappa_h,4)

            if kappa_h_input[0] != '-': kappa_h_input = ' '+kappa_h_input
            if kappa_v_input[0] != '-': kappa_v_input = fpformat.sci(OVC_windings*current*kappa_v,5)

            kappa_h_input=kappa_h_input[:-4]+kappa_h_input[-1]
            kappa_v_input=kappa_v_input[:-4]+kappa_v_input[-2:]

            print 'using template:', template_type
            if template_type == None:
                with open(template_filename,'r') as template_f:template=template_f.read()
            elif template_type == 'simple':
                template = heliac_simple
            elif template_type == 'jason':
                template = heliac_jason
            elif template_type == 'dave':
                template = heliac_comp
            
            template = template.replace(helical_current_placeholder, kappa_h_input)
            template = template.replace(outer_vertical_current_placeholder, kappa_v_input)
            template = template.replace(timestamp_placeholder, str(datetime.datetime.now()))
            template = template.replace(kappa_h_placeholder, kappa_h_string)
            template = template.replace(kappa_v_placeholder, kappa_v_string)
            template = template.replace(inner_launch_r_placeholder, launch_r0_str)
            template = template.replace(outer_launch_r_placeholder, launch_r1_str)
            template = template.replace(launch_z_placeholder, launch_z_str)
            template = template.replace(number_of_surfaces_placeholder, str(n_surfaces))
            template = template.replace(phi_placeholder, phi_values_str)
            template = template.replace(ibline_placeholder,ibline_value_str)
            template = template.replace('<<step_size>>','%.3f'%(trace_increment))
            template = template.replace('<<steps>>','%6d'%(trace_iterations))
            template = template.replace('<<save_inc>>','%2d'%(trace_save_increment))

            #print output_directory + output_filename
            output_file=open(output_directory + output_filename + '.hin','w')
            output_file.write(template)
            output_file.close()
    return output_filename_list


def do_heliac(input_data):
    '''Run heliac based on the data in input_data
    input_data is a list, [0] is the results directory, [1] is the input filename

    The output filenames will be input_filename.out, .plt, .trace
    SH: 2Apr2013
    '''
    results_dir, input_filename, kappa_h, kappa_v = input_data
    print os.getpid(), 'heliac run starting'
    #input_filename = input_data[1]
    
    output_filename = input_filename+'.out'
    plot_filename = input_filename+'.plt'
    trace_filename = input_filename+'.trace'
    run_dir = '/tmp/'
    temp_dir = tempfile.mkdtemp(prefix="heliac-", dir=run_dir)

    # copy heliac binary and input file to working dir
    heliac_location = '/home/srh112/bin/heliac2'
    executable_name = 'heliac2'
    shutil.copy(heliac_location, temp_dir)
    shutil.copy(results_dir + input_filename+'.hin', temp_dir)

    # move to working dir
    initial_dir = os.getcwd()
    os.chdir(temp_dir)

    # create and open HELIAC output file
    output_f = open(output_filename, 'w')
    input_f = open(input_filename+'.hin')

    # put output header in a temp file. Note that we can't just write to
    # output_f because the subprocess will overwrite it when we pass  it
    # to stdout.  Instead, we'll save it to a temp  file and merge it in
    # when we do the gzip compression
    header_filename = 'output_header.txt'
    header_f = open(header_filename,'w')
    header_f.writelines(['Quick interactive heliac run - no NAG routines',
                         str(datetime.datetime.now()),
                         'xheliactiny %(i)s %(o)s %(p)s' %{'i':input_filename,
                                                           'o':output_filename,
                                                           'p':plot_filename},
                     ])
    header_f.close()
    
    os.environ['fu08'] = plot_filename
    #print 'calling the run'
    nice = 19
    subprocess.call(['nice','-%d' %nice, './' + executable_name], stdin=input_f, stdout=output_f)

    output_f.close()
    input_f.close()
    os.rename('fort.8', plot_filename)
    os.rename('fort.10', trace_filename)
    for i in [input_filename + '.hin', output_filename, plot_filename, trace_filename]:
        try:
            shutil.move(i, results_dir)
        except:
            print "EXCEPTION was raised moving files"
            print i, results_dir
    # take the user back to the original working directory
    os.chdir(initial_dir)

    # remove temp folder
    shutil.rmtree(temp_dir)
    print os.getpid(), 'finished heliac'


heliac_simple = r''' &INPUT
  NEQ=2,       IDL=1,      ITERM=0,      IBAUD=110,
  RMIN=0.60,   RMAX=1.6,   ZMIN=-0.6,    ZMAX=0.6,
  RMAG=1.24,   MXITER=8,  DRMIN=1.0E-4,  LASTSF=0, NRDTH=5,
  IFPLOT=0,    IMAG=2,     IBLINE=<<ibline>>,     NAXIS=5,
  JFLAG=0,     MXSAVE=1,   IFSYM=0,      ZMAG=0.0,
  NSV =    0,  1,  3,  3, 3,  3,  6,  6 
  THPIC=<<phi values>>
  MSV =   -1, -1, -1, -2, 0, -3, -2, -3,
  IFFT=0,      IBETA = 0,
  ETA=4.0,     DW=0.1,
  MXITER=15,   DRMIN=1E-6,   NRDTH=15,  NAXIS=5, RMAG=1.2478,
 &END
3.000                                                                         03
New Vertical coil positions
(total) Main current = 500000, Hel curr = <<helical current>>, OV curr = <<outer vertical current>>
Kappa helical = <<kappa helical>>, Kappa (outer) vertical = <<kappa vertical>>
 3.000     160.0     1.000    0.2200                                          05
0.0000E-01 1.000    0.0000                                                    08
0.1389E+06                                                                    01
  1.212078  0.074496  0.039049-90.000000-86.482925  0.383000                  20
  1.162084  0.224804  0.113799-90.000000-79.051422  0.383000                  20
  1.051897  0.371787  0.176181-90.000000-70.534348  0.383000                  20
  0.878499  0.496336  0.205364-90.000000-60.534348  0.383000                  20
  0.669488  0.580923  0.173882-90.000000-49.051422  0.383000                  20
  0.475826  0.643443  0.070527-90.000000-36.482922  0.383000                  20
  0.319325  0.733799 -0.070527-90.000000-23.517076  0.383000                  20
  0.168350  0.870255 -0.173882-90.000000-10.948578  0.383000                  20
 -0.009410  1.008971 -0.205364-90.000000  0.534349  0.383000                  20
 -0.203972  1.096863 -0.176181-90.000000 10.534348  0.383000                  20
 -0.386356  1.118797 -0.113799-90.000000 19.051430  0.383000                  20
 -0.541523  1.086938 -0.039049-90.000000 26.482920  0.383000                  20
 -0.670555  1.012442  0.039049-90.000000 33.517071  0.383000                  20
 -0.775728  0.893993  0.113799-90.000000 40.948578  0.383000                  20
 -0.847925  0.725076  0.176181-90.000000 49.465645  0.383000                  20
 -0.869089  0.512635  0.205364-90.000000 59.465641  0.383000                  20
 -0.837838  0.289332  0.173882-90.000000 70.948586  0.383000                  20
 -0.795151  0.090356  0.070527-90.000000 83.517090  0.383000                  20
 -0.795151 -0.090356 -0.070527 90.000000 83.517067  0.383000                  20
 -0.837838 -0.289332 -0.173882 90.000000 70.948578  0.383000                  20
 -0.869089 -0.512635 -0.205364 90.000000 59.465649  0.383000                  20
 -0.847925 -0.725076 -0.176181 90.000000 49.465652  0.383000                  20
 -0.775728 -0.893993 -0.113799 90.000000 40.948570  0.383000                  20
 -0.670555 -1.012442 -0.039049 90.000000 33.517078  0.383000                  20
 -0.541523 -1.086938  0.039049 90.000000 26.482912  0.383000                  20
 -0.386355 -1.118797  0.113800 90.000000 19.051390  0.383000                  20
 -0.203971 -1.096863  0.176181 90.000000 10.534337  0.383000                  20
 -0.009410 -1.008971  0.205364 90.000000  0.534339  0.383000                  20
  0.168350 -0.870255  0.173882 90.000000-10.948575  0.383000                  20
  0.319325 -0.733798  0.070527 90.000000-23.517109  0.383000                  20
  0.475826 -0.643443 -0.070528 90.000000-36.482964  0.383000                  20
  0.669488 -0.580922 -0.173882 90.000000-49.051441  0.383000                  20
  0.878499 -0.496336 -0.205364 90.000000-60.534328  0.383000                  20
  1.051897 -0.371787 -0.176181 90.000000-70.534355  0.383000                  20
  1.162085 -0.224804 -0.113799 90.000000-79.051437  0.383000                  20
  1.212077 -0.074496 -0.039049 90.000000-86.482925  0.383000                  20
-.5000E+06                                                                    01
    0.0000    0.0000    0.0000  180.0000  180.0000    1.0030                  20
    1.0000    0.1000    0.0000    3.0000    0.0000    3.0000 <<helical current>>  400.0 06
   -0.0000    0.0000     0.000    0.0000      0.00     0.000    0.0000  999.0 06
    1.0000    0.1000    0.0000    3.0000    0.0000    3.0000 <<helical current>>  400.0 06
   -0.0000    0.0000     0.000    0.0000      0.00     0.000    0.0000  999.0 06
    1.0000    0.1000    0.0000    3.0000    0.0000    3.0000 <<helical current>>  400.0 06
   -0.0000    0.0000     0.000    0.0000      0.00     0.000    0.0000  999.0 06
    1.0000    0.1000    0.0000    3.0000    0.0000    3.0000 <<helical current>>  400.0 06
   -0.0000    0.0000     0.000    0.0000      0.00     0.000    0.0000  999.0 06
<<outer vertical current>>                                                                    01
    0.0000    0.0000   -0.7000  180.0000  180.0000    2.1300                  20
<<outer vertical current>>                                                                    01
    0.0000    0.0000    0.7000  180.0000  180.0000    2.1300                  20
2.2222E+05                                                                    01
    0.0000    0.0000   -1.0700  180.0000  180.0000    0.7200                  20
2.2222E+05                                                                    01
    0.0000    0.0000    1.0700  180.0000  180.0000    0.7200                  20
                                                                              09
1    1.264     0.000     0.02   200 1 1     1.35415
1    0.0050    0.000      1.00  3000 320    0.0500 5 V3 quick search
1    0.1300    0.000      0.50 80000 320    0.140010 V3 outer surfaces
1    1.2640    0.0000     0.0101000015    1.354 9
1    1.2640    0.0000     0.0101000015    1.354 9
1    0.0100    0.000      0.25160000 320    0.0700 6 V3 inner surface for iota
1    0.0800    0.000      .125320000 380    0.1200 4 V3 outer surface for iota
(I1,2F10.5,F10.5,I8,I2,I4,F10.5,I3)
0    0.0950    0.000       .125   96000 3  40     0.1200 0 V3 FORMAT 24/100 planes 10cm
'''
#0    <<r0>>    <<z>>     <<step_size>><<steps>> 3<<save_inc>>    <<r1>><<N>> gHg <<timestamp>>

heliac_comp = r''' &INPUT
  NEQ=2,       IDL=1,      ITERM=0,      IBAUD=110,
  RMIN=0.60,   RMAX=1.6,   ZMIN=-0.6,    ZMAX=0.6,
  RMAG=1.24,   MXITER=8,  DRMIN=1.0E-4,  LASTSF=0, NRDTH=5,
  IFPLOT=0,    IMAG=2,     IBLINE=<<ibline>>,     NAXIS=5,
  JFLAG=0,     MXSAVE=1,   IFSYM=0,      ZMAG=0.0,
  NSV =    0,  1,  3,  3, 3,  3,  6,  6, 
  THPIC=<<phi values>>
  MSV =   -1, -1, -1, -2, 0, -3, -2, -3,
  IFFT=0,      IBETA = 0,
  ETA=4.0,     DW=0.1,
  MXITER=4,   DRMIN=1E-9,   NRDTH=2,  NAXIS=5, 
 &END
3.000                                                                         03
New Vertical coil positions
(total) Main current = 500000, Hel curr = <<helical current>>, OV curr = <<outer vertical current>>
Kappa helical = <<kappa helical>>, Kappa (outer) vertical = <<kappa vertical>>
 1.000     480.0     1.000    0.2200                                          05
 9.500E-03 1.000    0.0000                                                    08
 -13888.89                                                                     2
  0.866680  0.611600  0.027750                                                99
  0.886540  0.558300  0.027750                                                99
  0.846180  0.615500  0.027750                                                99
  0.871430  0.547640  0.027750                                                99
  0.831070  0.604830  0.027750                                                99
  0.856310  0.536970  0.027750                                                99
  0.815950  0.594160  0.027750                                                99
  0.841200  0.526300  0.027750                                                99
  0.800840  0.583500  0.027750                                                99
  0.826080  0.515640  0.027750                                                99
  0.785720  0.572830  0.027750                                                99
  0.810970  0.504970  0.027750                                                99
  0.770610  0.562160  0.027750                                                99
  0.776020  0.547630  0.027750                                                99
  0.776020  0.547630  0.046250                                                99
  0.787900  0.537650  0.046250                                                99
  0.747540  0.594850  0.046250                                                99
  0.803020  0.548320  0.046250                                                99
  0.762660  0.605510  0.046250                                                99
  0.818140  0.558990  0.046250                                                99
  0.777770  0.616180  0.046250                                                99
  0.833250  0.569650  0.046250                                                99
  0.792890  0.626850  0.046250                                                99
  0.848370  0.580320  0.046250                                                99
  0.808010  0.637510  0.046250                                                99
  0.863480  0.590990  0.046250                                                99
  0.823120  0.648180  0.046250                                                99
  0.866680  0.611600  0.046250                                                99
  0.866680  0.611600  0.027750                                                99
 -13888.89                                                                     2
 -0.963000  0.444770 -0.009250                                                99
 -0.926780  0.488620 -0.009250                                                99
 -0.956130  0.425070 -0.009250                                                99
 -0.909980  0.480860 -0.009250                                                99
 -0.939330  0.417310 -0.009250                                                99
 -0.893190  0.473100 -0.009250                                                99
 -0.922540  0.409550 -0.009250                                                99
 -0.876390  0.465350 -0.009250                                                99
 -0.905740  0.401800 -0.009250                                                99
 -0.859600  0.457590 -0.009250                                                99
 -0.888950  0.394040 -0.009250                                                99
 -0.842800  0.449830 -0.009250                                                99
 -0.872150  0.386280 -0.009250                                                99
 -0.862270  0.398240 -0.009250                                                99
 -0.862270  0.398240  0.009250                                                99
 -0.859570  0.413520  0.009250                                                99
 -0.888920  0.349970  0.009250                                                99
 -0.876370  0.421280  0.009250                                                99
 -0.905720  0.357730  0.009250                                                99
 -0.893160  0.429030  0.009250                                                99
 -0.922510  0.365480  0.009250                                                99
 -0.909960  0.436790  0.009250                                                99
 -0.939310  0.373240  0.009250                                                99
 -0.926750  0.444550  0.009250                                                99
 -0.956100  0.381000  0.009250                                                99
 -0.943550  0.452300  0.009250                                                99
 -0.972900  0.388750  0.009250                                                99
 -0.963000  0.444770  0.009250                                                99
 -0.963000  0.444770 -0.009250                                                99
 -13888.89                                                                     2
  0.096320 -1.056370 -0.046250                                                99
  0.040230 -1.046920 -0.046250                                                99
  0.109940 -1.040570 -0.046250                                                99
  0.038550 -1.028500 -0.046250                                                99
  0.108260 -1.022140 -0.046250                                                99
  0.036870 -1.010070 -0.046250                                                99
  0.106580 -1.003720 -0.046250                                                99
  0.035190 -0.991650 -0.046250                                                99
  0.104900 -0.985290 -0.046250                                                99
  0.033510 -0.973230 -0.046250                                                99
  0.103220 -0.966870 -0.046250                                                99
  0.031830 -0.954800 -0.046250                                                99
  0.101540 -0.948450 -0.046250                                                99
  0.086250 -0.945870 -0.046250                                                99
  0.086250 -0.945870 -0.027750                                                99
  0.071670 -0.951170 -0.027750                                                99
  0.141380 -0.944820 -0.027750                                                99
  0.073350 -0.969600 -0.027750                                                99
  0.143060 -0.963240 -0.027750                                                99
  0.075030 -0.988020 -0.027750                                                99
  0.144740 -0.981660 -0.027750                                                99
  0.076710 -1.006440 -0.027750                                                99
  0.146420 -1.000090 -0.027750                                                99
  0.078390 -1.024870 -0.027750                                                99
  0.148100 -1.018510 -0.027750                                                99
  0.080070 -1.043290 -0.027750                                                99
  0.149780 -1.036930 -0.027750                                                99
  0.096320 -1.056370 -0.027750                                                99
  0.096320 -1.056370 -0.046250                                                99
0.1389E+06                                                                    01
  1.210616  0.074405  0.039700-90.000000-86.482994  0.383000                  20
  1.159895  0.224389  0.111500-90.000000-79.050995  0.383000                  20
  1.050229  0.371204  0.174000-90.000000-70.533997  0.383000                  20
  0.877091  0.495547  0.204500-90.000000-60.533997  0.383000                  20
  0.666271  0.588410  0.166500-90.000000-48.550999  0.383000                  20
  0.475905  0.643548  0.067200-90.000000-36.482998  0.383000                  20
  0.319576  0.734379 -0.074000-90.000000-23.517000  0.383000                  20
  0.167694  0.866828 -0.176000-90.000000-10.948996  0.383000                  20
 -0.009389  1.007356 -0.206500-90.000000  0.533993  0.383000                  20
 -0.203642  1.095127 -0.177000-90.000000 10.533997  0.383000                  20
 -0.386273  1.118583 -0.115000-90.000000 19.050997  0.383000                  20
 -0.525719  1.090269 -0.047000-90.000000 25.742990  0.383000                  20
 -0.669468  1.010805  0.038500-90.000000 33.516987  0.383000                  20
 -0.769832  0.902983  0.109500-90.000000 40.449001  0.383000                  20
 -0.845447  0.722948  0.174500-90.000000 49.466000  0.383000                  20
 -0.869596  0.512927  0.203500-90.000000 59.465992  0.383000                  20
 -0.841766  0.282484  0.164500-90.000000 71.448997  0.383000                  20
 -0.795282  0.090372  0.069000-90.000000 83.516991  0.383000                  20
 -0.798461 -0.090733 -0.071500 90.000000 83.516991  0.383000                  20
 -0.835960 -0.288677 -0.175000 90.000000 70.949013  0.383000                  20
 -0.865548 -0.510539 -0.208000 90.000000 59.466000  0.383000                  20
 -0.846207 -0.723598 -0.177000 90.000000 49.466000  0.383000                  20
 -0.773947 -0.891926 -0.115000 90.000000 40.949005  0.383000                  20
 -0.674262 -1.007011 -0.043000 90.000000 33.804981  0.383000                  20
 -0.540737 -1.085358  0.039000 90.000000 26.482990  0.383000                  20
 -0.385621 -1.116693  0.115500 90.000000 19.051008  0.383000                  20
 -0.203642 -1.095127  0.175500 90.000000 10.534008  0.383000                  20
 -0.009391 -1.007556  0.204500 90.000000  0.534011  0.383000                  20
  0.176344 -0.870722  0.163500 90.000000-11.449009  0.383000                  20
  0.320145 -0.733585  0.067500 90.000000-23.576992  0.383000                  20
  0.475311 -0.642744 -0.074000 90.000000-36.482990  0.383000                  20
  0.669492 -0.580935 -0.173500 90.000000-49.051003  0.383000                  20
  0.876394 -0.495154 -0.207500 90.000000-60.533997  0.383000                  20
  1.054943 -0.372871 -0.177000 90.000000-70.533997  0.383000                  20
  1.160091 -0.224427 -0.115500 90.000000-79.050995  0.383000                  20
  1.211115 -0.074435 -0.043000 90.000000-86.483002  0.383000                  20
-13.8888e3                                                                    01
   -0.0015    0.0020   -0.0458  180.0000  180.0000    0.9572                  20
   -0.0015    0.0020   -0.0275  180.0000  180.0000    0.9572                  20
   -0.0015    0.0020   -0.0092  180.0000  180.0000    0.9572                  20
   -0.0015    0.0020    0.0092  180.0000  180.0000    0.9572                  20
   -0.0015    0.0020    0.0275  180.0000  180.0000    0.9572                  20
   -0.0015    0.0020    0.0458  180.0000  180.0000    0.9572                  20
   -0.0015    0.0020   -0.0458  180.0000  180.0000    0.9755                  20
   -0.0015    0.0020   -0.0275  180.0000  180.0000    0.9755                  20
   -0.0015    0.0020   -0.0092  180.0000  180.0000    0.9755                  20
   -0.0015    0.0020    0.0092  180.0000  180.0000    0.9755                  20
   -0.0015    0.0020    0.0275  180.0000  180.0000    0.9755                  20
   -0.0015    0.0020    0.0458  180.0000  180.0000    0.9755                  20
   -0.0015    0.0020   -0.0458  180.0000  180.0000    0.9938                  20
   -0.0015    0.0020   -0.0275  180.0000  180.0000    0.9938                  20
   -0.0015    0.0020   -0.0092  180.0000  180.0000    0.9938                  20
   -0.0015    0.0020    0.0092  180.0000  180.0000    0.9938                  20
   -0.0015    0.0020    0.0275  180.0000  180.0000    0.9938                  20
   -0.0015    0.0020    0.0458  180.0000  180.0000    0.9938                  20
   -0.0015    0.0020   -0.0458  180.0000  180.0000    1.0122                  20
   -0.0015    0.0020   -0.0275  180.0000  180.0000    1.0122                  20
   -0.0015    0.0020   -0.0092  180.0000  180.0000    1.0122                  20
   -0.0015    0.0020    0.0092  180.0000  180.0000    1.0122                  20
   -0.0015    0.0020    0.0275  180.0000  180.0000    1.0122                  20
   -0.0015    0.0020    0.0458  180.0000  180.0000    1.0122                  20
   -0.0015    0.0020   -0.0458  180.0000  180.0000    1.0305                  20
   -0.0015    0.0020   -0.0275  180.0000  180.0000    1.0305                  20
   -0.0015    0.0020   -0.0092  180.0000  180.0000    1.0305                  20
   -0.0015    0.0020    0.0092  180.0000  180.0000    1.0305                  20
   -0.0015    0.0020    0.0275  180.0000  180.0000    1.0305                  20
   -0.0015    0.0020    0.0458  180.0000  180.0000    1.0305                  20
   -0.0015    0.0020   -0.0458  180.0000  180.0000    1.0488                  20
   -0.0015    0.0020   -0.0275  180.0000  180.0000    1.0488                  20
   -0.0015    0.0020   -0.0092  180.0000  180.0000    1.0488                  20
   -0.0015    0.0020    0.0092  180.0000  180.0000    1.0488                  20
   -0.0015    0.0020    0.0275  180.0000  180.0000    1.0488                  20
   -0.0015    0.0020    0.0458  180.0000  180.0000    1.0488                  20
<<outer vertical current>>                                                                    01
    0.0000    0.0000   -0.6949  180.0000  180.0000    2.1300                  20
    0.0000    0.0000    0.7622  180.0000  180.0000    2.1300                  20
2.2222E+05                                                                    01
    0.0000    0.0000   -1.0700  180.0000  180.0000    0.7200                  20
    0.0000    0.0000    1.0700  180.0000  180.0000    0.7200                  20
    1.0080    0.0910    -2.700    3.0000    0.0000    3.0000 <<helical current>>  400.0 06
   -0.0015    0.0020     0.000    0.0000      0.00     0.000    0.0000  999.0 06
    1.0080    0.0910    -0.900    3.0000    0.0000    3.0000 <<helical current>>  400.0 06
   -0.0015    0.0020     0.000    0.0000      0.00     0.000    0.0000  999.0 06
    1.0080    0.0910     0.900    3.0000    0.0000    3.0000 <<helical current>>  400.0 06
   -0.0015    0.0020     0.000    0.0000      0.00     0.000    0.0000  999.0 06
    1.0080    0.0910     2.700    3.0000    0.0000    3.0000 <<helical current>>  400.0 06
   -0.0015    0.0020     0.000    0.0000      0.00     0.000    0.0000  999.0 06
                                                                              09
1    1.2500    -.002     0.250360000 320    1.360010 original for wtover
1    1.0460    -.230     0.250360000 320    1.240010 for 85.0
1    1.1600    -.213     0.250360000 320    1.200010 for 85.0
1    0.9904    -.228     0.250360000 320    1.145005 8 turns in vertical full surfaces
1    1.1450    -.228     0.250360000 320    1.165005 for std matching 04060710.png
1    1.1350    -.213     0.250360000 320    1.155010 1st surf matching with 05120516.png
1    1.1650    -.213     0.250360000 320    1.185010 3rd surf matching with 05120501.png
1    1.1626    -.180     0.250720000 320    1.189010 island
1    1.1600    -.213     0.250360000 320    1.192010 for 85.0
1    1.1550    -.230     0.250360000 320    1.155001 for 85.0
1    1.3200    -.005     0.250360000 320    1.369010 original for wtover
1    1.1250    -.230     0.250180000 320    1.190010 for testing
1    1.0001    -.230     0.250360000 320    1.130050 tas for 85.0 launch inner surfaces
1    1.1320    -.236     0.250360000 320    1.148010 taken from old matching one
1    1.1350    -.236     0.250360000 320    1.138003 taken from old matching one
1    1.1380    -.236     0.250360000 320    1.138001 taken from old matching one
1    1.1370    -.236     0.250180000 320    1.137001 single surf
1    1.2000    -.236     0.250360000 320    1.460052 dgp 060403
(I1,2F10.5,F10.5,I8,I2,I4,F10.5,I3)
0    0.0950    0.000       .125   96000 3  40     0.1200 0 V3 FORMAT 24/100 planes 10cm
'''
#0    <<r0>>    <<z>>     <<step_size>><<steps>> 3<<save_inc>>    <<r1>><<N>> gHg <<timestamp>>
#0    <<r0>>    <<z>>     0.250360000 320    <<r1>><<N>> gHg <<timestamp>>


heliac_jason = r''' &INPUT
  NEQ=2,       IDL=1,      ITERM=0,      IBAUD=110,
  RMIN=0.60,   RMAX=1.6,   ZMIN=-0.6,    ZMAX=0.6,
  RMAG=1.24,   MXITER=8,  DRMIN=1.0E-4,  LASTSF=0, NRDTH=5,
  IFPLOT=0,    IMAG=2,     IBLINE=0,     NAXIS=5,
  JFLAG=0,     MXSAVE=1,   IFSYM=0,      ZMAG=0.0,
  NSV =    0,  1,  3,  3, 3,  3,  6,  6, THPIC=0, 30, 60, 120
  MSV =   -1, -1, -1, -2, 0, -3, -2, -3,
  IFFT=0,      IBETA = 0,
  ETA=4.0,     DW=0.1,
  MXITER=15,   DRMIN=1E-6,   NRDTH=15,  NAXIS=5, 
  IBLINE=1
 &END
3.000                                                                         03
H1 std case 0.43 helical sm steph1ase000c 1.003C,215,.2068 0% 3 0 3 May  93 16/8
  alpha=  0.3000 beta=  0.0000 gamma=  0.0000 delta=  0.0000 epsilon=  0.0000   
  Ro=  1.0030 Rho=  0.2150 FracOf =   0.5000 wedgea =   0.0000                  
 3.000     160.0     1.000    0.2200                                          05
0.0000E-01 1.000    0.0000                                                    08
0.1389E+06                                                                    01
  1.212078  0.074496  0.039049-90.000000-86.482925  0.383000                  20
  1.162084  0.224804  0.113799-90.000000-79.051422  0.383000                  20
  1.051897  0.371787  0.176181-90.000000-70.534348  0.383000                  20
  0.878499  0.496336  0.205364-90.000000-60.534348  0.383000                  20
  0.669488  0.580923  0.173882-90.000000-49.051422  0.383000                  20
  0.475826  0.643443  0.070527-90.000000-36.482922  0.383000                  20
  0.319325  0.733799 -0.070527-90.000000-23.517076  0.383000                  20
  0.168350  0.870255 -0.173882-90.000000-10.948578  0.383000                  20
 -0.009410  1.008971 -0.205364-90.000000  0.534349  0.383000                  20
 -0.203972  1.096863 -0.176181-90.000000 10.534348  0.383000                  20
 -0.386356  1.118797 -0.113799-90.000000 19.051430  0.383000                  20
 -0.541523  1.086938 -0.039049-90.000000 26.482920  0.383000                  20
 -0.670555  1.012442  0.039049-90.000000 33.517071  0.383000                  20
 -0.775728  0.893993  0.113799-90.000000 40.948578  0.383000                  20
 -0.847925  0.725076  0.176181-90.000000 49.465645  0.383000                  20
 -0.869089  0.512635  0.205364-90.000000 59.465641  0.383000                  20
 -0.837838  0.289332  0.173882-90.000000 70.948586  0.383000                  20
 -0.795151  0.090356  0.070527-90.000000 83.517090  0.383000                  20
 -0.795151 -0.090356 -0.070527 90.000000 83.517067  0.383000                  20
 -0.837838 -0.289332 -0.173882 90.000000 70.948578  0.383000                  20
 -0.869089 -0.512635 -0.205364 90.000000 59.465649  0.383000                  20
 -0.847925 -0.725076 -0.176181 90.000000 49.465652  0.383000                  20
 -0.775728 -0.893993 -0.113799 90.000000 40.948570  0.383000                  20
 -0.670555 -1.012442 -0.039049 90.000000 33.517078  0.383000                  20
 -0.541523 -1.086938  0.039049 90.000000 26.482912  0.383000                  20
 -0.386355 -1.118797  0.113800 90.000000 19.051390  0.383000                  20
 -0.203971 -1.096863  0.176181 90.000000 10.534337  0.383000                  20
 -0.009410 -1.008971  0.205364 90.000000  0.534339  0.383000                  20
  0.168350 -0.870255  0.173882 90.000000-10.948575  0.383000                  20
  0.319325 -0.733798  0.070527 90.000000-23.517109  0.383000                  20
  0.475826 -0.643443 -0.070528 90.000000-36.482964  0.383000                  20
  0.669488 -0.580922 -0.173882 90.000000-49.051441  0.383000                  20
  0.878499 -0.496336 -0.205364 90.000000-60.534328  0.383000                  20
  1.051897 -0.371787 -0.176181 90.000000-70.534355  0.383000                  20
  1.162085 -0.224804 -0.113799 90.000000-79.051437  0.383000                  20
  1.212077 -0.074496 -0.039049 90.000000-86.482925  0.383000                  20
-.5000E+06                                                                    01
    0.0000    0.0000    0.0000  180.0000  180.0000    1.0030                  20
    1.0000    0.1000    0.0000    3.0000    0.0000    3.0000 -17.222E3  400.0 06
1.1111E+05                                                                    01
    0.0000    0.0000   -0.7000  180.0000  180.0000    2.1300                  20
1.1111E+05                                                                    01
    0.0000    0.0000    0.7000  180.0000  180.0000    2.1300                  20
2.2222E+05                                                                    01
    0.0000    0.0000   -1.0700  180.0000  180.0000    0.7200                  20
2.2222E+05                                                                    01
    0.0000    0.0000    1.0700  180.0000  180.0000    0.7200                  20
                                                                              09
1    1.264     0.000     0.02   200 1 1     1.35415
1    0.0050    0.000      1.00  3000 320    0.0500 5 V3 quick search
1    0.1300    0.000      0.50 80000 320    0.140010 V3 outer surfaces
1    1.2640    0.0000     0.0101000015    1.354 9
1    1.2640    0.0000     0.0101000015    1.354 9
1    0.0100    0.000      0.25160000 320    0.0700 6 V3 inner surface for iota
1    0.0800    0.000      .125320000 380    0.1200 4 V3 outer surface for iota
(I1,2F10.5,F10.5,I8,I2,I4,F10.5,I3)
1    1.1159    0.150       0.25 4000000 3 360    1.1163  4 V3 microdetail 7/6, in it
1     1.148    0.1750      0.50 1600000 3  40     1.150 10 8/7 hi RES
1    0.0020    0.000       .25    48000 3  40     0.1200 0 V3 FORMAT check iota near axis
1    0.0200    0.000       .125   96000 3  80     0.1200 0 V3 FORMAT 12 planes 100punc inner
1    0.0950    0.000       .125   96000 3  80     0.1200 0 V3 FORMAT 12/100 planes 10cm
0    0.0950    0.000       .125   96000 3  40     0.1200 0 V3 FORMAT 24/100 planes 10cm
1    0.0950    0.000       .125  192000 3  40     0.1200 2 V3 FORMAT 24/200 planes 8cm
1    0.1000    0.000       .125  144000 3 240     0.1200 0 V3 Special FORMAT 12 planes outer
1    0.1000    0.000       .125  144000 3 120     0.1200 0 V3 FORMAT 800 sec IBM .5cm steps
1    0.1000    0.000       .125  144000 3 160     0.1200 0 V3 FORMAT 800 sec IBM .5cm steps
'''
