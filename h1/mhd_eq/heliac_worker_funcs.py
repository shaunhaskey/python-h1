from StringIO import StringIO
import matplotlib.patches as mpatches
import matplotlib.pyplot as pt
import sys, re,os, copy
import pylab as pl
import numpy as np
import matplotlib as mpl
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
        for i in range(0, len(self.plot_data)):
            if self.plot_data[i]['y_label']=='IOTA-BAR/N':
                tmp = np.array(self.plot_data[i]['data'])
                self.iota_bar_rad = tmp[:,0]
                self.iota_bar = -self.Nfp*tmp[:,1]
            if self.plot_data[i]['y_label']=='PSI (WB)':
                tmp = np.array(self.plot_data[i]['data'])
                self.psi_rad = tmp[1:,0]
                self.psi = tmp[1:,1]
                self.psi_max = np.max(self.psi)
                self.psi_norm = self.psi/np.max(self.psi)

        
        self.iota_poly =np.polyfit(self.psi_norm,self.iota_bar,4)

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

    def plot_heliac_iota_bar(self,ax, Nfp, x_axis='r',label='', plot_style='o-'):
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


