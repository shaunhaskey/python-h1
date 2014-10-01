#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as pt
import matplotlib.mlab as mlab
import MDSplus as MDS
import pyfusion as pf
import numpy as np
import time, os, copy
import cPickle as pickle
def extract_ne_kh(current_shot, ne_array = 0):
    count,success=(0,0)
    while count<6 and success==0:
        try:
            MDSTree=MDS.Tree('h1data', current_shot)
            main_current=MDSTree.getNode('.operations.magnetsupply.lcu.setup_main:I2').data()
            sec_current=MDSTree.getNode('.operations.magnetsupply.lcu.setup_sec:I2').data()
            try:
                heating_freq = MDSTree.getNode('.log.heating:snmp:t2:measured:frequency').data()
            except:
                print 'couldnt get heating frequency'
                heating_freq = None
            kh=float(sec_current)/float(main_current)
            electr_dens_tree = MDS.Tree('electr_dens',current_shot)
            ne_node=electr_dens_tree.getNode('NE_HET:NE_CENTRE')
            ne_value=ne_node.record.data()
            ne_time=ne_node.record.dim_of().data()
            ne_value_list = []; ne_time_list = []
            if ne_array == 1:
                for i in range(1,8):
                    ne_node_tmp = electr_dens_tree.getNode('NE_HET:NE_%d'%(i))
                    ne_value_list.append(ne_node_tmp.record.data())
                    ne_time_list.append(ne_node_tmp.record.dim_of().data())
            success = 1
        except:
            kh=None
            ne_value=None
            ne_time=None
            main_current = None
            heating_freq = None
            print 'Error obtaining kh and ne values'
            count=count+1
            time.sleep(0.5)
    return kh, ne_value, ne_time, main_current, heating_freq, ne_value_list, ne_time_list

def extract_data(current_shot,array):
    tries,success=(0,0)
    while tries<10 and success==0:
        try:
            data=pf.getDevice('H1').acq.getdata(current_shot,array)
            success=1
            #print 'Data extracted on Shot : %d'%(current_shot)
        except (MDS.TdiException, MDS.TreeException) as e:
            print current_shot, e
            tries=tries+1
            time.sleep(0.5)
            data=None
    return data

def extract_polarisation_data(current_shot):
    try:
        #MDSTree=MDS.Tree('mirnov',current_shot)
        coil_1x=pf.getDevice('H1').acq.getdata(current_shot,'H1ToroidalMirnov_1x')
        coil_1y=pf.getDevice('H1').acq.getdata(current_shot,'H1ToroidalMirnov_1y')
        coil_1z=pf.getDevice('H1').acq.getdata(current_shot,'H1ToroidalMirnov_1z')
        #Need to narrow the time down.... maybe extract the above data elsewhere
        #print 'successful extraction of pyrex coil'
        return coil_1x, coil_1y, coil_1z
    except:
        print 'Error getting polarisation data'
        return 0,0,0

def do_a_single2(include_amplifier_sig, include_ne, include_ne_spec, include_flucstrucs, kh, overall_shot_grouping, channel, fft_length, labels, plot_style, NFFT_ne, base_dir):
    nplots = include_amplifier_sig + include_ne + include_ne_spec + 1
    ncols = len(overall_shot_grouping)
    fig, ax = pt.subplots(nrows = nplots, ncols = ncols, sharex=1)#nrows = 2, sharex = 1)
    if nplots == 1 and ncols == 1:
        ax = [ax]
    xlim_list = []
    title_string = ''
    kh = kh/100.
    if include_flucstrucs:
        flucstruc_fig, flucstruc_ax = pt.subplots(nrows=2, sharex = 1)


def do_a_single(include_amplifier_sig, include_ne, include_ne_spec, include_flucstrucs, kh, overall_shot_grouping, channel, fft_length, labels, plot_style, NFFT_ne, base_dir):
    nplots = include_amplifier_sig + include_ne + include_ne_spec + 1
    ncols = len(overall_shot_grouping)
    fig, ax = pt.subplots(nrows = nplots, ncols = ncols, sharex=1)#nrows = 2, sharex = 1)
    if nplots == 1 and ncols == 1:
        ax = [ax]
    xlim_list = []
    title_string = ''
    kh = kh/100.
    if include_flucstrucs:
        flucstruc_fig, flucstruc_ax = pt.subplots(nrows=2, sharex = 1)
    for plot_loc, shot in enumerate(overall_shot_grouping):
        print shot
        #tree_name = MDS.Tree(tree, shot)
        tree_name = MDS.Tree('mirnov', shot)
        node = tree_name.getNode(channel)
        #node = h1tree.getNode(channel)
        data_raw = node.record.data()
        time_raw = node.dim_of().data()
        xlim_list.append([np.min(time_raw),np.max(time_raw)])
        h1tree = MDS.Tree('h1data', shot)
        main_current=h1tree.getNode('.operations.magnetsupply.lcu.setup_main:I2').data()
        sec_current=h1tree.getNode('.operations.magnetsupply.lcu.setup_sec:I2').data()
        if include_amplifier_sig:
            current_node = tree_name.getNode('ACQ132_8:input_32')
            current_raw = current_node.record.data()

        clr_fig = ax[0, plot_loc].specgram(data_raw.flatten(), NFFT=fft_length, Fs=2e6/1000, window=mlab.window_hanning, noverlap=int(fft_length*15./16.),cmap='jet',xextent=[np.min(time_raw),np.max(time_raw)])

        if include_amplifier_sig:
            clr_fig2 = ax[1, plot_loc].specgram(current_raw.flatten(), NFFT=fft_length, Fs=2e6, window=mlab.window_hanning, noverlap=int(fft_length*15./16.),cmap='jet',xextent=[np.min(time_raw),np.max(time_raw)])
        clr_fig[3].set_clim([-100,20])
        clr_fig[3].set_clim([-90,60])


        if include_ne or include_ne_spec:
            nenode = h1tree.getNode('.electr_dens:ne_het:ne_centre')
            ne = nenode.data()
            ne_time = nenode.dim_of().data()
            nep = ax[1,plot_loc]

        if include_ne or include_ne_spec:
            nep.plot(nenode.dim_of().data(), ne)
            try:
                ne_edge = h1tree.getNode('.electr_dens:ne_het:ne_7')
                nep.plot(ne_edge.dim_of().data(),ne_edge.data())
            except:
                pass
            nep.set_ylim([0, 2.5])

        if include_ne_spec:
            nesp = ax[2,plot_loc]
            nesp.specgram(ne.flatten(), NFFT=NFFT_ne, Fs=1e-3/(ne_time[1]-ne_time[0]), window=mlab.window_hanning, noverlap=int(NFFT_ne*15./16.),cmap='jet',xextent=[np.min(ne_time),np.max(ne_time)])
            nesp.set_ylim([0, 50])
            nesp.set_xlim([np.min(time_raw),np.max(time_raw)]) # use the mirnov - it is shorter

        if include_flucstrucs:
            serial_number = 0; 
            fs_dictionary, serial_number, success = single_shot_fluc_strucs(shot, 'H1ToroidalAxial', [0,0.08], 2048, power_cutoff = 0.05, n_svs = 2)#fft_length)
            tmp_time_list = []; tmp_freq_list = []; tmp_ne_list = []#ne_list = []
            for serial in fs_dictionary.keys():
                tmp_time_list.append(copy.deepcopy(fs_dictionary[serial]['time']/1000))
                tmp_freq_list.append(copy.deepcopy(fs_dictionary[serial]['freq']))
                tmp_ne_list.append(copy.deepcopy(fs_dictionary[serial]['ne']))
            del fs_dictionary
            ax[0, plot_loc].plot(tmp_time_list, np.array(tmp_freq_list)/1000., 'k,')
            flucstruc_ax[0].plot(tmp_time_list, np.array(tmp_freq_list)/1000.,plot_style[plot_loc],label=labels[plot_loc])
            flucstruc_ax[1].plot(tmp_time_list, np.array(tmp_freq_list)/1000.*np.sqrt(tmp_ne_list),plot_style[plot_loc],label=labels[plot_loc])
            del tmp_time_list, tmp_freq_list, tmp_ne_list
        
        title_string+='%s_'%(shot)

        ax[0, plot_loc].set_title(r'%d,$\kappa_H$:%.2f,$I_m$:%d'%(shot, kh, main_current), fontsize=6)
        ax[0, plot_loc].set_ylim([0, 150])
        #ax[0, plot_loc].set_xlabel('Time (ms)')
        #ax[0, plot_loc].set_ylabel('Frequency (kHz)')
    

    
    title_string = title_string.rstrip('_')+'.png'
    #file_name = '%s/%d_%d.png'%(base_dir, shot_1, shot_2)
    file_name = base_dir + '/' + title_string #'%s/%d_%d.png'%(base_dir, shot_1, shot_2)
    xaxis_min = np.min(xlim_list[0][0], xlim_list[1][0])
    xaxis_max = np.max(xlim_list[0][1], xlim_list[1][1])
    ax[0, plot_loc].set_xlim([xaxis_min, xaxis_max])
    flucstruc_ax[0].set_xlim([xaxis_min, xaxis_max])
    flucstruc_ax[0].set_ylim([0, 150])
    flucstruc_ax[0].set_title(r'$\kappa_H$:%.2f'%(kh), fontsize=8)
    #flucstruc_fig.canvas.draw(); flucstruc_fig.show()
    box1 = flucstruc_ax[0].get_position()
    box2 = flucstruc_ax[1].get_position()


    flucstruc_ax[0].set_position([box1.x0, box1.y0, box1.width * 0.8, box1.height])
    flucstruc_ax[1].set_position([box2.x0, box2.y0, box2.width * 0.8, box2.height])
    leg1 = flucstruc_ax[0].legend(loc='center left', bbox_to_anchor=(1, 0.5), fancybox = True)
    leg2 = flucstruc_ax[1].legend(loc='center left', bbox_to_anchor=(1, 0.5), fancybox = True)
    #leg1.get_frame().set_alpha(0.5)
    #pt.setp(leg.get_texts(), fontsize = 'small')
    flucstruc_fig.savefig(base_dir + '/' + 'Fluc_kh_%03d_'%(int(round(kh*100))) + title_string)
    print ' ', file_name, 
    fig.savefig(file_name)
    fig.clf()
    flucstruc_fig.clf()
    #pt.close('all')

    link_path = base_dir + '/kh_link/kh_%03d_shot_%s'%(int(round(kh*100)), title_string)
    command = 'ln -s -f ../%s %s'%(title_string, link_path)
    #print command
    os.system(command)

    print 'success'
    return 0

def polarisation_fft(data,start_location,samples,freq):
    fft_transform = np.fft.fft(data.signal[start_location:start_location+samples])
    #NEED TO THINK MORE ABOUT WHAT THIS IS???
    if len(fft_transform)==2048:
        pass
    else:
        print 'non 2048,' ,len(fft_transform), samples

    freq_array = np.fft.fftfreq(len(fft_transform),d=data.timebase[1]-data.timebase[0])[0:len(fft_transform)/2]
    location = np.argmin(np.abs(freq_array - freq))
    #print 'polarisation output ', fft_transform[location], freq_array[location], freq
    return fft_transform[location]

def single_shot_fluc_strucs_new(shot, array, time_bounds, samples, fs_dictionary = None, serial_number=0, power_cutoff = 0.1, n_svs = 2, ignore_ne_fail = 0, ne_array = 0, naked_coil = 0, overlap = 4):
    #This is supposed to replace the old single_shot_fluc struc with the new
    #instance_array, misc_details_dict approach....
    data = pf.getDevice('H1').acq.getdata(shot, array).reduce_time(time_bounds)
    data = data.subtract_mean(copy=False).normalise(method='v',separate=True,copy=False)
    data_segmented = data.segment(samples,overlap=overlap,datalist = 1)

    #Get the naked coil and interferometer array if required
    if ne_array:
        ne_data = pf.getDevice('H1').acq.getdata(shot, "ElectronDensity").change_time_base(data.timebase)
        ne_data_segmented = ne_data.segment(samples,overlap=overlap,datalist=1)
    else: ne_data_segmented = [None for i in data_segmented]
    if naked_coil:
        naked_coil = pf.getDevice('H1').acq.getdata(shot, "H1ToroidalNakedCoil").change_time_base(data.timebase)
        naked_coil_segmented = naked_coil.segment(samples,overlap=overlap,datalist=1)
    else: naked_coil_segmented = [None for i in data_segmented]

    meta_data = ['kh','heating_freq','main_current','sec_current', 'shot']
    #kh, heating_freq, main_current, sec_current = data.meta['kh'], data.meta['heating_freq'],data.meta['main_current'], data.meta['sec_current'] 
    #data.meta['shot'] = shot
    #remove the mean, normalise and segment the data....

    #Segment the data....
    instance_array_list = []
    misc_data_dict = {'RMS':[],'time':[], 'svs':[]}
    if naked_coil: misc_data_dict['naked_coil'] = []
    if ne_array: misc_data_dict['ne_mode'] = []; misc_data_dict['ne_static'] = []

    meta_values = ['kh','main_current','heating_freq','shot','sec_current']
    fs_values = ['p','a12','H','freq','E']
    for i in meta_values: misc_data_dict[i]=[]
    for i in fs_values: misc_data_dict[i]=[]

    for data_seg, ne_seg, naked_coil_seg in zip(data_segmented, ne_data_segmented, naked_coil_segmented):
        time_seg_average_time=np.mean([data_seg.timebase[0],data_seg.timebase[-1]])
        fs_set = data_seg.flucstruc()
        if ne_seg!=None:
            ne_fft = np.fft.rfft(ne_seg.signal)/samples
            if not np.allclose(ne_seg.timebase,data_seg.timebase): print "WARNING possible timebase mismatch between ne_data and data!!!"
        if naked_coil_seg!=None:
            naked_fft = np.fft.rfft(naked_coil_seg.signal)/samples
            if not np.allclose(naked_coil_seg.timebase,data_seg.timebase): print "WARNING possible timebase mismatch between ne_data and naked coil!!!"
        d = (data_seg.timebase[1] - data_seg.timebase[0])
        val = 1.0/(samples*d)
        N = samples//2 + 1
        frequency_base = np.round((np.arange(0, N, dtype=int)) * val,4)
        #get the valid flucstrucs
        valid_fs = []
        #make a list of the valid flucstrucs
        for fs in fs_set:
            if (fs.p>power_cutoff) and (len(fs.svs())>=n_svs): valid_fs.append(fs)
        #extract the useful information from the valid flucstrucs
        for fs in valid_fs:
            for i in fs_values: misc_data_dict[i].append(getattr(fs,i))
            misc_data_dict['svs'].append(fs.svs())
            #for i in meta_values: misc_data_dict[i].append(eval(i))
            for i in meta_values: #misc_data_dict[i].append(eval(i))
                try:
                    misc_data_dict[i].append(copy.deepcopy(data.meta[i]))
                except KeyError:
                    misc_data_dict[i].append(None)
                    
            phases = np.array([tmp_phase.delta for tmp_phase in fs.dphase])
            phases[np.abs(phases)<0.001]=0
            instance_array_list.append(phases)
            misc_data_dict['RMS'].append((np.mean(data.scales**2))**0.5)
            misc_data_dict['time'].append(time_seg_average_time)

            #ne_data and naked_coil_data
            tmp_loc = np.argmin(np.abs(misc_data_dict['freq'][-1]-frequency_base))
            if ne_seg!=None:
                misc_data_dict['ne_static'].append(np.abs(ne_fft[:,0]))
                misc_data_dict['ne_mode'].append(ne_fft[:,tmp_loc])
            if naked_coil_seg!=None:
                misc_data_dict['naked_coil'].append(naked_fft[:,tmp_loc])
    #convert lists to arrays....
    for i in misc_data_dict.keys():misc_data_dict[i]=np.array(misc_data_dict[i])
    return np.array(instance_array_list), misc_data_dict


def single_shot_fluc_strucs(current_shot, array, time_bounds, samples, fs_dictionary = None, serial_number=0, power_cutoff = 0.1, n_svs = 2, ignore_ne_fail = 0, ne_array = 0, overlap = 4):
    if fs_dictionary == None:
        fs_dictionary = {}
    #print ' in single shot fluc struc'
    kh, ne_value, ne_time, main_current, heating_freq, ne_value_list, ne_time_list = extract_ne_kh(current_shot, ne_array = ne_array)
    #print 'got ne, kh etc...'
    data=extract_data(current_shot, array)
    #print 'extracted data'
    coil_1x_data, coil_1y_data, coil_1z_data = extract_polarisation_data(current_shot) #NEW EDIT SHAUN
    #print 'got polarisation'
    data_reduced_time=data.reduce_time(time_bounds,copy=True).subtract_mean(copy=False).normalise(method='v',separate=True,copy=False)
    if (kh==None or data==None) and (ignore_ne_fail!=1):
        success = 0
    else:
        for t_seg in data_reduced_time.segment(samples,overlap=overlap):
            time_seg_average_time=np.mean([t_seg.timebase[0],t_seg.timebase[-1]])

            start_location = np.argmin(np.abs(coil_1x_data.timebase - t_seg.timebase[0]))

            fs_set=t_seg.flucstruc()
            ne_list = []
            if ne_value!=None:
                ne=ne_value[np.searchsorted(ne_time,time_seg_average_time)]
                if ne_array:
                    ne_list = []
                    for tmp_ne in range(0,len(ne_value_list)):
                        ne_list.append(ne_value_list[tmp_ne][np.searchsorted(ne_time_list[tmp_ne],time_seg_average_time)])
            else:
                ne = None
            for fs in fs_set:
                if fs.p>power_cutoff:
                    if len(fs.svs())>=n_svs:
                        fs_dictionary[serial_number]={}
                        phases=np.zeros(len(fs.dphase),dtype=float)
                        fs_dictionary[serial_number]['shot']=current_shot
                        fs_dictionary[serial_number]['coil_1x']=polarisation_fft(coil_1x_data,start_location,samples,fs.freq)
                        fs_dictionary[serial_number]['coil_1y']=polarisation_fft(coil_1y_data,start_location,samples,fs.freq)
                        fs_dictionary[serial_number]['coil_1z']=polarisation_fft(coil_1z_data,start_location,samples,fs.freq)
                        fs_dictionary[serial_number]['time']=time_seg_average_time*1000
                        fs_dictionary[serial_number]['kh']=kh
                        fs_dictionary[serial_number]['main_current']=main_current
                        fs_dictionary[serial_number]['heating_freq']=heating_freq
                        fs_dictionary[serial_number]['ne']=ne
                        fs_dictionary[serial_number]['p']=fs.p
                        fs_dictionary[serial_number]['a12']=fs.a12
                        fs_dictionary[serial_number]['H']=fs.H
                        fs_dictionary[serial_number]['freq']=fs.freq
                        fs_dictionary[serial_number]['RMS']=(np.mean(data_reduced_time.scales**2))**0.5
                        fs_dictionary[serial_number]['E']=fs.E
                        fs_dictionary[serial_number]['SVs']=(len(fs.svs()))
                        fs_dictionary[serial_number]['length_phases']=len(fs.dphase)
                        fs_dictionary[serial_number]['phases']=[]
                        fs_dictionary[serial_number]['ne_array']=ne_list
                        for fs_phase in range(0,len(fs.dphase)):
                            if np.abs(fs.dphase[fs_phase].delta)<0.001:
                                phases[fs_phase]=0
                            else:
                                phases[fs_phase] = fs.dphase[fs_phase].delta
                        fs_dictionary[serial_number]['phases']=phases                        
                        #fs_dictionary[serial_number]['fft']=fs.fft_values
                        serial_number+=1
        success = 1
    return fs_dictionary, serial_number, success

def flucstruc_step_through(shots,array='H1ToroidalAxial',file_name='/home/srh112/code/flucstrucs_pickle/outputlog_flucstruc_misc.txt',samples=4096,time_bounds=[0,0.08],worker_number=0, power_cutoff = 0.1, n_svs = 2, ne_array = 0):
    fs_dictionary={}; finished_count = 0; serial_number=1
    for current_shot in shots:
        start_shot_time=time.time()
        try:
            fs_dictionary, serial_number, success = single_shot_fluc_strucs(current_shot, array, time_bounds, samples, fs_dictionary = fs_dictionary, serial_number = serial_number, power_cutoff = power_cutoff, n_svs = n_svs, ne_array = ne_array)
        except:
            success = 0
        if success == 0:
            print 'Worker: %d, FAILED: %d, %d of %d'%(worker_number, current_shot, finished_count, len(shots))
            finished_count += 1
        else:
            finished_count += 1
            print 'Worker: %d, Finish: %d , time: %.2fs, %d of %d'%(worker_number, current_shot, time.time()-start_shot_time, finished_count, len(shots))
    temp1=time.time()
    print 'Worker_num : %d, dumping pickle data'%(worker_number)
    pickle.dump(fs_dictionary,open(file_name+'.pickle','a'))
    print 'Worker_num : %d, data dumped : time %.2f'%(worker_number,time.time()-temp1)


def flucstruc_step_through_old(shots,array='H1ToroidalAxial',file_name='/home/srh112/code/flucstrucs_pickle/outputlog_flucstruc_misc.txt',samples=4096,time_bounds=[0,0.08],worker_number=0):
    fs_dictionary={}
    finished_count = 0
    #f=csv.writer(open(file_name,'a'),delimiter=' ')
    #f=open(file_name,'a')
    serial_number=1
    for j in shots:
        start_shot_time=time.time()
        current_shot=j

        kh, ne_value, ne_time, main_current, heating_freq = extract_ne_kh(current_shot)

        data=extract_data(current_shot,array)


        coil_1x_data, coil_1y_data, coil_1z_data = extract_polarisation_data(current_shot) #NEW EDIT SHAUN

        data_extract_time=time.time()-start_shot_time
        if kh==None or data ==None:
            print 'Worker: %d, FAILED: %d, %d of %d'%(worker_number,current_shot,finished_count,len(shots))
            finished_count = finished_count +1
        else:
            data_reduced_time=data.reduce_time(time_bounds,copy=True).subtract_mean(copy=False).normalise(method='v',separate=True,copy=False)
            for t_seg in data_reduced_time.segment(samples,overlap=4):
                time_seg_average_time=np.mean([t_seg.timebase[0],t_seg.timebase[-1]])

                start_location = np.argmin(np.abs(coil_1x_data.timebase - t_seg.timebase[0]))

                fs_set=t_seg.flucstruc()
                ne=ne_value[np.searchsorted(ne_time,time_seg_average_time)]
                for fs in fs_set:
                    if fs.p>0.1:
                        if len(fs.svs())==2:

                            fs_dictionary[serial_number]={}
                            phases=np.zeros(len(fs.dphase),dtype=float)
                            fs_dictionary[serial_number]['shot']=current_shot
                            fs_dictionary[serial_number]['coil_1x']=polarisation_fft(coil_1x_data,start_location,samples,fs.freq)
                            fs_dictionary[serial_number]['coil_1y']=polarisation_fft(coil_1y_data,start_location,samples,fs.freq)
                            fs_dictionary[serial_number]['coil_1z']=polarisation_fft(coil_1z_data,start_location,samples,fs.freq)
                            fs_dictionary[serial_number]['time']=time_seg_average_time*1000
                            fs_dictionary[serial_number]['kh']=kh
                            fs_dictionary[serial_number]['main_current']=main_current
                            fs_dictionary[serial_number]['heating_freq']=heating_freq
                            fs_dictionary[serial_number]['ne']=ne
                            fs_dictionary[serial_number]['p']=fs.p
                            fs_dictionary[serial_number]['a12']=fs.a12
                            fs_dictionary[serial_number]['H']=fs.H
                            fs_dictionary[serial_number]['freq']=fs.freq
                            fs_dictionary[serial_number]['RMS']=(np.mean(data_reduced_time.scales**2))**0.5
                            fs_dictionary[serial_number]['E']=fs.E
                            fs_dictionary[serial_number]['SVs']=(len(fs.svs()))
                            fs_dictionary[serial_number]['length_phases']=len(fs.dphase)
                            fs_dictionary[serial_number]['phases']=[]

                            for fs_phase in range(0,len(fs.dphase)):
                                if np.abs(fs.dphase[fs_phase].delta)<0.001:
                                    phases[fs_phase]=0
                                else:
                                    phases[fs_phase] = fs.dphase[fs_phase].delta
                            fs_dictionary[serial_number]['phases']=phases                        
                            #fs_dictionary[serial_number]['fft']=fs.fft_values
                            serial_number+=1
            finished_count += 1
            print 'Worker: %d, Finish: %d , Extract: %.2fs, Comp: %.2fs, %d of %d'%(worker_number,current_shot,data_extract_time,time.time()-start_shot_time-data_extract_time,finished_count,len(shots))
    temp1=time.time()
    print 'Worker_num : %d, dumping pickle data'%(worker_number)
    pickle.dump(fs_dictionary,open(file_name+'.pickle','a'))
    print 'Worker_num : %d, data dumped : time %.2f'%(worker_number,time.time()-temp1)
