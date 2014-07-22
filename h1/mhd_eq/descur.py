'''
This runs DESCUR many times using different parameters to try and find
the best set of options that can provide a fit to the heliac output at
the furthest out trace possible.

The output from this can then be used to choose the trace, and DESCUR
settings to use for VMEC runs if you want the a surface that is very
far out.

SRH: 2Dec2013
'''

import h1.mhd_eq.heliac_vmec_utils as hv_utils
import h1.mhd_eq.heliac_worker_funcs as heliac
import os, tempfile,shutil,pickle, itertools
import subprocess, copy
import numpy as np

def run_many_descur(input_data, method_styles=['<','o','s','d'], r_styles = ['r','b','k','y'], surface_styles = ['-','--','-.'], points_per_surface_styles = ['k','b','r'],ax=None):
    kh_value, desc_points_per_surface_list, desc_surfaces, traces, method_list, mu_list, heliac_base_dir, nu,output_pickle,max_Ra = input_data
    print input_data
    results = dict.fromkeys(desc_points_per_surface_list, 
                            dict.fromkeys(desc_surfaces, {}))
    heliac_dependent_params = itertools.product(desc_points_per_surface_list, desc_surfaces)
    #print len(heliac_dependent_params), len(desc_points_per_surface_list)*len(desc_surfaces)
    #loop through all possibilities of heliac dependent parameters
    results_list = []
    for desc_points_per_surface, desc_surf in heliac_dependent_params:
        heliac_filename = heliac_base_dir + '/kh%.3f-kv1.000fixed/kh%.3f-kv1.000-desc%d_%d-fixed'%(kh_value,kh_value,desc_surf,desc_points_per_surface)
        temp_dir = tempfile.mkdtemp(prefix="descur-", dir='/tmp')
        shutil.copy(heliac_filename+'.trace', temp_dir)
        shutil.copy(heliac_filename+'.plt', temp_dir)
        os.chdir(temp_dir)
        curr_heliac = heliac.heliac(heliac_filename+'.plt',3)
        curr_heliac.read_trace_output2(heliac_filename+'.trace')
        if max_Ra!=None:
            tmp_loc = np.argmin(np.abs(curr_heliac.Ra_list-max_Ra))
            start_trace = curr_heliac.trace_list_order[tmp_loc]
            print 'Forcing smaller start trace because of max_Ra'
        else:
            start_trace = np.max(curr_heliac.trace_list_order)
        #start_trace = generate_descur_input(heliac_filename, 'in_descur',force_trace=None)
        trace_list = range(start_trace - traces, start_trace + 1)
        non_heliac_dependent_params = itertools.product(trace_list, method_list, mu_list)
        results[desc_points_per_surface][desc_surf]=dict.fromkeys(trace_list, dict.fromkeys(method_list,{}))

        #loop through all possibilities of non heliac dependent parameters
        for force_trace, method, mu in non_heliac_dependent_params:
            curr_heliac.format_for_descur_input(force_trace=force_trace, filename='in_descur')

            #generate_descur_input(heliac_filename, 'in_descur',force_trace=force_trace)
            nv=desc_points_per_surface
            executable = '/home/srh112/vmec_847/DESCUR/Release/xcurve_nu%d_mu%d_nv%d'%(nu,mu,nv)
            print 'PID: %6d, kh: %6.3f, pts_per_surf: %4d, n_surf:%3d,trace:%3d,method:%3d,mu:%3d'%(os.getpid(),kh_value, desc_points_per_surface, desc_surf, force_trace, method, mu)
            process=subprocess.Popen(executable,stdin = subprocess.PIPE,stdout=subprocess.PIPE)
            process.communicate('%d\nV\n0\nin_descur\nY\n'%(method))

            success = 0; read_in_success = 0
            try:
                with open('outcurve','r') as fin:
                    lines = fin.readlines();
                    read_in_success = 1
            except:
                print "Unable to open outcurve"
            if read_in_success:
                for i in lines:
                    if i.find('MB')>0 and i.find('NB')>0: success=1
            if success==1:
                x = hv_utils.DescurPlotout('plotout')
                RMS_error=[]
                for i in range(0,x.nplanes):
                    try:
                        RMS_error.append(x.RMS_error(plane=i))
                    except ValueError as e:
                        success=0
                if success:
                    results[desc_points_per_surface][desc_surf][force_trace][method][mu]={}
                    results[desc_points_per_surface][desc_surf][force_trace][method][mu]['rad']=np.max(x.rz[:,0])
                    results[desc_points_per_surface][desc_surf][force_trace][method][mu]['RMS']=copy.deepcopy(RMS_error)
                    #print np.max(RMS_error)
                    results_list.append([np.max(RMS_error)*1000.,desc_points_per_surface,desc_surf,int(force_trace),method,mu,np.max(x.rz[:,0])])

                    #overall_rms_errors.append(RMS_error)
                    print 'DESCUR Success'
                    #mu_success_list.append(mu)
            else:
                print 'DESCUR failed to finish'
    pickle.dump(results_list,file(output_pickle,'w'))
    return (copy.deepcopy(results), copy.deepcopy(results_list))

def print_results(results):
    RMS=0;pts_surf=1;n_surf=2;trace=3;method=4;mu=5;radius=6
    max_trace_val = int(np.max(results[:,trace]))
    min_trace_val = int(np.min(results[:,trace]))
    for i in range(min_trace_val,max_trace_val+1):
        truth = (np.abs(results[:,trace] - i))<0.01
        if np.sum(truth)>=1:
            tmp = results[truth,:]
            min_val = np.argmin(tmp[:,RMS])
            #print 'trace: %d, RMS: %.3fmm, pts/surf:%d,n_surf:%d,method %d,mu %d,radius:%.3f'%(i, tmp[min_val,RMS],tmp[min_val,pts_surf],tmp[min_val,n_surf],tmp[min_val,method],tmp[min_val,mu],tmp[min_val,radius])
            print 'force_trace= %d;desc_points_per_surface=%d;desc_surf=%d;method=%d;mu=%d;#radius:%.3f,RMS: %.3fmm,'%(i,tmp[min_val,pts_surf],tmp[min_val,n_surf],tmp[min_val,method],tmp[min_val,mu],tmp[min_val,radius],tmp[min_val,RMS])


def produce_ordered(results):
    print 'applying handicaps'
    #RMS=0;pts_surf=1;n_surf=2;trace=3;method=4;mu=5;radius=6
    working_copy = copy.deepcopy(results)
    max_trace = np.max(working_copy[:,3])
    min_trace = np.min(working_copy[:,3])
    handicap_list = [0,2,3.5,5,6.5,8,9.5]
    for i in range(0,6):
        curr_trace = max_trace - i
        handicap = handicap_list[i]
        truth = (np.abs(results[:,3] - curr_trace))<0.01
        print np.sum(truth), curr_trace, handicap
        #truth = working_copy[:,3] == curr_trace
        working_copy[truth,0] = working_copy[truth,0] + handicap
    print working_copy.shape
    tmp = working_copy.tolist()
    tmp.sort()
    return tmp

def run_vmec2(x):
    '''This function runs descur and vmec in a loop until a vmec output is obtained.
    There is no quality control on the DESCUR output other than it exists and
    VMEC can converge using it

    starting_dir, input_filename, kh_value, kv_value, descur_run_results, use_dave_iota_param = x
    SH: 2Apr2013
    '''
    #need to determine the list of best descur options here
    #need a way to get the kh value here also!
    starting_dir, input_filename, kh_value, kv_value, descur_run_results, use_dave_iota_param  = x
    os.chdir(starting_dir)
    #input_filename = 'input.'+x[1]
    input_filename = 'input.'+input_filename
    #kh_value = x[2]
    #kv_value = x[3]
    #descur_run_results = x[4]
    #current_trace_value,descur_success = run_descur(x,executable = 'xcurve',force_trace=None)
    vmec_count = 0; vmec_success = 0;
    #descur_results = pickle.load(file('/home/srh112/code/python/h1_eq_generation/descur_results3.pickle','r'))
    #descur_results = pickle.load(file('/home/srh112/code/python/h1_eq_generation/descur_tests2/tmp3_pt2.pickle','r'))
    #descur_results = pickle.load(file('/home/srh112/code/python/h1_eq_generation/descur_results.pickle','r'))
    descur_results = pickle.load(file(descur_run_results,'r'))
    best_ops = np.array(descur_results['%.3f'%(kh_value)]['handicapped'])
    #loop to make sure vmec finishes successfully
    heliac_base_dir = '/home/srh112/code/python/h1_eq_generation/results6/'
    heliac_base_dir = starting_dir
    while (not vmec_success) and (vmec_count<5):
        #Loop to make sure descur finishes successfully
        #These are the list items in the descur run through results
        RMS_key=0;pts_surf_key=1;n_surf_key=2;trace_key=3;method_key=4;mu_key=5;radius_key=6
        heliac_filename = heliac_base_dir + '/kh%.3f-kv1.000-desc%d_%d-fixed'%(kh_value,best_ops[vmec_count,n_surf_key],best_ops[vmec_count,pts_surf_key])

        #for i in ['.plt','.trace','.out']: shutil.copy(heliac_filename+i, starting_dir)
        curr_heliac = heliac.heliac(heliac_filename+'.plt',3)
        curr_heliac.read_trace_output2(heliac_filename+'.trace')
        nv=best_ops[vmec_count,pts_surf_key]#desc_points_per_surface
        nu=200
        mu=best_ops[vmec_count,mu_key]
        force_trace = int(best_ops[vmec_count,trace_key])

        #this should generate the descur input file
        curr_heliac.format_for_descur_input(force_trace=int(force_trace), filename='in_descur')

        #select the correct executable for these settings
        executable = '/home/srh112/vmec_847/DESCUR/Release/xcurve_nu%d_mu%d_nv%d'%(nu,mu,nv)
        process=subprocess.Popen(executable,stdin = subprocess.PIPE,stdout=subprocess.PIPE)

        #select the correct method when descur is running
        process.communicate('%d\nV\n0\nin_descur\nY\n'%(best_ops[vmec_count,method_key]))

        #check to see if descur completed properly
        with open('outcurve','r') as fin:lines = fin.readlines()
            
        descur_success = 0
        for i in lines:
            if i.find('MB')>0 and i.find('NB')>0: descur_success=1
        if descur_success==1:
            descur_plotout = hv_utils.DescurPlotout('plotout')
            RMS_error=[]
            for i in range(0,descur_plotout.nplanes):
                try:
                    RMS_error.append(descur_plotout.RMS_error(plane=i))
                except ValueError as e:
                    RMS_error.append(-0.01)
                    descur_success=0
            if descur_success:
                print 'DESCUR successful run max rms :%.3fmm, %.3f'%(np.max(RMS_error)*1000,best_ops[vmec_count,0])
        else:
            print 'DESCUR failed to finish'

            #current_trace_value, descur_success = run_descur(x,executable = 'xcurve',force_trace=current_trace_value)
        #if not descur_success:
        #    print os.getpid(), 'descur run failed', current_trace_value
        #    current_trace_value -= 1
        #    descur_count+=1
        #vmec_success = 1
        vmec_success = 0

        #Now we are ready to try VMEC
        vmec_boundary_details = hv_utils.make_vmec_input2('outcurve', 'header.test', delrbs=True, prefix='input.',suffix='.outcurve', fmt='%15.10f')
        #generate the iota profile

        final_loc = curr_heliac.trace_list_order.index(force_trace)
        if use_dave_iota_param:
            print 'Using Dave iota paramaterisation'
            final_ra = curr_heliac.Ra_list[final_loc]
            print 'Final_ra:{}'.format(final_ra)
            s = np.linspace(0,1,100)
            iota = [heliac.iota_dave_fit(i, kh_value) for i in np.sqrt(s)*final_ra/0.2]
            p = np.polyfit(s, iota,4)
        else:
            final_loc = curr_heliac.trace_list_order.index(force_trace)
            final_ra = curr_heliac.Ra_list[final_loc]
            iota = curr_heliac.iota_list[:final_loc]
            psi = curr_heliac.psi_list[:final_loc]
            psi_norm = psi/np.max(psi)
            p = np.polyfit(psi_norm, iota,4)
            #iota = heliac_run.iota_list[:final_loc]
            #p = np.polyfit(psi_norm, iota,4)
        psi_vmec_7MHz = np.max(curr_heliac.psi_list[:final_loc]) /1000. *6500./0.1389e5
        print 'PHIEDGE is being set to {}'.format(psi_vmec_7MHz)
        tmp = ['%.6f'%(i) for i in p]
        tmp.reverse()
        ai_string = ', '.join(tmp)
        print ai_string
        hv_utils.generate_VMEC_input2(x, vmec_boundary_details, ai = ai_string, phiedge='{:.4f}'.format(psi_vmec_7MHz))

        # #generate iota profile using the heliac run
        # final_loc = curr_heliac.trace_list_order.index(force_trace)
        # iota = curr_heliac.iota_list[:final_loc]
        # psi = curr_heliac.psi_list[:final_loc]
        # psi_norm = psi/np.max(psi)
        # p = np.polyfit(psi_norm, iota,4)
        # tmp = ['%.6f'%(i) for i in p]
        # tmp.reverse()
        # ai_string = ', '.join(tmp)
        # #ai_string = '  ai           = ' + ai_string +'\n'
        # print ai_string
        # print p
        # hv_utils.generate_VMEC_input2(x, vmec_boundary_details, ai = ai_string)


        os.chdir(starting_dir)
        if os.path.isfile('iAmFinishedVMEC'):
            vmec_success = 1
            print 'vmec already finshed successfully for this directory!!'
        else:
            #os.system('rm iAmFinishedVMEC')
            os.system('/home/srh112/bin/xvmec2000_nc ' + input_filename+' >VMEC.out 2>VMEC.err')
            VMEC_output_filename = 'wout_'+input_filename.lstrip('input.')+'.nc'
            vmec_count += 1
            try:
                with open(VMEC_output_filename): pass
                if os.stat(VMEC_output_filename).st_size>200000:
                    print os.getpid(), 'output file bigger than 200K, so assuming run finished ok'
                    with file('iAmFinishedVMEC','w') as f_tmp:
                        tmp_string = '%s\n%d\n'%(heliac_filename,force_trace)
                        tmp_string += '{:.4f}\n'.format(final_ra)
                        f_tmp.write(tmp_string)
                    #os.system('touch iAmFinishedVMEC')
                    vmec_success = 1
                else:
                    print os.getpid(), 'Oh dear, VMEC unsuccessful, moving to the next best descur fit'
            except IOError:
                #current_trace_value -= 1
                print os.getpid(), 'Oh dear, VMEC unsuccessful, moving to the next best descur fit'
    print os.getpid(), starting_dir, input_filename, 'finished'
    #run_xform(x)



def run_xform(x):
    '''This function runs booz_xform based on the directories in x
    SH: 2Apr2013
    '''
    starting_dir = x[0]
    input_filename = 'input.'+x[1]
    print os.getpid(), 'booz_xform',starting_dir, input_filename
    os.chdir(starting_dir)
    #make the booz_xform input file

    Ns=100
    xbooz_input = '70 70\n'+'wout_'+'%s\n'%(x[1])
    xbooz_input += ' '.join([str(i) for i in range(1,Ns+1)])+'\n'
    with file('in_boo.txt','w') as f:
        f.write(xbooz_input)

    #os.system('rm iAmFinishedBOOZ')
    if os.path.isfile('iAmFinishedBOOZ'):
        pass
    else:
        try:
            with open('iAmFinishedVMEC'): pass
            vmec_success=1
        except IOError:
            print 'VMEC finished file doesnt seem to be there breaking'
            vmec_success=0
        if vmec_success:
            os.system('/home/srh112/bin/xbooz_xform in_boo.txt >BOOZ_XFORM.out')
            BOOZ_output_filename = 'boozmn_wout_'+input_filename.lstrip('input.')+'.nc'
            try:
                with open(BOOZ_output_filename): pass
                if os.stat(BOOZ_output_filename).st_size>200000:
                    print 'output file bigger than 200K, so assuming run finished ok'
                    os.system('touch iAmFinishedBOOZ')
            except IOError:
                print 'Oh dear, something went wrong, output file probably doesnt exist'
    print os.getpid(), starting_dir, input_filename, 'finished '#, vmec_success


