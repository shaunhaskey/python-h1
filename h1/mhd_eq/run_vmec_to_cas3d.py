import commands
import os, time, copy
import subprocess
import numpy as np
import sys

vmec = False
mc3d = False
prerun = True
cas3d_matrix = True
cas3d_eig = True
conti = False
cas3d_single_eig = False
free_boundary = True
N = 1
ew = 2.6070e-4
ew = 4.7006e-4
ew = 3.6051e-4
ew = 4.9459e-4
ew = 2.6206e-4
ew = 3.2688e-4
ew = 2.7710e-4
ew = 2.8904e-4
ew = 2.1685e-4
ew = 3.0896e-4
ew = 6.4371e-5
#ew = 2.4481e-5
##ew = 1.7611e-4
#ew = 4.2934e-4
#ew = 2.8216e-4
#ew = 1.7010e-4
#ew = 1.1532e-3
#ew = 1.0249e-3
#ew = 2.7942e-4
mode_table_size = 'tiny' #'tiny', 'tiniest'
mode_table_size = 'tiniest' #'tiny', 'tiniest'
mode_table_size = 'full'
density_profile = 'lin' #'flat', 'lin', 'lin_off', 'expt'
density_profile_eig = 'expt' #'flat', 'lin', 'lin_off', 'expt', 'quad', 'cubic', 'quartic'
size_txt = '' if mode_table_size=='full' else '_{}'.format(mode_table_size)

try:
    kappain = np.float(sys.argv[1])
    print kappain
except:
    kappain = 0.45
    print 'problems'
    pass


input_vmec_file = 'input.kh{:.3f}-kv1.000fixed'.format(kappain)
#/short/y08/srh112/
print input_vmec_file
input_dir = '/short/y08/srh112/' + 'vmec_input_files/'
work_directory = '/short/y08/srh112/whale_tail/'+ input_vmec_file +'_dir'
vmec_surfaces = 353
prerun_surfaces = 300
vmec_dir = work_directory+'/vmec/'
mc3d_dir = work_directory+'/mc3d/'
prerun_dir = work_directory+'/prerun_{}/'.format(prerun_surfaces)
fort12_loc = r'../prerun_{}/fort.12'.format(prerun_surfaces)
os.system('mkdir '+work_directory)
os.system('mkdir '+vmec_dir)

#CAS3D
CAS3D_settings = {}
CAS3D_settings['extra_info'] = ',software=nag'
CAS3D_settings['ncpus'] = '32'
CAS3D_settings['mem'] = '80GB'
CAS3D_settings['walltime'] = '05:00:00'
CAS3D_settings['module_load'] = 'module load openmpi'
CAS3D_settings['executable_string'] = 'mpirun /home/112/srh112/test_source2/cas3d_2007_10_12/source/wkin_even/main_Linux_mpi'
CAS3D_settings['log_file'] = 'out_cas3d'
CAS3D_settings['final_commands'] = ''

CAS3D_m_settings = {}
CAS3D_m_settings['extra_info'] = ',software=nag'
CAS3D_m_settings['ncpus'] = '16'
CAS3D_m_settings['mem'] = '60GB'
CAS3D_m_settings['walltime'] = '05:00:00'
CAS3D_m_settings['module_load'] = 'module load openmpi'
CAS3D_m_settings['executable_string'] = 'mpirun /home/112/srh112/test_source2/cas3d_2007_10_12/source/wkin_even/main_Linux_mpi'
CAS3D_m_settings['log_file'] = 'out_cas3d_m'
CAS3D_m_settings['final_commands'] = ''

CAS3D_r_settings = {}
CAS3D_r_settings['extra_info'] = ',software=nag'
CAS3D_r_settings['ncpus'] = '32'
CAS3D_r_settings['mem'] = '120GB'
CAS3D_r_settings['walltime'] = '10:00:00'
CAS3D_r_settings['module_load'] = 'module load openmpi'
CAS3D_r_settings['executable_string'] = 'mpirun /home/112/srh112/test_source2/cas3d_2007_10_12/source/wkin_even/main_Linux_mpi'
CAS3D_r_settings['log_file'] = 'out_cas3d_r'
#CAS3D_r_settings['final_commands'] = 'tar cvzf out.tar.gz *.dat out_cas3d input2 in_* *.ps\nrm fort.2??\nrm *.ps\nrm *.dat'
CAS3D_r_settings['final_commands'] = 'rm fort.2??\n'

CAS3D_single_eig_settings = {}
CAS3D_single_eig_settings['extra_info'] = ',software=nag'
CAS3D_single_eig_settings['ncpus'] = '1'
CAS3D_single_eig_settings['mem'] = '16GB'
CAS3D_single_eig_settings['walltime'] = '1:00:00'
CAS3D_single_eig_settings['module_load'] = 'module load openmpi'
CAS3D_single_eig_settings['executable_string'] = 'mpirun /home/112/srh112/test_source2/cas3d_2007_10_12/source/wkin_even/main_Linux_mpi'
CAS3D_single_eig_settings['log_file'] = 'out_cas3d_single_eig'
#CAS3D_single_eig_settings['final_commands'] = 'tar cvzf out.tar.gz *.dat out_cas3d input2 in_* *.ps\nrm fort.2??\nrm *.ps\nrm *.dat'
CAS3D_single_eig_settings['final_commands'] = 'rm fort.51 fort.16\n'

VMEC_settings = {}
VMEC_settings['extra_info'] = ''
VMEC_settings['ncpus'] = '1'
VMEC_settings['mem'] = '4GB'
VMEC_settings['walltime'] = '10:00:00'
VMEC_settings['module_load'] = ''
VMEC_settings['executable_string'] = '~/bin/xvmec2000 ' + input_vmec_file
VMEC_settings['log_file'] = 'out_vmec'
VMEC_settings['final_commands'] = ''

MC3D_settings = {}
MC3D_settings['extra_info'] = ',software=nag'
MC3D_settings['ncpus'] = '1'
MC3D_settings['mem'] = '6GB'
MC3D_settings['walltime'] = '10:00:00'
MC3D_settings['module_load'] = ''
MC3D_settings['executable_string'] = '~/bin/mc3d < in_mc3d '
MC3D_settings['log_file'] = 'out_mc3d'
MC3D_settings['final_commands'] = ''

prerun_settings = {}
prerun_settings['extra_info'] = ',software=nag'
prerun_settings['ncpus'] = '1'
prerun_settings['mem'] = '6GB'
prerun_settings['walltime'] = '10:00:00'
prerun_settings['module_load'] = ''
prerun_settings['executable_string'] = '~/bin/prerun <in_prerun '
prerun_settings['log_file'] = 'out_prerun'
prerun_settings['final_commands'] = ''

CONTI_settings = {}
CONTI_settings['extra_info'] = ',software=nag'
CONTI_settings['ncpus'] = '16'
CONTI_settings['mem'] = '40GB'
CONTI_settings['walltime'] = '00:40:00'
CONTI_settings['module_load'] = 'module load openmpi\nmodule load nag'
CONTI_settings['executable_string'] = 'mpirun /home/112/srh112/conti/source/conti h1'
CONTI_settings['log_file'] = 'out_conti'
CONTI_settings['final_commands'] = ''


def submit_check_if_done():
    status, output = commands.getstatusoutput("qsub jobscript")
    job_id = output[:output.find('.')]
    job_output_filename = 'jobscript.o'+job_id
    job_error_filename = 'jobscript.e'+job_id
    print status
    print output, job_output_filename
    while not os.path.exists(job_output_filename):
        time.sleep(5)
        print 'not done'
    time.sleep(5)
    with file(job_error_filename,'r') as tmp: tmp_len = len(tmp.read())
    if tmp_len>0: raise(Exception)

#setenv PATH /opt/bin:/bin:/usr/bin:/usr/X11R6/bin:/opt/pbs/default/bin

job_string = '''#!/bin/csh 
#PBS -q normal 
#PBS -l walltime=<<walltime>>,mem=<<mem>>,ncpus=<<ncpus>><<extra_info>>
#PBS -l wd

limit stacksize unlimited
<<module_load>>
#for matrix assembly
<<executable_string>> > <<log_file>>

#mc3d < in_mc3d > out_mc3d_054

#xvmec2000 input.054 > out_vmec_054

#prerun < in_prerun > out_prerun_054

#for eval scan
#./cas3d_serial >  out_cas3d_h1ass027v1_p43_10_lowb_left
<<final_commands>>

'''
vmec_input_file_changes = {'ns_array':'ns_array = 3, 23, 87, 163, 207, 291, 353\n', 'ftol_array':'ftol_array = 1.e-05,1.e-08,1.0e-10,1.0e-11,8.0e-12,6.0e-12,5.0e-12\n'}

if vmec:
    with file(input_dir +'/'+input_vmec_file,'r') as tmp:
        vmec_lines = tmp.readlines()
    for i in vmec_input_file_changes.keys():
        for j in range(len(vmec_lines)):
            if vmec_lines[j].find(i)>=0:
                vmec_lines[j] = vmec_input_file_changes[i]
    os.chdir(vmec_dir)
    with file(input_vmec_file,'w') as tmp:
        tmp.writelines(vmec_lines)

    #os.system('cp '+input_vmec_file +' ' + vmec_dir)
    new_job_string = copy.deepcopy(job_string)
    for i in VMEC_settings.keys():
        print i
        new_job_string = new_job_string.replace('<<'+i+'>>', VMEC_settings[i])
    with file('jobscript','w') as tmp:
        tmp.writelines(new_job_string)
    submit_check_if_done()

# def submit_check_if_done():
#     status, output = commands.getstatusoutput("qsub jobscript")
#     job_id = output[:output.find('.')]
#     job_output_filename = 'jobscript.o'+job_id
#     job_error_filename = 'jobscript.e'+job_id
#     print status
#     print output, job_output_filename
#     while not os.path.exists(job_output_filename):
#         time.sleep(5)
#         print 'not done'
#     time.sleep(5)
#     with file(job_error_filename,'r') as tmp: tmp_len = len(tmp.read())
#     if tmp_len>0: raise(Exception)
#     submit_check_if_done()
#     #os.system('qsub jobscript')


mc3d_in_file = '''&mc3d_indata
    icode=2009,
    woutfile="wout_input",
    nsurf=7,
    indexsurf=2, 50, 75, 95
    m0=30,    n0=30,  nu=120,   nv=120,
    conv=0.005,
/
&end
'''
stub = input_vmec_file.lstrip('input.')
print stub
if mc3d:
    os.system('mkdir '+mc3d_dir)
    os.chdir(mc3d_dir)
    os.system('ln -sf ../vmec/wout_'+stub +'.txt wout_input')

    with file('fort.77','w') as tmp:
        tmp.writelines('mc3d')
        
    with file('in_mc3d','w') as tmp:
        tmp.writelines(mc3d_in_file)
    new_job_string = copy.deepcopy(job_string)
    for i in MC3D_settings.keys():
        print i
        new_job_string = new_job_string.replace('<<'+i+'>>', MC3D_settings[i])
    with file('jobscript','w') as tmp:
        tmp.writelines(new_job_string)
    submit_check_if_done()
    #os.system('qsub jobscript')

prerun_change_settings = {'final_surfaces':str(prerun_surfaces), 'initial_surfaces':str(vmec_surfaces-1)}

prerun_in_file='''BINARY FILE format (New=.true. OLD=.false.)           perfect equilibrium
.true.                                                .true.
SUBINTERVAL:   IA       IB      NS (MESH FOR INTERPOLATION)
                0       <<initial_surfaces>>    <<final_surfaces>>
subinterval: ileft  iright
              0        <<final_surfaces>>
original surface resolution: ipnu_old   ipnv_old
                              4           4
surface resolution: m0  n0  nu  nv (refined mesh, =0: no refining)
                    0   0   0   0
nper(periods of configuration)    ntopol(mode topological torus)
3                                 3
back-transformation data needed
.true.
'''

    
if prerun:
    os.system('mkdir '+prerun_dir)
    os.chdir(prerun_dir)
    os.system('ln -sf ../mc3d/for_cas3d_stp_14.dat fort.11')
    for i in prerun_change_settings.keys():
        print i
        prerun_in_file = prerun_in_file.replace('<<'+i+'>>', prerun_change_settings[i])
    with file('in_prerun','w') as tmp:
        tmp.writelines(prerun_in_file)
    new_job_string = copy.deepcopy(job_string)
    for i in prerun_settings.keys():
        print i
        new_job_string = new_job_string.replace('<<'+i+'>>', prerun_settings[i])
    with file('jobscript','w') as tmp:
        tmp.writelines(new_job_string)
    submit_check_if_done()
    # status, output = commands.getstatusoutput("qsub jobscript")
    # job_id = output[:output.find('.')]
    # job_output_filename = 'jobscript.o'+job_id
    # job_error_filename = 'jobscript.e'+job_id
    # print status
    # print output, job_output_filename
    # while not os.path.exists(job_output_filename):
    #     time.sleep(5)
    #     print 'not done'
    # time.sleep(5)
    # with file(job_error_filename,'r') as tmp: tmp_len = len(tmp.read())
    # if tmp_len>0: raise(Exception)
        
    print 'Finished!!!'

    #os.system('qsub jobscript')

conti_change_settings = {'m_min':'0',
                         'm_max':'18',
                         'n_tor':'18',
                         'tor_harms':'-23 -22 -20 -19 -17 -16 -14 -13 -11 -10 -8 -7 -5 -4 -2 -1 1 2 ',
                         'freq_low':'0.0',
                         'freq_high':'500.0',
                         'idens':'2 1 -1',
                         'rad_points':'500'
                         }


#'idens':'6 0.84827165 4.31756279 -35.54226959 76.85344266 -68.80442217 22.33691942',
#'idens':'2 1 -1'

CAS3D_m_change_settings = {'surfaces':str(prerun_surfaces),
                           'gamma':'1.6666',
                           'evp':'1',
                           'shift':'5.44182E-04',
                           'shiftmin':'0.0',
                           'shiftmax':'8.000e-03',
                           'number':'30000',
                           'number_dens':'0.06',
                           'idens':'2 1.  -1',
                           'IRAND':'0',
                           'IFIX':'1',
                           'use_ite':'.true.',
                           'ite':'4 1. 0. 15. -10.',
                           'matrix':'0',
                           'nxis':'95','neta':'95','nmu':'209',
                           'm_pot':'10', 'n_pot':'10','nu_vac':'100',
                           'nv_vac':'100'}
#'idens':'6 0.84827165 4.31756279 -35.54226959 76.85344266 -68.80442217 22.33691942',
#'idens':'4 1.  0. -3. 2.',
#Dave's fit to HELIAC iota

def iota(r,kappa):
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

#mode = 1
if mode_table_size == 'full':
    mmin = 0
    mmax = 18
    mummin = 0
    mummax = 20
    nmin = -21
    nmax = 2
    munmin = -25
    munmax = 5 
    etatol=5
    mutol = 10
elif mode_table_size == 'tiny':
    mmin = 1
    mmax = 7
    mummin = 0
    mummax = 10
    nmin = -13
    nmax = 2
    munmin = -16
    munmax = 5 
    etatol = 5
    mutol = 10
elif mode_table_size == 'tiniest':
    mmin = 3
    mmax = 5
    mummin = 2
    mummax = 6
    nmin = -13
    nmax = 2
    munmin = -16
    munmax = 5 
    etatol= 4
    mutol = 4
else:
    raise ValueError('mode_table_size must be full, tiny or tiniest')
if density_profile == 'flat':
    idens_txt = '1  1.'
elif density_profile == 'lin':
    idens_txt = '2 1.  -1'
elif density_profile == 'lin_off':
    idens_txt = '2 1  -0.8'
elif density_profile == 'expt':
    #idens_txt = '4 1. 0. -3. 2.'
    idens_txt = '6 0.84827165 4.31756279 -35.54226959 76.85344266 -68.80442217 22.33691942'
    idens_txt = '6 1.         2.01432225 -13.42055742  18.19353285  -9.1737854 1.38826713'
    #idens_txt = '6 1.        , -3.60869565,  4.82608696, -2.82608696,  0.60869565'

if density_profile_eig == 'flat':
    idens_eig_txt = '1  1.'
elif density_profile_eig == 'lin':
    idens_eig_txt = '2 1.  -1'
elif density_profile_eig == 'lin_off':
    idens_eig_txt = '2 1  -0.8'
elif density_profile_eig == 'quad':
    idens_eig_txt = '3 1. -2. 1.'
elif density_profile_eig == 'cubic':
    idens_eig_txt = '4 1. -3. 3. -1.'
elif density_profile_eig == 'quartic':
    idens_eig_txt = '5 1.        , -3.60869565,  4.82608696, -2.82608696,  0.60869565'
elif density_profile_eig == 'expt':
    #idens_txt = '4 1. 0. -3. 2.'
    #idens_txt = '6 0.84827165 4.31756279 -35.54226959 76.85344266 -68.80442217 22.33691942'
    idens_eig_txt = '5 1.        -3.60869565  4.82608696 -2.82608696  0.60869565'
    idens_eig_txt = '6 1.         2.01432225 -13.42055742  18.19353285  -9.1737854 1.38826713'
    #idens_txt = '6 1.         2.01432225 -13.42055742  18.19353285  -9.1737854 1.38826713'


CAS3D_m_change_settings['m_pot'] = str(np.max(np.abs([mmin,mmax])))
CAS3D_m_change_settings['title'] = 'h1 {:.2f}'.format(kappain)
CAS3D_m_change_settings['n_pot'] = str(np.max(np.abs([nmin,nmax])))
CAS3D_m_change_settings['nu_vac'] = str(int(CAS3D_m_change_settings['m_pot'])*5)
CAS3D_m_change_settings['nv_vac'] = str(int(CAS3D_m_change_settings['n_pot'])*5)
CAS3D_m_change_settings['idens'] = idens_txt
if free_boundary: 
    CAS3D_m_change_settings['IFIX'] = '0'
    CAS3D_m_change_settings['IRAND'] = '1'
sample = np.arange(0,1.01,0.01)
iotamin = min([iota(r,kappain) for r in sample])
iotamax = max([iota(r,kappain) for r in sample])
print 'iotamin=',iotamin,' iotamax=',iotamax

etatable = ''
count = 1
for n in range(nmin,nmax+1):
    etatable = etatable + repr(n).rjust(5)+'   '
    for m in range(mmin,mmax+1):
        #all pairs within a certain tolerance from resonance, no -ve n when m=0 and modefamily condition
        if m*(iotamin)-etatol<-n and -n<m*(iotamax)+etatol and (m!=0 or n<=0) and (N==2 or not(N^(n%3>0))):
            etatable = etatable + repr(count).rjust(4)
            count=count+1
        else:
            etatable = etatable + repr(0).rjust(4)
    etatable = etatable+'\n'
CAS3D_m_change_settings['nxis'] = str(count-1)
CAS3D_m_change_settings['neta'] = str(count-1)
#fout.write(repr(count-1).rjust(6))
#fout.write(repr(count-1).rjust(8))

mutable=''
count=1
for n in range(munmin,munmax+1):
    mutable = mutable + repr(n).rjust(5)+'   '
    for m in range(mummin,mummax+1):
        if m*(iotamin)-mutol<-n and -n<m*(iotamax)+mutol and (m!=0 or n<=0) and (N==2 or not(N^(n%3>0))):
            mutable = mutable + repr(count).rjust(4)
            count=count+1
        else:
            mutable = mutable + repr(0).rjust(4)
    mutable = mutable+'\n'

#fout.write(repr(count-1).rjust(8)+repr(0).rjust(8)+'\n')
CAS3D_m_change_settings['nmu'] = str(count-1)

#fin.readline()
#for i in range(5): fout.write(fin.readline())

eta_header = ''.rjust(8)
for i in range(mmin,mmax+1): 
    eta_header+=repr(i).rjust(4)
eta_header+='\n\n'
#eta_header+='\n'
etatable = ' '+repr(mmin)+' '+repr(mmax)+' '+repr(nmin)+' '+repr(nmax)+' '+'1\n\n' + eta_header + etatable
#fout.write(' '+repr(mmin)+' '+repr(mmax)+' '+repr(nmin)+' '+repr(nmax)+' '+'1\n\n')
print '==============='
print etatable


#fout.write(''.rjust(8))
#for i in range(mmin,mmax+1): 
#    fout.write(repr(i).rjust(4))
#fout.write('\n\n')
#fout.write(etatable)

#fout.write('c  eta\n '+repr(mmin)+' '+repr(mmax)+' '+repr(nmin)+' '+repr(nmax)+' '+'1\n\n')
#fout.write(''.rjust(8))
#for i in range(mmin,mmax+1): 
#    fout.write(repr(i).rjust(4))
#fout.write('\n\n')
#fout.write(etatable)

#fout.write('c  mu\n '+repr(mummin)+' '+repr(mummax)+' '+repr(munmin)+' '+repr(munmax)+' '+'1\n\n')
mu_header=''.rjust(8)
for i in range(mummin,mummax+1): 
    mu_header+=repr(i).rjust(4)
mu_header+='\n\n'
#mu_header+='\n'
mutable = ' '+repr(mummin)+' '+repr(mummax)+' '+repr(munmin)+' '+repr(munmax)+' '+'1\n\n'+ mu_header + mutable
print '==============='
print mutable

CAS3D_m_change_settings['mu'] = mutable.rstrip('\n')
CAS3D_m_change_settings['eta'] = etatable.rstrip('\n')
CAS3D_m_change_settings['xis'] = etatable.rstrip('\n')

# temp = fin.readline()
# while temp[:6] != '  m  n':
#     temp = fin.readline()

# while temp != '':	
#     fout.write(temp)
#     temp = fin.readline()

# fin.close()
# fout.close()
free_b_string = 'fixed'
if free_boundary: free_b_string='free'
suffix = '{}{}_{}_{}'.format(N, size_txt, free_b_string, density_profile)
suffix = '{}{}_{}_{}_{}'.format(N, size_txt, free_b_string, density_profile, prerun_surfaces)
cas3d_r_dir = work_directory+'/cas3d_r_n{}/'.format(suffix,)
cas3d_m_dir = work_directory+'/cas3d_m_n{}/'.format(suffix,)
#suffix_eig = '{}{}_{}_{}_{}'.format(N, size_txt, free_b_string, density_profile, density_profile_eig)
suffix_eig = '{}{}_{}_{}_{}_{}'.format(N, size_txt, free_b_string, density_profile, density_profile_eig, prerun_surfaces)
cas3d_single_eig_dir = work_directory+'/cas3d_m_n{}_ew_{:.4e}/'.format(suffix_eig, ew)

if cas3d_matrix:
    os.system('mkdir '+cas3d_m_dir)
    os.chdir(cas3d_m_dir)

    with file('/home/112/srh112/code_input_files/in_cas3d','r') as tmp:
        tmp_settings = tmp.read()
    for i in CAS3D_m_change_settings.keys():
        print i
        tmp_settings = tmp_settings.replace('<<'+i+'>>', CAS3D_m_change_settings[i])
    with file('in_cas3dm','w') as tmp:
        tmp.writelines(tmp_settings)
    os.system('ln -sf in_cas3dm input2')
    #os.system('ln -sf ../prerun/fort.12 fort.12')
    os.system('ln -sf {} fort.12'.format(fort12_loc))
    new_job_string = copy.deepcopy(job_string)
    for i in CAS3D_m_settings.keys():
        print i
        new_job_string = new_job_string.replace('<<'+i+'>>', CAS3D_m_settings[i])
    with file('jobscript','w') as tmp:
        tmp.writelines(new_job_string)
    with file('fort.77','w') as tmp:
        tmp.writelines('cas3d')
    submit_check_if_done()
    #proc = subprocess.Popen(["qsub", "jobscript"], stdout=subprocess.PIPE, shell=True)
    #(out, err) = proc.communicate()
    #print "program output:", out, err
    # status, output = commands.getstatusoutput("qsub jobscript")
    # job_id = output[:output.find('.')]
    # job_output_filename = 'jobscript.o'+job_id
    # job_error_filename = 'jobscript.e'+job_id
    # print status
    # print output
    # while not os.path.exists(job_output_filename):
    #     time.sleep(5)
    #     print 'not done'
    # time.sleep(5)
    # with file(job_error_filename,'r') as tmp: tmp_len = len(tmp.read())
    # if tmp_len>0: raise(Exception)

    #a = os.system('qsub jobscript')
    #print a

CAS3D_r_change_settings = CAS3D_m_change_settings
CAS3D_r_change_settings['evp'] = '2'
CAS3D_r_change_settings['matrix'] = '1'
CAS3D_r_change_settings['idens'] = idens_txt
#CAS3D_r_change_settings['idens'] = '4 1. 0. -3. 2.'
#CAS3D_r_change_settings['idens'] = '6 0.84827165 4.31756279 -35.54226959 76.85344266 -68.80442217 22.33691942'

# CAS3D_r_change_settings = {'surfaces':str(prerun_surfaces),
#                            'gamma':'1.6666',
#                            'evp':'2',
#                            'shift':'8.0329E-03',
#                            'shiftmin':'1.000e-06',
#                            'shiftmax':'2.000e-02',
#                            'number':'20000',
#                            'number_dens':'0.06',
#                            'IRAND':'2',
#                            'IFIX':'1',
#                            'idens':'2 1.  -1',
#                            'matrix':'0'}

if cas3d_eig:
    os.system('mkdir '+cas3d_r_dir)
    os.chdir(cas3d_r_dir)
    with file('/home/112/srh112/code_input_files/in_cas3d','r') as tmp:
        tmp_settings = tmp.read()
    for i in CAS3D_r_change_settings.keys():
        print i
        tmp_settings = tmp_settings.replace('<<'+i+'>>', CAS3D_r_change_settings[i])
    with file('in_cas3dr','w') as tmp:
        tmp.writelines(tmp_settings)
    os.system('ln -sf in_cas3dr input2')
    #os.system('ln -sf ../prerun/fort.12 fort.12')
    os.system('ln -sf {} fort.12'.format(fort12_loc,))
    os.system('ln -sf ../cas3d_m_n{}/fort.14 fort.14'.format(suffix,))
    for i in range(114,114+int(CAS3D_r_settings['ncpus'])):
        os.system('ln -sf fort.14 fort.'+str(i))
    new_job_string = copy.deepcopy(job_string)
    for i in CAS3D_r_settings.keys():
        print i
        new_job_string = new_job_string.replace('<<'+i+'>>', CAS3D_r_settings[i])
    with file('jobscript','w') as tmp:
        tmp.writelines(new_job_string)
    with file('fort.77','w') as tmp:
        tmp.writelines('cas3d')
    submit_check_if_done()

'''
new_job_string = copy.deepcopy(job_string)
for i in CAS3D_settings.keys():
    print i
    new_job_string = new_job_string.replace('<<'+i+'>>', CAS3D_settings[i])
print new_job_string

'''
conti_dir = work_directory+'/CONTI_n{}_{}_expt_ne/'.format(N, free_b_string)

if conti:
    os.system('mkdir '+conti_dir)
    os.chdir(conti_dir)
    with file('/home/112/srh112/code_input_files/in_conti','r') as tmp:
        tmp_settings = tmp.read()
    for i in conti_change_settings.keys():
        print i
        tmp_settings = tmp_settings.replace('<<'+i+'>>', conti_change_settings[i])
    with file('h1','w') as tmp:
        tmp.writelines(tmp_settings)
    #os.system('ln -sf ../prerun/fort.12 fort.12')
    os.system('ln -sf {} fort.12'.format(fort12_loc,))
    new_job_string = copy.deepcopy(job_string)
    for i in CONTI_settings.keys():
        print i
        new_job_string = new_job_string.replace('<<'+i+'>>', CONTI_settings[i])
    with file('jobscript','w') as tmp:
        tmp.writelines(new_job_string)
    submit_check_if_done()

CAS3D_r_change_settings = CAS3D_m_change_settings
CAS3D_r_change_settings['evp'] = '3'
CAS3D_r_change_settings['matrix'] = '1'
CAS3D_r_change_settings['idens'] = idens_eig_txt
CAS3D_r_change_settings['shift'] = '{:.4e}'.format(ew)

if cas3d_single_eig:
    os.system('mkdir '+cas3d_single_eig_dir)
    os.chdir(cas3d_single_eig_dir)
    with file('/home/112/srh112/code_input_files/in_cas3d','r') as tmp:
        tmp_settings = tmp.read()
    for i in CAS3D_r_change_settings.keys():
        print i
        tmp_settings = tmp_settings.replace('<<'+i+'>>', CAS3D_r_change_settings[i])
    with file('in_cas3dr','w') as tmp:
        tmp.writelines(tmp_settings)
    os.system('ln -sf in_cas3dr input2')
    os.system('ln -sf {} fort.12'.format(fort12_loc,))
    os.system('ln -sf ../cas3d_m_n{}/fort.14 fort.14'.format(suffix,))
    for i in range(114,114+int(CAS3D_single_eig_settings['ncpus'])):
        os.system('ln -sf fort.14 fort.'+str(i))
    new_job_string = copy.deepcopy(job_string)
    for i in CAS3D_single_eig_settings.keys():
        print i
        new_job_string = new_job_string.replace('<<'+i+'>>', CAS3D_single_eig_settings[i])
    with file('jobscript','w') as tmp:
        tmp.writelines(new_job_string)
    with file('fort.77','w') as tmp:
        tmp.writelines('cas3d')
    submit_check_if_done()
