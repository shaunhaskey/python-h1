#!/usr/bin/env python
#./run_script.py tree pause initial_trans_check manual_pause
#-- tree is the tree it 'runs' - usually h1data
#-- pause is a manual pause AFTER init - needed if there is no LAMwait. 1: include manual pause
#-- initial_trans_check, 1: include initial trans check, 0: exclude initial trans check
#-- manual_pause : 1:include pause before each init phase to allow settings to be changed etc...

#could set this up so that it creates its own jDispatcher.properties?
#read the text in the exception further below for a description of input arguments

import MDSplus
import time, re, datetime
import subprocess, os
import sys


#Cro is set to single sweep after data is saved
single_sweep = 0
verbose=1

include_transmitter_check = 0
store_rf_data_mdsplus = 1
transmitter_nums = [1,2]

try:
    tree = str(sys.argv[1])
    pause = int(sys.argv[2])
    initial_trans_check = int(sys.argv[3])
    manual_init = int(sys.argv[4])
except:
    print 'You need to specify a tree and whether or not you want a manual pause, if you want the initial transmitter check, and if you want a manual pause before init'
    print 'eg, run h1data without pause - i.e lamwait based, skip the initial transmitter check, no manual pause before init'
    print './run_script.py h1data 0 0 0'
    print 'eg, run sh_test3 with pause and include the initial transmitter check, include manual pause before init'
    print './run_script.py sh_test3 1 1 1'
    print 'exiting'
    sys.exit()

current_day = datetime.datetime.now().day

destination = 'prl40:8003'
local_comp = 'prl@prl40'
expected_shots = 300
computer_details = {}
computer_details['prl@prl24'] = {}
computer_details['prl@prl40'] = {}
computer_details['prl@prl45'] = {}
computer_details['prl@prl24']['servers'] = [('jServer',8010)]
computer_details['prl@prl40']['servers'] = [('jServer',8010),('jServer',8011),("/usr/local/mdsplus/bin/mdsip -s -p 8050 -h /usr/local/mdsplus/etc/mdsip.hosts",None)]
computer_details['prl@prl45']['servers'] = [('jServer',8010),('mdsip',8004),('jServer',8011)]
computer_details['prl@prl40']['programs'] = [('jDispatcherIp', tree),('jDispatchMonitor','prl40:8006'), ('run_script.py','python')]

def eval_tcl(cmd):
    '''wrapper for mdstcl
    '''
    MDSplus.Data.execute('tcl($)', cmd)

def dispatch_cmd(cmd, destination):
    '''create the mdstcl command, and then send it
    '''
    command_string = 'dispatch/command/server=%s %s\n' % (destination, cmd)
    now = datetime.datetime.now()
    print '%d-%d %d:%d:%d '%(now.day, now.month, now.hour, now.minute, now.second),
    print 'tcl:', command_string,
    eval_tcl(command_string)

def kill_servers(computer, prog_details, local=0):
    '''gets ps aux from 'computer', checks each line for any line that
    that has all the items in any of the touples in the prog_details list - if it does
    it kills that process as long as it isn't the process this script belongs to
    '''
    if local == 0:
        args = ['ssh', computer, 'ps aux']
    else:
        args = ['ps', 'aux']
    print computer, args
    output,error = subprocess.Popen(args,stdout = subprocess.PIPE, stderr= subprocess.PIPE,env={'LANG':'C'}).communicate()
    output = output.splitlines()
    for i in output:
        for tmp_prog_details in prog_details:
            truth = True
            for j in tmp_prog_details:
                truth = truth and i.find(str(j))>=0
            #if i.find(tmp_prog_details[0])>0 and i.find(str(tmp_prog_details[1]))>0:
            if truth:
                kill_pid = re.search('\d+',i).group()
                print computer, tmp_prog_details, kill_pid
                if local == 0 :
                    print 'remote kill called'
                    subprocess.call(['ssh', computer, 'kill -9 %s'%(kill_pid)])
                #make sure it doesnt try to kill itself
                elif local==1 and str(os.getpid())!=kill_pid:
                    print 'local kill called'
                    print os.getpid(), kill_pid
                    subprocess.call(['kill', '-9', kill_pid])
                else:
                    pass

def start_checking_transmitters(transmitter_numbers):
    '''starts the process on prl60 that checks the status of the transmitter
    numbers listed in transmitter_numbers
    '''
    a = []
    for i in transmitter_numbers:
        args = ['ssh', 'prl@prl60', '/home/prl/transmitter_snmp/check_trips2.py %d'%(i)]
        a.append(subprocess.Popen(args, stdout=subprocess.PIPE, 
                                  stderr=subprocess.PIPE, stdin=subprocess.PIPE))

    return a

def populate_snmp_rf(shot):
    '''starts the process on prl60 that populates MDSplus with transmitter settings etc...
    '''
    args = ['ssh', 'prl@prl60', '/home/prl/transmitter_snmp/rf_snmp_mdsplus.py %d'%shot]
    a = subprocess.Popen(args, stdout=subprocess.PIPE, 
                     stderr=subprocess.PIPE, stdin=subprocess.PIPE)
    return a


def finish_checking_transmitters(processes, transmitter_nums):
    '''Wait for the processes on prl60 that are checking the transmitters to finish
    and output the return, if it isn't zero, print out the stdout and stderr 
    which should include error messages
    '''
    return_list = []
    for i, a in enumerate(processes):
        b = None
        print 'Polling transmitter %d process :'%(transmitter_nums[i],),
        while b == None:
            b = a.poll()
            print '.',
            time.sleep(2)
        out, err = a.communicate()
        return_list.append(a.returncode)
        if return_list[-1] != 0:
            print out
            print err
    return return_list

def start_storing_cro(shot, single_sweep):
    '''start the process that stores data on the cro
    '''
    if verbose>0: 
        print('store CRO for shot {s}'.format(s=shot))

    args = ['ssh', 'prl@prl45', 'python /home/prl/python/save_shot.py %d %d'%(shot, single_sweep)]
    a = subprocess.Popen(args, stdout=subprocess.PIPE, 
                     stderr=subprocess.PIPE, stdin=subprocess.PIPE)
    return a

def finish_checking_cro(proc):
    '''currently not used, and it is blindly assumed the process finishes happily
    here for future use - needs to fit in with what Boyd does with the script on prl45
    '''
    # boyd's hack
    if proc is None: return(None, 0)
    else:
        if proc.poll():
            out, err = proc.communicate()
            print('***Error: Store CRO returns {err}'.format(err=err))
            if verbose>1: print('store CRO returns\n {std}'.format(std=out))
            return (None, proc.returncode)
        else:
            return(proc, 0)

    b = None
    print 'Polling cro store process:',
    while b == None:
        b = proc.poll()
        print '.',
        time.sleep(4)
    out, err = proc.communicate()
    return(proc.returncode)

def start_servers(computer, prog_details, local=0):
    '''If remote, make sure that the file jServer_XXX.sh or mdsip_XXX.sh where XXX is the port number exists
    because this essentially just runs that script. Weird work around for bash -l because ssh with command
    doesn't give a login shell - look into better ways of doing this later
    '''
    subproc_list = []
    for tmp_prog_details in prog_details:
        if local==0:
            args = ["xterm", "-T", "%s %s %d"%(computer, tmp_prog_details[0], tmp_prog_details[1]),"-e", "ssh", computer, "-t", "bash -l %s_%d.sh"%(tmp_prog_details[0], tmp_prog_details[1])]
        else:
            if tmp_prog_details[1]==None:
                args = ["xterm", "-T", "%s %s"%(computer, tmp_prog_details[0]),"-e", tmp_prog_details[0]]
            else:
                args = ["xterm", "-T", "%s %s %d"%(computer, tmp_prog_details[0], tmp_prog_details[1]),"-e", tmp_prog_details[0], str(tmp_prog_details[1])]
            #args = ["xterm", "-T", "%s jServer %d"%(computer,port),"-e", "jServer",str(port)]
        if tmp_prog_details[1]==None:
            print 'start %s on %s '%(tmp_prog_details[0], computer), args
        else:
            print 'start %s on %s : %d '%(tmp_prog_details[0], computer, tmp_prog_details[1]), args
        subproc_list.append(subprocess.Popen(args, stdout=subprocess.PIPE, 
                                             stderr=subprocess.PIPE, stdin=subprocess.PIPE))
    return subproc_list

def kill_everything(computer_details):
    '''kill all the processes that were created as part of this script
    so that there is no mess left behind for a future dispatcher
    '''
    print 'killing all subprocesses'
    for i in computer_details.keys():
        for i in computer_details[i]['processes']:
            i.kill()
    print 'Finished cleaning up'


##################################################################################
#Start of script
#kill all jServers
print '='*8, 'Checking for and kill all servers that could conflict', '='*8
for i in computer_details.keys():
    if i==local_comp:
        kill_servers(i, computer_details[i]['servers'], local=1)
    else:
        kill_servers(i, computer_details[i]['servers'], local=0)

#start all jServers ands mdsip servers that are required
print '='*8, 'Starting all jServers on the required computers', '='*8
for i in computer_details.keys():
    if i==local_comp:
        computer_details[i]['processes'] = start_servers(i, computer_details[i]['servers'], local=1)
    else:
        computer_details[i]['processes'] = start_servers(i, computer_details[i]['servers'], local=0)

time.sleep(2)
#see if any jDispatcherIp and jDispatchMonitors are running and kill them if they are
print '='*8, 'Killing jDispatcherIp and jDispatchMonitor', '='*8
kill_servers(local_comp, computer_details[local_comp]['programs'], local=1)

#Start a local jDispatcherIp and jDispatchMonitor
print '='*8, 'Starting jDispatcherIp and jDispatchMonitor', '='*8
args = ["xterm", "-T", "jDispatcherIp", "-e", "jDispatcherIp", tree]
computer_details[local_comp]['processes'].append(subprocess.Popen(args, stdout=subprocess.PIPE, 
                                                                  stderr=subprocess.PIPE, stdin=subprocess.PIPE))

time.sleep(3)
args = ["xterm", "-T", "jDispatchMonitor", "-e", "jDispatchMonitor", "prl40:8006"]
computer_details[local_comp]['processes'].append(subprocess.Popen(args, stdout=subprocess.PIPE, 
                                                                  stderr=subprocess.PIPE, stdin=subprocess.PIPE))


#Everything should be running now, give it a bit of time to get going
time.sleep(5)

#loop for running all the shots
print '='*8, 'Starting shot looper', '='*8
cur_shot =  MDSplus.Tree.getCurrent(tree)
print 'Current_shot:', cur_shot
executing_shot = cur_shot + 1
killed_all = 0 
cro_proc = None

try:
    for i in range(0, expected_shots):
        #check the transmitters if asked to... if this is the first shot, it needs to initialise the
        #transmitter check process and check its return, if it isn't hte first shot, it just needs to 
        #check the transmitter return because it was started during tehe store phase
        if include_transmitter_check:
            if (i == 0) and (initial_trans_check):
                print 'starting transmitter check'
                transmitter_process = start_checking_transmitters(transmitter_nums)
            if (i != 0) or (initial_trans_check):
                print 'Checking transmitter return code'
                return_list = finish_checking_transmitters(transmitter_process, transmitter_nums)
                print 'return codes:', return_list
            else:
                print 'skipping tramsitter check'
                return_list = [0]
            for return_code in return_list:
                if return_code != 0:
                    print '!!! WARNING transmitter error - FIX IT then press any key'
                    raw_input()
                    print 'continuing'

        #This is if you want to change settings before the init phase, i.e for tuning antennas etc...
        #requires intervention before the init phase starts 
        if manual_init:
            print 'Press enter to start init phase :',
            raw_input()

        #create pulse and perform init phase
        print '='*5,'Current_shot: %d, executing shot: %d'%(cur_shot, executing_shot),'='*5
        dispatch_cmd('set tree %s' % tree, destination)
        dispatch_cmd('create pulse %d' % executing_shot, destination)

        time.sleep(3)
        #print 'press enter.......'
        #raw_input()

        #maybe remove this???!
        if tree=='h1data':
            tmp_tree = MDSplus.Tree('h1data', executing_shot)
            #print 'opened tree'
            pulse_created_node = tmp_tree.getNode('.log.statistics:init_done')
            #print 'told tree it has been created'
            pulse_created_node.putData(MDSplus.makeData('1'))
            #print 'told tree it has been created'


        dispatch_cmd('dispatch /build', destination)
        tmp_time = time.time()
        print "sending event"
	# only send EXECUTING_SHOT as CURRENT_SHOT is EXECUTING_SHOT -1
        MDSplus.Event.setevent("EXECUTING_SHOT",executing_shot)
        dispatch_cmd('dispatch /phase INIT', destination)
        print '--> time to perform INIT %.2fs'%(time.time() - tmp_time)

        #SH 20May2013 - testing abort phase
        # print 'Init command has finished for some reason....'
        # print 'press any key to redo the INIT phase'
        # raw_input()
        # dispatch_cmd('set tree %s' % tree, destination)
        # dispatch_cmd('create pulse %d' % (executing_shot), destination)
        # time.sleep(3)
        # dispatch_cmd('dispatch /build', destination)
        # time.sleep(3)
        # dispatch_cmd('dispatch /phase INIT', destination)
        # print '--> time to perform INIT %.2fs'%(time.time() - tmp_time)
        # print 'Init command has finished for some reason....'
        # print 'press any key to redo the INIT phase'
        # raw_input()
        # dispatch_cmd('dispatch /phase INIT', destination)
        ##########


        #wait for user intervention before starting store phase
        #usually only useful if the init phase doesn't wait for a trigger - i.e lamwait
        (cro_proc, retcode) = finish_checking_cro(cro_proc)
        if pause==1:
            print 'Initialised, press enter when triggered to store :',
            raw_input()
        else:
            time.sleep(2)

        tmp_time = time.time()

        #start checking the transmitter status if required
        if include_transmitter_check:
            print 'starting transmitter check'
            transmitter_process = start_checking_transmitters(transmitter_nums)

        #start the store phase, send MDSEvent for LabVIEW Acq devices added MB 13/3/13
	MDSplus.Event.setevent("STORE_DATA_PHASE")
        #dispatch_cmd('close tree', destination)
        #dispatch_cmd('set tree %s' % tree, destination)


        dispatch_cmd('dispatch /phase STORE', destination)
        print '--> time to perform STORE %.2fs'%(time.time() - tmp_time)
        dispatch_cmd('close tree', destination)
        
        #finished store, increment current shot
        MDSplus.Tree.setCurrent(tree, executing_shot)
        #cur_shot =  MDSplus.Tree.getCurrent(tree)
        now = datetime.datetime.now()
        print '%d-%d %d:%d:%d '%(now.day, now.month, now.hour, now.minute, now.second),
        print '%s current shot set to : %d'%(tree, executing_shot)

        #start storing the cro data
        cro_proc = start_storing_cro(executing_shot, single_sweep)

        #start storing the rf data in MDSplus
        if store_rf_data_mdsplus:
            rf_pop_proc = populate_snmp_rf(executing_shot)

        # break to force an overnight pause so too many shots don't occur
        if (current_day - now.day) != 0:
            print '%d-%d-%d %d:%d : Midnight has passed! Press enter to continue:'%(now.day, now.month, now.year, now.hour, now.minute)
            raw_input()
            current_day = datetime.datetime.now().day
            print '--> continuing'
        executing_shot += 1

except Exception, details:
    #this catches ctrl-c or anything else, and starts killing all processes it started
    #so that its easy for the next script to start everything
    print '!!!! Some kind of exception !!!!!'
    print(details)
    kill_everything(computer_details)
    killed_all = 1
    print 'Finished with exception'

#tidy up, but check to see if it already happened first
#this one will probably happen because you have run out of shots
#might be better to reset 'shot clock' and go back to the beginning
if killed_all != 1:
    print 'killing all subprocesses'
    kill_everything(computer_details)
    print 'Program finished'

