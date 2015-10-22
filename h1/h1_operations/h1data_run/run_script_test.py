#!/usr/bin/env python
import MDSplus
import time, re, datetime
import subprocess, os
import sys

#could set this up so that it creates its own jDispatcher.properties?

try:
    tree = str(sys.argv[1])
    pause = int(sys.argv[2])
except:
    print 'You need to specify a tree and whether or not you want a manual pause'
    print 'eg, run h1data without pause - i.e lamwait based'
    print './run_script.py h1data 0'
    print 'eg, run sh_test3 with pause'
    print './run_script.py sh_test3 1'
    print 'exiting'
    sys.exit()

#pause = 0
#tree = 'h1data'
#tree = 'sh_test3'

current_day = datetime.datetime.now().day

destination = 'prl40:8003'
local_comp = 'prl@prl40'
expected_shots = 100
computer_details = {}
computer_details['prl@prl24'] = {}
computer_details['prl@prl40'] = {}
computer_details['prl@prl45'] = {}
computer_details['prl@prl24']['servers'] = [('jServer',8010)]
computer_details['prl@prl40']['servers'] = [('jServer',8010),('jServer',8011)]
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
            args = ["xterm", "-T", "%s %s %d"%(computer, tmp_prog_details[0], tmp_prog_details[1]),"-e", tmp_prog_details[0], str(tmp_prog_details[1])]
            #args = ["xterm", "-T", "%s jServer %d"%(computer,port),"-e", "jServer",str(port)]
        print 'start %s on %s : %d '%(tmp_prog_details[0], computer, tmp_prog_details[1]), args
        subproc_list.append(subprocess.Popen(args, stdout=subprocess.PIPE, 
                                             stderr=subprocess.PIPE, stdin=subprocess.PIPE))
    return subproc_list

def kill_everything(computer_details):
    '''kill all the processes that were created as part of this script
    '''
    print 'killing all subprocesses'
    for i in computer_details.keys():
        for i in computer_details[i]['processes']:
            i.kill()
    print 'Finished cleaning up'


#kill all jServers
print '='*8, 'Checking for and kill all servers that could conflict', '='*8
for i in computer_details.keys():
    if i==local_comp:
        kill_servers(i, computer_details[i]['servers'], local=1)
    else:
        kill_servers(i, computer_details[i]['servers'], local=0)

#start all jServers ands mdsip servers
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


#Everything should be running now, give it a bit of time
time.sleep(5)

#loop for running all the shots
print '='*8, 'Starting shot looper', '='*8
cur_shot =  MDSplus.Tree.getCurrent(tree)
print 'Current_shot:', cur_shot
executing_shot = cur_shot + 1
killed_all = 0 
try:
    for i in range(0, expected_shots):
        print '='*5,'Current_shot: %d, executing shot: %d'%(cur_shot, executing_shot),'='*5
        dispatch_cmd('set tree %s' % tree, destination)
        dispatch_cmd('create pulse %d' % executing_shot, destination)
        time.sleep(5)
        dispatch_cmd('dispatch /build', destination)
        tmp_time = time.time()
        dispatch_cmd('dispatch /phase INIT', destination)
        print '--> time to perform INIT %.2fs'%(time.time() - tmp_time)
        if pause==1:
            print 'Initialised, press enter when triggered to store :',
            raw_input()
        else:
            time.sleep(5)
        tmp_time = time.time()
        dispatch_cmd('dispatch /phase STORE', destination)
        print '--> time to perform STORE %.2fs'%(time.time() - tmp_time)
        dispatch_cmd('close tree', destination)
        MDSplus.Tree.setCurrent(tree, executing_shot)
        cur_shot =  MDSplus.Tree.getCurrent(tree)
        now = datetime.datetime.now()
        print '%d-%d %d:%d:%d '%(now.day, now.month, now.hour, now.minute, now.second),
        print '%s current shot set to : %d'%(tree, cur_shot)

        # break to force an overnight pause
        if (current_day - now.day) != 0:
            print '%d-%d-%d %d:%d : Midnight has passed! Press enter to continue:'%(now.day, now.month, now.year, now.hour, now.minute)
            raw_input()
            current_day = datetime.datetime.now().day
            print '--> continuing'
        executing_shot += 1

except:
    print '!!!! Some kind of exception !!!!!'
    kill_everything(computer_details)
    killed_all = 1
    print 'Finished with exception'

#tidy up, but check to see if it already happened first
if killed_all != 1:
    print 'killing all subprocesses'
    kill_everything(computer_details)
    print 'Program finished'

