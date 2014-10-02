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
import threading
import SocketServer

import pygtk
import pango
import gobject
pygtk.require('2.0')
import gtk
import atexit

#Cro is set to single sweep after data is saved
single_sweep = 0
verbose= 1

#key should be the same as the label next to the checkbox
#key must also exist in disable_buttons below to set the default value
#item is a list of nodes that are to be disabled/enabled for that particular key
disable_nodes = {'mirnov' : ['.mirnov.ACQ132_7', 
                             '.mirnov.ACQ132_8', 
                             '.mirnov.ACQ132_9', 
                             '.mirnov.HMA_AMPS'],
                 'electron density': ['.electr_dens.camac.TR612_7',
                                      '.electr_dens.camac.TR612_8',
                                      '.electr_dens.camac.TR612_9',
                                      '.electr_dens.camac.TR612_10'],
                 'Reentrant Camera':['.mirnov.cam_alex'],
                 'ISOPlane spectr':['.mirnov.pimax']}

#Default values for the disable_nodes above
disable_buttons = {'mirnov': True,
                   'electron density': True,
                   'Reentrant Camera':True,
                   'ISOPlane spectr':True}

class GUI_control():
    def __init__(self):
        '''Supposed to be GUI control of the run script
        SRH: 19July2013 
        '''

        #setup the GUI with a button and 3 labels
        self.window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        self.window.connect("delete_event", self.delete_event)
        self.window.connect("destroy", self.destroy)
        self.window.set_border_width(10) 
        self.vbox = gtk.VBox(homogeneous=False, spacing=0)
        self.l_phase = gtk.Label('Phase :')
        self.l_executing = gtk.Label('Executing :')
        self.l_current = gtk.Label('Current :')

        self.l_phase.modify_font(pango.FontDescription("sans 48"))
        self.l_executing.modify_font(pango.FontDescription("sans 48"))
        self.l_current.modify_font(pango.FontDescription("sans 48"))

        self.check_buttons = {}
        for i in disable_nodes:self.check_buttons[i] = gtk.CheckButton(i)
        tmp_label = gtk.Label('Ticked is enabled')
        self.button = gtk.Button("Reinitialise Shot")
        self.vbox.pack_start(self.button, expand=True, fill=True, padding=0)
        self.vbox.pack_start(self.l_phase, expand=True, fill=True, padding=0)
        self.vbox.pack_start(self.l_executing, expand=True, fill=True, padding=0)
        self.vbox.pack_start(self.l_current, expand=True, fill=True, padding=0)

        self.vbox.pack_start(tmp_label, expand=True, fill=True, padding=0)
        for i in self.check_buttons.keys(): self.vbox.pack_start(self.check_buttons[i], expand=True, fill=True, padding=0)

        self.button.connect("clicked", self.reinit, None)
        self.window.add(self.vbox)
        self.vbox.show()
        self.button.show()
        for i in self.check_buttons.keys(): self.check_buttons[i].set_active(True)
        for i in self.check_buttons.keys(): self.check_buttons[i].show()
        self.l_phase.show()
        self.l_executing.show()
        self.l_current.show()
        tmp_label.set_justify(gtk.JUSTIFY_LEFT)
        tmp_label.show()
        self.window.show()

    def reinit(self, widget, data=None):
        print "Reinitialise called"
        self.shot_cycle.abort_shot()

    def delete_event(self, widget, event, data=None):
        print "delete event occurred"
        return False

    def destroy(self, widget, data=None):
        print "destroy signal occurred"
        gtk.main_quit()

    def main(self):
        gtk.main()



class ThreadedTCPRequestHandler(SocketServer.BaseRequestHandler):
    def handle(self):
        termination_character = '\r\n'
        sep_char = chr(0x09) 
        data = ''
        count = 0
        while (not data.endswith(termination_character)) or count>10:
            tmp = self.request.recv(1024)
            data += tmp
            count+=1
        #data = self.request.recv(1024)
        print('hello')
        print(data)
        data = data.rstrip(termination_character)
        if data == "LOOPER_STATUS":
            response = '{}{}{}{}'.format(self.server.shot_cycle.cycle, sep_char, self.server.shot_cycle.executing_shot, termination_character)
        else:
            cur_thread = threading.current_thread()
            response = "{}: {} : {} :{}".format(cur_thread.name, data, self.server.shot_cycle.cycle, self.server.shot_cycle.executing_shot)
        self.request.sendall(response)
        if data.find('ABORT')>=0: self.server.shot_cycle.abort_shot()
class ThreadedTCPServer(SocketServer.ThreadingMixIn, SocketServer.TCPServer):
    pass

def start_TCP_server(PORT):
    '''Starts up the TCP server in a separate thread
    SRH: 19July2013 
    '''
    HOST = ""
    server = ThreadedTCPServer((HOST, PORT), ThreadedTCPRequestHandler, False)
    server.allow_reuse_address = True # Prevent 'cannot bind to address' errors on restart
    server.server_bind()     # Manually bind, to support allow_reuse_address
    server.server_activate() # (see above comment)
    # Start a thread with the server -- that thread will then start one
    # more thread for each request
    server_thread = threading.Thread(target=server.serve_forever)
    # Exit the server thread when the main thread terminates
    server_thread.daemon = True
    server_thread.start()
    print "status server loop running in thread:", server_thread.name
    return server_thread, server



def eval_tcl(cmd):
    '''wrapper for mdstcl
    '''
    MDSplus.Data.execute('tcl($)', cmd)

def dispatch_cmd(cmd, destination):
    '''create the mdstcl command, and then send it using MDSplus module
    '''
    command_string = 'dispatch/command/server=%s %s\n' % (destination, cmd)
    now = datetime.datetime.now()
    print '%d-%d %d:%d:%d '%(now.day, now.month, now.hour, now.minute, now.second),
    print 'tcl:', command_string,
    eval_tcl(command_string)

def dispatch_cmd2(cmd, destination):
    '''create the mdstcl command, and then send it using os.system
    '''
    string_is = "mdstcl 'dispatch/command/server={} {}'".format(destination, cmd)
    print(string_is)
    a = os.system(string_is)
    print('finished...', a)
    # command_string = 'dispatch/command/server=%s %s\n' % (destination, cmd)
    # now = datetime.datetime.now()
    # print '%d-%d %d:%d:%d '%(now.day, now.month, now.hour, now.minute, now.second),
    # print 'tcl:', command_string,
    # eval_tcl(command_string)


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

    kill_list = []
    for i in output:
        for tmp_prog_details in prog_details:
            truth = True
            for j in tmp_prog_details:
                truth = truth and i.find(str(j))>=0 and i.find('tee')<0 and i.find('grep')<0
            #if i.find(tmp_prog_details[0])>0 and i.find(str(tmp_prog_details[1]))>0:
            if truth:
                kill_pid = i.split()[1]
                #kill_pid = re.search('\d+',i).group()
                print computer, tmp_prog_details, kill_pid
                if local==1 and str(os.getpid()) == kill_pid:
                    pass
                else:
                    kill_list.append(kill_pid)
    if local == 0 and len(kill_list) > 0:
        print 'remote kill called'
        subprocess.call(['ssh', computer, 'kill -9 {}'.format(' '.join(kill_list))])
    #make sure it doesnt try to kill itself
    elif local==1 and len(kill_list) > 0:
        print 'local kill called'
        #print os.getpid(), kill_pid
        #print kill_list
        args_list = ['kill', '-9']
        args_list.extend(kill_list)
        print args_list
        subprocess.call(args_list)
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

def start_camera(shot,):
    '''start the process that stores data on the cro
    '''
    if verbose>0: 
        print('Storing breakdown webcam for shot {s}'.format(s=shot))
    args = ['ssh', 'prl@prl45', '/home/prl/python/webcam_capture_burst shot={shot}'.format(shot=shot)]
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

def start_servers(computer, prog_details, local=0, execute = True):
    '''If remote, make sure that the file jServer_XXX.sh or mdsip_XXX.sh where XXX is the port number exists
    because this essentially just runs that script. Weird work around for bash -l because ssh with command
    doesn't give a login shell - look into better ways of doing this later
    '''
    subproc_list = []
    gnome_term_list = []
    for tmp_prog_details in prog_details:
        if local==0:
            args = ["xterm", "-T", "%s %s %d"%(computer, tmp_prog_details[0], tmp_prog_details[1]),"-e", "ssh", computer, "-t", "bash -l %s_%d.sh"%(tmp_prog_details[0], tmp_prog_details[1])]
            args = ["ssh", computer, "-t", "bash -l '%s_%d.sh'"%(tmp_prog_details[0], tmp_prog_details[1])]
            term_title = "%s %s %d"%(computer[computer.find('@')+1:], tmp_prog_details[0], tmp_prog_details[1])
        else:
            if tmp_prog_details[1]==None:
                args = ["xterm", "-T", "%s %s"%(computer, tmp_prog_details[0]),"-e", tmp_prog_details[0]]
                args = [tmp_prog_details[0]]
                term_title = "%s %s"%(computer[computer.find('@')+1:], tmp_prog_details[0]),"-e", tmp_prog_details[0]
            else:
                args = ["xterm", "-T", "%s %s %d"%(computer, tmp_prog_details[0], tmp_prog_details[1]),"-e", tmp_prog_details[0], str(tmp_prog_details[1])]
                args = [tmp_prog_details[0], str(tmp_prog_details[1])]
                term_title = "%s %s %d"%(computer[computer.find('@')+1:], tmp_prog_details[0], tmp_prog_details[1])
            #args = ["xterm", "-T", "%s jServer %d"%(computer,port),"-e", "jServer",str(port)]
        if tmp_prog_details[1]==None:
            print 'start %s on %s '%(tmp_prog_details[0], computer), args
        else:
            print 'start %s on %s : %d '%(tmp_prog_details[0], computer, tmp_prog_details[1]), args
        gnome_term_list.extend(['--tab','-t',term_title,'-e'])
        gnome_term_list.extend(['''{}'''.format(' '.join(args))])
        if execute:
            subproc_list.append(subprocess.Popen(args, stdout=subprocess.PIPE, 
                                                 stderr=subprocess.PIPE, stdin=subprocess.PIPE))
    if execute:
        return subproc_list, gnome_term_list
    else:
        return [], gnome_term_list


def kill_everything(computer_details):
    '''kill all the processes that were created as part of this script
    so that there is no mess left behind for a future dispatcher
    '''
    print 'killing all subprocesses'
    for i in computer_details.keys():
        for i in computer_details[i]['processes']:
            i.kill()
    print 'Finished cleaning up'


class shot_cycle():
    def __init__(self, cur_shot, executing_shot, include_transmitter_check, manual_init, tree, destination, transmitter_nums, single_sweep, store_rf_data_mdsplus, initial_trans_check, cro_proc):
        '''This class contains all the information about the shot, and
        the various phases methods for each of the phases. It needs
        another object to call its init_shot and store_shot, and abort
        phases

        SRH: 19July2013 
        '''
        print('Hello world, the shot cycle has started')
        self.cur_shot = cur_shot
        self.executing_shot = executing_shot
        self.include_transmitter_check = include_transmitter_check
        self.manual_init = manual_init
        self.tree = tree
        self.destination = destination
        self.shot_count = 0
        self.transmitter_nums = transmitter_nums
        self.single_sweep = single_sweep
        self.current_day = datetime.datetime.now().day
        self.initial_trans_check = initial_trans_check
        self.abort = 0
        self.store_rf_data_mdsplus = store_rf_data_mdsplus
        self.cro_proc = cro_proc
    def update_label(self,):
        '''update the text on the GUI
        SRH: 19July2013
        '''
        self.GUI.l_phase.set_text('Phase :' + self.cycle)
        self.GUI.l_current.set_text('Current :' + str(self.cur_shot))
        self.GUI.l_executing.set_text('Executing :' + str(self.executing_shot))

    def init_shot(self,):
        '''Perform all the required tasks to initialise the shot
        SRH: 19July2013
        '''
        self.cycle = 'INIT'
        gobject.idle_add(self.update_label,)

        #Disable nodes in the tree as requested with the check boxes
        T = MDSplus.Tree('h1data', -1)
        for i in disable_nodes.keys():
            for j in disable_nodes[i]:
                N = T.getNode(j)
                #Need to check the status of the relevant checkbutton and set it accordingly
                #print i, j, disable_buttons[i]
                N.setOn(self.GUI.check_buttons[i].get_active())

        if self.include_transmitter_check:
            if (self.shot_count == 0) and (self.initial_trans_check):
                print 'starting transmitter check'
                self.transmitter_process = start_checking_transmitters(self.transmitter_nums)
            if (self.shot_count != 0) or (self.initial_trans_check):
                print 'Checking transmitter return code'
                return_list = finish_checking_transmitters(self.transmitter_process, self.transmitter_nums)
                print 'return codes:', return_list
            else:
                print 'skipping tramsitter check'
                return_list = [0]
            for return_code in return_list:
                if return_code != 0:
                    print '!!! WARNING transmitter error - FIX IT then press any key'
                    raw_input()
                    print 'continuing'
        print '='*5,'Current_shot: %d, executing shot: %d'%(self.cur_shot, self.executing_shot),'='*5
        dispatch_cmd('set tree %s' % self.tree, self.destination)
        dispatch_cmd('create pulse %d' % self.executing_shot, self.destination)
        (self.cro_proc, retcode) = finish_checking_cro(self.cro_proc)
        time.sleep(3)
        start_camera(self.executing_shot)
        #Need to comment this out - for some reason the getNode part of this throws an
        #Error after an abort... not sure why
        # if self.tree=='h1data':
        #     tmp_tree = MDSplus.Tree('h1data', self.executing_shot)
        #     pulse_created_node = tmp_tree.getNode('.log.statistics:init_done')
        #     pulse_created_node.putData(MDSplus.makeData('1'))
        dispatch_cmd('dispatch /build', self.destination)
        tmp_time = time.time()
        print "sending event"
	# only send EXECUTING_SHOT as CURRENT_SHOT is EXECUTING_SHOT -1
        MDSplus.Event.setevent("EXECUTING_SHOT",self.executing_shot)
        dispatch_cmd('dispatch /phase INIT', self.destination)
        print '--> time to perform INIT %.2fs'%(time.time() - tmp_time)
        if self.abort!=0:print('Shot was aborted')
        return self.abort

    def store_shot(self,):
        '''Perform all the required tasks to store a shot
        SRH: 19July2013
        '''
        #SRH MINISTER VISIT
        #time.sleep(18)
        self.cycle = 'STORE'
        gobject.idle_add(self.update_label,)
        if self.include_transmitter_check:
            print 'starting transmitter check'
            self.transmitter_process = start_checking_transmitters(self.transmitter_nums)
	MDSplus.Event.setevent("STORE_DATA_PHASE")
        tmp_time = time.time()

        #start storing the cro data
        self.cro_proc = start_storing_cro(self.executing_shot, self.single_sweep)

        dispatch_cmd2('ABORT PHASE', self.destination)
        dispatch_cmd('close tree', self.destination)
        dispatch_cmd('set tree %s /shot=%d'%(self.tree,self.executing_shot), self.destination)
        dispatch_cmd('dispatch /build', self.destination)

        dispatch_cmd('dispatch /phase STORE', self.destination)
        print '--> time to perform STORE %.2fs'%(time.time() - tmp_time)
        dispatch_cmd2('ABORT PHASE', self.destination)
        dispatch_cmd('close tree', self.destination)
        
        #finished store, increment current shot
        MDSplus.Tree.setCurrent(self.tree, self.executing_shot)

        #cur_shot =  MDSplus.Tree.getCurrent(tree)
        self.cur_shot =  MDSplus.Tree.getCurrent(self.tree)
        now = datetime.datetime.now()
        print '%d-%d %d:%d:%d '%(now.day, now.month, now.hour, now.minute, now.second),
        print '%s current shot set to : %d'%(self.tree, self.executing_shot)

        #Update the summary database
        #args = ['curl','-X', 'POST', '--data', '''"shot={}"'''.format(self.cur_shot), 'http://h1ds.anu.edu.au/data/h1/latest/']
        #print args
        #sum_update_proc = subprocess.Popen(args, stdout=subprocess.PIPE,stderr=subprocess.PIPE, stdin=subprocess.PIPE)
        #sum_update_out, sum_update_err = sum_update_proc.communicate()
        #print sum_update_err
        time.sleep(3)
        command = """curl -X POST --data "shot={}" http://h1ds.anu.edu.au/data/h1/latest/ > curl_log""".format(self.cur_shot)
        print('doing curl ', command)
        os.system(command)
#        #start storing the cro data
#        cro_proc = start_storing_cro(self.executing_shot, self.single_sweep)

        #start storing the rf data in MDSplus
        if self.store_rf_data_mdsplus:
            rf_pop_proc = populate_snmp_rf(self.executing_shot)

        # break to force an overnight pause so too many shots don't occur
        if (self.current_day - now.day) != 0:
            print '%d-%d-%d %d:%d : Midnight has passed! Press enter to continue:'%(now.day, now.month, now.year, now.hour, now.minute)
            raw_input()
            self.current_day = datetime.datetime.now().day
            print '--> continuing'
        self.executing_shot += 1
        self.shot_count += 1
        return 1

    def abort_shot(self,):
        '''Perform all the required tasks to abort a shot
        SRH: 19July2013
        '''
        if self.cycle == 'INIT':
            print('aborting method called')
            self.abort = 1
            dispatch_cmd2('ABORT PHASE', self.destination)
            dispatch_cmd2('close tree', self.destination)
        else:
            print('Cannot abort if not in the INIT phase')
        
class shot_engine():
    def __init__(self, shot_cycle, expected_shots):
        '''This cycles shot_cycle through the shots

        SRH: 19July2013
        '''
        self.shot_cycle = shot_cycle
        self.expected_shots = expected_shots
        print 'shot_engine initialised'

    def start(self,):
        print 'Hello world, the shot engine has started'
        for i in range(0, self.expected_shots):
            a = 1
            while a==1:
                print('cycle',a)
                a = self.shot_cycle.init_shot()
                self.shot_cycle.abort = 0
            b = self.shot_cycle.store_shot()


if __name__=='__main__':
    #Cro is set to single sweep after data is saved
    print 'using python-h1 version'
    single_sweep = 1
    include_transmitter_check = 0
    store_rf_data_mdsplus = 1
    transmitter_nums = [1,2]
    cro_proc = None

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
    destination = 'h1svr:8003'
    #local_comp = 'prl@prl40'
    local_comp = 'h1operator@h1svr'
    expected_shots = 300
    computer_details = {}
    computer_details['prl@prl24'] = {}
    computer_details['prl@prl40'] = {}
    computer_details['h1operator@h1svr'] = {}
    computer_details['prl@prl45'] = {}
    computer_details['datasys@h1svr'] = {}
    computer_details['prl@prl24']['servers'] = [('jServer',8010)]
    computer_details['datasys@h1svr']['servers'] = [('jServer',8010)]
    computer_details['prl@prl40']['servers'] = [('jServer',8010),('jServer',8011)]#,("/usr/local/mdsplus/bin/mdsip -s -p 8050 -h /usr/local/mdsplus/etc/mdsip.hosts",None)]
    computer_details['prl@prl45']['servers'] = [('jServer',8010),('mdsip',8004),('jServer',8011), ('jServer',8012)]
    computer_details['h1operator@h1svr']['programs'] = [('/usr/local/mdsplus/bin/jDispatcherIp', tree),('/usr/local/mdsplus/bin/jDispatchMonitor','prl40:8006'), ('run_script.py','python')]

    computer_details['h1operator@h1svr']['programs'] = [('jDispatcherIp', tree),('jDispatchMonitor','h1svr:8006')]#, ('run_script.py','python')]

    computer_details['h1operator@h1svr']['servers'] = []

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

    gnome_list = ['gnome-terminal']
    for i in computer_details.keys():
        if i==local_comp:
            computer_details[i]['processes'], tmp_list = start_servers(i, computer_details[i]['servers'], local=1, execute = False)
        else:
            computer_details[i]['processes'], tmp_list = start_servers(i, computer_details[i]['servers'], local=0, execute = False)
        gnome_list.extend(tmp_list)


    #time.sleep(2)
    #see if any jDispatcherIp and jDispatchMonitors are running and kill them if they are
    print '='*8, 'Killing jDispatcherIp and jDispatchMonitor', '='*8
    kill_servers(local_comp, computer_details[local_comp]['programs'], local=1)
    time.sleep(3)

    #Start a local jDispatcherIp and jDispatchMonitor
    print '='*8, 'Starting jDispatcherIp and jDispatchMonitor', '='*8
    args = ["xterm", "-T", "jDispatcherIp", "-e", "jDispatcherIp", tree]
    args = ["--tab", "-t","jDispatcherIp", "-e", "bash -c 'sleep 3;jDispatcherIp {}'".format(tree)]
    #computer_details[local_comp]['processes'].append(subprocess.Popen(args, stdout=subprocess.PIPE, 
    #                                                                  stderr=subprocess.PIPE, stdin=subprocess.PIPE))
    #subprocess.call(gnome_string)
    #time.sleep()
    gnome_list.extend(args)

    args = ["xterm", "-T", "jDispatchMonitor", "-e", "jDispatchMonitor", "h1svr:8006"]
    args = ["--tab", '-t','jDispatchMonitor', "-e", "bash -c 'sleep 7;jDispatchMonitor {}'".format('h1svr:8006')]
    #computer_details[local_comp]['processes'].append(subprocess.Popen(args, stdout=subprocess.PIPE, 
    #                                                                  stderr=subprocess.PIPE, stdin=subprocess.PIPE))

    gnome_list.extend(args)

    print gnome_list
    subproc = subprocess.Popen(gnome_list, stdout=subprocess.PIPE, 
                               stderr=subprocess.PIPE, stdin=subprocess.PIPE)

    for i in computer_details.keys():
        computer_details[i]['processes'] = [subproc]
    computer_details[local_comp]['processes'] = [subproc]

    #time.sleep(6)


    #Everything should be running now, give it a bit of time to get going
    time.sleep(7)

    #loop for running all the shots
    print '='*8, 'Starting shot looper', '='*8
    cur_shot =  MDSplus.Tree.getCurrent(tree)
    print 'Current_shot:', cur_shot
    executing_shot = cur_shot + 1
    killed_all = 0
    shot_cycle = shot_cycle(cur_shot, executing_shot, include_transmitter_check, manual_init, tree, destination, transmitter_nums, single_sweep, store_rf_data_mdsplus, initial_trans_check, cro_proc)

    TCP_PORT = 8111
    TCP_server_thread, TCP_server = start_TCP_server(TCP_PORT)

    #attach the shot_cycle object to the TCP_server object so it can modify it
    TCP_server.shot_cycle = shot_cycle

    #Need this to allow threads in addition to gtk GUI
    gobject.threads_init()

    #Create the GUI
    GUI = GUI_control()

    #Attach various objects to the GUI so it can modify them - i.e abort
    GUI.shot_cycle = shot_cycle
    GUI.TCP_server = TCP_server
    GUI.shot_cycle = shot_cycle
    GUI.expected_shots = expected_shots

    #Attach the GUI to shot_cycle so it can modify the GUI
    shot_cycle.GUI = GUI

    #Create the shot engine and start it in a separate daemon thread
    shot_engine_instance = shot_engine(shot_cycle, expected_shots)
    shot_thread = threading.Thread(target=shot_engine_instance.start)
    shot_thread.daemon = True
    shot_thread.start()

    # Start a thread with the server -- that thread will then start one
    # more thread for each request
    try:
        GUI.main()
    except KeyboardInterrupt:
        print 'killing all subprocesses'

    if killed_all!= 1:
        TCP_server.shutdown()
        kill_everything(computer_details)
        for i in computer_details.keys():
            if i==local_comp:
                kill_servers(i, computer_details[i]['servers'], local=1)
            else:
                kill_servers(i, computer_details[i]['servers'], local=0)
        kill_servers(local_comp, computer_details[local_comp]['programs'], local=1)
        print 'Program finished'
        
    #tidy up, but check to see if it already happened first
    #this one will probably happen because you have run out of shots
    #might be better to reset 'shot clock' and go back to the beginning
    # if killed_all != 1:
    #     print 'killing all subprocesses'
    #     kill_everything(computer_details)
    #     for i in computer_details.keys():
    #         if i==local_comp:
    #             kill_servers(i, computer_details[i]['servers'], local=1)
    #         else:
    #             kill_servers(i, computer_details[i]['servers'], local=0)
    #     TCP_server.shutdown()
    #     print 'Program finished'

