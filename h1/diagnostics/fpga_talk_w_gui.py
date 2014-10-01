import h1.diagnostics.SPE_reader as SPE_reader
import sys
import serial, copy
import time
import numpy as np
import threading
import SocketServer
import os
import gtk
import gobject
import shutil
import MDSplus as MDS

class GUI_control():
    def __init__(self):
        '''
        SRH: 17Aug2013 
        '''
        #setup the GUI with a button and 3 labels
        self.builder = gtk.Builder()
        gladefile = 'PLL_GUI.glade'
        self.builder.add_from_file(gladefile)
        self.window = self.builder.get_object('TopLevel')

        self.min_phase_text_entry = self.builder.get_object('min_phase_entry')
        self.automatic_toggle = self.builder.get_object('automatic_scan')
        self.max_phase_text_entry = self.builder.get_object('max_phase_entry')
        self.increment_text_entry = self.builder.get_object('n_phases_entry')
        self.active_phases_label = self.builder.get_object('actual_phases')
        self.next_automatic_label = self.builder.get_object('next_automatic_label')
        self.increment_automatic_button = self.builder.get_object('increment_automatic_button')
        self.decrement_automatic_button = self.builder.get_object('decrement_automatic_button')
        self.increment_automatic_button.set_sensitive(False)
        self.decrement_automatic_button.set_sensitive(False)
        self.manual_trigger_button = self.builder.get_object('manual_trigger_button')

        self.PLL_setting = self.builder.get_object('actual_PLL')
        self.active_phases_entry = self.builder.get_object('actual_phases_entry')
        self.MDSplus_shot = self.builder.get_object('MDSplus_shot')
        self.MDSplus_phase = self.builder.get_object('MDSplus_phase')
        self.fpga_status_label = self.builder.get_object('status_fpga_label')
        self.MDSplus_textview = self.builder.get_object('MDSplus_textview')
        self.MDSplus_text_buffer = gtk.TextBuffer()
        #self.PLL_phase_entry = self.builder.get_object('PLL_phase_entry')
        combo_box_parent = self.builder.get_object('combo_box_parent')

        enditer = self.MDSplus_text_buffer.get_end_iter()
        self.MDSplus_text_buffer.insert(enditer, 'hello\nworld!')

        self.active_phases_label.set_label("WARNING!! NO PHASES SET")
        self.PLL_setting.set_label("WARNING!! PLL NOT SET")
        self.MDSplus_textview.set_buffer(self.MDSplus_text_buffer)
        self.handlers = {
            "on_min_phase_entry_changed": self.update_phases,
            "on_max_phase_entry_changed": self.update_phases,
            "on_n_phases_entry_changed": self.update_phases,
            "on_update_fpga_clicked": self.fpga_update_clicked,
            "on_reset_fpga_clicked": self.fpga_reset_clicked,
            "on_status_fpga_clicked": self.fpga_status_clicked,
            "on_PLL_update_phase_clicked": self.PLL_phase_update_clicked,
            "on_automatic_scan_toggled": self.automatic_scan_toggled,
            "on_increment_automatic_button_clicked": self.automatic_scan_increment,
            "on_decrement_automatic_button_clicked": self.automatic_scan_decrement,
            "on_manual_trigger_button_clicked": self.manual_trigger_clicked
            }

        self.combobox = gtk.combo_box_new_text()
        phase_items = ['0 User cap',
                       '1 130kHz-180kHz, 90kHz-240kHz*',
                       '2 72kHz-110kHz, 50kHz-160kHz*',
                       '3 44kHz-75kHz, 35kHz-98kHz*',
                       '4 45kHz-82kHz, 33kHz-110kHz*',
                       '5 31kHz-59kHz, 24kHz-76kHz*',
                       '6 23kHz-48kHz, 19kHz-57kHz*',
                       '7 19kHz-39kHz, 16kHz-45kHz*',
                       '8 45kHz-82kHz, 33kHz-110kHz-',
                       '9 31kHz-59kHz, 24kHz-76kHz-',
                       '10 23kHz-48kHz, 19kHz-57kHz-',
                       '11 19kHz-39kHz, 16kHz-45kHz-',
                       '12 17kHz-39kHz, 17kHz-45kHz*',
                       '13 15kHz-33kHz, 13kHz-38kHz*',
                       '14 13kHz-29kHz, 12kHz-32kHz',
                       '15 11kHz-25kHz, 11kHz-25kHz*']
        for i in phase_items:
            self.combobox.append_text(i)

        combo_box_parent.pack_start(self.combobox)
        self.combobox.set_active(13)
        self.builder.connect_signals(self.handlers)
        self.window.connect("delete_event", self.delete_event)
        self.window.connect("destroy", self.destroy)
        self.automatic_scan = self.automatic_toggle.get_active()

        #self.update_phases('')
        self.window.show_all()
        #gtk.main()
        self.update_phases('')
        self.window.show()

    def manual_trigger_clicked(self,tmp):
        self.fpga.fpga_manual_trigger()

    def automatic_scan_increment(self,tmp):
        self.automatic_scan_item += 1
        self.automatic_scan_item = self.automatic_scan_item%self.n_automatic_items
        self.next_automatic_label.set_label('Next:'+self.automatic_items[self.automatic_scan_item])


    def automatic_scan_decrement(self,tmp):
        self.automatic_scan_item -= 1
        self.automatic_scan_item = self.automatic_scan_item%self.n_automatic_items
        self.next_automatic_label.set_label('Next:'+self.automatic_items[self.automatic_scan_item])


    def fpga_update_clicked(self,tmp):
        #self.fpga.phase_range = self.new_phase_range
        if self.automatic_scan:
            cur_item = self.automatic_scan_item
            self.automatic_scan_item += 1
            self.automatic_scan_item = self.automatic_scan_item%self.n_automatic_items
            self.active_phases_entry.set_text(self.automatic_items[cur_item])
            self.next_automatic_label.set_label('Next:' + self.automatic_items[self.automatic_scan_item])
        self.fpga.phase_range = [int(i) for i in self.active_phases_entry.get_text().split(',')]
        #self.fpga.phase_range = self.GUI.new_phase_range
        #self.fpga.initialise()

        #self.fpga.phase_range = [int(i) for i in self.active_phases_entry.get_text().split(',')]
        #print self.fpga.phase_range
        self.fpga.initialise()

    def fpga_reset_clicked(self,tmp):
        self.fpga.fpga_reset()
        #print 'fpga reset clicked'

    def automatic_scan_toggled(self,tmp):
        self.automatic_scan = self.automatic_toggle.get_active()
        print self.automatic_scan, not self.automatic_scan
        self.min_phase_text_entry.set_sensitive(not self.automatic_scan)
        self.max_phase_text_entry.set_sensitive(not self.automatic_scan)
        self.increment_text_entry.set_sensitive(not self.automatic_scan)
        self.active_phases_entry.set_sensitive(not self.automatic_scan)
        self.increment_automatic_button.set_sensitive(self.automatic_scan)
        self.decrement_automatic_button.set_sensitive(self.automatic_scan)
        if self.automatic_scan:
            with file('automatic_scan.txt','r') as input_file:
                self.automatic_items = input_file.readlines()
                self.n_automatic_items = len(self.automatic_items)
                for i in range(self.n_automatic_items):
                    self.automatic_items[i] = self.automatic_items[i].rstrip('\r\n')
                    self.automatic_items[i] = self.automatic_items[i].rstrip('\n')
                print self.automatic_items
            self.automatic_scan_item =0
        self.next_automatic_label.set_label('Next:'+self.automatic_items[self.automatic_scan_item])

    def PLL_phase_update_clicked(self,tmp):
        #print 'PLL_update_phase_clicked'
        PLL_item = self.combobox.get_active_text()
        PLL_phase = PLL_item.split(' ')[0]
        transmit_code = int(PLL_phase)
        print 'Updating PLL to :{}, setting:{}'.format(PLL_item, PLL_phase)
        fpga.PLL_update(transmit_code)
        self.PLL_setting.set_label(PLL_item)

    def fpga_status_clicked(self,tmp):
        #print 'fpga status clicked'
        self.fpga.get_state()
        self.fpga_status_label.set_text('FPGA_status:' + self.fpga.fpga_state_txt)
        #gobject.idle_add(self.update_fpga_status_label,)

    def update_phases(self,tmp):
        #print 'hello'
        min_phase = self.min_phase_text_entry.get_text()
        max_phase = self.max_phase_text_entry.get_text()
        n_phases = self.increment_text_entry.get_text()
        if min_phase=='' or max_phase =='' or n_phases =='':
            self.new_phase_range = []
        else:
            min_phase = int(min_phase)
            max_phase = int(max_phase)
            n_phases = int(n_phases)
            #increment = (max_phase - min_phase) / n_phases
            increment = n_phases
            self.new_phase_range = [min_phase]
            while ((self.new_phase_range[-1]+increment)<max_phase):
                self.new_phase_range.append(self.new_phase_range[-1]+increment)
            #if max_phase>min_phase and increment >0:
            #    print min_phase, max_phase, n_phases, increment
            #    self.new_phase_range = range(min_phase, max_phase, increment)
            #else:
            #    self.new_phase_range = []
        phase_range_txt = ''
        for i in self.new_phase_range:
            phase_range_txt += str(i)+','
        #print phase_range_txt
        phase_range_txt = phase_range_txt.rstrip(',')
        #print phase_range_txt
        #self.active_phases.set_label(phase_range_txt)
        self.active_phases_entry.set_text(phase_range_txt)
        #print min_phase, max_phase, n_phases

    def delete_event(self, widget, event, data=None):
        print "!!! delete event occurred"
        return False

    def destroy(self, widget, data=None):
        print "!!! destroy signal occurred"
        gtk.main_quit()

    def main(self):
        gtk.main()

class shot_object():
    def __init__(self,):
        self.state = 'IDLE'
        self.shot_number = None

        if sys.platform == 'win32' or sys.platform == 'win64':
            self.SPE_dir = 'C:/LightField/shaun_stuff/'
            self.SPE_filename = 'C:/LightField/test.SPE'
        else:
            self.SPE_filename = '/home/srh112/test.SPE'
            self.SPE_dir = '/home/srh112/'

        self.SPE_filename_begins_with = 'test'
        self.rename_prefix = 'shaun_'
        #last time the SPE file was stored
        self.SPE_store_time = 0
        self.SPE_filename = 0
        self.SPE_store_time, self.SPE_filename = self.last_modified_shot()
        self.first_shot = 1
        # try:
        #     tmp = os.stat(self.SPE_filename)
        #     self.SPE_store_time = tmp.st_mtime
        #     print self.SPE_store_time
        # except OSError:
        #     print 'SPE file probably doesnt exist yet, setting last store time to zero'
        #     self.SPE_store_time = 0

    def last_modified_shot(self,):
        dir_contents = os.listdir(self.SPE_dir)
        most_recent_name = ''
        SPE_store_time = 0
        min_filesize = 400000
        for i in dir_contents:
            #check the file is an SPE and has the test part in it
            if (i.find(self.SPE_filename_begins_with)>=0) and ((i.find('SPE')>0) or i.find('spe')>0):
                tmp = os.stat(self.SPE_dir + i)
                store_time = tmp.st_mtime
                if (store_time>SPE_store_time) and (tmp.st_size>min_filesize):
                    SPE_store_time = store_time
                    most_recent_name = self.SPE_dir + i
        print 'latest spe file,',SPE_store_time, most_recent_name
        return SPE_store_time, most_recent_name

    def update_phase(self, phase, shot, tree):
        if phase == 'INIT':
            #Moving from idle to init
            if self.state == 'IDLE':
                self.state = 'INIT'
                self.shot_number = shot
                self.tree_path = tree
                self.init_phase_tasks()
                print '##move from IDLE to INIT'
            #Moving to the next shot
            elif self.state == 'STORE':
                self.state = 'INIT'
                self.tree_path = tree
                self.shot_number = shot
                self.init_phase_tasks()
                print '##move from store to INIT - next shot'
            #reinitialise
            elif self.state == 'INIT' and self.shot_number == shot:
                self.state = 'INIT'
                self.shot_number = shot
                self.tree_path = tree
                print '##Reinitialise'
                self.init_phase_tasks()
            #missed store cycle.... not sure what to do here....
            elif self.state == 'INIT' and self.shot_number != shot:
                self.state = 'INIT'
                self.tree_path = tree
                self.shot_number = shot
                self.init_phase_tasks()
                print '!!!! Something strange has happened, not reinit, or move to store'
            #Give up and go back to idle
            else: 
                self.state = 'IDLE'
        elif phase == 'STORE':
            #correct store
            if self.state == 'INIT' and self.shot_number == shot:
                self.shot_number = shot
                self.tree_path = tree
                self.state = 'STORE'
                self.store_phase_tasks()
                print '##Performing STORE'
            #initialise and store numbers don't match up... what to do here?
            elif self.state == 'INIT' and self.shot_number != shot:
                self.state = 'INIT'
                self.tree_path = tree
                self.shot_number = shot
                print '!!!! INIT and store numbers dont match up.... not doing anything'
            #received two stores in a row?? what to do here?
            elif self.state == 'STORE':
                self.state = 'INIT'
                self.tree_path = tree
                self.shot_number = shot
                print '!!!! multiple stores in a row??? not doing anything'
            #Give up and go back to idle
            else: 
                self.state = 'IDLE'
                print '!!!! Going to IDLE (probably started after INIT)'
        else:
            pass
        gobject.idle_add(self.update_label,)

    def update_label(self,):
        #print 'updating label'
        self.GUI.MDSplus_phase.set_text('MDSplus phase:' + self.state)
        self.GUI.MDSplus_shot.set_text('MDSplus shot:' + str(self.shot_number))

    def init_phase_tasks(self,):
        print 'Performing init ', self.state, self.shot_number
        self.GUI.fpga_update_clicked(None)
        #Leave a semaphore for the autoarm
        success = 0; count = 0
        self.init_filename = self.rename_prefix + str(self.shot_number)
        if sys.platform == 'win32' or sys.platform == 'win64':
            while (success==0) and (count<10):
                try:
                    f = file("C:/LightField/INIT",'w')
                    f.write(self.init_filename)
                    f.close()
                    f = file("C:/LightField/ABORT",'w')
                    f.close()
                    success = 1
                except:
                    print 'exception leaving INIT file'
                    time.sleep(0.1)
                    count += 1


    def store_phase_tasks(self,):
        print 'performing store ', self.state, self.shot_number
        fpga_success = self.fpga.store()
        print(' fpga success : {}'.format(fpga_success))
        old_storage_method = 0
        if old_storage_method:
            new_SPE_store_time, new_SPE_filename = self.last_modified_shot()
            new_file_check = (new_SPE_store_time>self.SPE_store_time)
            new_filename_check = (new_SPE_filename!=self.SPE_filename)
        else:
            new_SPE_filename = self.SPE_dir + self.init_filename + '.spe'
            new_filename_check = 1
            new_file_check = 1
        try:
            camera_data, xml_string = SPE_reader.extract_SPE_data(new_SPE_filename)
            able_to_read_data = 1
        except Exception, e:
            able_to_read_data = 0
            print '!!!!!! cant read data in SPE FILE'
            print e

        if new_file_check and new_filename_check and able_to_read_data:
            if old_storage_method:
                new_filename = self.SPE_dir + self.rename_prefix + str(self.shot_number) + '.spe'
                print ' copying {} to new :{}'.format(new_SPE_filename, new_filename)
                shutil.copy(new_SPE_filename, new_filename)
                self.SPE_store_time = new_SPE_store_time
                self.SPE_filename = new_SPE_filename 
            self.store_MDSplus_data = 1
            if self.store_MDSplus_data:
                #camera_data, xml_string = SPE_reader.extract_SPE_data(new_SPE_filename)
                os.environ['h1data_path']='h1data.anu.edu.au::'
                Tree = MDS.Tree('h1data', self.shot_number)

                Tree.getNode(self.tree_path+'.PIMAX:IMAGES').putData(camera_data)
                Tree.getNode(self.tree_path+'.PIMAX:SETTINGS').putData(xml_string)
                #tmp = [str(i) for i in self.fpga.current_stored_phases]
                
                #tmp = self.GUI.active_phases_entry.get_text())
                Tree.getNode(self.tree_path+'.scan:phases').putData(self.GUI.active_phases_label.get_label())
                Tree.getNode(self.tree_path+'.PLL:lockrange').putData(self.GUI.PLL_setting.get_label())
        else:
            print '!!!!! SPE file has not been updated'
            print '!!!!! Leaving ABORT file for LightField'
            #send abort to LightField
            if sys.platform == 'win32' or sys.platform == 'win64':
                success = 0; count = 0
                while (success==0) and (count<10):
                    try:
                        f = file("C:/LightField/ABORT",'w')
                        f.close()
                        success = 1
                        print ' wrote ABORT file'
                    except:
                        print '!!!!! exception leaving ABORT file'
                        time.sleep(0.1)
                        count += 1

class ThreadedTCPRequestHandler(SocketServer.BaseRequestHandler):
    def handle(self):
        data = self.request.recv(1024)
        #print('hello - connection received')
        print 'TCP connection, data received : {}'.format(data)
        termination_character = '\r\n'
        sep_char = chr(0x09)
        data = data.rstrip(termination_character)
        phase, tree, shot = data.split(sep_char)
        shot = int(shot)
        self.server.shot_object.update_phase(phase, shot, tree)
        cur_thread = threading.current_thread()
        response = '0'+termination_character
        #response = "{}: {} : {} :{}".format(cur_thread.name, data, self.server.shot_cycle.cycle, self.server.shot_cycle.executing_shot)
        self.request.sendall(response)
        #if data.find('ABORT')>=0: self.server.shot_cycle.abort_shot()

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
    print "Status server is running in thread:", server_thread.name
    return server_thread, server

####################################################
class fpga_PLL():
    def __init__(self,):
        if sys.platform == 'win32' or sys.platform == 'win64':
            self.serial_port = 'COM1'
        else:
            self.serial_port = '/dev/ttyUSB0'
        self.baud_rate = 9600
        self.parity = serial.PARITY_NONE
        self.stopbits = 2
        self.get_status = 63
        self.disarm_command = 64
        self.arm_command = 65
        self.write_to_memory_code = 66
        self.read_from_memory_code = 67
        self.finish_data_to_memory_code = 68
        self.PLL_cap_update_code = 69
        self.manual_trigger_code = 70
        self.reset_command = 255
        self.phase_range = []
        self.old_way = 0
    def connect_serial(self,):
        try:
            ser=serial.Serial(self.serial_port, self.baud_rate, timeout=1, parity = self.parity, stopbits = self.stopbits)
        except serial.SerialException:
            raise RuntimeError('!!! EXCEPTION - Unable to open serial port - check no other program is using it')
        #print 'Serial Port Open::'
        return ser

    def get_state(self,):
        ser = self.connect_serial()
        ser.flushInput()
        self.fpga_state = self.sendMicroCommand(ser, self.get_status, echo = 0, expect_return = 1, give_return_value = 1)
        tmp = ord(self.fpga_state)
        if tmp == 0:
            self.fpga_state_txt = 'IDLE'
        elif tmp ==1:
            self.fpga_state_txt = 'CHANGE_PHASE'
        elif tmp ==2:
            self.fpga_state_txt = 'RUNNING'
        print 'fpga state:{}'.format(ord(self.fpga_state))
        #gobject.idle_add(self.GUI.update_fpga_status_label,)
        ser.close()

    def fpga_reset(self,):
        ser = self.connect_serial()
        ser.flushInput()
        print 'Disarm FPGA'
        self.sendMicroCommand(ser, self.disarm_command)
        ser.close()

    def fpga_manual_trigger(self,):
        ser = self.connect_serial()
        ser.flushInput()
        print 'Manual trigger FPGA'
        self.sendMicroCommand(ser, self.manual_trigger_code)
        ser.close()

    def PLL_update(self, phase_setting):
        if self.old_way:
            phase_setting = phase_setting + 32
        ser = self.connect_serial()
        ser.flushInput()
        print 'Update PLL lock range'
        if self.old_way:
            self.sendMicroCommand(ser, phase_setting)
        else:
            self.sendMicroCommand(ser, self.PLL_cap_update_code)
            self.sendMicroCommand(ser, phase_setting)
        ser.close()

    def initialise(self,  arm = 1):
        ser = self.connect_serial()
        ser.flushInput()
        print 'Initialise FPGA'
        #print ' disarm fpga'
        self.sendMicroCommand(ser, self.disarm_command)
        #print ' reset fpga active phases'
        self.sendMicroCommand(ser, self.reset_command)
        if not self.old_way:
            self.sendMicroCommand(ser, self.write_to_memory_code)
        print ' program phases:', self.phase_range
        for i in self.phase_range:
            self.sendMicroCommand(ser, i)
        #print ' arm fpga'
        if not self.old_way:
            self.sendMicroCommand(ser, self.finish_data_to_memory_code)
        self.sendMicroCommand(ser, self.arm_command)
        #print(' close serial port')
        ser.close()
        self.current_stored_phases = copy.deepcopy(self.phase_range)
        tmp = ''
        for i in self.current_stored_phases: tmp+='{},'.format(i)
        self.GUI.active_phases_label.set_label(tmp.rstrip(','))

    def store(self, arm = 1):
        print 'Store FPGA'
        ser = self.connect_serial()
        ser.flushInput()
        self.fpga_state = self.sendMicroCommand(ser, self.get_status, echo = 0, expect_return = 1, give_return_value = 1)
        tmp = ord(self.fpga_state)
        return_value = 0
        if tmp == 0:
            self.fpga_state_txt = 'IDLE'
            return_value = 1
        elif tmp ==1:
            self.fpga_state_txt = 'CHANGE_PHASE'
        elif tmp ==2:
            self.fpga_state_txt = 'RUNNING'
        elif tmp ==3:
            self.fpga_state_txt = 'WAIT'
        print ' fpga state:{}'.format(ord(self.fpga_state))
        #gobject.idle_add(self.GUI.update_fpga_status_label,)
        ser.close()
        return return_value

    def waitForData(self,ser, numberOfBytes):
        delayCount=0
        while ((ser.inWaiting()<numberOfBytes)&(delayCount<50)):
            delayCount=delayCount+1
            time.sleep(0.05)
        if delayCount == 50:
            ser.close()
            raise RuntimeError('!!! EXCEPTION - Waited too long for data')

    def sendMicroCommand(self,ser, command, echo = 1, expect_return = 1, give_return_value = 0):
        command_ascii = chr(int(np.binary_repr(command),2))
        ser.write(command_ascii)
        #print('sending : {}, ascii : {}'.format(command, command_ascii))
        if expect_return:
            self.waitForData(ser, 1)
            tmp = ser.read()
            #print tmp
            if echo == 1:
                if tmp != command_ascii:
                    ser.close()
                    raise RuntimeError('!!! EXCEPTION - Incorrect Echoed Value')
        if give_return_value ==1:
            return tmp

if __name__ == '__main__':
    gobject.threads_init()

    fpga = fpga_PLL()
    shot_object_cycle = shot_object()
    shot_object_cycle.fpga = fpga

    TCP_PORT = 8111
    TCP_server_thread, TCP_server = start_TCP_server(TCP_PORT)
    TCP_server.shot_object = shot_object_cycle

    GUI = GUI_control()
    GUI.fpga = fpga
    fpga.GUI = GUI
    shot_object_cycle.GUI = GUI
    GUI.main()

