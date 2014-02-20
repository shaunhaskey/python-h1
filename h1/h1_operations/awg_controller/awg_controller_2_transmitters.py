#!/usr/bin/env python
'''
Notes:
-It is important for the output to go into 50Ohms, otherwise you don't get what you think you will...
-The voltage range should go from -32767-> 32767, then use amplitude and offset to set the meaning of the max and min
- currently a problem with capping at the top!!
'''

from matplotlib.figure import Figure
from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas
from matplotlib.backends.backend_gtkagg import NavigationToolbar2GTKAgg as NavigationToolbar
import gtk
import numpy as np
import socket, time
import datetime

import matplotlib.pyplot as pt
class awg():
    def __init__(self,host, PORT):
        self.host = host
        self.PORT = PORT
        self.connect_to_awg()

    def connect_to_awg(self):
        '''
        Open a socket to the AWG
        SH : 11Mar2013
        '''
        self.s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        print 'connecting to host: %s:%d'%(self.host,self.PORT)
        self.s.connect((self.host, self.PORT))

    def reset_awg(self):
        '''
        Clear the settings on the AWG
        SH : 11Mar2013
        '''
        print ' Resetting AWG'
        self.s.send("*RST\n")
        print ' Clearing AWG'
        self.s.send("*CLS\n")

    def download_data_to_awg(self, data_to_send, filename, chan=1):
        '''
        Download waveform to the AWG, and load it into memory
        SH : 11Mar2013
        '''
        bytes_to_send = len(data_to_send)
        len_bytes_to_send = len(str(bytes_to_send))
        print ' Sending waveform to AWG'
        new_string = 'MMEM:DOWN:DATA #%d%d%s\n'%(len_bytes_to_send, bytes_to_send, data_to_send)
        print ' Save waveform to AWG'
        self.s.send('''MMEM:DOWN:FNAM "%s"\n'''%(filename+'.arb'))
        self.s.send(new_string)
        print ' Load waveform into memory'
        self.s.send('MMEM:LOAD:DATA%d "%s"\n'%(chan,filename+'.arb'))
        #self.send_command('MMEM:LOAD:DATA{}'.format(chan),'"{}"'.format(filename+'.arb'))

    def setup_arb_mode(self):
        '''
        Setup up for ARB waveform on the AWG and load file we sent
        SH : 11Mar2013
        '''
        #self.s.send('FUNCtion ARB\n')
        success = 1
        for i in [1,2]:
            ret_val = self.send_command('SOUR{}:FUNCtion'.format(i), 'ARB')
            success*=ret_val
        ret_val = self.send_command('SOUR{}:FUNC:ARB'.format(1), '"'+self.remote_file_location+'.arb'+'"')
        success*=ret_val
        ret_val = self.send_command('SOUR{}:FUNC:ARB'.format(2), '"'+self.remote_file_location2+'.arb'+'"')
        success*=ret_val
        return success

    def send_command(self, stump, arg):
        command = '{} {}'.format(stump, arg)
        self.s.send('{} {}\n'.format(stump, arg))
        self.s.send('{}?\n'.format(stump,))
        reply = self.s.recv(1024).rstrip('\n')
        try:
            success = np.allclose(float(reply),float(arg))
        except ValueError as e:
            success = reply==arg
        print(' command:{}, sucessfully set:{}'.format(command, success))
        if not success:print " !!!!!!!!!"
        return success

    def setup_trigger(self):
        '''
        Setup the trigger settings for the AWG
        SH : 11Mar2013
        '''
        success = True
        for i in [1,2]:
            ret_val = self.send_command('SOUR{}:BURS:MODE'.format(i), 'TRIG')
            success*=ret_val
            ret_val = self.send_command('SOUR{}:BURS:NCYC'.format(i), '1')
            success*=ret_val
            ret_val = self.send_command('TRIG{}:SOUR'.format(i), 'EXT')
            success*=ret_val
            ret_val = self.send_command('SOUR{}:BURS:STAT'.format(i), '1')
            success*=ret_val
        return success

    def set_status_on(self):
        '''
        Start the AWG waiting for a trigger
        SH : 11Mar2013
        '''
        success = True
        for i in [1,2]:
            ret_value = self.send_command('OUTPut{}'.format(i), '1')
            success*=ret_value
        return success

    def close_socket(self):
        '''
        close the socket
        SH : 11Mar2013
        '''
        print 'Closing connection to AWG'
        self.s.close()


class DialogTest:
    def __init__(self):
        self.dia = gtk.Dialog('TEST DIALOG', self.window, 
           gtk.DIALOG_MODAL  | gtk.DIALOG_DESTROY_WITH_PARENT)
        self.dia.vbox.pack_start(gtk.Label('This is just a Test'))
        
        #self.button.connect("clicked", self.rundialog, None)
        self.rundialog()
        #self.window.add(self.button)
        #self.button.show()
        #self.window.show()
        self.dia.show_all()
        result = self.dia.run() 
        self.dia.hide()

class nice_gui():
    def __init__(self, host, PORT):
        try:
            self.AWG = awg(host, PORT)
        except:
            print 'failed to connect to AWG'
        self.builder = gtk.Builder()
        #gladefile = 'SCPI_GUI_V2.glade'
        gladefile = 'awg_gui_2_transmitters.glade'
        self.builder.add_from_file(gladefile)

        #extract the useful objects from the gui
        self.win = self.builder.get_object('Charter')
        self.vbox1 = self.builder.get_object('vbox2')
        self.end_time_inp= self.builder.get_object('end_time_input')
        self.amplitude_inp = self.builder.get_object('amplitude_input')
        self.sample_rate_inp = self.builder.get_object('sample_rate_input')
        #self.offset_inp = self.builder.get_object('offset_input')
        #self.offset_inp.set_text('0.0')
        self.waveform_inp = self.builder.get_object('waveform_input')
        self.waveform_inp_T2 = self.builder.get_object('waveform_input1')
        self.execute_button = self.builder.get_object('execute_button')
        self.execute_button.set_sensitive(False)
        self.check_button_50ohm = self.builder.get_object('check_button_50ohm')
        self.same_check_toggle = self.builder.get_object('same_check')
        self.status_label = self.builder.get_object('status_label')
        self.status_label.set_use_markup(True)
        self.status_label.set_alignment(0,0.5)

        #Create a matplotlib figure to show everything
        self.fig = Figure(figsize=(5,4), dpi=100)
        self.ax = self.fig.add_subplot(211)
        self.ax2 = self.fig.add_subplot(212, sharex=self.ax, sharey = self.ax)
        self.ax.grid()
        self.ax2.grid()
        self.canvas = FigureCanvas(self.fig)  # a gtk.DrawingArea
        self.vbox1.pack_start(self.canvas)
        self.toolbar = NavigationToolbar(self.canvas, self.win)
        self.vbox1.pack_start(self.toolbar, False, False)

        #Connect handles to functions here
        self.handlers = {
            "gtk_main_quit": gtk.main_quit,
            "on_validate_button_clicked": self.validate,
            "on_execute_button_clicked": self.execute,
            "on_retrieve_button_clicked": self.retreive,
            "on_store_button_clicked": self.store,
            "on_end_time_input_changed": self.text_changed,
            "on_amplitude_input_changed": self.text_changed,
            "on_sample_rate_input_changed": self.text_changed,
            "on_waveform_input_changed": self.text_changed,
            "on_check_button_50ohm_toggled": self.text_changed,
            "on_check_button_50ohm_clicked": self.text_changed,
            "on_same_check_toggled": self.same_check_toggled,
            "on_default_settings_button_clicked": self.default_settings,
            }
        self.builder.connect_signals(self.handlers)
        self.duplicate_transmitters = 0
        self.default_settings(None)
        self.same_check_toggled(None)
        
        #show the gui and give control to gtk
        self.win.show_all()
        gtk.main()

    def default_settings(self,dum):
        self.waveform_inp.set_text('2*step(0.001,0.1)')
        self.waveform_inp_T2.set_text('2*step(0.001,0.1)')
        self.end_time_inp.set_text('0.15')
        self.amplitude_inp.set_text('4.5')
        self.sample_rate_inp.set_text('20000')
        self.same_check_toggle.set_active(True)
        self.check_button_50ohm.set_active(False)
        self.same_check_toggled(None)


    def same_check_toggled(self,toggle):
        self.duplicate_transmitters = self.same_check_toggle.get_active()
        print self.duplicate_transmitters
        if self.duplicate_transmitters:
            self.waveform_inp_T2.set_editable(False)
            self.waveform_inp_T2.set_text(self.waveform_inp.get_text())
        else:
            self.waveform_inp_T2.set_editable(True)

    def text_changed(self, widget):
        '''
        Some of the settings have changed - disable the execute button
        SH : 11Mar2013
        '''
        self.execute_button.set_sensitive(False)
        self.status_label.set_markup('<span color="red">STATUS : WAVEFORM CHANGED, NEED TO VALIDATE AND PROGRAM AWG</span>')
        if self.duplicate_transmitters:
            self.waveform_inp_T2.set_text(self.waveform_inp.get_text())
            

    def validate(self, button):
        '''
        Validate the settings, and create the arb file to send to the AWG
        Checks to implement : minimum voltage, maximum voltage, maximum samples
        maximum sample rate, output impedance
        '''
        print('###################################################')
        print('######## {} ############'.format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
        print('######## Validating waveforms ########')
        self.status_label.set_markup('<span color="red">STATUS : VALIDATING WAVEFORM</span>')
        self.execute_button.set_sensitive(False)

        self.end_time = float(self.end_time_inp.get_text())
        self.amplitude = float(self.amplitude_inp.get_text())
        self.sample_rate = int(self.sample_rate_inp.get_text())
        #self.offset = float(self.offset_inp.get_text())

        self.waveform = self.waveform_inp.get_text()
        self.waveform2 = self.waveform_inp_T2.get_text()
        
        print ' waveform T1:', self.waveform
        print ' waveform T2:', self.waveform2

        def step(start, stop):
            sig_time = self.t
            start_point = np.argmin(np.abs(sig_time - start))
            end_point = np.argmin(np.abs(sig_time - stop))
            output = np.zeros(len(sig_time))
            output[start_point:end_point] = 1
            return output
        
        t = np.arange(0,self.end_time,1./self.sample_rate)
        self.t = t
        def check_waveform(waveform, t, step):
            if self.waveform[0]=='[':
                data = t*0
                tmp = waveform.lstrip('[').rstrip(']')
                tmp_list = tmp.split(';')
                for i in range(len(tmp_list)-1):
                    tmp1 = tmp_list[i].lstrip('(').rstrip(')')
                    tmp1 = tmp1.split(',')
                    tmp2 = tmp_list[i+1].lstrip('(').rstrip(')')
                    tmp2 = tmp2.split(',')
                    start_pt = np.argmin(np.abs(t-float(tmp1[0])))
                    end_pt = np.argmin(np.abs(t-float(tmp2[0])))
                    tmp_poly = np.polyfit([float(tmp1[0]),float(tmp2[0])],[float(tmp1[1]),float(tmp2[1])], 1)
                    data[start_pt:end_pt] = tmp_poly[0]*t[start_pt:end_pt] +tmp_poly[1]
                pass
            else:
                data = eval(waveform)
            error = 0
            if np.min(data) <0:
                print '!!! Error : data below minimim'
                error = 1
            if data[0] != 0:
                print '!!! Error : first value isnt zero'
                error = 1
            if data[-1] != 0:
                print '!!! Error : last value isnt zero'
                error = 1
            if np.max(data) >4.5 :
                print '!!! Error : asking for a drive greater than 4.5'
                error = 1
            return error, data
        error1 = 0; error2 = 0
        try:
            error1, self.data = check_waveform(self.waveform, self.t, step)
        except Exception as e:
            print(' !!!! Error validating channel 1') 
            error1 = 1
            print e
        try:
            error2, self.data2 = check_waveform(self.waveform2, self.t, step)
        except Exception as e:
            print(' !!!! Error validating channel 1') 
            error2 = 1
            print e
        error = (error1 + error2)>=1
        if not error:
            self.amplitude = np.max(self.data)
            self.offset = np.min(self.data)
            self.amplitude2 = np.max(self.data2)
            self.offset2 = np.min(self.data2)
            print ' end time : %.2f, amp : %.2f, sample : %d, offset : %.2f'%(self.end_time, self.amplitude, self.sample_rate, self.offset)
            print ' end time : %.2f, amp : %.2f, sample : %d, offset : %.2f'%(self.end_time, self.amplitude2, self.sample_rate, self.offset2)

            if not self.check_button_50ohm.get_active():
                self.amplitude = np.max(self.data)*0.5
                self.offset = np.min(self.data)*0.5
                self.amplitude2 = np.max(self.data2)*0.5
                self.offset2 = np.min(self.data2)*0.5
            
            self.output_string = self.create_arb_file(1, self.sample_rate, self.amplitude, self.offset, 60 , 'off', self.data)
            self.output_string2 = self.create_arb_file(1, self.sample_rate, self.amplitude2, self.offset2, 60 , 'off', self.data2)
            if not error:
                self.execute_button.set_sensitive(True)
                self.status_label.set_markup('<span color="red">STATUS : SUCCESSFUL VALIDATION, NEED TO SEND TO AWG!!</span>')
            else:
                self.execute_button.set_sensitive(False)
                print '!!! Unable to allow execution due to validation problems'
                self.status_label.set_markup('<span color="red">STATUS : FAILED TO VALIDATE!!</span>')
        else:
            self.status_label.set_markup('<span color="red">STATUS : FAILED TO VALIDATE!!</span>')
        
        #Update the figure with the new waveform
        self.ax.cla()
        self.ax2.cla()
        if not error:
            self.ax.plot(t, self.data, '.-')
            self.ax2.plot(t, self.data2, '.-')
        self.ax.grid()
        self.ax2.grid()
        self.ax2.set_xlabel('Time (s)')
        self.ax.set_ylabel('Output (V)')
        self.ax2.set_ylabel('Output (V)')
        self.canvas.draw()

    def create_arb_file(self, chan_count, samp_rate, amp, offset, mark_point, filt, data):
        '''
        Create the arb file based on the settings given
        SH : 11Mar2013
        '''
        output = 'File Format:1.10\n'
        output += 'Checksum:0\n'
        output += 'Channel Count:%d\n'%(chan_count)
        output += 'Sample Rate:%.4f\n'%(samp_rate)
        output += 'High Level:%.5f\n'%(amp)
        output += 'Low Level:%.5f\n'%(offset)
        output += 'Marker Point:%d\n'%(mark_point)
        output += 'Data Type:"short"\n'
        output += 'Filter:"%s"\n'%(filt)
        output += 'Data Points:%d\n'%(len(data))
        output += 'Data:\n'
        max_value = np.max(data)
        min_value = np.min(data)
        offset = (min_value+max_value)/2
        multiplier = 32767/(max_value-offset)
        print multiplier, offset
        raw_data = []
        cropped_data = []
        for i in data:
            tmp = int(np.round((i-offset)*multiplier))
            raw_data.append(tmp)
            #tmp = int(np.round(32767*i))
            if tmp>32767:
                tmp=32767
                print 'too large'
            if tmp<-32767:
                tmp = -32767
                print 'too low'
            cropped_data.append(tmp)
            output+=str(tmp)+'\n'
        show_dac_values = 0
        if show_dac_values:
            fig, ax = pt.subplots()
            ax.plot(raw_data,'.')
            ax.plot(cropped_data,'x')
            fig.canvas.draw(); fig.show()
        
        tmp = file('tmp_output.txt','w')
        tmp.write(output)
        tmp.close()
        return output


    def execute(self,button):
        '''
        Download the data to the AWG, and set everything up so its ready to trigger
        SH : 11Mar2013
        '''
        success = True
        print('###################################################')
        print('######## {} ############'.format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
        print('########Sending waveforms and setting up AWG########')
        self.AWG.reset_awg()
        self.AWG.remote_file_location = 'INT:\MYFILE4'
        self.AWG.remote_file_location2 = 'INT:\MYFILE5'

        self.AWG.download_data_to_awg(self.output_string, self.AWG.remote_file_location)
        self.AWG.download_data_to_awg(self.output_string2, self.AWG.remote_file_location2, chan=2)

        ret_value = self.AWG.setup_arb_mode()
        success *= ret_value
        ret_value = self.AWG.setup_trigger()
        success *= ret_value
        ret_value = self.AWG.set_status_on()
        success *= ret_value
        print('######## {} ############'.format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
        if success:
            self.status_label.set_markup('<span color="black">STATUS : WAVEFORM SUCCESSFULLY LOADED</span>')
            print '######## SUCCESSFULLY SETUP AWG  ########'
        else:
            self.status_label.set_markup('<span color="red">STATUS : ERROR DOWNLOADING WAVEFORM</span>')
            print '########!!!!! ERROR SETTING UP AWG ########'
        print '#######################################'

    def retreive(self,button):
        '''
        Going to implement something with MDSplus here
        SH : 11Mar2013
        '''
        
        print 'Retreive!!'

    def store(self,button):
        '''
        Going to implement something with MDSplus here for the future
        SH : 11Mar2013
        '''
        print 'Store!!'
        self.error('BLAHHHHHH')

    def error(self, error_msg):
        self.dia = gtk.Dialog('TEST DIALOG', self.win, 
           gtk.DIALOG_MODAL  | gtk.DIALOG_DESTROY_WITH_PARENT)
        self.dia.vbox.pack_start(gtk.Label(error_msg))
        
        self.dia.show_all()
        result = self.dia.run() 
        self.dia.hide()

if __name__ == '__main__':
    PORT = 5025
    host = "192.168.1.153"
    a = nice_gui(host, PORT)
