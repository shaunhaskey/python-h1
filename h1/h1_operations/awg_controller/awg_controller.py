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
        print 'Resetting AWG'
        tmp = self.s.send("*RST\n")
        print tmp
        print 'Clearing AWG'
        tmp = self.s.send("*CLS\n")

    def download_data_to_awg(self, data_to_send):
        '''
        Download waveform to the AWG, and load it into memory
        SH : 11Mar2013
        '''
        self.remote_file_location = 'INT:\MYFILE4'
        bytes_to_send = len(data_to_send)
        len_bytes_to_send = len(str(bytes_to_send))
        print 'Sending waveform to AWG'
        new_string = 'MMEM:DOWN:DATA #%d%d%s\n'%(len_bytes_to_send, bytes_to_send, data_to_send)
        print 'Save waveform to AWG'
        self.s.send('''MMEM:DOWN:FNAM "%s"\n'''%(self.remote_file_location+'.arb'))
        self.s.send(new_string)
        print 'Load waveform into memory'
        self.s.send('MMEM:LOAD:DATA "%s"\n'%(self.remote_file_location+'.arb'))

    def setup_arb_mode(self):
        '''
        Setup up for ARB waveform on the AWG and load file we sent
        SH : 11Mar2013
        '''
        self.s.send('FUNCtion ARB\n')
        self.s.send('FUNC:ARB "%s"\n'%(self.remote_file_location+'.arb'))

    def setup_trigger(self):
        '''
        Setup the trigger settings for the AWG
        SH : 11Mar2013
        '''
        self.s.send('BURS:MODE TRIG\n')
        self.s.send('BURS:NCYC 1\n')
        self.s.send('TRIG:SOUR EXT\n')
        self.s.send('BURS:STAT ON\n')


    def set_status_on(self):
        '''
        Start the AWG waiting for a trigger
        SH : 11Mar2013
        '''
        tmp = self.s.send('OUTPut 1\n')


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
        gladefile = 'SCPI_GUI_V2.glade'
        self.builder.add_from_file(gladefile)

        #extract the useful objects from the gui
        self.win = self.builder.get_object('Charter')
        self.vbox1 = self.builder.get_object('vbox2')
        self.end_time_inp= self.builder.get_object('end_time_input')
        self.end_time_inp.set_text('0.15')
        self.amplitude_inp = self.builder.get_object('amplitude_input')
        self.amplitude_inp.set_text('4.5')
        self.sample_rate_inp = self.builder.get_object('sample_rate_input')
        self.sample_rate_inp.set_text('20000')
        #self.offset_inp = self.builder.get_object('offset_input')
        #self.offset_inp.set_text('0.0')
        self.waveform_inp = self.builder.get_object('waveform_input')
        self.waveform_inp.set_text('np.sin(100*t*np.pi*2.)*step(0.01,0.1)')
        self.execute_button = self.builder.get_object('execute_button')
        self.execute_button.set_sensitive(False)
        self.check_button_50ohm = self.builder.get_object('check_button_50ohm')


        #Create a matplotlib figure to show everything
        self.fig = Figure(figsize=(5,4), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.ax.grid()
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
            }
        self.builder.connect_signals(self.handlers)
        
        #show the gui and give control to gtk
        self.win.show_all()
        gtk.main()
        
    def text_changed(self, widget):
        '''
        Some of the settings have changed - disable the execute button
        SH : 11Mar2013
        '''
        self.execute_button.set_sensitive(False)
        

    def validate(self, button):
        '''
        Validate the settings, and create the arb file to send to the AWG
        Checks to implement : minimum voltage, maximum voltage, maximum samples
        maximum sample rate, output impedance
        '''
        self.end_time = float(self.end_time_inp.get_text())
        self.amplitude = float(self.amplitude_inp.get_text())
        self.sample_rate = int(self.sample_rate_inp.get_text())
        #self.offset = float(self.offset_inp.get_text())

        self.waveform = self.waveform_inp.get_text()
        
        print 'waveform :', self.waveform

        def step(start, stop):
            sig_time = self.t
            start_point = np.argmin(np.abs(sig_time - start))
            end_point = np.argmin(np.abs(sig_time - stop))
            output = np.zeros(len(sig_time))
            output[start_point:end_point] = 1
            return output
        
        t = np.arange(0,self.end_time,1./self.sample_rate)
        self.t = t
        if self.waveform[0]=='[':
            self.data = self.t*0
            tmp = self.waveform.lstrip('[').rstrip(']')
            tmp_list = tmp.split(';')
            for i in range(len(tmp_list)-1):
                tmp1 = tmp_list[i].lstrip('(').rstrip(')')
                tmp1 = tmp1.split(',')
                tmp2 = tmp_list[i+1].lstrip('(').rstrip(')')
                tmp2 = tmp2.split(',')
                start_pt = np.argmin(np.abs(self.t-float(tmp1[0])))
                end_pt = np.argmin(np.abs(self.t-float(tmp2[0])))
                tmp_poly = np.polyfit([float(tmp1[0]),float(tmp2[0])],[float(tmp1[1]),float(tmp2[1])], 1)
                self.data[start_pt:end_pt] = tmp_poly[0]*self.t[start_pt:end_pt] +tmp_poly[1]
            pass
        else:
            self.data = eval(self.waveform)
        error = 0
        if np.min(self.data) <0:
            print '!!! Error : data below minimim'
            error = 1
        if self.data[0] != 0:
            print '!!! Error : first value isnt zero'
            error = 1
        if self.data[-1] != 0:
            print '!!! Error : last value isnt zero'
            error = 1
        if np.max(self.data) >4.5 :
            print '!!! Error : asking for a drive greater than 4.5'
            error = 1

        self.amplitude = np.max(self.data)
        self.offset = np.min(self.data)
        print 'Validate!! end time : %.2f, amp : %.2f, sample : %d, offset : %.2f'%(self.end_time, self.amplitude, self.sample_rate, self.offset)

        if not self.check_button_50ohm.get_active():
            self.amplitude = np.max(self.data)*0.5
            self.offset = np.min(self.data)*0.5
            
        self.output_string = self.create_arb_file(1, self.sample_rate, self.amplitude, self.offset, 60 , 'off', self.data)
        if error==0:
            self.execute_button.set_sensitive(True)
        else:
            print '!!! Unable to allow execution due to validation problems'
        #Update the figure with the new waveform
        self.ax.cla()
        self.ax.plot(t, self.data, '.-')
        self.ax.grid()
        self.ax.set_xlabel('Time (s)')
        self.ax.set_ylabel('Output (V)')
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
        self.AWG.reset_awg()
        self.AWG.download_data_to_awg(self.output_string)
        self.AWG.setup_arb_mode()
        self.AWG.setup_trigger()
        self.AWG.set_status_on()

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
    host = "192.168.1.150"
    a = nice_gui(host, PORT)
