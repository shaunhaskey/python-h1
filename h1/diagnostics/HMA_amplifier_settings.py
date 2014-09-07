#!/usr/bin/env python
'''
This widget is for changing the setting for the HMA amplifiers. The MDSplus write has not been implemented yet - still need to check it works properly first.

#f1 : 1.45MHz 2 pole LP
#f2 : 1.66MHz 2 pole LP
#f3 : 927Hz 2 pole HP
#f4 : 380kHz 2 pole LP
#AA : 1MHz 4 pole LP (1.45MHz 2 pole and 1.66MHz 2 pole)
#gain : 1=125x, 2=350x, 3=625x, 4=1500x


SRH : 4Aug2014
'''


import pygtk
pygtk.require('2.0')
import gtk
import MDSplus
class MyProgram:
    def __init__(self):
        # create a new window


        labels = []
        self.filt1 = []; self.filt2 = []; self.filt3 = []; self.filt4 = []
        self.aa = []; self.gain = []


        app_window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        app_window.set_size_request(1100, 800)
        app_window.set_border_width(10)
        app_window.set_title("HMA amplifier settings")

        #Add a parent vbox
        vbox_app = gtk.VBox(False, 0)
        app_window.add(vbox_app)
        vbox_app.show()

        #Add a hbox for the MDSplus buttons
        hbox_app = gtk.HBox(False, 0)
        hbox_app.show()
        vbox_app.pack_start(hbox_app, expand = False)

        #Add the MDSplus buttons
        button_get_mds = gtk.Button(label = 'get from MDSplus')
        button_get_mds.connect("clicked", self.get_mds)
        hbox_app.pack_start(button_get_mds,)
        button_get_mds.show()
        button_save_mds = gtk.Button(label = 'save to MDSplus')
        button_save_mds.connect("clicked", self.save_mds)
        hbox_app.pack_start(button_save_mds,)
        button_save_mds.show()

        #Add a table for the header section
        table_layout_hdr = gtk.Table(rows=2, columns=7, homogeneous=True)
        table_layout_hdr.show()

        count = 0
        overall_filts = []
        lab_list = ['', '1.45MHz LP 2pole', '1.66MHzLP 2pole', '927Hz HP 2pole', '380kHz LP 2pole', '1MHz 4 pole LP']
        for i, (filt_list, lab) in enumerate(zip([self.filt1, self.filt1, self.filt2, self.filt3, self.filt4, self.aa],lab_list)):
            #overall_filts.append(gtk.CheckButton("Filt{}".format(i+1)))
            tmp = (gtk.Label(lab))
            tmp.show()
            table_layout_hdr.attach(tmp, i, i+1, count, count + 1, 0,0,0,0)
        count += 1
        labels.append(gtk.Label("All"))
        labels[-1].show()
        table_layout_hdr.attach(labels[-1], 0, 1, count, count+1, 0,0,0,0)
        for i, (filt_list, lab) in enumerate(zip([self.filt1, self.filt2, self.filt3, self.filt4, self.aa],lab_list[1:])):
            #overall_filts.append(gtk.CheckButton("Filt{}".format(i+1)))
            overall_filts.append(gtk.CheckButton(lab))
            overall_filts[-1].connect("toggled", self.all_checkbox_changed, filt_list)
            overall_filts[-1].set_active(True)  # Set the defaut
            overall_filts[-1].show()
            table_layout_hdr.attach(overall_filts[-1], i+1, i+2, count, count + 1, 0,0,0,0)
        overall_filts.append(gtk.combo_box_new_text())

        #Overall gain settings
        gains = ['125x','350x', '625x', '1500x']
        for gain_lab in gains:overall_filts[-1].append_text(gain_lab)
        overall_filts[-1].set_active(0)
        overall_filts[-1].show()
        overall_filts[-1].connect("changed", self.all_gain_changed, self.gain)

        vbox_app.pack_start(table_layout_hdr,expand = False)

        #The main table
        table_layout = gtk.Table(rows=48+2, columns=7, homogeneous=True)
        table_layout.attach(overall_filts[-1], i+2, i+3, count, count+1, 0,0,0,0)

        #Everything goes in a scroll window
        scrolled_window = gtk.ScrolledWindow(hadjustment=None, vadjustment=None)
        scrolled_window.show()
        vbox_app.pack_start(scrolled_window)

        count = 0
        for i in range(16):
            for j in ['x','y','z']:
                labels.append(gtk.Label("coil {}{}".format(i+1,j)))
                labels[-1].show()
                #table_layout.attach(labels[-1], count, 0, 0, 1, 0,0,0,0)
                table_layout.attach(labels[-1], 0, 1, count, count+1, 0,0,0,0)
                for ii, (filt_list, lab) in enumerate(zip([self.filt1, self.filt2, self.filt3, self.filt4, self.aa],lab_list[1:])):
                    #filt_list.append(gtk.CheckButton("Filt{}".format(ii+1)))
                    filt_list.append(gtk.CheckButton(lab))
                    #filt_list.append(gtk.CheckButton(''))
                    #filt_list[-1].connect("toggled", self.entry_checkbox, check_box)
                    filt_list[-1].set_active(True)  # Set the defaut
                    filt_list[-1].show()
                    table_layout.attach(filt_list[-1], ii+1, ii+2, count, count+1, 0,0,0,0)
                self.gain.append(gtk.combo_box_new_text())
                for iii in gains:self.gain[-1].append_text(iii)
                self.gain[-1].set_active(0)
                self.gain[-1].show()
                table_layout.attach(self.gain[-1], ii+2, ii+3, count, count+1, 0,0,0,0)
                count +=1
        table_layout.show()
        scrolled_window.add_with_viewport(table_layout)
        app_window.connect("delete_event", lambda w,e: gtk.main_quit())
        app_window.show()
        return

    def all_checkbox_changed(self, widget, data):
        act = widget.get_active()
        #print act
        for i in data:
            i.set_active(act)
    def all_gain_changed(self, widget, data):
        model = widget.get_active()
        #print model
        for i in data:
            i.set_active(model)
        

    def get_mds(self, widget):

        t=MDSplus.Tree('h1data',-1)#, 'READONLY')
        node=t.getNode('.mirnov:HMA_AMPS:SETTINGS')
        self.data2=[]; self.pastSettings=[]; self.newSettings=[]

        try:
            fileData=node.data()
            for line in fileData.split('\n'):
                formatted_line = line.strip('\n').strip('\r').strip(' ').strip('<').strip('>').strip(' ')
                if formatted_line != '':
                    formatted_line = map(int,formatted_line.split(' '))
                    self.pastSettings.append(formatted_line)
                    self.data2.extend(formatted_line)
            #print out the existing settings that are to be modified
            print self.pastSettings
            print 'Coil Settings read in Successfully'
        except:
            print 'exception reading data in from amp_setup'

        count = 0
        for i in range(16):
            for j, jj in enumerate(['x','y','z']):
                switchChip1 = self.pastSettings[i][j*2]
                switchChip2 = self.pastSettings[i][j*2+1]
                if switchChip1>=128:
                    self.filt4[count].set_active(True)
                    switchChip1 -= 128
                else:
                    self.filt4[count].set_active(False)
                    switchChip1 -= 64
                if switchChip1>=32:
                    self.filt3[count].set_active(True)
                    switchChip1 -= 32
                else:
                    self.filt3[count].set_active(False)
                    switchChip1 -= 16
                if switchChip1>=8:
                    self.filt2[count].set_active(True)
                    switchChip1 -= 8
                else:
                    self.filt2[count].set_active(False)
                    switchChip1 -= 4
                if switchChip1>=2:
                    self.filt1[count].set_active(True)
                    switchChip1 -= 2
                else:
                    self.filt1[count].set_active(False)
                    switchChip1 -= 1
                #print switchChip1,

                if switchChip2>=2:
                    self.filt1[count].set_active(True)
                    switchChip1 -= 2
                else:
                    self.filt1[count].set_active(False)
                    switchChip1 -= 1

                tmp_val = '{0:08b}'.format(switchChip2)[::-1]
                #print tmp_val
                if tmp_val[2]=='1':
                    self.aa[count].set_active(True)
                elif tmp_val[3]=='1':
                    self.aa[count].set_active(False)
                else:
                    raise ValueError('AA filter setting weird')

                if tmp_val[7]=='0' and tmp_val[6]=='1' and tmp_val[5]=='0':
                    self.gain[count].set_active(0)
                elif tmp_val[7]=='1' and tmp_val[6]=='1' and tmp_val[5]=='0':
                    self.gain[count].set_active(1)
                elif tmp_val[7]=='0' and tmp_val[6]=='0' and tmp_val[5]=='0':
                    self.gain[count].set_active(2)
                elif tmp_val[7]=='1' and tmp_val[6]=='0' and tmp_val[5]=='0':
                    self.gain[count].set_active(3)
                else:
                    raise ValueError("Gain setting weird")
                # if self.aa[count].get_active():
                #     switchChip2=4
                # else:
                #     switchChip2=8
                # val = self.gain[count].get_active()+1
                # if val==1:
                #     switchChip2=switchChip2+64
                # elif val==2:
                #     switchChip2=switchChip2+128+64
                # elif val==3:
                #     switchChip2=switchChip2+0
                # elif val==4:
                #     switchChip2=switchChip2+128
                # elif val==0:
                #     switchChip1=0
                #     switchChip2=0


                count += 1
        #print ''
        #print 'get mds data'
    def save_mds(self, widget):
        count = 0
        self.newSettings = []
        for i in range(16):
            cur_former = []
            for j in ['x','y','z']:
                switchChip1_new = int('{}{}{}{}{}{}{}{}'.format(int(self.filt4[count].get_active()), int(not self.filt4[count].get_active()),
                                                                int(self.filt3[count].get_active()), int(not self.filt3[count].get_active()),
                                                                int(self.filt2[count].get_active()), int(not self.filt2[count].get_active()),
                                                                int(self.filt1[count].get_active()), int(not self.filt1[count].get_active()),), 2)
                switchChip1 = 2 if self.filt1[count].get_active() else 1
                switchChip1 = switchChip1 + 8 if self.filt2[count].get_active() else switchChip1 + 4
                switchChip1 = switchChip1 + 32 if self.filt3[count].get_active() else switchChip1 + 16
                switchChip1 = switchChip1 + 128 if self.filt4[count].get_active() else switchChip1 + 64
                #print i, j, switchChip1, switchChip1_new
                cur_former.append(switchChip1)
                #Chip 2
                
                if self.aa[count].get_active():
                    switchChip2=4
                else:
                    switchChip2=8
                val = self.gain[count].get_active()+1
                if val==1:
                    switchChip2=switchChip2+64
                elif val==2:
                    switchChip2=switchChip2+128+64
                elif val==3:
                    switchChip2=switchChip2+0
                elif val==4:
                    switchChip2=switchChip2+128
                elif val==0:
                    switchChip1=0
                    switchChip2=0
                #print switchChip2

                cur_former.append(switchChip2)
                count +=1
            self.newSettings.append(cur_former)
        print self.newSettings

        #Format the data to be output to MDSplus
        anotherOutput=''
        for i in range(16):
            tmp_line = self.newSettings[i][0:6]
            line_str = ' '.join('%d' %i for i in tmp_line)
            output_str = '< %s >\n' %line_str
            anotherOutput+=output_str
        print anotherOutput


        # print 'save mds data'

def main():
    gtk.main()
    return 0

if __name__ == "__main__":
    MyProgram()
    main()
