import numpy as np
import matplotlib.pyplot as pt
import re, copy, sys
import ctypes as c
import time as time_mod
import pb_class
import MDSplus

#include a delta option
#maybe specify a particular shot for it to read in instead of a file?
#multipliers must have the form ;(12msdL;15msdH)x30;

dict_multiplier = {'ns':1, 'us':1000, 'ms':1000000, 's':1000000000}
dict_string = {'l':'0', 'h':'1'}
clock_freq = 120; board = 0
file_loc = sys.argv[1]
try:
    model_shot = int(file_loc)
    print 'have model shot', model_shot
except ValueError:
    model_shot = None
    print 'model shot not given, using cfg file'

#define commands
CONTINUE = 0; STOP = 1; LOOP = 2; END_LOOP = 3; JSR = 4
RTS = 5; BRANCH = 6; LONG_DELAY = 7; WAIT = 8; RTI = 9
debug = 1
PULSE_PROGRAM = 0
produce_plots = 0
initialise_board = 1
clock_list = [1e6, 1e6, 5e3, 10e6] #Hz, set to None if not wanted...
write_MDSplus_data = 1
read_MDSplus_data = 0
min_time_command = 1000

if len(sys.argv)>2: produce_plots=int(sys.argv[2])

def write_MDSplus(put_data, t):
    node = t.getNode('.log.machine:timing')
    data = MDSplus.makeData(put_data)
    node.putData(data)

def read_MDSplus(t):
    node = t.getNode('.log.machine:timing')
    return node.record

def plot_waveforms(total_output):
    fig, ax = pt.subplots(nrows = 20, sharex = 1, sharey = 1)
    for i in total_output:
        #print i[1], int(i[2])
        ax[i[0]-1].plot(i[1], int(i[2]),'o')
    ax[0].set_ylim(-0.1,1.5)
    fig.canvas.draw(); fig.show()
                        
def extract_data(example, debug = 0):
    ex_split = example.split(';')
    try:
        channel = int(ex_split[0])
    except:
        raise ValueError("Error reading in channel number : %s"%(example))
    if channel>20:
        raise ValueError("Error reading line : %s, channel number too high"%(example))
    else:
        output = []
        for i in range(1, len(ex_split)):
            tmp = ex_split[i]
            if tmp[0]=='d':
                delta_time = 1
                tmp = tmp[1:]
            else:
                delta_time = 0
            a = re.findall(r'\d+', tmp)
            level = tmp[-1]
            if ((level!='l') and (level!='h')):
                raise ValueError("Error reading line : %s, part : %s level high or low not specified properly"%(example, ex_split[i]))
            level = dict_string[level]
            
            if len(a)>1:
                raise ValueError("Error reading line : %s, part : %s multiple time specification??"%(example, ex_split[i]))
            try:
                time = int(a[0])
            except:
                raise ValueError("Error reading line : %s, part : %s time part not integer??"%(example, ex_split[i]))
            time_unit = tmp[:-1].replace(a[0],'')
            if time_unit not in dict_multiplier:
                raise ValueError("Error reading line : %s, part : %s, time unit : %s,  not one of s, ns, us, or ms??"%(example, ex_split[i], time_unit))
            time =  time * dict_multiplier[time_unit]
            if (delta_time == 1) and (i==1):
                raise ValueError("Error reading line : %s, part : %s, can't have delta time as first command - delta to what?"%(example, ex_split[i]))
            if delta_time == 1:
                time = time + output[i-2][1]
            if debug:
                print channel, time, time_unit, level
            output.append((channel, time, level))
    return output

def plot_lists(output, ax_list):
    time_list = [0]
    level_list = [0]
    for i in output:
        time_list.append(i[1])
        level_list.append(level_list[-1])
        time_list.append(i[1])
        level_list.append(int(i[2]))
    ax_list[output[0][0]-1].plot(time_list, level_list, 'o-')

def look_for_multipliers(line):
    '''
    must have the form ;(12msdL;15msdH)x30;
    '''

    while line.find('(')>=0:
        start_loc = line.find('(')
        beginning = line[:start_loc].rstrip(';')+';'
        end_loc = line.find(')')
        if line[end_loc:].find(';')>=0:
            end_mult = end_loc + line[end_loc:].find(';')
            end = line[end_mult:].lstrip(';')
            last_item = False
        else:
            end_mult = end_loc + line[end_loc:].find('\n')
            end = line[end_mult:].lstrip(';')
            last_item = True
        mult_text = line[start_loc+1:end_loc].lstrip(';').rstrip(';') + ';'
        mult_num = int(line[end_loc+1:end_mult].lstrip('x'))
        print mult_text, mult_num
        line = beginning + mult_text*mult_num
        if last_item:
            line = line.rstrip(';') + end
        else:
            line = line + end
        print line

#READ INPUT FILE AND MAKE LIST OF COMMANDS

#mds_tree = MDSplus.Tree('h1data', -1)

if model_shot!=None:
    print 'reading data from model shot'
    mds_tree = MDSplus.Tree('h1data', model_shot)
    lines = read_MDSplus(mds_tree).rstrip('\n').split('\n')
    print lines
else:
    print 'reading data from cfg file'
    lines = file(file_loc, 'r').readlines()
    
#if read_MDSplus_data:
#    lines = read_MDSplus(mds_tree).rstrip('\n').split('\n')
#    print lines
#else:
#    lines = file(file_loc, 'r').readlines()

total_output = []
if produce_plots:
    fig, ax = pt.subplots(nrows = 20, sharex = 1, sharey = 1)

mdsplus_output_text = ''
for i in lines:
    #remove end of lines, putting ; at the end using spaces
    i = i.rstrip('\n').rstrip(';').replace(' ','').replace(':',';').replace(',',';').lower()
    if ((i[0]!='#') and (len(i)>2)):
        mdsplus_output_text += i + '\n'
        if debug:
            print i
        output = extract_data(i, debug = debug)
        if len(output)<=1:
            raise ValueError("Error not enough commands : %s"%(i))
        if output[-1][2]!='0':
            print 'Final logic state is not low for line : %s'%(i)
            print 'Are you sure you want to continue? (y,n)'
            logic_not_low_end = raw_input()
            if logic_not_low_end!='y':
                raise ValueError("Error reading line : %s , final logic state is not low"%(i))
        for j in total_output:
            if output[0][0] == j[0]:
                raise ValueError("Error reading line : %s Duplicate channel listing in config file"%(i))
        for j in range(1,len(output)):
            if (output[j][1] - output[j-1][1])<=0:
                raise ValueError("Error reading line : %s times not listed in increasing order"%(i))
        if produce_plots:
            plot_lists(output, ax)
        if debug:
            print output
        for j in output:
            total_output.append(j)


if produce_plots:
    ax[0].set_ylim(-0.1,1.5)
    for i, j in enumerate(ax):
        j.set_ylabel('%d'%(i+1))
        j.set_yticklabels([])
    ax[-1].set_xlabel('time ns')
    fig.canvas.draw(); fig.show()
    pt.show()

#sort the list by the times that the commands happen
total_output.sort(key=lambda x: int(x[1]))
if produce_plots:
    pass
    #plot_waveforms(total_output)
current_list =  ['0' for i in range(0,24)]
command_list = []
for j, i in enumerate(total_output):
    tmp_chan = i[0]; tmp_time = i[1]; tmp_str = i[2]
    current_list[24-tmp_chan] = tmp_str
    if j == 0:
        #make sure there is a first command in there
        command_list.append((tmp_time, copy.deepcopy(current_list), CONTINUE))
    else:
        if (tmp_time == command_list[-1][0]) or (j == 0):
            #if the time is hte same as the last command, replace it
            command_list[-1] = (tmp_time, copy.deepcopy(current_list), CONTINUE)
        else:
            #unique time, so gets a new command for it
            command_list.append((tmp_time, copy.deepcopy(current_list), CONTINUE))

#run through the command list to check that none of the temp times are too low
#check the list in output is in increasing order of time
for j in range(1,len(command_list)):
    if (command_list[j][0] - command_list[j-1][0])<min_time_command:
        raise ValueError("Command list error, some times are not far enough clock cycles apart")


################################################
#BUILD COMMAND LIST
#create a new command list, this time with dt instead
command_list_new = []
for i in range(0, len(command_list)):
    if i == (len(command_list)-1):
        #This is the final command, its not clear how long it should last for
        #Give it 3 seconds??
        command_list_new.append((10*10**9, command_list[i][1], command_list[i][2]))
        pass
    else:
        dt = command_list[i+1][0]-command_list[i][0]
        command_list_new.append((dt, command_list[i][1], command_list[i][2]))


if initialise_board!=1:
    for i in command_list_new:
        print "%-15d %s %s"%(i[0], ''.join(i[1]), i[2])
    exit("initialise_board not set to 1 therefore, not initialising boards")

if write_MDSplus_data:
    print mdsplus_output_text        
    if model_shot==None:
        mds_tree = MDSplus.Tree('h1data', -1)
        print 'no model_shot opening shot -1 to write data into'
    elif model_shot!=(-1):
        mds_tree = MDSplus.Tree('h1data', -1)
        print 'model_shot is not equal to -1, opening shot -1 to write data into'
    write_MDSplus(mdsplus_output_text, mds_tree)


print '################################################'
print '#INITIALISE THE BOARD AND SOFTWARE LIBRARY'
print 'load spinapi library'
PB = pb_class.pb() #this loads library and executs pb_init()
print 'select board number %d return :'%(board,), 
tmp = PB.pb_select_board(board)
print tmp
if tmp <0:
    raise ValueError("Error selecting board number")
print 'initialise board return :',
tmp = PB.pb_init()
print tmp
if tmp < 0:
    raise ValueError("Error initialising board")
print 'set board clock frequency  %d MHz'%(clock_freq,)
PB.pb_set_clock(clock_freq)
#print 'BOARD STATUS  : ', PB.pb_read_status()
print 'sending pb_stop command return :', 
tmp = PB.pb_reset()
#tmp = PB.pb_stop()
print tmp
if tmp < 0:
    raise ValueError("Error stopping board")
time_mod.sleep(0.1)
#print 'BOARD STATUS  : ', PB.pb_read_status()
################################################


print '#################################################'
print '#SET THE CLOCKS'
for i in range(0,len(clock_list)):
    freq = clock_list[i]
    if freq != None:
        print 'set clock %d, frequency %.2fHz'%(3-i,freq)
        period = int(float(clock_freq)/100.*10**9/freq) #60MHz correction included
        print '@@@@',period, period/2
        tmp = PB.pb_set_pulse_regs(3-i, period, int(period/2.0), 0)
        if tmp < 0:
            raise ValueError("Error setting output clock")
#################################################



print '###################################################'
print '#PROGRAM THE BOARD'
print 'pb_start_programming return :', 
tmp = PB.pb_start_programming(PULSE_PROGRAM)
print tmp
if tmp < 0:
    raise ValueError("Error start programming board")

#insert a dummy command to begin with, this leads to the wait command which is stalling for a trigger
#This current one flashes channel 20 so you know a new program has been received and started
start_command = PB.pb_inst_pbonly('0'*24, CONTINUE, 0, dict_multiplier['ms']*500)
if start_command != 0:
    raise ValueError("Error programming board - first command return not 0")
curr_command = 0
tmp = PB.pb_inst_pbonly(4*'0'+'1'+'0'*19, CONTINUE, 0, dict_multiplier['ms']*500)
curr_command += 1
if tmp != curr_command:
    raise ValueError("Error programming board - program return not incrementing")
tmp = PB.pb_inst_pbonly('0'*24, CONTINUE, 0, dict_multiplier['ms']*500)
curr_command += 1
if tmp != curr_command:
    raise ValueError("Error programming board - program return not incrementing")

wait_command = PB.pb_inst_pbonly('0'*24, WAIT, 0, command_list[0][0])
curr_command += 1
if wait_command != curr_command:
    raise ValueError("Error programming board - program return not incrementing")

#send off the commands
for i in command_list_new:
    word = ''.join(i[1])
    #print bin(int(word,2)), i[2],  0, i[0]/1000000
    tmp = PB.pb_inst_pbonly(word, i[2], 0, i[0])
    curr_command += 1
    if tmp != curr_command:
        raise ValueError("Error programming board - program return not incrementing")
    time_mod.sleep(0.1)

#final command to set all channels to '0'
tmp = PB.pb_inst_pbonly('0'*24, CONTINUE, 0, 500*dict_multiplier['ms'])
curr_command += 1
if tmp != curr_command:
    raise ValueError("Error programming board - program return not incrementing")

#branch back to the WAIT for trigger command
tmp = PB.pb_inst_pbonly('0'*24, BRANCH, wait_command, 500*dict_multiplier['ms'])
curr_command += 1
if tmp != curr_command:
    raise ValueError("Error programming board - program return not incrementing")

#FINISH PROGRAMMING AND SOFTWARE TRIGGER TO START EVERYTHING
print 'pb_stop_programming return:', 
tmp = PB.pb_stop_programming()
print tmp
if tmp != 0:
    raise ValueError("Error finish programming board")

#print 'BOARD STATUS  : ', PB.pb_read_status()

time_mod.sleep(0.1)
tmp = PB.pb_reset()
print 'PB_reset return :', tmp
time_mod.sleep(0.1)
tmp = PB.pb_start()
print 'pb_start return : ', tmp
if tmp != 0:
    raise ValueError("Error starting the program")

#print 'BOARD STATUS  : ', PB.pb_read_status()
#for i in range(0,15):
#    print 'BOARD STATUS  : ', PB.pb_read_status()
#    time_mod.sleep(2)
tmp = PB.pb_close()
print 'pb_close return :', tmp
if tmp != 0:
    raise ValueError("Error closing the board")
print 'Program finished'
