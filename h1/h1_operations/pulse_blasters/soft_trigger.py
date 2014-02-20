import pb_class

#include a delta option
#maybe specify a particular shot for it to read in instead of a file?

dict_multiplier = {'ns':1, 'us':1000, 'ms':1000000, 's':1000000000}
dict_string = {'l':'0', 'h':'1'}
clock_freq = 120; board = 0

print '################################################'
print '#TRIGGER THE PULSE BLASTERS'
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
################################################
print 'triggering board return :',
tmp = PB.pb_start()
print tmp
if tmp != 0:
    raise ValueError("Error starting the program")
tmp = PB.pb_close()
print 'pb_close return :', tmp
if tmp != 0:
    raise ValueError("Error closing the board")
print 'Program finished'
