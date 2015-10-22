import socket, time, sys

class antenna():
    def __init__(self, hostip, tree_node, shot):
        self.shot = shot
        self.ip, self.port = hostip.split(':')
        self.tree_loc = tree_node
    def init(self,string_to_send=None):
        '''initialise the labview device over a socket
        SH:20Mar2013
        '''
        try:
            error = None
            print '####################################################'
            print 'Antenna init has started : MDSplus Device, SH 1May2013'

            #build the string to send, and send it
            if string_to_send == None:
                self.build_string('INIT')
            else:
                self.send_string = string_to_send
                print self.send_string.rstrip(self.termination_character)
            self.send_receive()
            return 1
        except Exception,e:
            #Catch exceptions, then throw another exception so the dispatcher fails on this item
            if error is not None:
                e=error
            print "%s" % (str(e),)
            raise RuntimeError('!!! EXCEPTION - %s'%(str(e)))


    def build_string(self,phase):
        '''Put together a string to send to the labview server
        SH:20Mar2013
        '''
        self.termination_character = '\r\n'
        self.sep_char = chr(0x09) 
        print 'Building string to send to Labview:'
        self.send_string = phase + self.sep_char + str(self.tree_loc) + self.sep_char + str(self.shot) + self.termination_character
        print self.send_string.rstrip(self.termination_character)

    def send_receive(self):
        '''Send the command, and get the return status
        SH:20Mar2013
        '''
        ip = self.ip
        port = self.port
        #print 'Connecting to host: %s:%d'%(self.hostip.record,self.port.record)
        print 'Connecting to host: %s:%s'%(ip, port)
        self.s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.s.settimeout(5) #10 second timeout
        #self.s.connect((str(self.hostip.record), int(self.port.record)))
        self.s.connect((str(ip), int(port)))
        print 'Connected, sending string:%s'%(self.send_string.rstrip(self.termination_character))
        self.s.send(self.send_string)
        print '  Waiting for return message'
        data = ''
        count=0
        while (not data.endswith(self.termination_character)) or count>10:
            tmp = self.s.recv(100)
            #print '  received some data : %s'%(tmp)
            data += tmp
            count+=1
        print 'Finished receiving data, returned string :%s'%(data.rstrip(self.termination_character))
        self.data = data.rstrip(self.termination_character)
        self.return_value = self.data[0]
        if len(self.data)>2:
            self.return_message = self.data[1:]
        else:
            self.return_message = ''
        time.sleep(0.5)
        print 'Closing socket'
        self.s.close()
        print 'Socket closed'
        print 'Check return status from Labview'
        #check the returned value for success (0) or fail (1)
        if self.return_value=='0':
            print ' Initialisation success, returned message:%s'%(self.return_message)
        elif self.return_value=='1':
            raise RuntimeError(' Labview failed, returned message:%s'%(self.return_message))
        else:
            raise RuntimeError(' Labview failed unknown return value - not 1 or 0, returned message:%s'%(self.data))


    def store(self):
        '''tell the labview device it is time to store
        SH : 20Mar2013
        '''
        try:
            error = None
            print '####################################################'
            print 'Camera Store has started : MDSplus Device, SH 19Mar2013'
            self.build_string('STORE')
            self.send_receive()
            return 1
        except Exception,e:
            if error is not None:
                e=error
            print "%s\n" % (str(e),)
            raise RuntimeError('!!! EXCEPTION - %s'%(str(e)))
            return 0


