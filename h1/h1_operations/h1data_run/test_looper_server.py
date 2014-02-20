import MDSplus
import time, re, datetime
import subprocess, os
import sys
import threading
import SocketServer

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

class shot_cycle():
    def __init__(self, cycle, executing_shot):
        '''
        print('Hello world, the shot cycle has started')
        '''
        self.cycle = cycle
        self.executing_shot = executing_shot

def start_up(cycle, executing_shot, TCP_PORT=8111):
    shot_cycle_obj = shot_cycle(cycle, executing_shot)

    #TCP_PORT = 8111
    TCP_server_thread, TCP_server = start_TCP_server(TCP_PORT)

    TCP_server.shot_cycle = shot_cycle_obj
    return shot_cycle_obj
