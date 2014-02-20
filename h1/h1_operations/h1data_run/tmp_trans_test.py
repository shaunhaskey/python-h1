import subprocess
import time
def start_checking_transmitters():
    args = ['ssh', 'prl@prl60', '/home/prl/transmitter_snmp/check_trips2.py']
    a = subprocess.Popen(args, stdout=subprocess.PIPE, 
                     stderr=subprocess.PIPE, stdin=subprocess.PIPE)
    return a

def finish_checking_transmitters(a):
    b = None
    while b == None:
        b = a.poll()
        print 'b:', b
        time.sleep(2)
    out, err = a.communicate()
    return a.returncode

print 'starting transmitter check'
a = start_checking_transmitters()
print 'checking transmitter return code'
return_code = finish_checking_transmitters(a)
print 'finished'
print 'return code:',return_code
