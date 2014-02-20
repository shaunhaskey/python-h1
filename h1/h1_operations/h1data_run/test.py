import subprocess, time


subproc = subprocess.Popen(["xterm", "-e","ssh", "prl@prl24", "-t", "bash -l server1.sh"], 
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE, 
                                     stdin=subprocess.PIPE)
#subprocess.call(['ssh', 'prl@prl24','-t','bash -l server1.sh'])
time.sleep(10)
subproc.kill()
#check to see if it is still alive?
