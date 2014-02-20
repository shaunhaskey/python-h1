import subprocess, re
def kill_rem_jServer(computer):
    args = ['ssh',computer, 'ps aux |grep jServer|grep java']
    output,error = subprocess.Popen(args,stdout = subprocess.PIPE, stderr= subprocess.PIPE).communicate()
    output = output.split('\n')
    kill_pid_list = []
    for i in output:
        print i
        tmp = re.search('\d+',i)
        if tmp!=0 and i.find('grep')<0 and len(i)>3:
            kill_pid_list.append(re.search('\d+',i).group())
        else:
            print 'passing'
    for i in kill_pid_list:
        subprocess.call(['ssh', computer, 'kill -9 %s'%(i)])

computer = 'prl@prl24'
kill_rem_jServer(computer)
