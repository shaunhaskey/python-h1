import numpy as np

def parse_bin(s):
    t = s.split('.')
    return int(t[0], 2) + int(t[1], 2) / 2.**len(t[1])

def parse_bin_signed(s):
    orig = s
    sign = int(s[0])
    #print s
    s = ''.join([s[i] for i in range(1,len(s))])
    #print s
    new = 2**(s.find('.')) * sign * (-1)
    #print new
    t = s.split('.')
    if t[0]=='':t[0]='0'
    return int(t[0], 2) + int(t[1], 2) / 2.**len(t[1]) + new

# def parse_bin_signed2(s,):
#     orig = s
#     sign = int(s[0])
#     #print s
#     s = ''.join([s[i] for i in range(1,len(s))])
#     #print s
#     new = 2**(s.find('.')) * sign * (-1)
#     #print new
#     t = s.split('.')
#     return int(t[0], 2) + int(t[1], 2) / 2.**len(t[1]) + new
