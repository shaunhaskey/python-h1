import MDSplus as MDS
import matplotlib.pyplot as pt
import numpy as np
T = MDS.Tree('electr_dens',81443)
devices = [1,2,4,7,8,9,10]
nodes = [1,1,1,4]
#nodes = [1,1,1,1,1,1,4]
fig, ax = pt.subplots()
fig2, ax2 = pt.subplots()
devices = range(7,11)
for i,node in zip(devices, nodes):
    node_name = '.camac.TR612_{}:input_{}'.format(i,node)
    n = T.getNode(node_name)
    print node_name, n.data().shape
    time_base = n.dim_of().data()
    period = (time_base[1]-time_base[0])/3.
    print 1./period/1000000.
    fft = np.fft.fft(n.data())
    loc = np.argmax(np.abs(fft))
    fft_freq = np.fft.fftfreq(len(fft),period)
    print fft_freq[loc]
    #ax.plot(time_base, n.data(),label='{}'.format(i))
    #ax2.plot(fft_freq, np.abs(fft))
ax.legend()
fig.canvas.draw(); fig.show()
fig2.canvas.draw(); fig2.show()
