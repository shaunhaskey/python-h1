import ctypes, os
import numpy as np
import struct
import xml.dom.minidom as dom
import MDSplus as MDS
#import matplotlib.pyplot as pt

def extract_SPE_data(filename):
    in_file = file(filename,'rb')
    in_file.seek(1992)
    version, = struct.unpack('f',in_file.read(4))
    if version==3.0:
        print '3.0 SPE file,', 
        in_file.seek(678)
        xml_offset, = struct.unpack('Q',in_file.read(8))
        in_file.seek(xml_offset)
        xml_string = in_file.read()
        xml_parsed = dom.parseString(xml_string)
        #DataFormat, MetaFormat, Calibrations, DataHistories, GeneralInformation = xml_parsed.childNodes[0].childNodes
        tmp = xml_parsed.childNodes[0].childNodes
        DataFormat = tmp[0]
        frame_information = DataFormat.childNodes[0]
        atts = ['type','count','pixelFormat', 'size', 'stride','metaFormat']
        atts = ['type','count','pixelFormat', 'size', 'stride']#,'metaFormat']
        frame_info = {}
        for i in atts:
            print i, frame_information.attributes[i].value, 
            frame_info[i] = frame_information.attributes[i].value
        print ''
        data_info = frame_information.childNodes[0]
        region_info = ['count', 'stride', 'height', 'width', 'calibrations', 'type', 'size']
        frame_info['reg_info'] = {}
        for i in region_info:
            frame_info['reg_info'][i] = data_info.attributes[i].value

    #pixels -> rows -> regions -> frames
    n_frames = int(frame_info['count'])
    width = int(frame_info['reg_info']['width'])
    height = int(frame_info['reg_info']['height'])
    camera_data = np.ones((n_frames, width, height), dtype = np.uint16)
    starting_point = 4100
    #fig, ax = pt.subplots(ncols = n_frames)
    for i in  range(n_frames):
        in_file.seek(starting_point)
        c = in_file.read(int(frame_info['reg_info']['size']))
        b = np.fromstring(c,dtype=np.uint16)
        b.resize([int(frame_info['reg_info']['width']), int(frame_info['reg_info']['height'])])
        camera_data[i,:,:] = b.copy()
        starting_point += int(frame_info['reg_info']['size'])
        #cax = ax[i].imshow(b, aspect='auto')
        #cax.set_clim([0,20000])
    #fig.canvas.draw(); fig.show()
    #DataHistory = DataHistories.childNodes[0]
    #camera = DataHistory.childNodes[0].childNodes[0].childNodes[0].childNodes[0].childNodes[0]
    return camera_data, xml_string

    #Saving the data into MDSplus
def write_data_MDSplus(tree_name, shot, camera_data, xml_string, camera_data_node, xml_string_node):
    T = MDS.Tree(tree_name, shot)
    T.getNode(camera_data_node).putData(camera_data)
    T.getNode(xml_string_node).putData(xml_string)

    #T = MDS.Tree('imax_test',1)
    #print np.allclose(T.getNode('IMAX_DATA').data(), camera_data)
    #print T.getNode('XML_DATA').data() == xml_string

#filename = '/home/srh112/test.SPE'
#filename = '/home/srh112/imax_new/8.spe'
#camera_data, xml_string
