"""
Supplementary codes for calibration
"""
import numpy as np
import pyfits

def display_header(filelist,option='M81'):
    """Display the necessary information of each file"""
    for i in filelist:
        if 'list' in i[0]: # this is the list...
            continue
        if option not in i[0]:
            continue
        R = pyfits.open(i[0])
        hdr = R[0].header
        print i[0]
        print 'starttime ',hdr[8]
        print 'duration  ',hdr[9]
        print 'ccd temp  ',hdr[12]
        print 'filter    ',hdr[19]
        print ''
        
def image_range(filelist,option='bias'):
    """Print the max, med, min values in the images"""
    for i in filelist:
        if 'list' in i[0]: # this is the list...
            continue
        if option not in i[0]:
            continue
        R = pyfits.open(i[0])
        img = R[0].data
        print i[0]
        print 'Max = {0}, Median = {1}, Min = {2}'.format(np.max(img),np.median(img),np.min(img))