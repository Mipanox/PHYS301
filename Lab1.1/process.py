import pyfits
import numpy as np
import sys,os

def AverageBias(biasfiles):
    ''' 
    AverageBias produces a master bias image from a list of individual
    bias exposures.
           biasfiles  - list of file names
           returns a 2D numpy array
    '''
    biasdata = [pyfits.open(i)[0].data for i in biasfiles] 
    biascube = np.array(biasdata)
    
    ## use median, as suggested by Chromey in the bok
    medianimage = np.median(biascube,axis=0)
    masterbias  = medianimage
    return masterbias

def AverageDark(darkfiles,masterbias):
    darkdata = [pyfits.open(i)[0].data for i in darkfiles]
    darkexpo = [pyfits.open(i)[0].header['exposure'] for i in darkfiles]

    darkcube     = np.array(darkdata)
    darkexpocube = np.array(darkexpo)
    
    darklist=[]
    for (image,time) in zip(darkcube,darkexpocube):
        ## subtraict bias first
        cleandark = image - masterbias
        normdark  = cleandark / time
        darklist.append(normdark)
        
    cleandark2 = np.array(darklist) 
    ## median again. Let's not bother doing more complicated computations here
    masterdark = np.median(cleandark2,axis=0)
    return masterdark

def AverageFlat(flatfiles,masterbias,masterdark):
    flatdata = [pyfits.open(i)[0].data for i in flatfiles]
    flatexpo = [pyfits.open(i)[0].header['exposure'] for i in flatfiles]
    
    ## let's exploit numpy's broadcasting...
    flatcube  = np.array(flatdata)
    flatexpo_ = np.array(flatexpo).reshape(1,1,len(flatexpo)).transpose()
    flatexpocube = np.broadcast_to(flatexpo_,(len(flatexpo),
                                              flatcube.shape[1],
                                              flatcube.shape[2]))
    
    darkcorr  = masterdark * flatexpocube
    cleanflat = flatcube - darkcorr - masterbias
    
    normflat = cleanflat / np.median(cleanflat)
    masterflat = np.median(normflat,axis=0)
    return masterflat

def ScienceExposure(rawscidata,masterbias,masterdark,masterflat):
    rawimage = rawscidata.data
    expotime = rawscidata.header['exposure']
    scienceimage = (rawimage - expotime * masterdark - masterbias) / masterflat
    return scienceimage