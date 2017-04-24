"""
Functions for calibrating CCD images:
- Bias : selective rejection (default 3.5 sigma clip)
- Dark frame : median
- Flat field : mild selective rejection mainly to remove 
               residual hot pixels (default to 1.5 sigma clip)

(Last modified: 04/23/17)
"""

import pyfits
import numpy as np
import sys,os


def AverageBias(biasfiles, k=3.5, method='rej'):
    """
    AverageBias produces a master bias image from a list of individual
    bias exposures.
           biasfiles  - list of file names
           returns a 2D numpy array
           
    Algorithm: selective rejection or median
    
    Inputs
    - biasfiles: list of strings
    - k: level of clipping (in sigma)
    
    Options
    - method: string
      'rej' or 'med'. Choice for algorithm.
      Some observations may be better calibrated by one or the other.
    """
    # opens each bias image file and stores the 2d images in a list
    biasdata = [pyfits.open(i)[0].data for i in biasfiles] 
    biascube = np.array(biasdata)

    meanimage  = np.mean(biascube,axis=0)
    sigmaimage = np.std(biascube,axis=0)
    
    if method=='rej':
        ## clipping
        mask = np.where(biascube[:] > meanimage + k*sigmaimage)
        biascube[mask] = np.nan
        masterbias = np.nanmedian(biascube,axis=0)
        
    elif method=="med":
        mastebias = np.median(biascude,axis=0)
        
    else:
        raise NameError('No such algorithm! Use "rej" or "med"')

    return masterbias

def AverageDark(darkfiles,masterbias):
    """
    AverageDark produces a master dark frame image from a list of individual
    dark frame exposures.
           darkfiles  - list of file names
           returns a 2D numpy array

    Inputs:
    - darkfiles: list of strings
    - masterbias: the master bias image
    """
    darkdata = [pyfits.open(i)[0].data for i in darkfiles]
    darkexpo = [pyfits.open(i)[0].header['exposure'] for i in darkfiles]

    darkcube = np.array(darkdata)
    darkexpocube = np.array(darkexpo)
    
    darklist=[]
    for (image,time) in zip(darkcube,darkexpocube): 
        cleandark = (image - masterbias)
        normdark  = cleandark / time 
        darklist.append(normdark)
        
    cleandark2=np.array(darklist) 
    ## median
    masterdark = np.median(cleandark2,axis=0) 
    
    return masterdark


def AverageFlat(flatfiles,masterbias,masterdark,k=1.5):
    """
    AverageFlat produces a master flat field image from a list of individual
    flat field exposures.
           flatfiles  - list of file names
           returns a 2D numpy array

    Inputs:
    - flatfiles: list of strings
    - masterbias: the master bias image
    - masterbias: the master dark frame image
    
    Options:
    - k: float
      The level of sigma clipping to remove any residual 'hot pixels'
      Defaults to 1.5
    """
    
    flatdata = [pyfits.open(i)[0].data for i in flatfiles]
    flatexpo = [pyfits.open(i)[0].header['exposure'] for i in flatfiles]
    flatcube = np.array(flatdata)
    flatexpocube = np.array(flatexpo)
    
    cleanflatlist=[]
    for (image,time) in zip(flatcube,flatexpocube):
        cleanflat = image - time*masterdark - masterbias
        normflat = cleanflat / np.nanmedian(cleanflat,axis=None)
        cleanflatlist.append(normflat)
        
    cleancube  = np.array(cleanflatlist)
    meanimage  = np.mean(cleancube,axis=0)
    sigmaimage = np.std(cleancube,axis=0)
    
    ## sigma clipping
    mask = np.where(cleancube[:]>meanimage+k*sigmaimage)
    cleancube[mask] = np.nan
    masterflat = np.nanmedian(cleancube,axis=0)
    
    return masterflat

def ScienceExposure(rawscidata,masterbias,masterdark,masterflat):
    """
    Do the calibration based on the master bias, dark frame,
    and flat images.
    """
    rawimage = rawscidata.data
    expotime = rawscidata.header['exposure']
    
    scienceimage = (rawimage-expotime*masterdark-masterbias)/masterflat
    
    return scienceimage

def batchScienceExposure(imageList, masterbias,masterdark,masterflat, outDir=None):
    """
    This function helps you to reduce a bunch of images. 
    
    Parameters:    
    - imageList: list of string
        a list of path for images needed to be reduced. 
    - masterbias: np.ndarray
        bias field
    - masterdark: np.array
        dark field
    - masterflat: np.ndarray
        flat field
    
    Returns:
    - outList: list of np.ndarray
        list of reduced images. 
    """
    outList=[]
    for imageName in imageList:
        rawData = pyfits.open(imageName)[0]
        newData = ScienceExposure(rawData,masterbias,masterdark,masterflat)
        if outDir is not None:
            name = imageName.split("/")[-1]
            sciencehdu = pyfits.PrimaryHDU(newData,header=rawData.header)
            sciencehdu.writeto(outDir+"Processed_"+name, clobber=True)
            del sciencehdu
        outList.append(newData)
    return outList