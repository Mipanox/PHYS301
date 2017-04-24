"""
Photometry stuff.
Compute fluxes for sources/background defined by aperture
and creates source catalog.

No error propagation yet.

(Last modified: 04/23/17)
"""

###############################################

import numpy as np
import fitmodel
import pyfits, ldac

################################################
def backgroundEstimator(image):
    """
    Estimate background flux
    """
    ## maximum 5 iterations
    numberofTimes=5
    
    mean   = np.nanmean(image)
    median = np.nanmedian(image)
    if median > mean:
        return np.median(image)
    
    modeold=0.0
    ## consider cases w/ large "noises"
    while(numberofTimes>0):
        mean   = np.nanmean(image)
        median = np.nanmedian(image)
        
        mode = 3*median - 2*mean
        std = np.nanstd(image)
        
        if np.abs(modeold-mode)<=100.:
            return mode
        
        ## omit weird-behaving pixels
        image[np.where(image>median+3*std)] = np.nan
        modeold = mode
        numberofTimes -= 1
        
    return modeold

def aperFlux(image, x, y, source_radius = 10, background_width = 10, 
             gain = 3.0, nimages = 1, errflag=False):
    """ 
    Measures the flux (in adc counts), gaussian flux
    uncertainty, and the probability that the background is
    described by a constant.
    
    Parameters
    - image: 
      A 2D numpy array containing the image, from hdu.data
    - x: 
      The x position of the object, in pixels
    - y: 
      The y position of the object, in pixels
    - source_radius:
      The radius of a circular region centered on the object 
      representing the source region
    - background_width: 
      Width of an annulus, with
      source_radius < R <= source_radius + background_width,
      representing the background region
    - gain: The gain of the detector
    - nimages: 
      The number of images that were coadded to form this image
    - errflaf: 
      This bool determines whether or not errors
      should be determined for the flux. Used later.
    """
    
    ysize, xsize = image.shape
    X,Y = np.meshgrid(np.arange(xsize), np.arange(ysize))

    ## distance. '-1' accounts for convention difference in indexing
    dR = np.sqrt((X - x - 1)**2 + (Y - y - 1)**2)
    
    ## in source or in background masks
    inSource = dR <= source_radius
    inBackground = np.logical_and(dR > source_radius, 
                                  dR <= (background_width + source_radius))

    ## counting pixels
    nsourcepix = len(image[inSource])
    nbackgroundpix = len(image[inBackground])

    
    A     = np.sum(image[inSource])
    B_bar = backgroundEstimator(image[inBackground])

    ## correct for background
    flux = (A - nsourcepix*B_bar)
    
    ## error handling
    fluxerr=0.
    background_prob=0.
                   
    #if (errFlag):
        #fluxerr=?
        #background_prob=?

    return flux, fluxerr, background_prob

#############################################

def aperMag(flux, fluxerr, errflag=False):
    """
    Converts a measured flux and fluxerr into a magnitude
    and magnitude error. Assumes gaussian error propagation.
    
    Parameters
    
    - flux: 
      the flux of the object (may be a 1D array)
    - fluxerr: 
      the gaussian flux uncertainty of the object (may be a 1D array)
      
    Returns
    - mag, magerr: 
      the magnitude and gaussian uncertainty in the magnitude
    """
    
    ## flux to magnitude
    mag = -2.5*np.log10(flux)
    magerr = 0.0

    return mag, magerr
    
    
###############################################

def createPhotCat(imageArray, detectcat, filters, **otherparams):
    """
    Given an image and a detection catalog, this
    function will calculate the flux at every x,y position in the
    detection catalog and return the results as a new catalog
    
    Parameters
    - imageArray:
      A list/array for pyfits HDU containing coadded images in different filters
    - detectcat
      An ldac catalog containing an xcol and a ycol
    - otherparams: 
      Any additional function arguments passed to this function with 
      key=val pairs will be passed to aperFlux.

    Returns
    - photcat:
      an LDAC catalog with a flux, fluxerr, mag, magerr, and
      background_prob column
    """

    ## check colors
    if len(imageArray)!=len(filters):
        print "You need to spectify the name for each image in imageArray"
        return None
    
    HDU = detectcat
    detectcat = detectcat.data
    newColumnsArray = []
    
    ## iterate through each filter
    for index, image in enumerate(imageArray):
        nimages = image.header['NCOMBINE']

        ## variables
        fluxs = np.zeros(len(detectcat['X_IMAGE']))
        fluxerrs = np.zeros(len(detectcat['X_IMAGE']))
        backprobs = np.zeros(len(detectcat['X_IMAGE']))
        mags = np.zeros(len(detectcat['X_IMAGE']))
        magerrs = np.zeros(len(detectcat['X_IMAGE']))

        ## loop over the objects
        for i, x,y in zip(np.arange(len(detectcat['X_IMAGE'])),
                          detectcat['X_IMAGE'], 
                          detectcat['Y_IMAGE']):

            ## flux
            fluxs[i], fluxerrs[i], backprobs[i] = \
                     aperFlux(image.data, x,y,nimages = nimages, **otherparams)
            ## magnitude
            mags[i], magerrs[i] = aperMag(fluxs[i], fluxerrs[i])

        ## create columns for a new catalog
        newcols=[]
        newcols.extend([pyfits.Column(name = 'flux_'+filters[index],
                                  format = 'E', array = fluxs),
                        pyfits.Column(name = 'fluxerr_'+filters[index],
                                  format = 'E', array = fluxerrs),
                        pyfits.Column(name = 'backprob_'+filters[index],
                                  format = 'E', array = backprobs),
                        pyfits.Column(name = 'mag_'+filters[index],
                                  format = 'E', array = mags),
                        pyfits.Column(name = 'magerr_'+filters[index],
                                  format = 'E', array = magerrs)])
        
        newColumnsArray.append(newcols) 
    
    ## convert the columns into an LDAC catalog, copying the existing
    ## columns from the detection cat
    totalColumns=HDU.columns
    
    for column in newColumnsArray:
        totalColumns+=pyfits.ColDefs(column)
        
    photcat = ldac.LDACCat(pyfits.new_table(totalColumns,header=HDU.header))
    return photcat