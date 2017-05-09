"""
Photometry stuff.
Compute fluxes for sources/background defined by aperture
and creates source catalog.

Including error and probability of background being a constant

(Last modified: 05/08/17)
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
             gain = 3.0, nimages = 1, errflag=False,
             RN2 = 1., DarkMed = 1., expTime = 1.):
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
    - gain: The gain of the detector. 
      For 24 inch: 6
      For 16 inch: 1.5 (estimated by TAs)
    - nimages: 
      The number of images that were coadded to form this image
      
    - errflag: 
      This bool determines whether or not errors
      should be determined for the flux. Used later.
    - RN2:
      readout noise squared determined by variance of bias images
    - DarkMed:
      Median of dark frame. Multiplied by gain equals variance in counts
      due to Poissonality.
    - expTime:
      Total exposure time
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
    flux = (A - nsourcepix*B_bar )/ expTime * gain
    
    ## error handling
    #-- initial
    fluxerr=0.
    background_prob=0.
                   
    if (errflag):
        ## Poisson. Total noise squared = F_obj*t + F_sky(per pix)*t*n +
        #                                 D*gain*t*n + RN2*n
        #  where F_obj and F_sky are given by the fluxes time gain
        fluxErrSq = A*gain
        skyErrSq  = B_bar*gain*nsourcepix
        darkErrSq = DarkMed*gain*expTime*nsourcepix
        RNErrSq   = RN2*gain*nsourcepix
        CountErr  = np.sqrt(fluxErrSq+skyErrSq+darkErrSq+RNErrSq)
        fluxerr   = (CountErr/expTime)/np.sqrt(nimages)
        
        ## background probability- assume every background pixel follows Poisson
        #  add up the chi square
        #  (where 1.25 represents the diff. betweem the error for median and for mean)
        b_std = np.sqrt(gain*image[inBackground].flatten()+DarkMed*gain*expTime+RN2*gain)/np.sqrt(nimages)*1.25
        b_res = (image[inBackground].flatten() - B_bar)*gain
        chi2 = np.sum((b_res/b_std)**2)
        
        ## probability of background being constant ( = B-bar)
        background_prob = sp.stats.chi2.sf(chi2, nbackgroundpix-1)

    return flux, fluxerr, background_prob

#############################################

def aperMag(flux, fluxerr):
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
    mag    = -2.5*np.log10(flux)
    magerr = 2.5/(np.log(10)*flux)*fluxerr

    return mag, magerr
    
    
###############################################

def createPhotCat(imageArray, detectcat, filters, refStarInf, RN2DarkMedInf,**otherparams):
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
    for index, images in enumerate(imageArray):
        refIdx, refTrue = refStarInf[index] 
        RN2DarkMeds = RN2DarkMedInf[index]
        
        ## loop through images, for reference star
        ReferenceFlux=[]
        ReferenceFluxErr=[]
        for i, image in enumerate(images):
            nimages = image.header['NCOMBINE']
            expTime = image.header['EXPTIME']
            
            ## define some arrays to hold our calculations
            x = detectcat['X_IMAGE'][refIdx]
            y = detectcat['Y_IMAGE'][refIdx]
            
            RN2, DarkMed = RN2DarkMeds[i]
            flux, fluxerr, background_prob = aperFlux(image.data,x,y,
                                                      nimages=nimages,
                                                      expTime=expTime,
                                                      RN2=RN2,
                                                      DarkMed=DarkMed,**otherparams)
            ReferenceFlux.append(flux)
            ReferenceFluxErr.append(fluxerr)
         
        ## loop through images, for targets
        rawFlux=[]
        rawFluxErr=[]
        for iRN, image in enumerate(images):
            nimages = image.header['NCOMBINE']
            expTime = image.header['EXPTIME']
            
            fluxs     = np.zeros(len(detectcat['X_IMAGE']))
            fluxerrs  = np.zeros(len(detectcat['X_IMAGE']))
            backprobs = np.zeros(len(detectcat['X_IMAGE']))
            
            RN2, DarkMed = RN2DarkMeds[iRN]
            ## over objects
            for i, x,y in zip(np.arange(len(detectcat['X_IMAGE'])),
                              detectcat['X_IMAGE'], 
                              detectcat['Y_IMAGE']):

                ## calculate the flux
                fluxs[i], fluxerrs[i], backprobs[i] = aperFlux(image.data,x,y,
                                                               nimages=nimages,
                                                               expTime=expTime, 
                                                               RN2=RN2, 
                                                               DarkMed=DarkMed,**otherparams)

            rawFlux.append(fluxs)
            rawFluxErr.append(fluxerrs)
            
        ReferenceFlux    = np.array(ReferenceFlux)
        ReferenceFluxErr = np.array(ReferenceFluxErr)
        rawFlux    = np.array(rawFlux)
        rawFluxErr = np.array(rawFluxErr)
        
        refmag, refmagerr = aperMag(ReferenceFlux, ReferenceFluxErr)
        rawmag, rawmagerr = aperMag(rawFlux, rawFluxErr)
        truemag    = rawmag-refmag[:,np.newaxis]+refTrue
        #-- error propagations
        truemagerr = np.sqrt(rawmagerr**2+refmagerr[:,np.newaxis]**2)
        
        weightedTrueMag    = np.average(truemag, weights=1/truemagerr**2., axis=0)
        weightedTrueMagErr = np.sqrt(1/np.sum(1/truemagerr**2., axis=0))
        
        ## create columns for a new catalog
        newcols=[]
        newcols.extend([pyfits.Column(name = 'truemag_'+filters[index],
                                  format = 'E', array = weightedTrueMag),
                        pyfits.Column(name = 'truemagerr_'+filters[index],
                                  format = 'E', array = weightedTrueMagErr)])
        newColumnsArray.append(newcols) 
        
    ## convert the columns into an LDAC catalog, 
    #  copying the existing columns from the detection cat
    totalColumns=HDU.columns
    
    for column in newColumnsArray:
        totalColumns+=pyfits.ColDefs(column)
        
    photcat = ldac.LDACCat(pyfits.new_table(totalColumns,header=HDU.header))
    
    return photcat    