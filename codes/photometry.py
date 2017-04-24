'''
photometry.py

This is a stub file provided for lab 1.3 part II. 

COPY THIS FILE TO A WORK DIRECTORY OF YOUR CHOOSING FIRST!!!!!!

Fill in the missing blanks in the following functions.
'''

###############################################

import numpy as np
import fitmodel
import pyfits, ldac


################################################
def backgroundEstimator(image):
    numberofTimes=5
    mean=np.nanmean(image)
    median=np.nanmedian(image)
    if median>mean:
       return np.median(image)
    modeold=0.0
    while(numberofTimes>0):
       mean=np.nanmean(image)
       median=np.nanmedian(image)
       mode=3*median-2*mean
       std=np.nanstd(image)
       if np.abs(modeold-mode)<=100.:
          return mode
       image[np.where(image>median+3*std)]=np.nan
       modeold=mode
       numberofTimes-=1
    return modeold

def aperFlux(image, x, y, source_radius = 10, background_width = 10, 
             gain = 3.0, nimages = 1, errflag=False):
    ''' aperFlux - Measures the flux (in adc counts), gaussian flux
             uncertainty, and the probability that the background is
             described by a constant.
         @parameter image - A 2D numpy array containing the image, from
             hdu.data
         @parameter x - The x position of the object, in pixels
         @parameter y - The y position of the object, in pixels
         @parameter source_radius - The radius of a circular region
             centered on the object representing the source region
         @parameter background_width - Width of an annulus, with
             source_radius < R <= source_radius + background_width,
             representing the background region
         @parameter gain - The gain of the detector
         @parameter nimages - The number of images that were coadded to
             form this image
         @parameter errflag- This bool determines whether or not errors
             should be determined for the flux. Only used starting with
             Lab 1.5, once we figure out the calculation of uncertainties

         @returns flux, fluxerr, background_prob
         '''
    ysize, xsize = image.shape

    #create 2D arrays with the X,Y coordinates of every pixel
    X,Y = np.meshgrid(np.arange(xsize), np.arange(ysize))

    #calculate the distance of every pixel from the object
    #note: the x-1 accounts for the offset between sextractor and array indicies
    dR = np.sqrt((X - x - 1)**2 + (Y - y - 1)**2)

    #2D boolean arrays selecting the background and source regions
    inSource = dR <= source_radius
    #if (x>285-50)&(x<285+50)&(y<275-50)&(y>275+50):
    #    dR=np.sqrt((X - x - 1-100)**2.)
    inBackground = np.logical_and(dR > source_radius, 
                                  dR <= (background_width + source_radius))

    #counting pixels
    nsourcepix = len(image[inSource])
    nbackgroundpix = len(image[inBackground])

    ### TODO: FINISH THESE LINES
    # calculate the flux of the source, the flux uncertainty, and the
    # probability that the background is described by a constant
    #
    # Feel free to add additional intermediate steps as needed. We will
    # need to calculate
    # the flux with Lab 1.3, but the uncertainties on the flux will wait
    # until Lab 1.5.  Set the fluxerr and background_prob to zero until
    # then. 

#    flux = 
    A=np.sum(image[inSource])
    B_bar=backgroundEstimator(image[inBackground])

    flux=(A-nsourcepix*B_bar)
    

    fluxerr=0.
    background_prob=0.
                      

    #if (errFlag):
        #fluxerr=?
        #background_prob=?

    
    return flux, fluxerr, background_prob

#############################################

def aperMag(flux, fluxerr, errflag=False):
    ''' aperMag - converts a measured flux and fluxerr into a magnitude
    and magnitude error. Assumes gaussian error propagation.
         @parameter flux  the flux of the object (may be a 1D array)
         @parameter fluxerr the gaussian flux uncertainty of the object
                          (may be a 1D array)
         @returns mag, magerr  - the magnitude and gaussian uncertainty in
                                    the magnitude
                                    z
    '''

    ### TODO: FINISH THESE LINES
    # calculate the magnitude and the gaussian magnitude uncertainty.
    #
    # Feel free to add additional intermediate steps as needed. Raw
    # magnitudes will be calculated starting in Lab 1.3, and errors will
    # be propagated in Lab 1.5. Set the magerr term to zero until then. 

    mag = -2.5*np.log10(flux)
    magerr = 0.0

    return mag, magerr
    
    
###############################################

def createPhotCat(imageArray, detectcat, filters, **otherparams):
    ''' createPhotCat - Given an image and a detection catalog, this
    function will calculate the flux at every x,y position in the
    detection catalog and return the results as a new catalog
         @parameter image - A pyfits HDU containing a coadded image
         @parameter detectcat - An ldac catalog containing an xcol and a
                                ycol
         @parameter otherparams - Any additional function arguments passed
    to this function with key=val pairs will be passed to aperFlux.

         @returns an LDAC catalog with a flux, fluxerr, mag, magerr, and
         background_prob column
         '''

##### TODO:  Either 1.) Use this function and write your own to combine
##### the photcats from an R, G, & B exposure into one catalog.
##### or
##### 2.) Rewrite this function to do photometry on RGB at put it into one catalog

    if len(imageArray)!=len(filters):
        print "you need to spectify the name for each image in imageArray"
        return None
    HDU =detectcat
    detectcat=detectcat.data
    newColumnsArray=[]
    for index, image in enumerate(imageArray):
        nimages = image.header['NCOMBINE']

    #define some arrays to hold our calculations
        fluxs = np.zeros(len(detectcat['X_IMAGE']))
        fluxerrs = np.zeros(len(detectcat['X_IMAGE']))
        backprobs = np.zeros(len(detectcat['X_IMAGE']))
        mags = np.zeros(len(detectcat['X_IMAGE']))
        magerrs = np.zeros(len(detectcat['X_IMAGE']))

    #loop over the objects to measure
        for i, x,y in zip(np.arange(len(detectcat['X_IMAGE'])),
                          detectcat['X_IMAGE'], 
                          detectcat['Y_IMAGE']):

            #calculate the flux; propagate parameters to the aperFlux function
            fluxs[i], fluxerrs[i], backprobs[i] = \
                     aperFlux(image.data, x,y,nimages = nimages, **otherparams)

            mags[i], magerrs[i] = aperMag(fluxs[i], fluxerrs[i])

    #create columns for a new catalog
        newcols=[]
        newcols.extend([pyfits.Column(name = 'flux_'+filters[index],
                                  format = 'E',
                                  array = fluxs),
                        pyfits.Column(name = 'fluxerr_'+filters[index],
                                  format = 'E',
                                  array = fluxerrs),
                        pyfits.Column(name = 'backprob_'+filters[index],
                                  format = 'E',
                                  array = backprobs),
                        pyfits.Column(name = 'mag_'+filters[index],
                                  format = 'E',
                                  array = mags),
                        pyfits.Column(name = 'magerr_'+filters[index],
                                  format = 'E',
                                  array = magerrs)])
        newColumnsArray.append(newcols) 
    #convert the columns into an LDAC catalog, copying the existing
    #columns from the detection cat
    totalColumns=HDU.columns
    for column in newColumnsArray:
        totalColumns+=pyfits.ColDefs(column)
    photcat = ldac.LDACCat(pyfits.new_table(totalColumns,header=HDU.header))
    return photcat

        
##############
