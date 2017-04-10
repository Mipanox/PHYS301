import pyfits
import numpy as np
import sys,os

#This function statistically combines bias images into a masterbias image
def AverageBias(biasfiles):
    ''' 
    AverageBias produces a master bias image from a list of individual
    bias exposures.
           biasfiles  - list of file names
           returns a 2D numpy array
    '''
    biasdata= [pyfits.open(i)[0].data for i in open(biasfiles)] 
    
    biascube=np.array(biasdata)
    medianimage=np.median(biascube,axis=0)
    meanimage=np.mean(biascube,axis=0)
    sigmaimage=np.std(biascube,axis=0)
    
    masterbias=0#What goes in here? Mean, median, or sigma? Or something else?

    return masterbias #This is the end of the function 




#This function does the combining of dark currents
def AverageDark(darkfiles,masterbias):

    darkdata=[pyfits.open(i)[0].data for i in open(darkfiles)]
    
    #What is happening here? What is being put into this list?
    darkexpo=[pyfits.open(i)[0].header['exposure'] for i in
              open(darkfiles)]

    darkcube=np.array(darkdata) #This is an array of 2-D images
    darkexpocube=np.array(darkexpo) #This is an array of scalars.
    darklist=[]
    for (image,time) in zip(darkcube,darkexpocube): #The zip command here
        #loops over two lists simultaneously. The iterator of darkcube is
        #image, the iterator of darkexpocube is time

        ### !!! TODO FINISH THIS FUNCTION !!!

        cleandark=(image-masterbias)
        normdark=cleandark/time #(Dark Image- Bias Image) / Time, gives
#image/s which we need, since dark current for any observation depends on
#the exposure time.
        darklist.append(normdark)
    cleandark2=np.array(darklist) 
    masterdark=0#What goes in here? 
    return masterdark


#This function creates a combined flat field image
def AverageFlat(flatfiles,masterbias,masterdark):
    flatdata=[pyfits.open(i)[0].data for i in open(flatfiles)]
    flatexpo=[]#We need an array of exposure times. How do we get these?
    flatcube=np.array(flatdata)
    flatexpocube=np.array(flatexpo)
    cleanflatlist=[]
    for (image,time) in zip(flatcube,flatexpocube):
        cleanflat=image-0#bias and darkframe. we don't actually need 'flatexpocude', do we?
        normflat=cleanflat/np.median(cleanflat,axis=None)
        cleanflatlist.append(normflat)
    cleancube=np.array(cleanflatlist)
    masterflat=np.median(cleancube,axis=0)
    return masterflat

#This function creates the processed science image after combined bias,
#dark, and flat images have been created.  
def ScienceExposure(rawscidata,masterbias,masterdark,masterflat):
    rawimage=rawscidata.data
    expotime=rawscidata.header['exposure']
    scienceimage=(rawimage-expotime*masterdark-masterbias)/masterflat
    return scienceimage


# This is the end of the functions. 

