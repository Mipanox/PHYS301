"""
Functions for HR diagram with flux correction and error propagation
See `clusterHRdiag.py` for some details omitted here
"""

from __future__ import division
import sys
sys.path.append("/afs/ir.stanford.edu/class/physics100/workdir/g2/Jason/codes/")

import sextractor as se
import pyfits
import numpy as np
import matplotlib.pyplot as plt
import photometry as pho
import glob
import ldac
from matplotlib.colors import LogNorm
from util import *

plt.rcParams['axes.titlesize'] = 20
plt.rcParams['axes.labelsize'] = 20
plt.rcParams['xtick.labelsize'] = 20
plt.rcParams['ytick.labelsize'] = 20

def disp_catalog(workdir,cadddir,batches,target,flagthr=8,
                 edge=15,vmin=50,vmax=1000,stupidthing=0):
    
    ## number of batches
    number_batches = int(len([item for item in glob.glob(cadddir+"*") \
                              if "coadd_batch" in item.split("/")[-1]])/3/2)
    batchImages = [batches.format(stupidthing+i) for i in range(number_batches)]
    
    ## use this one for display
    imageFits_G = pyfits.open(batchImages[0])[0]
    image = imageFits_G.data
    
    ## read-in pre-created catalog
    cat = pyfits.open(workdir+"Processed_{0}_G.cat".format(target))
    catData = cat[2].data
    
    mask = np.where((catData['FLAGS']<flagthr)&
                  (catData['X_IMAGE']>edge)&(catData['X_IMAGE']<(image.shape[1]-edge))&
                  (catData['Y_IMAGE']>edge)&(catData['Y_IMAGE']<(image.shape[0]-edge)))
    
    catData_mask=catData[mask]
    
    for c in cat[2].columns:
         c.array = c.array[mask]
    cat[2].data = catData_mask
    
    plt.figure(figsize=(20,20))
    
    med = np.median(image[np.where(image>0)])
    plt.imshow(image, 
               vmin=med-vmin,
               vmax=med+vmax,
               norm=LogNorm(),origin='lower', cmap="Greys")
    plt.scatter(catData_mask['X_IMAGE']-1,catData_mask['Y_IMAGE']-1,s=240, 
                facecolors='none', edgecolors='r')
    
    ## flags
    for i in range(len(catData['FLAGS'])):
        plt.annotate(str(catData['FLAGS'][i]),
                     ((catData['X_IMAGE'][i]-1,catData['Y_IMAGE'][i]-1)),
                     size=10,color='blue')
        
    fwhmMedian = np.median(catData_mask['FWHM_IMAGE'])
    print 'Median of FWHM is {0:.3f}'.format(fwhmMedian)
    
    return catData_mask, fwhmMedian, cat
    
def bgd_prob_plot(catData_mask,workdir,cadddir,batches,
                  radlim_steps,gain,fwhmMedian,stupidthing=0):
    """
    Returns the radius of maximum probability
    """
    ## number of batches
    number_batches = int(len([item for item in glob.glob(cadddir+"*") \
                              if "coadd_batch" in item.split("/")[-1]])/3/2)
    batchImages = [batches.format(stupidthing+i) for i in range(number_batches)]
    
    ##
    image = pyfits.open(batchImages[0])[0].data
    imageFits_G = pyfits.open(batchImages[0])[0]
    
    radius_background=[]
    radius = np.arange(radlim_steps[0],radlim_steps[1],radlim_steps[2])
    
    ## read in readout noise and median of dark arrays
    RN2, DarkMed = np.load(batchImages[0][:-5]+"_RN2_DarkMed.npy")
    
    ## loop over possible radii
    for outer in radius:
        background_prob_save=[]
        for i in range(len(catData_mask['X_IMAGE'])):
            expTime = imageFits_G.header['EXPTIME']
            flux, fluxerr, background_prob=pho.aperFlux(image, 
                                                        catData_mask['X_IMAGE'][i], 
                                                        catData_mask['Y_IMAGE'][i],
                                                        source_radius = 1.3*fwhmMedian, 
                                                        background_width = outer*fwhmMedian,
                                                        gain = gain, nimages =10., 
                                                        errflag=True, 
                                                        RN2=RN2, DarkMed=DarkMed, expTime=expTime)
            background_prob_save.append(background_prob)
        background_prob_save = np.array(background_prob_save)
        radius_background.append(np.mean(background_prob_save)) 
        
        radius_bgd = np.array(radius_background)
    
    plt.figure(figsize=(16,9))
    plt.plot(radius, radius_background)
    plt.xlabel("radius_background-radius_source (FWHM) ")
    plt.ylabel('Prob. of being constant')
    
    return radius[radius_bgd.argmax()]

def flux_of_ref_star(catData_mask,workdir,cadddir,batches,
                     ref_idx,gain,max_r,fwhmMedian,stupidthing=0):
    
    ## number of batches
    number_batches = int(len([item for item in glob.glob(cadddir+"*") \
                              if "coadd_batch" in item.split("/")[-1]])/3/2)
    batchImages = [batches.format(stupidthing+i) for i in range(number_batches)]
    
    ReferenceFlux=[]
    ReferenceFluxErr=[]
    for name in batchImages:
        imageFits_G = pyfits.open(name)[0]
        image = imageFits_G.data
        
        RN2, DarkMed = np.load(name[:-5]+"_RN2_DarkMed.npy")
        expTime = imageFits_G.header['EXPTIME']
        
        flux, fluxerr, background_prob = pho.aperFlux(image, catData_mask['X_IMAGE'][599], catData_mask['Y_IMAGE'][599],
                                                      source_radius = 1.3*fwhmMedian, 
                                                      background_width= max_r*fwhmMedian,
                                                      gain = gain, 
                                                      nimages = imageFits_G.header['NCOMBINE'], 
                                                      errflag=True, 
                                                      RN2=RN2, DarkMed=DarkMed, expTime=expTime)
        ReferenceFlux.append(flux)
        ReferenceFluxErr.append(fluxerr)
        
    plt.figure(figsize=(16,9))
    x = np.arange(len(batchImages))
    plt.errorbar(x, ReferenceFlux, yerr = ReferenceFluxErr)
    plt.xlabel("batch number"); plt.ylabel("flux")
    
def corr_catalog(cat,workdir,cadddir,batches,target,refStarInf,
                 gain,max_r,fwhmMedian,stupidthing=0):
    """
    Generate catalog for corrected magnitudes 
    after differential photometry and error propagation
    """
    # number of batches
    number_batches = int(len([item for item in glob.glob(cadddir+"*") \
                              if "coadd_batch" in item.split("/")[-1]])/3/2)
    batchImages = [batches.format(stupidthing+i) for i in range(number_batches)]
    
    ##
    batches_B = cadddir+target+"_B_coadd_batch{0}.fits"
    batches_G = cadddir+target+"_G_coadd_batch{0}.fits"
    
    batchImages_B = [batches_B.format(stupidthing+i) for i in range(number_batches)]
    batchImages_G = [batches_G.format(stupidthing+i) for i in range(number_batches)]
    
    imageFits_G = [pyfits.open(name)[0] for name in batchImages_G]
    imageFits_B = [pyfits.open(name)[0] for name in batchImages_B]
    
    ## reference star information
    refStarInf = refStarInf
    
    RN2DarkMedInf_G = [np.load(name[:-5]+"_RN2_DarkMed.npy") for name in batchImages]
    RN2DarkMedInf_B = [np.load(name[:-5]+"_RN2_DarkMed.npy") for name in batchImages_B]
    phoCat = pho.createPhotCat([imageFits_G, imageFits_B], cat[2], 
                               ['G','B'],  
                               refStarInf  =refStarInf,
                               RN2DarkMedInf = [RN2DarkMedInf_G,RN2DarkMedInf_B],
                               source_radius = 1.3*fwhmMedian,
                               background_width = max_r*fwhmMedian, 
                               gain = gain, errflag=True)
    
    phoCat.saveas(target+"_correctedMag.cat",clobber=True)
    
    return phoCat
    
def color_cat(target):
    """ Make a new catalog with 'color' """
    
    ## readin magnitude catalog generated from `corr_catalog`
    hdu = pyfits.open(target+'_correctedMag.cat')
    phoCat = hdu[1].data
    
    ## add color information
    colorB_G    = phoCat['truemag_B'] - phoCat['truemag_G']
    colorB_GErr = np.sqrt(phoCat['truemagerr_B']**2.+phoCat['truemagerr_G']**2.)
    
    cols = [] 
    cols.extend([pyfits.Column(name='colorB_G',    format = 'E', array= colorB_G),
                 pyfits.Column(name='colorB_GErr', format = 'E', array= colorB_GErr)]
               )
    
    orig_cols = hdu[1].columns
    new_cols = pyfits.ColDefs(cols)
    
    newhdu = pyfits.new_table(orig_cols+new_cols,header=hdu[1].header)
    newhdu.writeto(target+'_correctedMag_color.cat')
    
def HR_diagram(target,dist=None,dist_ref=None):
    """ 
    Plot full HR diagram with error bars and correct magnitudes 
    
    Options
    - dist: float
      The distance to the cluster, in parcsec
    """
    
    ## colored catalog
    newhdu = pyfits.open(target+'_correctedMag_color.cat')
    phoCat = newhdu[1].data
    
    mag = phoCat['truemag_G']
    
    ## absolute magnitude?
    if dist is not None:
        #-- correct back the difference in distance
        #   for target and reference star
        mag -= 5*(np.log10(dist)-np.log10(dist_ref))
        
    plt.figure(figsize=(10,10))
    plt.errorbar(newhdu[1].data['colorB_G'], mag, 
                 yerr = phoCat['truemagerr_G'], 
                 xerr = phoCat['colorB_GErr'], fmt='.')
    #plt.ylim(14,8)
    plt.xlabel('mag(B)-mag(G)'); plt.ylabel('mag(G)')
    plt.title(target)