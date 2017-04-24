"""
Make HR (color-magnitude) diagram for a given cluster
Color is defined as Mag(B)-Mag(G) (B-V)
and magnitude is represented by Mag(G)
"""
import sys
sys.path.append("/afs/ir.stanford.edu/class/physics100/workdir/g2/Jason/codes/")

import sextractor as se
import pyfits
import numpy as np
import matplotlib.pyplot as plt
import photometry as pho
from matplotlib.colors import LogNorm
from util import *

plt.rcParams['axes.titlesize'] = 20
plt.rcParams['axes.labelsize'] = 20
plt.rcParams['xtick.labelsize'] = 20
plt.rcParams['ytick.labelsize'] = 20

class clusterHRdiag(object):
    def __init__(self,workdir,name,basename):
        ## directory of source coadded files
        self.name = name
        
        ## working directory for output catalogs
        self.workdir = workdir
        
        ## check files
        self.check_image_name = self.workdir+"check/check.fits"

        ## source name
        self.basename = basename
        
    def readfits(self):
        name = self.name
        
        ## use G as magnitude
        self.imageFits_B = pyfits.open(name.format("B"))[0]
        self.imageFits_G = pyfits.open(name.format("G"))[0]
        
        self.image = self.imageFits_G.data
        
    def runSE(self,flagthr=4,edge=15,**kwargs):
        """
        Parameters
        - flagthr: integer (float)
          Threshold for selecting good sources based on SExtractor flags.
          Defaults to 4 for contrasting good/bad sources in `dispCatalog`
          Should change to 1 before plotting HR diagram.

        SExtractor flags:
        ----------------

        1     The object has neighbors, bright and close enough to 
              significantly bias the photometry, or bad pixels 
              (more than 10% of the integrated area affected).
        
        2     The object was originally blended with another one.
        
        4     At least one pixel of the object is saturated 
              (or very close to).
        
        8     The object is truncates (to close to an image boundary).
        
        16    Object's aperture data are incomplete or corrupted.
        
        32    Object's isophotal data are incomplete or corrupted.
        
        64    A memory overflow occurred during deblending.
        
        128   A memory overflow occurred during extraction.
        
        
        For instance, to isolate blended and/or multiple sources,
        we would like to reject sources with FLAG < 4.
        
        =====================
        column of a catalog:
        - 'NUMBER': 
           Running object number 
        - 'MAG_ISO', 'MAGERR_ISO': 
           Isophotal magnitude
        - 'FLUX_APER', 'FLUXERR_APER':  
           Flux vector within fixed circular aperture(s) 
        - 'MAG_APER', 'MAGERR_APER': 
           Fixed aperture magnitude vector
        - 'MAG_AUTO', 'MAGERR_AUTO': 
           Kron-like elliptical aperture magnitude
        - 'X_IMAGE', 'Y_IMAGE':
           pixel location along the axes
        - 'THETA_IMAGE': 
           Position angle (CCW/x)[deg]
        - 'ELONGATION': 
           A_IMAGE/B_IMAGE  
        - 'ELLIPTICITY':
           1 - B_IMAGE/A_IMAGE  
        - 'FWHM_IMAGE': 
           FWHM assuming a gaussian core
        - 'FLAGS': 
           Extraction flags (see above)
        """
        self.readfits()
        
        se.sextractor(self.name.format("G"), \
                      "Processed_{0}_G.cat".format(self.basename),\
                      backphoto_type='LOCAL',
                      detection_type='CCD',
                      checkimage_type='SEGMENTATION',
                      checkimage_name=self.check_image_name,
                      **kwargs)
        
        ## catalog
        cat = pyfits.open(self.workdir+"Processed_{0}_G.cat".format(self.basename))
        catData = cat[2].data
        
        ## get rid of edges and inappropriate sources
        mask = np.where((catData['FLAGS']<flagthr)&
                        (catData['X_IMAGE']>edge)&(catData['X_IMAGE']<(self.image.shape[1]-edge))&
                        (catData['Y_IMAGE']>edge)&(catData['Y_IMAGE']<(self.image.shape[0]-edge)))
        
        self.catData_mask = catData[mask]
        cat[2].data = self.catData_mask
        self.cat = cat
        self.catData = cat[2].data
        
        ## check that FWHM is reasonable
        self.fwhmMedian = np.median(self.catData_mask['FWHM_IMAGE'])
        print 'The median of FWHM of the cataloged sources is {0:.3f} pixels'.format(self.fwhmMedian)
        
    
    def fluxArray(self,src_width=1.5,bgd_width=2,low=3e2,high=15e2):
        """ Flux subtraction and cataloging """
        
        fluxArray=[]
        for i in range(len(self.catData_mask['X_IMAGE'])):
            ## x, y positions
            x = self.catData_mask['X_IMAGE'][i]
            y = self.catData_mask['Y_IMAGE'][i]
            
            fwhm = self.catData_mask['FWHM_IMAGE'][i]
            flux, fluxerr, bgd_prob = pho.aperFlux(self.image, x, y, 
                                                   source_radius = src_width*self.fwhmMedian , 
                                                   background_width = bgd_width*self.fwhmMedian,
                                                   gain = 3.0, nimages = 1, errflag=False)
            fluxArray.append([x,y,fwhm,flux])
        
        ## Is there nan, zero, or negative fluxes?
        nanArrray = [[i[0],i[1],i[2]] for i in fluxArray if np.isnan(i[3])]
        print 'Nan fluxes: [X, Y, FWHM]\n {0}'.format(nanArrray)
        zeroArray = np.array([[i[0],i[1],i[2]] for i in fluxArray if i[3]==0])
        print 'Zero fluxes: [X, Y, FWHM]\n {0}'.format(zeroArray)
        negArray = [[i[0],i[1],i[2]] for i in fluxArray if i[3]<0]
        print 'Negative fluxes: [X, Y, FWHM]\n {0}'.format(negArray)
        
        ## create catalog
        phoCat = pho.createPhotCat([self.imageFits_G,self.imageFits_B], 
                                   self.cat[2], ['G','B'], 
                                   source_radius = src_width*self.fwhmMedian,
                                   background_width = bgd_width*self.fwhmMedian)
        
        self.color = phoCat['mag_B']-phoCat['mag_G']
        self.magni = phoCat['mag_G']
        
        ## examine if the "radii" are well-chosen
        plt.figure(figsize=(48,48))
        plt.imshow(self.image, vmin=low, vmax=high, 
                   origin="lowerleft",norm=LogNorm(), 
                   cmap="Greys")
        circles(self.catData_mask['X_IMAGE']-1,
                self.catData_mask['Y_IMAGE']-1,
                src_width*self.fwhmMedian,lw=0.5,c='r',fc='none')
        circles(self.catData_mask['X_IMAGE']-1,
                self.catData_mask['Y_IMAGE']-1,
                bgd_width*self.fwhmMedian,lw=0.5,c='r',ls='--',fc='none')

        
    def dispCatalog(self,low=3e2,high=15e2):
        """ Display extracted catalog with flags """
        
        plt.figure(figsize=(96,48))
        plt.subplot(121); plt.title('G')
        plt.imshow(self.image, vmin=low, vmax=high, 
                   origin="lowerleft",norm=LogNorm(), 
                   cmap="Greys")
        plt.scatter(self.catData_mask['X_IMAGE']-1,
                    self.catData_mask['Y_IMAGE']-1,s=240,
                    facecolors='none',edgecolors='r')
        for i in range(len(self.catData['FLAGS'])):
            plt.annotate(str(self.catData['FLAGS'][i]),
                         ((self.catData['X_IMAGE'][i]-1,self.catData['Y_IMAGE'][i]-1)),
                         size=50,color='yellow')
        
        plt.subplot(122); plt.title('B')
        plt.imshow(self.imageFits_B.data, vmin=low, vmax=high,
                   origin="lowerleft",norm=LogNorm(),
                   cmap="Greys")
        plt.scatter(self.catData_mask['X_IMAGE']-1,
                    self.catData_mask['Y_IMAGE']-1,s=240, 
                    facecolors='none', edgecolors='r')
        for i in range(len(self.catData['FLAGS'])):
            plt.annotate(str(self.catData['FLAGS'][i]),
                         ((self.catData['X_IMAGE'][i]-1,self.catData['Y_IMAGE'][i]-1)),
                         size=50,color='yellow')
            
    def dispHR(self):
        """ Plot the HR diagram """
        plt.figure(figsize=(24,18))
        plt.scatter(self.color,self.magni,s=10,facecolors='r',edgecolors='k')
        plt.ylim(-4,-12); plt.xlim(-2,2)
        plt.xlabel('Color Mmag(B)-Mag(G))')
        plt.ylabel('Magnitude (Mag(G))')