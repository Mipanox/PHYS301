"""
Calibration for spectral data

(Last modified: 05/26/17)
"""

import sys
sys.path.append("/afs/ir.stanford.edu/class/physics100/workdir/g2/Jason/codes/")

import pyfits
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import calibration as cp
from matplotlib.colors import LogNorm
from batch import *
import glob
import os

plt.rcParams['axes.titlesize'] = 60
plt.rcParams['axes.labelsize'] = 40
plt.rcParams['xtick.labelsize'] = 40
plt.rcParams['ytick.labelsize'] = 40

class CaliProcess(object):
    def __init__(self,datadir,basename,procdir,arcname):
        self.datadir  = datadir
        self.basename = basename
        
        ## where to store the processed files
        self.procdir  = procdir
        
        ## arc
        self.arcname = arcname
        
    def __call__(self,**kwargs):
        """
        **kwargs: slices for calibration file indices.
          Namely selecting appropriate files.
          - biasidx
          - darkidx
          - flatidx
        """
        ## initialize, in case of empty entries
        self.biasfiles, \
        self.darkfiles, \
        self.flatfiles  = [],[],[]
        
        ## lists of calibration files
        biaslist = [f for f in glob.glob(self.datadir+'/*') if 'bias' in f or 'Bias' in f]
        darklist = [f for f in glob.glob(self.datadir+'/*') if 'dark' in f or 'Dark' in f]
        flatlist = [f for f in glob.glob(self.datadir+'/*') if 'flat' in f or 'Flat' in f]
        
        if len(kwargs)>0:
            for key, value in kwargs.iteritems():
                setattr(self, key, value)
            
            self.biasfiles = biaslist[self.biasidx[0]:self.biasidx[1]]
            self.darkfiles = darklist[self.darkidx[0]:self.darkidx[1]]
            self.flatfiles = flatlist[self.flatidx[0]:self.flatidx[1]]
            
        else: ## all files
            self.biasfiles = biaslist
            self.darkfiles = darklist
            self.flatfiles = flatlist
            
    def makeCali(self,method,**kwargs):
        """
        Make calibration master images
        
        Parameters
        - method: 'rej' or 'med'
          Algorithm for master bias
        """
        
        print "Have you checked the calibration images?"
        print "If not, call `dispCali`               \n"
        
        ## 
        self.finalbias = cp.AverageBias(self.biasfiles,method=method,**kwargs) 
        self.finaldark = cp.AverageDark(self.darkfiles,self.finalbias,**kwargs) 
        self.finalflat = cp.AverageFlat(self.flatfiles,self.finalbias,self.finaldark,**kwargs) 
        
        ## write fits files
        newbias = self.procdir+'MasterBias.fits'
        newdark = self.procdir+'MasterDark.fits'
        newflat = self.procdir+'MasterFlat.fits'
        
        biashdu = pyfits.PrimaryHDU(self.finalbias)
        biashdu.writeto(newbias, clobber=True)

        darkhdu = pyfits.PrimaryHDU(self.finaldark)
        darkhdu.writeto(newdark, clobber=True)

        flathdu = pyfits.PrimaryHDU(self.finalflat)
        flathdu.writeto(newflat, clobber=True)
        
    def runCali(self,specflat=None,**kwargs):
        """
        Run calibration
        
        Options
        - Can calibrate only for bias, or bias+dark, or bias+flat.
          Use 'nodark' or 'noflat' accordingly
        - specflat: boolean
          If specified, use this SpecMasterFlat.fits to calibrate
        """
        datalist = [f for f in glob.glob(self.datadir+'/*') if self.basename in f]
        arclist  = [f for f in glob.glob(self.datadir+'/*') if self.arcname in f]
        
        if specflat is None:
            self.data_cln = cp.batchScienceExposure(datalist, self.finalbias, \
                                                    self.finaldark, self.finalflat, \
                                                    outDir=self.procdir,**kwargs)
        
            self.arc_cln  = cp.batchScienceExposure(arclist, self.finalbias, \
                                                    self.finaldark, self.finalflat, \
                                                    outDir=self.procdir,**kwargs)
        else:
            specflat = pyfits.open(specflat)[0].data
            
            self.data_cln = cp.batchScienceExposure(datalist, self.finalbias, \
                                                    self.finaldark, specflat, \
                                                    outDir=self.procdir,**kwargs)
        
            self.arc_cln  = cp.batchScienceExposure(arclist, self.finalbias, \
                                                    self.finaldark, specflat, \
                                                    outDir=self.procdir,**kwargs)
    #####################
    def dispCali(self,which,low=None,high=None):
        
        try:
            if which=='bias':
                filelist = self.biasfiles
            elif which=='dark':
                filelist = self.darkfiles
            elif which=='flat':
                filelist = self.flatfiles
                
        except:
            raise NameError('No such calibration file(s)!\n \
                             Use "bias", "dark", "flat"')
        
        plt.figure(figsize=(60,30))
        for i,f in enumerate(filelist):
            temp = pyfits.open(f)
            data = temp[0].data
            plt.subplot(int(len(filelist)/4)+1,4,i+1)
            plt.imshow(data, origin='lower')
            plt.axis('off'); plt.title(str(i))
            
            if low is not None:
                plt.clim(low,high)