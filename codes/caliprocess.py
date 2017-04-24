"""
Calibration and coaddition.

(Last modified: 04/23/17)
"""

import sys
sys.path.append("/afs/ir.stanford.edu/class/physics100/workdir/g2/Jason/codes/")

import coaddition
import chto_coadd
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
    def __init__(self,datadir,basename):
        self.datadir  = datadir
        self.basename = basename
        
    def makeCali(self,calilist,method,**kwargs):
        """
        Make calibration master images
        
        Inputs
        - calilist: np.ndarray
          List of calibration files
        
        Parameters
        - method: 'rej' or 'med'
          Algorithm for master bias
          
        **kwargs: slices for calibration file indices.
          Namely selecting appropriate files.
          - biasidx
          - darkidx
          - flatRidx, flatGidx, flatBidx
        """
        
        print "Have you checked the calibration images?"
        print "If not, call `dispCali`               \n"
        
        if len(kwargs)>0:
            for key, value in kwargs.iteritems():
                self.key = value
            
            self.biasfiles = [i for i in calilist if 'bias' in i or 'Bias' in i][self.biasidx[0],self.biasidx[1]]
            self.darkfiles = [i for i in calilist if 'dark' in i or 'Dark' in i][self.darkidx[0],self.darkidx[1]]
                
            flatfiles = [i for i in calilist if 'flat' in i or 'Flat' in i][self.flatidx[0],self.flatidx[1]]
            
        else: ## all files
            self.biasfiles = [i for i in calilist if 'bias' in i or 'Bias' in i]
            self.darkfiles = [i for i in calilist if 'dark' in i or 'Dark' in i]
            flatfiles = [i for i in calilist if 'flat' in i or 'Flat' in i]
            
        self.flatfiles_R=[i for i in flatfiles if 'R' in i]
        self.flatfiles_G=[i for i in flatfiles if 'G' in i]
        self.flatfiles_B=[i for i in flatfiles if 'B' in i]
        
        ## 
        self.finalbias = cp.AverageBias(self.biasfiles,method=method) 
        self.finaldark = cp.AverageDark(self.darkfiles,self.finalbias) 
        self.finalflat_R = cp.AverageFlat(self.flatfiles_R,self.finalbias,self.finaldark) 
        self.finalflat_G = cp.AverageFlat(self.flatfiles_G,self.finalbias,self.finaldark) 
        self.finalflat_B = cp.AverageFlat(self.flatfiles_B,self.finalbias,self.finaldark) 
                        
        ## write fits files
        newbias = self.datadir+'Cali/'+self.basename+'_Bias.fits'
        newdark = self.datadir+'Cali/'+self.basename+'_Master_Dark.fits'
        newflat_R = self.datadir+'Cali/'+self.basename+'_Master_Flat_R.fits'
        newflat_G = self.datadir+'Cali/'+self.basename+'_Master_Flat_G.fits'
        newflat_B = self.datadir+'Cali/'+self.basename+'_Master_Flat_B.fits'
        
        biashdu = pyfits.PrimaryHDU(self.finalbias)
        biashdu.writeto(newbias, clobber=True)

        darkhdu = pyfits.PrimaryHDU(self.finaldark)
        darkhdu.writeto(newdark, clobber=True)

        flathdu_R = pyfits.PrimaryHDU(self.finalflat_R)
        flathdu_R.writeto(newflat_R, clobber=True)
        flathdu_G = pyfits.PrimaryHDU(self.finalflat_G)
        flathdu_G.writeto(newflat_G, clobber=True)
        flathdu_B = pyfits.PrimaryHDU(self.finalflat_B)
        flathdu_B.writeto(newflat_B, clobber=True)
    
    def runCali(self,dataRlist,dataGlist,dataBlist):
        """
        Run calibration
        
        Inputs
        - dataRlist, dataGlist, dataBlist:
          Lists of R, G, B data images
        """
        self.data_cln_R = cp.batchScienceExposure(dataRlist, self.finalbias, \
                                                  self.finaldark, self.finalflat_R, \
                                                  outDir=self.datadir+"R/")
        self.data_cln_G = cp.batchScienceExposure(dataGlist, self.finalbias, \
                                                  self.finaldark, self.finalflat_G, \
                                                  outDir=self.datadir+"G/")
        self.data_cln_B = cp.batchScienceExposure(dataBlist, self.finalbias, \
                                                  self.finaldark, self.finalflat_B, \
                                                  outDir=self.datadir+"B/")
        
    def coaddIndiv(self):
        pass
        
    def coaddWhole(self):
        """
        Coaddition, 30 files as a group.
        """
        print "Will use 'R' to determine reference catalog"
        
        self.data_coadd_R = RunCoadd("R", self.datadir, referenceCat=True)
        self.data_coadd_G = RunCoadd("G", self.datadir)
        self.data_coadd_B = RunCoadd("B", self.datadir)
        
    def colorImage(self):
        """ Generate color image of coadded data """
        fin_Rdir = datadir+"{0}/coadd/coadd.fits".format("R")
        fin_Gdir = datadir+"{0}/coadd/coadd.fits".format("G")
        fin_Bdir = datadir+"{0}/coadd/coadd.fits".format("B")
        
        os.system("/afs/ir.stanford.edu/class/physics100/stiff/bin/stiff {0} {1} {2} \
                   -OUTFILE_NAME {3}_colored.tiff".format(fin_Rdir,fin_Gdir,fin_Bdir,self.basename))
    
    #####################
    def dispCali(self,which,low=None,high=None):
        
        try:
            if which=='bias':
                filelist = self.biasfiles
            elif which=='dark':
                filelist = self.darkfiles
            elif which=='flat_R':
                filelist = self.flatfiles_R
            elif which=='flat_G':
                filelist = self.flatfiles_G
            elif which=='flat_B':
                filelist = self.flatfiles_B
                
        except:
            raise NameError('No such calibration file(s)!\n \
                             Use "bias", "dark", "flat_R", "flat_G" or "flat_B"')
        
        plt.figure(figsize=(60,30))
        for i,f in enumerate(filelist):
            temp = pyfits.open(f)
            data = temp[0].data
            plt.subplot(len(filelist),1,i+1)
            plt.imshow(data, origin='lower')
            plt.axis('off'); plt.title(str(i))
            
            if low is not None:
                plt.clim(low,high)
    
    def dispCoadd(self,Rlim,Glim,Blim):
        """ Before and after """
        
        plt.figure(figsize=(48,60))
        plt.subplot(231)
        plt.imshow(self.data_coadd_R,origin='lower',cmap='Reds',norm=LogNorm()); plt.clim(Rlim[0],Rlim[1])
        plt.subplot(232); plt.title('After')
        plt.imshow(self.data_coadd_G,origin='lower',cmap='Greens',norm=LogNorm()); plt.clim(Glim[0],Glim[1])
        plt.subplot(233)
        plt.imshow(self.data_coadd_B,origin='lower',cmap='Blues',norm=LogNorm()); plt.clim(Blim[0],Blim[1])

        ## representative first images
        plt.subplot(234)
        plt.imshow(self.data_cln_R[0],origin='lower',cmap='Reds',norm=LogNorm()); plt.clim(Rlim[0],Rlim[1])
        plt.subplot(235); plt.title('Before')
        plt.imshow(self.data_cln_G[0],origin='lower',cmap='Greens',norm=LogNorm()); plt.clim(Glim[0],Glim[1])
        plt.subplot(236)
        plt.imshow(self.data_cln_B[0],origin='lower',cmap='Blues',norm=LogNorm()); plt.clim(Blim[0],Blim[1])