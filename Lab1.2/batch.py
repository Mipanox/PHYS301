"""
Run the coaddition script.

(Issue: I cannot use the `visualizeAlignmentStars` method
        in the `coaddition` module.
        Have to invoke the modified one.)

Writen and modified by Chun-Hao and Jason (04/18/17)
"""
import sys
sys.path.append("/afs/ir.stanford.edu/class/physics100/workdir/g2/")

import glob
import chto_coadd

def RunCoadd(filterName,datadir):
    processedImage = [file for file in glob.glob(datadir+"{0}/*".format(filterName)) if "Processed" in file]
    
    ## working directory
    coaddDir = datadir+"{0}/coadd/".format(filterName)
    RefDir   = datadir+"reference.cat"
    
    ## coaddition
    coaddition.setupWorkdir(processedImage,coaddDir)
    coaddition.detectStars(coaddDir)
    coaddition.measureSeeing(coaddDir, showseeingplots=True)
    coaddition.createRefCatalog(coaddDir, outputcatpath = RefDir,
                                handpickedcat = None,
                                minfwhm = 1.0, maxfwhm = 8.5, maxell = 1.25)
    chto_coadd.visualizeAlignmentStars(coaddDir)
    coaddition.matchAlignmentStars(coaddDir, refcat =RefDir , tolerance = 5.0,
                        minnbrstars = 4, maxdist = 15)
    coaddition.createTransforms(coaddDir)
    coaddition.transformImages(coaddDir)
    coaddition.transformMasks(coaddDir, maskext='.mask')
    coaddition.coaddImages(coaddDir, coaddfile = None, method='median')
    
    product  = pyfits.open(coaddDir+"coadd.fits")[0].data
    return product