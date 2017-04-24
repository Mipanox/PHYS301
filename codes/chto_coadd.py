"""
Improving the original coaddition script.

Fixing the deprecated syntax in `f2n.py`
and getting around the number limit of images.
Currently set to "30", and then the 'coadded'
groups are coadded thereafter.

Adding rejection based on unphysical alignment.

(Created and maintained by Chun-Hao.
 Modified by Jason.
 Last change: 04/23/17)
"""

import sys, os, shutil, math
import numpy as np
import pyfits
from numpy import *
from pyraf import iraf
from glob import glob
from modules.variousfct import *
from modules.star import *
from datetime import datetime, timedelta
import imp
f2n = imp.load_source('f2n', '/afs/ir.stanford.edu/class/physics100/workdir/g2/f2n.py')

##################################################

def visualizeAlignmentStars(workdir):
    '''
    Creates an image with the alignment reference stars highlighted,
    for inspection purposes.

    @parameter workdir : directory where workfiles are stored for
                          alignment
    '''

    print "Reading database..."
    db = readpickle(os.path.join(workdir, "db.pkl"))
    print "Ok, we have %i images." % len(db)
    print "Reading alignment star catalog..."
    alistars = readpickle(os.path.join(workdir, "alistars.pkl"))

    alistarsdicts = [{"name":star.name, "x":star.x, "y":star.y} for star in alistars]

    # Get the ref image (similar to previous script)
    print "Automatic reference image selection :"
    sorteddb = sorted(db, key=lambda k: k['nbrsources'])
    refimg = sorteddb[-1]
    refimgname = refimg['rawimgname']
    print refimgname

    pnginfostring = "Alignment star selection (%i stars)" % len(alistarsdicts)

    # We are ready to proceed with the png creation

    f2nimage = f2n.fromfits(refimg['rawimgpath'], hdu=0, verbose=True)
    f2nimage.setzscale("auto", "auto")
    f2nimage.rebin(3)


    f2nimage.makepilimage(scale="log", negative=False)
    f2nimage.drawstarslist(alistarsdicts, r=20)
    f2nimage.writeinfo([pnginfostring])
    f2nimage.writetitle(os.path.basename("Ref : " + refimgname))

    pngfilepath = os.path.join(workdir, "alistars_refimage.png")
    f2nimage.tonet(pngfilepath)

    print "Done."
    print "Now have a look at :"
    print pngfilepath
    print "... and decide if you are happy with this selection."


def coaddImages(workdir, coaddfile = None, method='median'):
    '''
    Combine aligned images, masking out bad pixels using transformed
    masks

    @parameter coaddfile path to output coadded file to be produced by
    this method. None defaults to the workdir.
    
    @parameter method Method for combining images. See IRAF's imcombine
    for additional options
    '''
    
    print "Reading database..."
    db = readpickle(os.path.join(workdir, "db.pkl"))
    rejectimages=[]
    for image in db:
        ## reject images with unphysical transformations.
        ##  e.g., those rotated too much
        if image['geomapangle'] > 15. and image['geomapangle'] < 345:
            image['okforalignment'] = False
            rejectimages.append(image['rawimgname']) 
            
    db = [image for image in db if image['okforalignment'] == True]
    if len(db) == 0:
        print "Cannot continue, no matching images!"
        return

    print "Ok, we have %i images to combine." % len(db)

    curdir = os.getcwd()
    os.chdir(workdir)

    if coaddfile is None:
        coaddfile = 'coadd.fits'

    starttime = datetime.now()

    iraf.images()
    iraf.immatch()
    inputimages=[]
    outnames=[]
    
    ## sub-grouping the files
    numberofGroups=int(len(db)/30.)+1
    for i in range(numberofGroups-1):
        inputimages.append(','.join([image['aliimgname'] for image in db[i*30:(i+1)*30]]))
        outnames.append(coaddfile+'_{0}.fits'.format(i))
        
    if len(db)%30!=0:
        ## handling the remaining not-grouped files
        inputimages.append(','.join([image['aliimgname'] for image in db[(numberofGroups-1)*30:]]))
        outnames.append(coaddfile+'_{0}.fits'.format(numberofGroups-1))
    print "divied the images in {} groups".format(numberofGroups)
    
    for i in range(len(inputimages)):
        iraf.flprcache()
        iraf.imcombine(input = inputimages[i], output = outnames[i],
                       combine=method, masktype='goodvalue')
    iraf.flprcache()
    iraf.imcombine(input = outnames[i], output = coaddfile,
                   combine=method, masktype='goodvalue')
    
    print "- " * 40
    endtime = datetime.now()
    timetaken = nicetimediff(endtime - starttime)

    os.chdir(curdir)

    print "Dear user, I'm done with the alignment. I did it in %s." % timetaken
    print "The following are rejected: \n{0}".format(rejectimages)
    ## Update the db :
    writepickle(db, os.path.join(workdir, "db.pkl"))