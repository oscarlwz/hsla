#! /usr/bin/env python

"""
 this code will open the datapile catalog (cos_exposures_Oct26_2015.txt)
 and sort the fits files in a particular datapile into 
 directories for each target
"""

import os
import sys 
from astropy.io import ascii 
from astropy.io import fits
from astropy.table import Table 
import glob
import argparse
#

def sort_targets(catalog): 

    #### first thing we do is open the MAST-provided data catalog
    targets = ascii.read(catalog) 

    #del targets['Ref', 'Start Time', 'Stop Time', 'Apertures', 'Release Date', 'Preview Name', 'High-Level Science Products'] 

    #### then create target name directories if they don't exist 
    for target in targets['Target Name']: 
        if not (os.path.exists(target)): 
            os.system('mkdir ./'+ target) 
            print 'SORT_TARGETS: created directory for target: ', target

    #### now get list of all x1d files in this directory and move them to the appropriate directories 
    filelist = glob.glob(os.path.join('.', '*x1d.fits'))

    for thisfile in filelist: 
        targ = fits.getheader(thisfile, 0)['TARGNAME']
        print(thisfile, targ) 
        move = 'mv '+thisfile+'  ./'+ targ
        print(move)
        ret = os.system(move)
        print(ret)
       
    return "targets sorted!"


#-----------------------------------------------------------------------------------------------------

 
def parse_args():
    '''
    Parse command line arguments.  Returns args object.
    '''
    parser = argparse.ArgumentParser(description="sorts targets into directories")
    parser.add_argument('catalog', metavar='catalog', type=str, action='store',
                        help='datapile catalog grabbed from MAST query, csv')

    args = parser.parse_args()
    return args
       

#-----------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    args = parse_args()

    message = sort_targets(args.catalog)
    sys.exit(message)

     
 
