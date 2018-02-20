#! /usr/bin/env python

from astropy.io import ascii
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import os
import multiprocessing as mp
import argparse

import quick_look
import add
import plot_interactive 
import sys


def parse_args():
    '''
    Parse command line arguments.  Returns args object.
    '''
    parser = argparse.ArgumentParser(description="makes co-adds and quicklooks for everything in 'targets.list' file")
    parser.add_argument('targets', metavar='list', type=str, action='store',
                        help="""targets.list is the file to be read in;
                              first column = flag (0,1) if target is to be used,
                              second column = target/directory name""")

    parser.add_argument('--clobber', dest='clobber', action='store_true')
    parser.add_argument('--no-clobber', dest='clobber', action='store_false')
    parser.set_defaults(clobber=False)

    args = parser.parse_args()
    return args

#-----------------------------------------------------------------------------------------------------

def drive_targets(targets):
    clobber = targets[1]
    canonical = ascii.read(targets[0]+'.list')
    print canonical
    dirlist = canonical['targname'][np.where(canonical['flag'] == 1)]
    print "dirlist in driver"
    print dirlist
    print 'clobber ', clobber
    dirc = []
    for d in dirlist:
        dirc.append((d,clobber))

    # put this back in to get parallel performance 
    mp_drive_dirlist(dirc)

#-----------------------------------------------------------------------------------------------------

def drive_dirlist(dirc):
    dirname = dirc[0]
    clobber = dirc[1]
    print "Driving target:  ", dirname, "with clobber = ",clobber
    if os.path.isdir(dirname):
        os.chdir(dirname)

        if os.path.isfile('DONE_'+dirname): 
            print 'found a DONE file for ', dirname, ' so will skip it' 
        if not os.path.isfile('DONE_'+dirname): 
 
            print 'did not file a DONE file for ', dirname 

            os.system('touch WORKING_'+dirname)
            os.system('rm DONE_'+dirname)
            # PUT ANYTHING THAT WILL BE DONE IN TARGET DIRECTORY HERE
            filelist = os.listdir('.')
                #if(clobber != 1 and os.path.exists(dirname+'_FUV_M_coadd.dat') or os.path.exists(dirname+'_FUV_L_coadd.dat')):
                #    print dirname+ ":  Coadd already exists, skipping  "
                #else:
                #    print dirname+ ":  Running coadd_x1d"
                #    os.system('/Applications/exelis/idl82/bin/idl -e "@../../code/coadd_script.pro"')
    
            if clobber != 1 and (os.path.exists(dirname+'_coadd_G130M_final_all.fits.gz') or os.path.exists(dirname+'_coadd_G140L_final_all.fits.gz')):
                print dirname+ ":  JT Coadd already exists, skipping  "
            else:
                print dirname+ ":  Running add"
                coadd = add.main('FUVM')
                coadd = add.main('FUVL')
    
            if clobber != 1 and os.path.exists(dirname+'_quicklook.html'):
                print dirname+ ":  Quicklook already exists, skipping  "
            else:
                print dirname+":  Creating Quick Look for ", dirname
                quick_look.get_quick_look()
        
            tarstring = 'tar -cvf '+dirname+'_target.tar *png *fits.gz *txt *html' 
            print tarstring 
            os.system(tarstring) 
            gzipstring = 'gzip -f *_target*tar' 
            print gzipstring 
            os.system(gzipstring) 
    
            if clobber != 1 and os.path.exists(dirname+'_quicklook.html'):
                print dirname+ ":  Quicklook already exists, skipping  "
            else:
                print dirname+":  Creating Quick Look for ", dirname
                quick_look.get_quick_look()
    
    
            os.system('touch DONE_'+dirname)
            os.system('rm WORKING_'+dirname)

        os.chdir('..')  		# go back to "datapile"


#-----------------------------------------------------------------------------------------------------

def mp_drive_targets(targets):
    pool = mp.Pool(processes=10)
    pool.map(drive_targets, targets)

#-----------------------------------------------------------------------------------------------------

def mp_drive_dirlist(dirc):
    pool = mp.Pool(processes=10)
    pool.map(drive_dirlist, dirc)


#-----------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    args = parse_args()
    targets = (args.targets, args.clobber)
    
    drive_targets(targets)
    sys.exit("""
    
    ~~~~~~~*~*~*~*~
    ~~~~~~~*~*~*~*~  all done!!!! spectra are fun!
    ~~~~~~~*~*~*~*~
    """)
