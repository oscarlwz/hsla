#! /usr/bin/env python
from astropy.io import ascii
from astropy.io import fits
import glob
import numpy as np
import os

def ban_programs(targets):

    banned_programs_list = ascii.read('/grp/hst/HST_spectro/hsla_releases/code/banned_programs')
    banned_visits = ascii.read(glob.glob('banned_visits*.txt')[0]) # works for both COS and STIS datapiles
    canonical = ascii.read(targets+'.list')
    dirlist = canonical['targname']

    for n,dirname in enumerate(dirlist):
        print("Banning Programs for Target: ", dirname)
        if os.path.isdir(dirname):
            os.chdir(dirname)
 
            if os.path.exists("all_exposures.txt"):
                exposure_cat = ascii.read("all_exposures.txt")
                print("read exposures for ", dirname)
                
                # Ban bad proposals
                for pid_to_ban in banned_programs_list['pid']:
                    mask = np.where(exposure_cat['PropID'] == pid_to_ban)
                    exposure_cat['Flag'][mask] = 0

                # Ban failed visits and bad exposures
                for i,root in enumerate(exposure_cat['Rootname']):
                    # Ban failed visits
                    select = ((banned_visits['Proposal'] == exposure_cat['PropID'][i]) & 
                              (banned_visits['Visit'] == root[4:6].upper()))
                    if len(banned_visits[select]) != 0:
                        exposure_cat['Flag'][i] = 0  
                    # Ban bad exposures
                    fname = root + '_x1d.fits'
                    if fits.getheader(fname,1)['EXPTIME'] / fits.getheader(fname,1)['PLANTIME'] == 0.0:
                        exposure_cat['Flag'][i] = 0

                ascii.write(exposure_cat, "all_exposures.txt", overwrite=True)

            os.chdir('..')          # go back to "datapile"

        if n%100==0:
            print('{}/{} targets completed.'.format(n, len(dirlist)))
            