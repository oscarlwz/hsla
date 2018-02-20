#! /usr/bin/env python
from astropy.io import ascii
import numpy as np
import os

def ban_programs(targets):

    banned_programs_list = ascii.read('/grp/hst/HST_spectro/samples/banned_programs')

    canonical = ascii.read(targets+'.list')

    dirlist = canonical['targname'] # [np.where(canonical['flag'] == 1)]
			# programs shouldl be banned whther or not target is used 

    for dirname in dirlist:

        print "Banning Programs for Target: ", dirname
        if os.path.isdir(dirname):

            os.chdir(dirname)
 
            if os.path.exists("all_exposures.txt"):
                exposure_cat = ascii.read("all_exposures.txt")
                print "read exposures for ", dirname

                for pid_to_ban in banned_programs_list['pid']:
                    mask = np.where(exposure_cat['PropID'] == pid_to_ban)
                    print pid_to_ban, mask
                    exposure_cat['Flag'][mask] = 0

                ascii.write(exposure_cat, "all_exposures.txt")

            os.chdir('..')          # go back to "datapile"
