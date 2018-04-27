#! /usr/bin/env python
from astropy.io import ascii
import numpy as np
import os

def ban_programs(targets):

    #JRD changed path
    banned_programs_list = ascii.read('banned_programs')

    canonical = ascii.read(targets+'.list')

    dirlist = canonical['targname'] # [np.where(canonical['flag'] == 1)]
			# programs shouldl be banned whther or not target is used 

    for dirname in dirlist:

        #JRD python 3
        print("Banning Programs for Target: ", dirname)
        if os.path.isdir(dirname):

            os.chdir(dirname)
 
            if os.path.exists("all_exposures.txt"):
                exposure_cat = ascii.read("all_exposures.txt")
                #JRD python3
                print("read exposures for ", dirname)

                for pid_to_ban in banned_programs_list['pid']:
                    mask = np.where(exposure_cat['PropID'] == pid_to_ban)
                    #JRD python 3
                    print(pid_to_ban, mask)
                    exposure_cat['Flag'][mask] = 0

                #JRD Deprecated automatic writing. Needs overwrite=True
                ascii.write(exposure_cat, "all_exposures.txt", overwrite=True)

            os.chdir('..')          # go back to "datapile"
