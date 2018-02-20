from astropy.io import fits
from astropy.table import Table, Column, vstack 
from astropy.io import ascii
import numpy as np
import glob
import os 
import sys  
import copy 

def check_db(): 

    canonical_files = ['EXOPLANETS', 'STAR_EARLY', 'STAR_LATE', 'STAR_NOVAE', 'STAR_TTAURI', 'STAR_POST_AGB', 'STAR_WD', 'STAR_LMXB','STAR_OTHER',\
	'SOLAR_SYSTEM', 'QSOALS', \
         'GALAXY_STARBURST','GALAXY_SPIRAL','GALAXY_STARFORMING', 'GALAXY_DWARF_COMPACT','GALAXY_EMISSION_LINE', 'GALAXY_IRREGULAR', 'GALAXY_CLUSTER',\
        'SUPERNOVAE'] 

    canonical_sum = ascii.read(canonical_files[0]+'.list') 

    for cfile in canonical_files[1:]: 
	print cfile+'.list' 
        c = ascii.read(cfile+'.list') 
        canonical_sum = vstack([canonical_sum, c])  

    all_canonical = ascii.read('/grp/hst/HST_spectro/samples/all_targets.list') 
    copy_canonical = copy.deepcopy(all_canonical) 

    for ii in np.arange(np.size(all_canonical)): 
        if all_canonical['targname'][ii] in canonical_sum['targname']: 
	    print 'yes', all_canonical['targname'][ii] 
            jj = np.where(copy_canonical['targname'] == all_canonical['targname'][ii])[0] 
            copy_canonical.remove_row(jj[0]) 
	else: 
	    print 'no', all_canonical[ii] 

    copy_canonical.write('UNALLOCATED.html', format='jsviewer') 
    copy_canonical.write('UNALLOCATED.temp',format='ascii.fixed_width', delimiter=',') 

    os.system("cat UNALLOCATED.temp | awk '{print $2, $4}' > UNALLOCATED.list") 
    os.system('rm UNALLOCATED.temp')

    return copy_canonical 
