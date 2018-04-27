
#### this code will open the datapile catalog (cos_exposures_Oct26_2015.txt) 
#### and sort the fits files in a particular datapile into the 

import os
from astropy.io import ascii 
from astropy.table import Table 

def defeat_aliases(aliasfile): 

    #### in the alias file, the first column is the default MAST name for each target. 
    #### the second column is the name we WANT to use, so we'll move stuff from the 
    #### directory for the first name with the directory with the second name. 
    #### this is done non-destructively, so the original directory stays in place 
    #### even though it should remain empty. 

    aliases = ascii.read(aliasfile) 

    i = 0 

    for original_targname in aliases['Target Name']:

        if aliases['alias0'][i] == original_targname:  
            print("The target name and alias are identical, so I'll do nothing!")
        else: 
            move_command = 'mv ./'+original_targname+'/*.fits.gz '+'  ./'+aliases['alias0'][i] 
            print(move_command)
            os.system(move_command) 

        i += 1 

