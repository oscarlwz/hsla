import os
from astropy.io.votable import parse_single_table
from astropy.io.votable import parse
from astropy.table import unique as tabunique
from astropy.io import ascii
import numpy as np

def get_simbad_data(catalog):

    #### THIS SHOWS BASIC USAGE FOR SIMBAD TO GET OUT THE CATALOG DATA    
    #### ON A NEAR-POSITION SEARCH THAT OUTPUTS AND THEN PARSES AN XML
    #### VOTABLE (!)

    targets = ascii.read(catalog)
    if ('Ref' in targets.keys()) : del t['Ref']
    if ('Start Time' in targets.keys()) : del t['Start Time']
    if ('Stop Time' in targets.keys()) : del t['Stop Time']
    if ('Release Date' in targets.keys()) : del t['Release Date']
    if ('Preview Name' in targets.keys()) : del t['Preview Name']
    if ('High-Level Science Products' in targets.keys()) : del t['High-Level Science Products']
    if ('Exp Time' in targets.keys()): del targets['Exp Time'] 
    if ('Apertures' in targets.keys()): del targets['Apertures'] 

    unique_targets = tabunique(targets, keys='Target Name') # now cull out those with no aliases
    unique_targets['SIMBAD Name'] = 'No match found in SIMBAD        '
    unique_targets['Redshift'] = '  . . .  '
    unique_targets['SIMBAD Angsep'] = '  . . . '
    unique_targets['SIMBAD Class'] = ' . . .           '

    for i in np.arange(np.size(unique_targets['RA (J2000)'])):

        found_a_name = False 

        print "First we are going to try to get a match to "+str(unique_targets['Target Name'][i])+" based on the identifier" 

        command1 = 'wget -q -O simbad.xml "http://simbad.u-strasbg.fr/simbad/sim-id?Ident='
        command2 = '&NbIdent=1&Radius=2&Radius.unit=arcmin&output.format=VOTable&submit=submit+id"' 
        command = command1 + str(unique_targets['Target Name'][i]) + command2 
        print str(unique_targets['Target Name'][i]), command 
        os.system(command)

        votable = parse("simbad.xml")  # otabin the output of SIMBAD as a VOTABLE
        if np.size(votable.resources) > 0: # if the number of resources is > 0, use it
            votable = parse_single_table("simbad.xml").to_table()
    
            print votable 

            i_close = np.where(votable['TYPED_ID'] == str(unique_targets['Target Name'][0]))
    
            if np.size(i_close) > 0:
                canonical_name = votable['MAIN_ID'][i_close[0]]

                ned_class = votable['OTYPE_S'][i_close[0]]

                redshift = votable['Z_redshift'][i_close[0]]

                angsep = votable['ANG_DIST'][i_close[0]]

                unique_targets['SIMBAD Name'][i] = canonical_name[0]
                unique_targets['Redshift'][i] = redshift[0]
                unique_targets['SIMBAD Class'][i] = ned_class[0]
                unique_targets['SIMBAD Angsep'][i] = angsep[0]

                found_a_name = True 

        if not found_a_name: 
            print "Now we are going to try to get a match to "+str(unique_targets['Target Name'][i])+" based on the coordinates. " 
            command1 = 'wget -q -O simbad.xml "http://simbad.u-strasbg.fr/simbad/sim-coo?Coord='
            command2 = 'd+'
            command3 = 'd&CooFrame=FK5&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=none&Radius=2&Radius.unit=arcmin&output.format=VOTable&submit=submit+query"'
    
            command = command1 + str(unique_targets['RA (J2000)'][i]) + \
                        command2 + str(unique_targets['Dec (J2000)'][i]) + command3
            #print unique_targets['Target Name'][i], unique_targets['RA (J2000)'][i], unique_targets['Dec (J2000)'][i] 
    
            os.system(command)
    
            votable = parse("simbad.xml")  # otabin the output of SIMBAD as a VOTABLE
            if np.size(votable.resources) > 0: # if the number of resources is > 0, use it
                votable = parse_single_table("simbad.xml").to_table()
        
                #print votable 
    
                i_close = np.where(votable['ANG_DIST'] < 3.0)
    
                if np.size(i_close) > 0:
                    canonical_name = votable['MAIN_ID'][i_close[0]]
    
                    ned_class = votable['OTYPE_S'][i_close[0]]
    
                    redshift = votable['Z_redshift'][i_close[0]]
    
                    angsep = votable['ANG_DIST'][i_close[0]]
    
                    unique_targets['SIMBAD Name'][i] = canonical_name[0]
                    unique_targets['Redshift'][i] = redshift[0]
                    unique_targets['SIMBAD Class'][i] = ned_class[0]
                    unique_targets['SIMBAD Angsep'][i] = angsep[0]
 
                found_a_name = True 

        print found_a_name, unique_targets['Target Name'][i], unique_targets['SIMBAD Name'][i]  

        if not found_a_name: 
            print "No object found by SIMBAD, sorry!"

    unique_targets.write('unique_targets_simbad_test.html', format='jsviewer')
    unique_targets.write('unique_targets_simbad_test.fits', format='fits', overwrite=True)

    return unique_targets
