#! /usr/bin/env python

from astropy.io import ascii 
from astropy.table import unique as tabunique
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
import copy 


def target_alias(catalog): 

    t = ascii.read(catalog) 
    if ('Status' in t.keys()) : del t['Status'] 
    if ('Dataset' in t.keys()) : del t['Dataset'] 
    if ('Apertures' in t.keys()) : del t['Apertures'] 

    rasort = t.argsort(keys='RA (J2000)') 
    hold = t[rasort] 
    unique_targs = tabunique(hold, keys='Target Name')

    targets = unique_targs
    targets.write('unique_target_names.html', format='jsviewer') 
    targets.write('unique_target_names.list', format='ascii') 

    if ('Target Descrip' in targets.keys()) : del targets['Target Descrip'] 
    if ('Broad Category' in targets.keys()) : del targets['Broad Category'] 

    aliases = targets
    aliases['alias0'] = '          . . .                   ' 
    aliases['alias1'] = '          . . .                   ' 
    aliases['alias2'] = '          . . .                   ' 
    aliases['alias3'] = '          . . .                   ' 
    aliases['alias4'] = '          . . .                   ' 
    aliases['alias5'] = '          . . .                   ' 
    aliases['alias6'] = '          . . .                   ' 
    aliases['alias7'] = '          . . .                   ' 
    aliases['alias8'] = '          . . .                   ' 
    aliases['alias9'] = '          . . .                   ' 
    aliases['alias10']= '          . . .                   ' 

    c = SkyCoord(ra=targets['RA (J2000)']*u.degree, dec=targets['Dec (J2000)']*u.degree)

    for index in np.arange(np.size(targets['RA (J2000)'])): 
        sep = c.separation(c[index]) 
        i_where = np.where(sep.arcsecond <= 3) 
        unique_by_name = tabunique(targets[:][i_where], keys='Target Name')
        print 
        print index, targets['Target Name'][index]  
        print unique_by_name['Target Name'] 
 	aliases['alias0'][index] = unique_by_name['Target Name'][0] 

        #print index, targets['Target Name'][index], targets['Target Name'][i_where] 
        unique_by_name = tabunique(targets[:][i_where], keys='Target Name')
        print index, targets['Target Name'][index],unique_by_name['Target Name'] 
        aliases['alias0'][index] = unique_by_name['Target Name'][0] 
        if (np.size(unique_by_name['Target Name']) > 1): aliases['alias1'][index] = unique_by_name['Target Name'][1] 
        if (np.size(unique_by_name['Target Name']) > 2): aliases['alias2'][index] = unique_by_name['Target Name'][2] 
        if (np.size(unique_by_name['Target Name']) > 3): aliases['alias3'][index] = unique_by_name['Target Name'][3] 
        if (np.size(unique_by_name['Target Name']) > 4): aliases['alias4'][index] = unique_by_name['Target Name'][4] 
        if (np.size(unique_by_name['Target Name']) > 5): aliases['alias5'][index] = unique_by_name['Target Name'][5] 
        if (np.size(unique_by_name['Target Name']) > 6): aliases['alias6'][index] = unique_by_name['Target Name'][6] 
        if (np.size(unique_by_name['Target Name']) > 7): aliases['alias7'][index] = unique_by_name['Target Name'][7] 
        if (np.size(unique_by_name['Target Name']) > 8): aliases['alias8'][index] = unique_by_name['Target Name'][8] 
        if (np.size(unique_by_name['Target Name']) > 9): aliases['alias9'][index] = unique_by_name['Target Name'][9] 
        if (np.size(unique_by_name['Target Name']) > 10): aliases['alias10'][index] = unique_by_name['Target Name'][10] 
    
    del aliases['Exp Time', 'Instrument', 'Filters/Gratings', 'Central Wavelength', 'Proposal ID'] 

    unique_names = tabunique(aliases, keys='alias1') # now try to cull out those with no aliases 
    del unique_names['alias0'] 

    rasort = aliases.argsort('RA (J2000)') 
    aliases = aliases[rasort] 
    del aliases['RA (J2000)', 'Dec (J2000)', 'alias6', 'alias7', 'alias8', 'alias9', 'alias10'] 
    aliases.write(catalog[0:-4]+'.alias', format='ascii.fixed_width_two_line') 

    return aliases 
