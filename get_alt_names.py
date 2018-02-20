import os 
from astropy.io.votable import parse, is_votable, validate, parse_single_table, exceptions 
from astropy.table import unique as tabunique
from astropy.table import Table
from astropy.io import ascii
import numpy as np 
import copy 
import targetname_munge as tname

def get_alt_names(catalog, database, outfile): 

    targets = ascii.read(catalog) 
    if ('Ref' in targets.keys()) : del t['Ref']
    if ('Start Time' in targets.keys()) : del t['Start Time']
    if ('Stop Time' in targets.keys()) : del t['Stop Time']
    if ('Release Date' in targets.keys()) : del t['Release Date']
    if ('Preview Name' in targets.keys()) : del t['Preview Name']
    if ('High-Level Science Products' in targets.keys()) : del t['High-Level Science Products']
    del targets['Exp Time','Apertures'] 

    unique_targets = tabunique(targets, keys='Target Name') # now try to cull out those with no aliases 
    unique_targets['AltName'] = 'No match                                          ' 
    unique_targets['AltNameRedshift'] = '  . . .  ' 
    unique_targets['AltNameAngsep'] = '. . .   ' 
    unique_targets['AltNameClass'] = ' . . .  ' 
    unique_targets['AltNameSource'] = '. . . ' 
    unique_targets['AltNameLink'] = ' ' * 500    #### 500 spaces -  this will contain the alternate name as a link to the external database page, scrape_headers will use it 

    for i in np.arange(np.size(unique_targets['RA (J2000)'])): 
	ra = unique_targets['RA (J2000)'][i] 
	dec = unique_targets['Dec (J2000)'][i] 

        if database == 'NED':

            found_a_name = False 
            target_to_search = tname.targetname_munge(str(unique_targets['Target Name'][i])) 			# first, munge the targetname with rules we have implemented in targetname_munge 

            print 
            print 
            print "First we are going to try to get a NED match to "+target_to_search+" based on the coordinates" 

            command1 = 'wget -q -O altnames.xml "http://ned.ipac.caltech.edu/cgi-bin/objsearch?in_csys=Equatorial&in_equinox=J2000.0&lon=' 
            command2 = 'd&lat='
            command3 = 'd&corr_z=1&z_constraint=Unconstrained&z_value1=&z_value2=&z_unit=z&ot_include=ANY&nmp_op=ANY&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=Distance+to+search+center&of=xml_main&radius=1&zv_breaker=30000.0&list_limit=5&img_stamp=YES&search_type=Near+Position+Search"' 
            command = command1 + str(ra) + command2 + str(dec) + command3 
            print 'Searching NED for target: ', unique_targets['Target Name'][i], ra, dec 
	    print command 
            os.system(command) 
            votable = parse_single_table("altnames.xml", pedantic=False)
	    if (np.size(votable.fields) > 2):
                table = votable.to_table() 

                i_close = np.where((table['main_col10'] < 0.05) & (table['main_col5'] != 'AbLS') & (table['main_col5'] != 'UvS')) 

                if (np.size(i_close) > 0): 
                    canonical_name = table['main_col2'][i_close[0]] 
                    ned_class = table['main_col5'][i_close[0]]
                    redshift = table['main_col7'][i_close[0]] 
 	            angsep = table['main_col10'][i_close[0]] 

                    unique_targets['AltName'][i] = canonical_name[0] 
                    unique_targets['AltNameRedshift'][i] = redshift[0] 
                    unique_targets['AltNameClass'][i] = ned_class[0] 

                    unique_targets['AltNameAngsep'][i] = angsep[0] 
                    unique_targets['AltNameSource'][i] = 'NED' 
                    unique_targets['AltNameLink'][i] = '<a href="http://ned.ipac.caltech.edu/cgi-bin/objsearch?in_csys=Equatorial&in_equinox=J2000.0&lon='+str(ra)[0:8]+'d&lat='+str(dec)[0:8]+'d&radius=2.0&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&search_type=Near+Position+Search&z_constraint=Unconstrained&z_value1=&z_value2=&z_unit=z&ot_include=ANY&nmp_op=ANY&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=Distance+to+search+center&of=pre_text&zv_breaker=30000.0&list_limit=5&img_stamp=YES">'+unique_targets['AltName'][i]+'</a>' 
                    found_a_name = True 

                print 
                print 

# can I find a better object labeled QSO, but with a redshift? 
                if ('QSO' in table[i_close]['main_col5']): 
                    print 'I found a QSO that meets the match criterion!' 
                    q = table[i_close]  
                    qso_i_found = q[q['main_col5'] == 'QSO'] 

                    print 'Closest object redshift:  ' , table['main_col7'][i_close[0]][0], unique_targets['AltNameRedshift'][i]
                    print 'Found QSO redshift: ' , qso_i_found['main_col7'][0] 
                    if ('0.0' in unique_targets['AltNameRedshift'][i]): 
                        print 'I would rather switch the match to this QSO here:'  
                        print qso_i_found 

                        canonical_name = qso_i_found['main_col2'][0] 
                        ned_class = qso_i_found['main_col5'][0]
                        redshift = qso_i_found['main_col7'][0] 
                        angsep = qso_i_found['main_col10'][0]  
       
                        print redshift 
    
                        unique_targets['AltName'][i] = canonical_name
                        unique_targets['AltNameRedshift'][i] = redshift
                        unique_targets['AltNameClass'][i] = ned_class
    
                        unique_targets['AltNameAngsep'][i] = angsep
                        unique_targets['AltNameSource'][i] = 'NED' 
                        unique_targets['AltNameLink'][i] = '<a href="http://ned.ipac.caltech.edu/cgi-bin/objsearch?in_csys=Equatorial&in_equinox=J2000.0&lon='+str(ra)[0:8]+'d&lat='+str(dec)[0:8]+'d&radius=2.0&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&search_type=Near+Position+Search&z_constraint=Unconstrained&z_value1=&z_value2=&z_unit=z&ot_include=ANY&nmp_op=ANY&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=Distance+to+search+center&of=pre_text&zv_breaker=30000.0&list_limit=5&img_stamp=YES">'+unique_targets['AltName'][i]+'</a>' 
                        found_a_name = True 
         
                    
                if ('G' in table[i_close]['main_col5']):
                    print 'I found a G that meets the match criterion!' 
                    g = table[i_close]  
                    g_i_found = g[g['main_col5'] == 'G'] 
                    print g_i_found.sort('main_col7')
                    

                    print 'Closest object redshift:  ' , table['main_col7'][i_close[0]][0], unique_targets['AltNameRedshift'][i]
                    print 'Found G redshift: ' , g_i_found['main_col7'][0] 
                    if ('0.0' in unique_targets['AltNameRedshift'][i]): 
                        print 'I would rather switch the match to this G here:'  

                        canonical_name = g_i_found['main_col2'][0] 
                        ned_class = g_i_found['main_col5'][0]
                        redshift = g_i_found['main_col7'][0] 
 	                angsep = g_i_found['main_col10'][0]  
       
                        print redshift 

                        unique_targets['AltName'][i] = canonical_name
                        unique_targets['AltNameRedshift'][i] = redshift
                        unique_targets['AltNameClass'][i] = ned_class

                        unique_targets['AltNameAngsep'][i] = angsep
                        unique_targets['AltNameSource'][i] = 'NED' 
                        unique_targets['AltNameLink'][i] = '<a href="http://ned.ipac.caltech.edu/cgi-bin/objsearch?in_csys=Equatorial&in_equinox=J2000.0&lon='+str(ra)[0:8]+'d&lat='+str(dec)[0:8]+'d&radius=2.0&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&search_type=Near+Position+Search&z_constraint=Unconstrained&z_value1=&z_value2=&z_unit=z&ot_include=ANY&nmp_op=ANY&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=Distance+to+search+center&of=pre_text&zv_breaker=30000.0&list_limit=5&img_stamp=YES">'+unique_targets['AltName'][i]+'</a>' 
                        found_a_name = True 

      

            if found_a_name: print "Found a name match in NED using the coordinates, that name is . . . ", unique_targets['AltName'][i], ' and the redshift is z = ', unique_targets['AltNameRedshift'][i] 

# THIS WAS ALL COMMENTED OUT BY JT ON 011916 BECAUSE THE NECESSARY SHELL SCRIPT ALT_NAME_COUNT.SH HAS BEEN LOST!!!! 
            if not found_a_name: 
                print "Now we are going to try to get a NED match to "+target_to_search+" based on the identifier" 

                command1 = 'wget -q -O altnames.xml --no-check-certificate "https://ned.ipac.caltech.edu/cgi-bin/objsearch?objname='
                command2 = '&extend=no&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=RA+or+Longitude&of=xml_main&zv_breaker=30000.0&list_limit=5&img_stamp=YES"' 
                command = command1 + target_to_search + command2 
                print command 
                os.system(command) 
#                os.system("alt_name_count.sh") 							# this is shameless hackery to parse the XML file in the shell but I have no alternative - VO tables are too hard to deal with 
#                counts = ascii.read('altnames_object_count') 
#
#                if counts['number'][0] > 0: 							# if ths XML is OK, use it 
#                    votable = parse_single_table("altnames.xml",pedantic=False).to_table()
#                    canonical_name = table['main_col2'][0] 
#                    ned_class = table['main_col5'][0] 
#                    redshift = table['main_col7'][0]
# 	            angsep = table['main_col10'][0] 
#
#                    unique_targets['AltName'][i] = canonical_name
#                    unique_targets['AltNameRedshift'][i] = redshift
#                    unique_targets['AltNameClass'][i] = ned_class
#                    unique_targets['AltNameAngsep'][i] = angsep
#                    unique_targets['AltNameSource'][i] = 'NED' 
#                    unique_targets['AltNameLink'][i] = '<a href="https://ned.ipac.caltech.edu/cgi-bin/objsearch?objname='+unique_targets['AltName'][i]+'&extend=no&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=RA+or+Longitude&of=pre_text&zv_breaker=30000.0&list_limit=5&img_stamp=YES">'+unique_targets['AltName'][i]+'</a>' 

#                    found_a_name = True 
#                    if found_a_name: print "Found a name match in NED using the identifier, that name is . . . ", unique_targets['AltName'][i]  

            if not found_a_name: print "No object found by NED, sorry!" 
            print 
            print 

        elif database == 'SIMBAD': 
 
            found_a_name = False 
            target_to_search = tname.targetname_munge(str(unique_targets['Target Name'][i])) 			# first, munge the targetname with rules we have implemented in targetname_munge 
            target_name_string = copy.deepcopy(target_to_search) 
            if '+' in target_to_search: target_name_string = str.split(target_to_search,'+')[0] +'%2B'+ str.split(target_to_search,'+')[1]
  
            print 
            print 
            print "First we are going to try to get a match to "+target_to_search+" based on the identifier" 

            command1 = 'wget -q -O altnames.xml "http://simbad.u-strasbg.fr/simbad/sim-id?Ident='
            command2 = '&NbIdent=1&Radius=2&Radius.unit=arcmin&output.format=VOTable&submit=submit+id"' 
            command = command1 + target_name_string + command2 
            os.system(command)

	    votable = parse("altnames.xml")  # otabin the output of SIMBAD as a VOTABLE
            if np.size(votable.resources) > 0: # if the number of resources is > 0, use it
                votable = parse_single_table("altnames.xml",pedantic=False).to_table()
    
                if target_to_search in votable['TYPED_ID'][0]: 
                    canonical_name = votable['MAIN_ID'][0] 

                    simbad_class = votable['OTYPE_S'][0] 

                    redshift = votable['Z_redshift'][0] 

                    angsep = votable['ANG_DIST'][0] 

                    unique_targets['AltName'][i] = canonical_name
                    unique_targets['AltNameRedshift'][i] = redshift
                    unique_targets['AltNameClass'][i] = simbad_class
                    unique_targets['AltNameAngsep'][i] = angsep
                    unique_targets['AltNameSource'][i] = 'SIMBAD' 
                    #unique_targets['AltNameLink'][i] = '<a href="http://simbad.u-strasbg.fr/simbad/sim-coo?CooDefinedFrames=none&CooEpoch=2000&Coord='+str(ra)[0:8]+'d'+str(dec)[0:8]+'d&submit=submit%20query&Radius.unit=arcsec&CooEqui=2000&CooFrame=FK5&Radius=4">'+unique_targets['AltName'][i]+' </a>'  
                    unique_targets['AltNameLink'][i] = '<a href="http://simbad.u-strasbg.fr/simbad/sim-id?Ident='+target_name_string+'&NbIdent=1&Radius=2&Radius.unit=arcmin&submit=submit+id">'+unique_targets['AltName'][i]+' </a>'  
                    found_a_name = True 
 
            if found_a_name: print "Found a name match in SIMBAD using the identifier, that name is . . . ", unique_targets['AltName'][i]  



            if not found_a_name: 
                print "Now we are going to try to get a match to "+str(unique_targets['Target Name'][i])+" based on the coordinates. " 

                command1 = 'wget -q -O altnames.xml "http://simbad.u-strasbg.fr/simbad/sim-coo?Coord='
                command2 = 'd+'
                command3 = 'd&CooFrame=FK5&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=none&Radius=2&Radius.unit=arcmin&output.format=VOTable&submit=submit+query"'
                command = command1 + str(unique_targets['RA (J2000)'][i]) + command2 + str(dec) + command3
                os.system(command) 
                votable = parse("altnames.xml")
                if np.size(votable.resources) > 0: # if the number of resources is > 0, use it
                    votable = parse_single_table("altnames.xml",pedantic=False).to_table()
    
                    i_close = np.where(votable['ANG_DIST'] < 3.0)
    
                    if np.size(i_close) > 0:
                        canonical_name = votable['MAIN_ID'][i_close[0]]
    
                        simbad_class = votable['OTYPE_S'][i_close[0]]
    
                        redshift = votable['Z_redshift'][i_close[0]]
    
                        angsep = votable['ANG_DIST'][i_close[0]]
    
                        unique_targets['AltName'][i] = canonical_name[0]
                        unique_targets['AltNameRedshift'][i] = redshift[0]
                        unique_targets['AltNameClass'][i] = simbad_class[0]
                        unique_targets['AltNameAngsep'][i] = angsep[0]
                        unique_targets['AltNameSource'][i] = 'SIMBAD' 
                        unique_targets['AltNameLink'][i] = '<a href="http://simbad.u-strasbg.fr/simbad/sim-coo?CooDefinedFrames=none&CooEpoch=2000&Coord='+str(ra)[0:8]+'d'+str(dec)[0:8]+'d&submit=submit%20query&Radius.unit=arcsec&CooEqui=2000&CooFrame=FK5&Radius=4">'+unique_targets['AltName'][i]+' </a>'  

            print found_a_name, unique_targets['Target Name'][i], unique_targets['AltName'][i]  
    
            if not found_a_name: 
                print "No object found by SIMBAD, sorry!"
            print 
            print 

        else: 
            print 'cannot do external database '+database+' yet!!!!!' 

    #print unique_targets[i] 
    #os.system('rm altnames.xml') 

    unique_targets.write(outfile+'_'+database+'.html', format='jsviewer') 
    unique_targets.write(outfile+'_'+database+'.fits', format='fits', overwrite=True) 

    return unique_targets 

 
