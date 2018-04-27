#! /usr/bin/env python

# JRD Astropy
#import fitsio
from astropy.io import fits
from astropy.table import Table, vstack
from astropy.io import ascii
import numpy as np
import argparse
import glob
import os 
import sys  

def parse_args():
    '''
    Parse command line arguments.  Returns args object.
    '''
    parser = argparse.ArgumentParser(description="scrapes headers for everything in 'filename.list' file")
    parser.add_argument('filename', metavar='filename', type=str, action='store',
                        help='filename.list is the file to be read in')
    parser.add_argument('--altnames', metavar='altnames', type=str, default='blank', action='store',
                        help='altnames is the file containing the NED names list')
    parser.add_argument('--redshifts', metavar='redshifts', type=bool, default=False, action='store',
                        help='do or do not include redshifts, default if not')
    #JRD adding one input parameter
    parser.add_argument('--dir', metavar='data directory name', type=str, default='blank', action='store',
                        help='name of the directory, default = "datapile"')
    args = parser.parse_args()
    return args

#-----------------------------------------------------------------------------------------------------

def scrape_headers(targets,altnames,redshifts,dir = 'datapile'): 
    canonical_filename = targets
    canonical = ascii.read(canonical_filename+'.list')

    have_alt_names = False 
    if altnames != 'blank':
       alt_namefile = altnames
       #JRD altname_table = fitsio.read(alt_namefile, ext=1, header=True)[0]
       altname_table = fits.open(alt_namefile)[1]
       print('Successfully opened '+alt_namefile+' with the alt names in it.')
       have_alt_names = True 

    #JRD python 3
    print("CANONICAL ", canonical)
    print()
    dirlist = canonical['targname'].data[np.where(canonical['flag'] == 1)]
    print()
    print("dirlist ", dirlist)
    print()

    sample_fitstable = Table(names=('Number','Target Name', 'RA','DEC','Nexp','Target Category', \
        'Target Description','Alt Name','median S/N'), dtype=('i4','S200','f4','f4','i8','S20','S20','S30','f4')) 

    sample_webtable = Table(names=('Number','Target Name', 'RA','DEC','Nexp','Target Category', \
                                   'Target Description','AltName','AltClass', 'Redshift','MAST', 'Median S/N', 'FUV M', 'FUV L', 'NUV M', 'NUV L', 'Download', 'L Download'), \
        dtype=('i4','S200','f4','f4','S240','S20','S20','S600','S30', 'S5','S240','f4','S350','S350','S350', 'S350', 'S350','S350')) 

    targets = Table(names=('Flag','Target Name', 'Target Category', 'Target Description'), dtype=('i4','S200','S25', 'S350')) 

    ## exposures = Table(names=('Rootname','Target Name', 'RA','DEC','PropID','PI Name','Detector','Segment',\
        # 'LP','Grating', 'Cenwave','FPPOS','Exptime','Nevents','Extended','Date','Target Description'),   
       # dtype=('S20','S35','f4','f4','I5','S20','S4','S5','S2','S10','S10','I2','f10','f8','S4','S12','S200'))
    exposures = []
       
    #### set up the master "header table"
    ###### ---->>>>>> why does this have to be done each time ????? <<<<<---------
    # generic = "generic_x1d.fits"
    # hdr0 = fitsio.read_header(generic, 0) 
    # hdr1 = fitsio.read_header(generic, 1)
    # for key in ["HISTORY", "COMMENT", ""]:
    #     del hdr0[key]
    # names = sorted(hdr0.keys())
    # rows = [hdr0[k] for k in names]
    # header_table0 = Table(rows=[rows], names=names)

    # for key in ["HISTORY", "COMMENT", ""]:
    #     del hdr1[key]
    # names = sorted(hdr1.keys())
    # rows = [hdr1[k] for k in names]
    # header_table1 = Table(rows=[rows], names=names)

    redshift_string = ' . . . ' 
    altname_string = ' skldjf' 

    counter = 1 
    total_number_of_headers = 0 

    for dirname in dirlist:
        #JRD python3
        print("Driving target:  ", dirname)
        if (os.path.isdir(dirname)):
            os.chdir(dirname)

            filelist = glob.glob(os.path.join('.', '*x1d.fits.gz'))
            #filelist = filelist + glob.glob(os.path.join('.', '*x1d.fits'))

            nfiles = np.size(filelist)
            #JRD python 3
            print("There are ", nfiles, " exposures for target ", dirname)

            if nfiles > 0:			### if there are no files, then this target was aliased or something else happened and you won't be using it. 
                webtable_row, fitstable_row, targetstable_row = get_webtable_info(filelist[0], nfiles, counter, dir = dir)

                if have_alt_names: 
                    i_alt = np.where(altname_table[:]['Target Name'] == dirname)
                    #print 'Redshift = ', altname_table[int(i_alt[0])][13] # not sure 
                    #altname_string = altname_table[int(i_alt[0])][12] 
                    #if 'No match' in altname_table[int(i_alt[0])][11]: altname_string = ' . . . ' # not sure 
                    #altname_class = altname_table[int(i_alt[0])][15] 
                    #redshift_string = str(altname_table[int(i_alt[0])][13]) 

                    webtable_row[7] = altname_table['AltName'][int(i_alt[0])]
                    if 'No match' in webtable_row[7]: webtable_row[7] = ' . . . ' 
                    webtable_row[8] = altname_table['AltNameClass'][int(i_alt[0])]
                    webtable_row[9] = altname_table['AltNameRedshift'][int(i_alt[0])]


                sample_webtable.add_row(webtable_row) 
                sample_fitstable.add_row(fitstable_row)
                targets.add_row(targetstable_row) 

                dataset_list = glob.glob(os.path.join('.', '*x1d.fits.gz'))
                #JRD python 3
                print("SCRAPE_HEADERS: Making Exposure Catalog: " , filelist)
  
                exposure_cat = make_exposure_catalog(filelist)
                print(len(exposures))
                print("this is exposure cat:")
                print(exposure_cat)
                if len(exposures) == 0:
                    exposures = exposure_cat
                else:
                    exposure_tmp = exposures
                    exposures = vstack([exposure_tmp, exposure_cat])
                    # exposures.add_row(exposure_cat)

                counter = counter + 1 

            os.chdir('..')          # go back to "datapile" 

    if not have_alt_names:  
      if ('AltName' in sample_webtable.keys()): del sample_webtable['AltName'] 
      if ('AltClass' in sample_webtable.keys()): del sample_webtable['AltClass'] 

    if not redshifts:  
      del sample_webtable['Redshift'] 

    sample_fitstable.write(canonical_filename+'_sample.fits', format='fits', overwrite=True) 

    del sample_webtable['Target Category','L Download']				# just omit this for now 
    sample_webtable.write(canonical_filename+'_websample.fits', format='fits', overwrite=True) 
    #sample_webtable.write(canonical_filename+'_sample_webtable.txt' ,format='ascii') 
    sample_webtable.write('sample_webtable.temp', format='jsviewer') 
    os.system('sed "s/&lt;/</g" sample_webtable.temp | sed "s/&gt;/>/g" > '+canonical_filename+'_sample.html') 
    os.system('rm sample_webtable.temp')

    targets.write(canonical_filename+'.info',format='ascii.fixed_width', delimiter=',') 

    exposures.write(canonical_filename+'_exposures.fits', format='fits', overwrite=True) 
    exposures.write(canonical_filename+'_exposures.html',format='jsviewer') 
 
    print()
    print("END OF SCRAPE_HEADERS")
    print()

    # only do this if you have created the master header table in the commented out bits above 
    # header_table0.write(canonical_filename+'_headers0.fits', format='fits', overwrite=True) 
    # header_table0.write(canonical_filename+'_headers0.html', format='jsviewer')
    
    # header_table1.write(canonical_filename+'_headers1.fits', format='fits', overwrite=True) 
    # header_table1.write(canonical_filename+'_headers1.html', format='jsviewer')



#-----------------------------------------------------------------------------------------------------

def get_webtable_info(filename, nfiles, counter, dir = 'datapile'):

    #JRD astropy.io fits
    observation = fits.open(filename)
    hdr0 = observation[0].header
    hdr1 = observation[1].header
    #hdr0 = fitsio.read_header(filename, 0) 
    #hdr1 = fitsio.read_header(filename, 1)

    targname = str(hdr0['TARGNAME']).strip()  
    targdesc = str(hdr0['TARDESCR']).strip()  
    ra = hdr0['RA_TARG']
    dec = hdr0['DEC_TARG']

    median_sn = 0.0     

    #REPLACE WITH CORRECT PATH FOR RELEASE
    #targname_urlstring = '<a href="../datapile/'+targname+'/'+targname+'_quicklook.html">'+targname+'</a>'
    targname_urlstring = '<a href="../' + dir + '/' + targname  + '/' +targname+'_quicklook.html">'+targname+'</a>'

    altname_string = ' . . . ' 
    altname_class = ' . . . ' 
    redshift_string = ' . . . ' 

    mast_string = '<a href="https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html?searchQuery='+str(ra)+','+str(dec)+'"> MAST  </a>'  

    #DUPATE PATH
    #n_exp_string = '<a href="../datapile/'+targname+'/all_exposures.html">'+str(nfiles)+'</a>' 
    n_exp_string = '<a href="../' + dir + '/'+targname+'/all_exposures.html">'+str(nfiles)+'</a>' 

    fuv_m_quicklook_urlstring = ' . . . ' 
    fuv_l_quicklook_urlstring = ' . . . ' 
    nuv_l_quicklook_urlstring = ' . . . ' 
    nuv_m_quicklook_urlstring = ' . . . ' 

    download_string = ' . . . ' 

    G130M_coadd_exists = False 
    G160M_coadd_exists = False 
    G140L_coadd_exists = False
    G185M_coadd_exists = False
    G225M_coadd_exists = False
    G285M_coadd_exists = False
    G230L_coadd_exists = False

    #JRD updating filenames
    if os.path.exists(str(hdr0['targname']).strip()+'_coadd_G130M_final_lpALL.fits.gz'): G130M_coadd_exists = True 
    if os.path.exists(str(hdr0['targname']).strip()+'_coadd_G160M_final_lpALL.fits.gz'): G160M_coadd_exists = True 
    if os.path.exists(str(hdr0['targname']).strip()+'_coadd_G140L_final_lpALL.fits.gz'): G140L_coadd_exists = True 
    if os.path.exists(str(hdr0['targname']).strip()+'_coadd_G185M_final_lp1.fits.gz'): G185M_coadd_exists = True 
    if os.path.exists(str(hdr0['targname']).strip()+'_coadd_G225M_final_lp1.fits.gz'): G225M_coadd_exists = True 
    if os.path.exists(str(hdr0['targname']).strip()+'_coadd_G285M_final_lp1.fits.gz'): G285M_coadd_exists = True 
    if os.path.exists(str(hdr0['targname']).strip()+'_coadd_G230L_final_lp1.fits.gz'): G230L_coadd_exists = True 

    if G130M_coadd_exists or G160M_coadd_exists:
        #DUPATE PATH
        #fuv_m_quicklook_urlstring = '<a href="../datapile/'+targname+'/'+targname+'_coadd_FUVM_final_all.png"><img height="40" src="../datapile/'+targname+'/'+targname+'_coadd_FUVM_final_all.png"></a>'
        fuv_m_quicklook_urlstring = '<a href="../' + dir + '/'+targname+'/'+targname+'_coadd_FUVM_final_all.png"><img height="40" src="../test_nuv_data/'+targname+'/'+targname+'_coadd_FUVM_final_all.png"></a>'

    if (os.path.exists(str(hdr0['targname']).strip()+'_coadd_G140L_final_all.png')):
        #UPDATE PATH
        #fuv_l_quicklook_urlstring = '<a href="../datapile/'+targname+'/'+targname+'_coadd_G140L_final_all.png"><img height="40" src="../datapile/'+targname+'/'+targname+'_coadd_G140L_final_all.png"></a>'
        #l_download_string = '<a href="../datapile/'+targname+'/'+targname+'_coadd_G140L_final_all.fits.gz">ALL</a> |'
        fuv_l_quicklook_urlstring = '<a href="../' + dir + '/'+targname+'/'+targname+'_coadd_G140L_final_all.png"><img height="40" src="../test_nuv_data/'+targname+'/'+targname+'_coadd_G140L_final_all.png"></a>'
        l_download_string = '<a href="../' + dir + '/'+targname+'/'+targname+'_coadd_G140L_final_all.fits.gz">ALL</a> |'

        this_coadd = Table.read(targname+'_coadd_G140L_final_all.fits.gz') 
        i_good = np.where((this_coadd['FLUX'] > 0) & (this_coadd['WAVE'] > 1100) & (this_coadd['WAVE'] < 1900))  ## take regions that are not lousy S/N 
        median_sn_140 = np.median(this_coadd['SN'][i_good]) 

        if (os.path.exists(str(hdr0['targname']).strip()+'_coadd_G140L_final_lp1.fits.gz')):
            print('I found a coadd for G140L LP=1')
            #l_download_string = l_download_string  + \
            #'  '+'<a href="../datapile/'+targname+'/'+targname+'_coadd_G140L_final_lp1.fits.gz">LP1</a> | '
            #DUPATE PATH
            l_download_string = l_download_string  + \
            '  '+'<a href="../' + dir + '/'+targname+'/'+targname+'_coadd_G140L_final_lp1.fits.gz">LP1</a> | '
        else:
            l_download_string = l_download_string+'  ' + \
            '. . . .  | '
        if (os.path.exists(str(hdr0['targname']).strip()+'_coadd_G140L_final_lp2.fits.gz')):
            print('I found a coadd for G140L LP=2')
            #DUPATE PATH
            #l_download_string = l_download_string + \
            #'  '+'<a href="../datapile/'+targname+'/'+targname+'_coadd_G140L_final_lp2.fits.gz">LP2</a> | '
            l_download_string = l_download_string + \
            '  '+'<a href="../' + dir + '/'+targname+'/'+targname+'_coadd_G140L_final_lp2.fits.gz">LP2</a> | '
        else:
            l_download_string = l_download_string + \
            '  '+'. . . .  | '
        if (os.path.exists(str(hdr0['targname']).strip()+'_coadd_G140L_final_lp3.fits.gz')):
            print('I found a coadd for G140L LP=3')
            #UPDATE PATH
            #l_download_string = l_download_string + \
            #'  '+'<a href="../datapile/'+targname+'/'+targname+'_coadd_G140L_final_lp3.fits.gz">LP3</a>   '
            l_download_string = l_download_string + \
            '  '+'<a href="../' + dir + '/'+targname+'/'+targname+'_coadd_G140L_final_lp3.fits.gz">LP3</a>   '
        else:
            l_download_string = l_download_string+'  '+'. . . .  '
    else:
        l_download_string = '. . . | . . . | . . . | . . . '

    if G185M_coadd_exists or G185M_coadd_exists or G285M_coadd_exists:
        #UPDAETE PATH
        #nuv_m_quicklook_urlstring = '<a href="../datapile/'+targname+'/'+targname+'_coadd_NUVM_final_lp1.png"><img height="40" src="../datapile/'+targname+'/'+targname+'_coadd_NUVM_final_lp1.png"></a>'
        nuv_m_quicklook_urlstring = '<a href="../' + dir + '/'+targname+'/'+targname+'_coadd_NUVM_final_lp1.png"><img height="40" src="../test_nuv_data/'+targname+'/'+targname+'_coadd_NUVM_final_lp1.png"></a>'

    if G230L_coadd_exists: 
        #nuv_l_quicklook_urlstring = '<a href="../datapile/'+targname+'/'+targname+'_coadd_NUVL_final_lp1.png"><img height="40" src="../datapile/'+targname+'/'+targname+'_coadd_NUVL_final_lp1.png"></a>'
        nuv_l_quicklook_urlstring = '<a href="../' + dir + '/'+targname+'/'+targname+'_coadd_NUVL_final_lp1.png"><img height="40" src="../' + dir + '/'+targname+'/'+targname+'_coadd_NUVL_final_lp1.png"></a>'
#    if (os.path.exists(hdr0['targname'].strip()+'_coadd_G130M_final_all.fits')): 
#        download_string = '<a href="../datapile/'+targname+'/'+targname+'_coadd_G130M_final_all.fits">ALL</a> |'
#
#        this_coadd = Table.read(targname+'_coadd_G130M_final_all.fits') 
#        i_good = np.where(this_coadd['FLUX'] > 0) 
#        median_sn_130 = np.median(this_coadd['SN'][i_good]) 
#        print 'Median G130M SN for ', targname, ' = ', str(median_sn)[0:6] 
#    
#        if (os.path.exists(hdr0['targname'].strip()+'_coadd_G130M_final_lp1.fits')): 
#            print 'I found a coadd for G130M LP=1' 
#            download_string = download_string  + \
#            '  '+'<a href="../datapile/'+targname+'/'+targname+'_coadd_G130M_final_lp1.fits">LP1</a> | '
#        else: 
#            download_string = download_string+'  ' + \
#            '. . . .  | ' 
#        if (os.path.exists(hdr0['targname'].strip()+'_coadd_G130M_final_lp2.fits')): 
#            print 'I found a coadd for G130M LP=2' 
#            download_string = download_string + \
#            '  '+'<a href="../datapile/'+targname+'/'+targname+'_coadd_G130M_final_lp2.fits">LP2</a> | '
#        else: 
#            download_string = download_string + \
#            '  '+'. . . .  | ' 
#        if (os.path.exists(hdr0['targname'].strip()+'_coadd_G130M_final_lp3.fits')): 
#            print 'I found a coadd for G130M LP=3' 
#            download_string = download_string + \
#            '  '+'<a href="../datapile/'+targname+'/'+targname+'_coadd_G130M_final_lp3.fits">LP3</a>   '
#        else: 
#            download_string = download_string+'  '+'. . . .  ' 
#    else: 
#        download_string = '. . . | . . . | . . . | . . . ' 

### why is SN a string here ??? ####

    if G140L_coadd_exists:
        this_coadd = Table.read(targname+'_coadd_G140L_final_lpALL.fits.gz') 
        i_good = np.where((this_coadd['FLUX'] > 0) & (this_coadd['WAVE'] > 1100) & (this_coadd['WAVE'] < 1900))  ## take regions that are not lousy S/N 
        median_sn = np.median(this_coadd['SN'][i_good]) * np.sqrt(7.0)
        print('Median G130M SN for ', targname, ' = ', str(median_sn)[0:6])
        #median_sn = str(median_sn)[0:6] 

    if G160M_coadd_exists:
        this_coadd = Table.read(targname+'_coadd_G160M_final_lpALL.fits.gz') 
        i_good = np.where(this_coadd['FLUX'] > 0) 
        median_sn = np.median(this_coadd['SN'][i_good])  * np.sqrt(7.0)
        print('Median G160M SN for ', targname, ' = ', str(median_sn)[0:6])
        #median_sn = str(median_sn)[0:6] 

    if G130M_coadd_exists:
        this_coadd = Table.read(targname+'_coadd_G130M_final_lpALL.fits.gz') 
        i_good = np.where(this_coadd['FLUX'] > 0) 
        median_sn = np.median(this_coadd['SN'][i_good])  * np.sqrt(7.0)
        print('Median G130M SN for ', targname, ' = ', str(median_sn)[0:6])
        #median_sn = str(median_sn)[0:6]

    if G185M_coadd_exists:
        this_coadd = Table.read(targname+'_coadd_G185M_final_lp1.fits.gz') 
        i_good = np.where(this_coadd['FLUX'] > 0) 
        median_sn = np.median(this_coadd['SN'][i_good])  * np.sqrt(3.0)
        print('Median G185M SN for ', targname, ' = ', str(median_sn)[0:6])
        #median_sn = str(median_sn)[0:6] 

    if G225M_coadd_exists:
        this_coadd = Table.read(targname+'_coadd_G225M_final_lp1.fits.gz') 
        i_good = np.where(this_coadd['FLUX'] > 0) 
        median_sn = np.median(this_coadd['SN'][i_good])  * np.sqrt(3.0)
        print('Median G225M SN for ', targname, ' = ', str(median_sn)[0:6])
        #median_sn = str(median_sn)[0:6] 

    if G285M_coadd_exists:
        this_coadd = Table.read(targname+'_coadd_G285M_final_lp1.fits.gz') 
        i_good = np.where(this_coadd['FLUX'] > 0) 
        median_sn = np.median(this_coadd['SN'][i_good])  * np.sqrt(3.0)
        print('Median G285M SN for ', targname, ' = ', str(median_sn)[0:6])
        #median_sn = str(median_sn)[0:6] 

    if G230L_coadd_exists:
        this_coadd = Table.read(targname+'_coadd_G230L_final_lp1.fits.gz') 
        i_good = np.where(this_coadd['FLUX'] > 0) 
        median_sn = np.median(this_coadd['SN'][i_good])  * np.sqrt(3.0)
        print('Median G230L SN for ', targname, ' = ', str(median_sn)[0:6])
        #median_sn = str(median_sn)[0:6] 
        
    if (os.path.exists(str(hdr0['targname']).strip()+'_target.tar.gz')): 
        #download_string = '<a href="../datapile/'+targname+'/'+targname+'_target.tar.gz">ALL</a>'
        download_string = '<a href="../' + dir + '/'+targname+'/'+targname+'_target.tar.gz">ALL</a>'
        

    webtable_row = [counter, targname_urlstring, ra, dec, n_exp_string, str.split(targdesc,';')[0],
        targdesc, altname_string, altname_class, redshift_string, mast_string, median_sn, 
                    fuv_m_quicklook_urlstring, fuv_l_quicklook_urlstring, nuv_m_quicklook_urlstring, nuv_l_quicklook_urlstring, download_string, l_download_string]

    print('WEBTABLE_ROW', webtable_row)
    fitstable_row = [counter, targname, ra, dec, nfiles, str.split(targdesc,';')[0], targdesc, altname_string, median_sn]
    print(fitstable_row)
    
    targetstable_row = [1,targname,str.split(targdesc,';')[0], targdesc]
    return webtable_row, fitstable_row, targetstable_row


#-----------------------------------------------------------------------------------------------------

def make_exposure_catalog(filelist):
    # exposure_cat contains database of all exposures for this target
    exposure_cat = Table(\
    names=('Flag', 'Rootname', 'Target Name', 'RA', 'DEC', 'PropID',\
                'PI Name', 'Detector', 'Segment', 'LP', 'Grating', 'Cenwave', 'FPPOS',\
                'Exptime', 'Nevents', 'Mean Flux', 'Median Flux', 'Date', 'Target Description'),   
            dtype=('i2', 'S20', 'S35', 'f4', 'f4', 'i',\
                   'S20', 'S4', 'S5', 'S2', 'S10', 'i', 'i',\
                    'f', 'f', 'f', 'f', 'S12', 'S200'))
#            dtype=('i2', 'S20', 'S35', 'f4', 'f4', 'i',\
#                   'S20', 'S4', 'S5', 'S2', 'S10', 'i', 'i',\
#                    'f10', 'f8', 'f8', 'f8', 'S4', 'S12', 'S200'))

    if len(filelist) == 0:
        return []

                    
    ## does alias.txt exist?
    alias_file = "alias.txt"
    if (os.path.exists(alias_file)):
        af = open(alias_file, 'r')
        targname = af.readline()
        af.close()
    else:
        print("---->>>> make_expsoure_catalog can't find ",alias_file,"!!!!! making one instead.....  <<<<----------")
        #JRD astropy
        observation = fits.open(filelist[0])
        hdr0 = observation[0].header
        #hdr0 = fitsio.read_header(filelist[0], 0) 
        targname = str(hdr0['TARGNAME']).strip()
        af = open(alias_file, 'w')
        af.write(targname)
        af.close()
                    
    LYA_MIN = 1206 ## should this depend on M vs L grating?
    LYA_MAX = 1226

    for filename in filelist:
        #JRD astropy
        observation =  fits.open(filename)
        hdr0 = observation[0].header
        #hdr0 = fitsio.read_header(filename, 0) 
        #data, hdr1 = fitsio.read(filename, ext=1, header=True)
        data = observation[1].data
        hdr1 = observation[1].header
        
        print("Obtaining headers for :", filename)
        if (np.shape(data)[0] < 1):
            print("no data:",filename)
        else:
            hdr0['TARGNAME'] = str(hdr0['TARGNAME']).strip()
            second_order = np.zeros(np.shape(data["FLUX"]))
            ### Need to exclude Stripe C on G230L because it's all 2nd Order Light
            if hdr0['OPT_ELEM'].strip() == "G230L":
                second_order[2] = 1


            segment_string = ' ' 
            if ('COS' in hdr0['PRIMESI']): 
                 indices = np.where((data["DQ_WGT"] > 0) & (data["DQ"] == 0) & ((data["WAVELENGTH"] > LYA_MAX) | (data["WAVELENGTH"] < LYA_MIN)) & (second_order == 0))
                 segment_string = hdr0['SEGMENT'].strip() 
                 lp_string = hdr0['LIFE_ADJ'] 
                 fppos_string = hdr0['FPPOS'] 
                 n_events = hdr1['NEVENTS']
            if ('STIS' in hdr0['PRIMESI']): 
                 indices = np.where((data["DQ"] == 0) & ((data["WAVELENGTH"] > LYA_MAX) | (data["WAVELENGTH"] < LYA_MIN)) & (second_order == 0))
                 segment_string = hdr0['DETECTOR'].strip() 
                 lp_string = 'N/A'
                 fppos_string = 0 
                 n_events = 0 
                
            flag = 1 
            
           
            exposure_cat.add_row([flag, hdr0['ROOTNAME'].strip(), hdr0['TARGNAME'].strip(), hdr0['RA_TARG'], hdr0['DEC_TARG'], \
                hdr0['PROPOSID'], hdr0['PR_INV_L'].strip(), hdr0['DETECTOR'].strip(), segment_string, lp_string, \
                hdr0['OPT_ELEM'].strip(), hdr0['CENWAVE'], fppos_string, hdr1['EXPTIME'], n_events, \
                np.mean(data["FLUX"][indices]), np.median(data['FLUX'][indices]), \
                hdr1['DATE-OBS'].strip(), hdr0['TARDESCR'].strip()] )  

            ## want a way to consolidate every header keyword for every exposure into single table / file,
            ## but this is really slow. method other than vstack?
            if (False): 
                for key in ["HISTORY", "COMMENT",""]: 
                    del hdr0[key] 
                names = sorted(hdr0.keys()) 
                rows = [hdr0[k] for k in names]
                header_table0 = vstack([header_table0, Table(rows=[rows], names=names)])
    
                for key in ["HISTORY", "COMMENT",""]: 
                    del hdr1[key] 
                names = sorted(hdr1.keys()) 
                rows = [hdr1[k] for k in names]
                header_table1 = vstack([header_table1, Table(rows=[rows], names=names)])
                total_number_of_headers = total_number_of_headers + 1 
                print('TOTAL NUMBER OF HEADERS : ', total_number_of_headers)

    ascii.write(exposure_cat, 'all_exposures.txt')  # write out the exposures for this target by itself 
    exposure_cat.write('all_exposures.html', format='jsviewer') # write out the exposures for this target by itself

    print("ALL EXPOSURES")   
    print(exposure_cat)
    return exposure_cat



#-----------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    args = parse_args()
    targets = args.filename
    
    scrape_headers(targets,args.altnames,args.redshifts)
    sys.exit("""
    
    ~~~~~~~*~*~*~*~
    ~~~~~~~*~*~*~*~  all done!!!! spectra are fun!
    ~~~~~~~*~*~*~*~
    """)
