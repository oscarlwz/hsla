from astropy.io import fits as f
import fitsio
from astropy.table import Table, Column, vstack 
from astropy.io import ascii
from datetime import datetime
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import copy
import sys
import os

import counts_error 

def main(screen, lifetime, Linda=False):

    #######################################################################################################
    #               STEP 1: obtain all exposures in the "all_exposures.txt" file                          #
    #######################################################################################################

    tt, h0, h1  = get_dict_of_exposures()             # obtain the headers and data for ALL exposures in the target directory

    if tt == 0: 
        print "It appears there is no all_exposures.txt file in this directory, attemping to exit gracefully."  
        return 0  

    targname = h0[h0.keys()[0]]['TARGNAME']
    print 'ADD.MAIN: We will now perform a coadd Target ', targname

    #######################################################################################################
    #               STEP 2: screen out the exposures that we don't want to have                           #
    #######################################################################################################

    #lifetime = 'ALL'
    print 'We will perform coadds for lifetime position : ', lifetime 
    tt, h0, h1 = screen_dict_of_exposures(tt, h0, h1, screen=screen, lifetime=lifetime)
    number_of_remaining_exposures = np.size(h0.keys())
    print 'ADD.MAIN: Exposure cuts are: ', screen, ' and ', lifetime
    print 'ADD.MAIN: Number of exposures surviving the cuts: ', number_of_remaining_exposures

    #######################################################################################################
    #               STEP 3: create the output dictionary and then populate it with a things we have       #
    #######################################################################################################

    out_dict = {}             ##### start with an empty dictionary to contain the output stuff

    out_dict['lifetime'] = lifetime 

    out_dict['exptable'] = tt  ##### these keys contain the *input* header and data file contents
    out_dict['h0table'] = h0
    out_dict['h1table'] = h1
    out_dict['targname'] = targname

    #######################################################################################################
    #               STEP 4: Add DQ_WGT entries to exposure tables AND do DQ screening                     #
    #######################################################################################################

    for nn in tt.keys():                                    #### Add a few columns to the exposure tables to include per-pixel exptime and DQ weights
        tt[nn]['EXP_PIX'] = tt[nn]['FLUX'] * 0.0 + tt[nn][0]['EXPTIME']            ####    every pixel in the output array gets the same exposure time, for now.
        tt[nn]['EXP_PIX'].unit = 's'
        tt[nn]['DQ_WEIGHT'] = tt[nn]['FLUX'] * 0.0 + 1.0                #### DQ_WGT is 1 by default
        tt[nn]['DQ_WEIGHT'].unit = '   '
        tt[nn]['FLUXFACTOR'] = np.nan_to_num(tt[nn]['FLUX'] / tt[nn]['NET'])      #### use nan_to_num to eat the nonsensical zeros in gross counts (denominator)
  

        tt[nn]['GRATING'] = h0[nn]['OPT_ELEM']					  #### add these useful bits of information to each exposure in the table 
        tt[nn]['CENWAVE'] = h0[nn]['CENWAVE'] 

        number_of_segments = np.size(tt[nn]['SEGMENT'])
        
        for nseg in np.arange(number_of_segments):
            test = (tt[nn]['DQ'][nseg] & 2**0) | (tt[nn]['DQ'][nseg] & 2**3)  | (tt[nn]['DQ'][nseg] & 2**6) | (tt[nn]['DQ'][nseg] & 2**7) | \
                (tt[nn]['DQ'][nseg] & 2**8) | (tt[nn]['DQ'][nseg] & 2**9) | (tt[nn]['DQ'][nseg] & 2**11) | (tt[nn]['DQ'][nseg] & 2**13)
            i_omit = np.where(test > 0)
            tt[nn]['FLUX'][nseg][i_omit] = 0.
            tt[nn]['ERROR'][nseg][i_omit] = 0.
            tt[nn]['GROSS'][nseg][i_omit] = 0.
            tt[nn]['GCOUNTS'][nseg][i_omit] = 0.
            tt[nn]['NET'][nseg][i_omit] = 0.
            tt[nn]['BACKGROUND'][nseg][i_omit] = 0.
            tt[nn]['EXP_PIX'][nseg][i_omit] = 0.
            tt[nn]['DQ_WEIGHT'][nseg][i_omit] = 0.
            tt[nn]['DQ'][nseg][i_omit] = 0.
            # set the DQ flags to zero 

            test = (tt[nn]['DQ'][nseg] & 2**2) | (tt[nn]['DQ'][nseg] & 2**4) | (tt[nn]['DQ'][nseg] & 2**5) | (tt[nn]['DQ'][nseg] & 2**10) | (tt[nn]['DQ'][nseg] & 2**12)
            i_deweight = np.where(test > 0)
            tt[nn]['EXP_PIX'][nseg][i_omit] = tt[nn]['EXP_PIX'][nseg][i_omit] / 2.             #### set exposure time down by 2x to deweight
            tt[nn]['DQ_WEIGHT'][nseg][i_omit] = 0.5
            # leave the DQ flags unmodified 

            tt[nn][nseg]['GRATING'] = h0[nn]['OPT_ELEM'] 

    #######################################################################################################
    #               STEP 5: Create the main dictionary entries for final spectra                          #
    #######################################################################################################

    count_exposures(out_dict)

    if (out_dict['number_of_G130M_exposures'] > 0):
        final_wave_130 = get_wavelength_grid('G130M')
        out_dict['G130M_final_all'] = create_coadd_format(final_wave_130)

    if (out_dict['number_of_G160M_exposures'] > 0):
        final_wave_160 = get_wavelength_grid('G160M')
        out_dict['G160M_final_all'] = create_coadd_format(final_wave_160)

    if (out_dict['number_of_G140L_exposures'] > 0):
        final_wave_140 = get_wavelength_grid('G140L')
        out_dict['G140L_final_all'] = create_coadd_format(final_wave_140)

    if (out_dict['number_of_G185M_exposures'] > 0):
        final_wave_185 = get_wavelength_grid('G185M')
        out_dict['G185M_final_all'] = create_coadd_format(final_wave_185)

    if (out_dict['number_of_G225M_exposures'] > 0):
        final_wave_225 = get_wavelength_grid('G225M')
        out_dict['G225M_final_all'] = create_coadd_format(final_wave_225)

    #######################################################################################################
    #               STEP 6: obtain the new pixel coordinates of the pixels in each x1d epxosure          #
    #######################################################################################################

    for nn in tt.keys():                                    ##### Add a few column to the exposure table for the "newpix"

        if ('FUV' in screen): tt[nn]['NEWPIX'] = tt[nn]['FLUX'] * 0.0 +  np.arange(16384)
        if ('NUV' in screen): tt[nn]['NEWPIX'] = tt[nn]['FLUX'] * 0.0 +  np.arange(1274)

        tt[nn]['NEWPIX'].unit = ' '

    a = get_pix_to_wave(out_dict, screen)


    #######################################################################################################
    #               STEP 7: replicate the structure of the exposure table and interpolate to it.          # 
    #######################################################################################################

    get_interp_exptable(out_dict) 			     #### this creates an exposure table paralleling "exptable" 
							     #### but with the new wavelength grid instead of the old one 

    #############################################################################################################
    #               STEP 8: "interpolate" the input exposures to the output grid, populate exp_interp_table    # 
    #############################################################################################################

    exposure_interpolate(out_dict)			     #### this interpolates the input exposures (exptable) 
							     #### onto the new wavelength grid (exp_interp_table) 

    #######################################################################################################
    #               STEPS 9: go ahead and do the final coadd                                              #
    #######################################################################################################

    print "Passing ", np.size(tt.keys()), ' exposures to cos_counts_coadd at LP = ALL '
    number_of_remaining_exposures = np.size(h0.keys())
    total_number_of_exposures = number_of_remaining_exposures
    print 'Prior to calling cos_counts_coadd, the exposure counts are: '
    print '       G130M: ', out_dict['number_of_G130M_exposures']
    print '       G160M: ', out_dict['number_of_G160M_exposures']
    print '       G140L: ', out_dict['number_of_G140L_exposures']
    print '       G185M: ', out_dict['number_of_G185M_exposures']
    print '       G225M: ', out_dict['number_of_G225M_exposures']

    print "POTENTIAL PROBLEM HERE PRIOR TO RUNNING COS_COUNTS_COADD ON ALL LPs: " 
    print 'number of keys in h0table', np.size(out_dict['h0table'].keys()), out_dict['h0table'].keys() 
    print 'number of keys in exp_interp_table', np.size(out_dict['exp_interp_table'].keys()), out_dict['exp_interp_table'].keys() 

    if (out_dict['number_of_G130M_exposures'] > 0):
        out_dict = cos_counts_coadd(out_dict, 'G130M_final_all')
    if (out_dict['number_of_G160M_exposures'] > 0):
        out_dict = cos_counts_coadd(out_dict, 'G160M_final_all')
    if (out_dict['number_of_G140L_exposures'] > 0):
        out_dict = cos_counts_coadd(out_dict, 'G140L_final_all')
    if (out_dict['number_of_G185M_exposures'] > 0):
        out_dict = cos_counts_coadd(out_dict, 'G185M_final_all')
    if (out_dict['number_of_G225M_exposures'] > 0):
        out_dict = cos_counts_coadd(out_dict, 'G225M_final_all')

#    lp1_dict = copy.deepcopy(out_dict)                             #### now we're going to do LP1 with cos_counts_coadd
#    count_exposures(lp1_dict)
#    print "Screening ", np.size(lp1_dict['exptable'].keys()), ' exposures to get LP = 1'
#    tt_lp1, h0_lp1, h1_lp1 = screen_dict_of_exposures(lp1_dict['exptable'], lp1_dict['h0table'], lp1_dict['h1table'], screen=screen, lifetime=1)
#    get_interp_exptable(lp1_dict) 			     
#    exposure_interpolate(lp1_dict)			     
#    count_exposures(lp1_dict)
#    print "There are ", number_of_remaining_exposures, ' exposures at LP = 1'
#    if (lp1_dict['number_of_G130M_exposures'] > 0):
#        print "We will now coadd ", lp1_dict['number_of_G130M_exposures'], " at LP = 1"
#        lp1_dict['G130M_final_lp1'] = create_coadd_format(final_wave_130)
#        lp1_dict = cos_counts_coadd(lp1_dict, 'G130M_final_lp1')
#        if (lp1_dict['number_of_G130M_exposures'] > 0): out_dict['G130M_final_lp1'] = lp1_dict['G130M_final_all']
#    if (lp1_dict['number_of_G160M_exposures'] > 0):
#        print "We will now coadd ", lp1_dict['number_of_G160M_exposures'], " at LP = 1"
#        lp1_dict['G160M_final_lp1'] = create_coadd_format(final_wave_160)
#        lp1_dict = cos_counts_coadd(lp1_dict, 'G160M_final_lp1')
#        if (lp1_dict['number_of_G160M_exposures'] > 0): out_dict['G160M_final_lp1'] = lp1_dict['G160M_final_all']
#    if (lp1_dict['number_of_G140L_exposures'] > 0):
#        print "We will now coadd ", lp1_dict['number_of_G140L_exposures'], " at LP = 1"
#        lp1_dict['G140L_final_lp1'] = create_coadd_format(final_wave_140)
#        lp1_dict = cos_counts_coadd(lp1_dict, 'G140L_final_lp1')
#        if (lp1_dict['number_of_G140L_exposures'] > 0): out_dict['G140L_final_lp1'] = lp1_dict['G140L_final_all']
#
#
#    lp2_dict = copy.deepcopy(out_dict)                             #### now we're going to do LP2 with cos_counts_coadd
#    count_exposures(lp2_dict)
#    print "Screening ", np.size(lp2_dict['exptable'].keys()), ' exposures to get LP = 2'
#    tt_lp2, h0_lp2, h1_lp2 = screen_dict_of_exposures(lp2_dict['exptable'], lp2_dict['h0table'], lp2_dict['h1table'], screen=screen, lifetime=2)
#    get_interp_exptable(lp2_dict) 			     
#    exposure_interpolate(lp2_dict)			     
#    count_exposures(lp2_dict)
#    number_of_remaining_exposures = np.size(h0_lp2.keys())
#    print "There are ", number_of_remaining_exposures, ' exposures at LP = 2'
#    if (lp2_dict['number_of_G130M_exposures'] > 0):
#        print "We will now coadd ", lp2_dict['number_of_G130M_exposures'], " at LP = 2"
#        lp2_dict['G130M_final_lp2'] = create_coadd_format(final_wave_130)
#        lp2_dict = cos_counts_coadd(lp2_dict, 'G130M_final_lp2')
#        if (lp2_dict['number_of_G130M_exposures'] > 0): out_dict['G130M_final_lp2'] = lp2_dict['G130M_final_all']
#    if (lp2_dict['number_of_G160M_exposures'] > 0):
#        print "We will now coadd ", lp2_dict['number_of_G160M_exposures'], " at LP = 2"
#        lp2_dict['G160M_final_lp2'] = create_coadd_format(final_wave_160)
#        lp2_dict = cos_counts_coadd(lp2_dict, 'G160M_final_lp2')
#        if (lp2_dict['number_of_G160M_exposures'] > 0): out_dict['G160M_final_lp2'] = lp2_dict['G160M_final_all']
#    if (lp2_dict['number_of_G140L_exposures'] > 0):
#        print "We will now coadd ", lp2_dict['number_of_G140L_exposures'], " at LP = 2"
#        lp2_dict['G140L_final_lp2'] = create_coadd_format(final_wave_140)
#        lp2_dict = cos_counts_coadd(lp2_dict, 'G140L_final_lp2')
#        if (lp2_dict['number_of_G140L_exposures'] > 0): out_dict['G140L_final_lp2'] = lp2_dict['G140L_final_all']
#
#    lp3_dict = copy.deepcopy(out_dict)                             #### now we're going to do LP3 with cos_counts_coadd
#    count_exposures(lp3_dict)
#    print "Screening ", np.size(lp3_dict['exptable'].keys()), ' exposures to get LP = 3'
#    tt_lp3, h0_lp3, h1_lp3 = screen_dict_of_exposures(lp3_dict['exptable'], lp3_dict['h0table'], lp3_dict['h1table'], screen=screen, lifetime=3)
#    get_interp_exptable(lp3_dict) 			     
#    exposure_interpolate(lp3_dict)			     
#    count_exposures(lp3_dict)
#    number_of_remaining_exposures = np.size(h0_lp3.keys())
#    print "There are ", number_of_remaining_exposures, ' exposures at LP = 3'
#    if (lp3_dict['number_of_G130M_exposures'] > 0):
#        print "We will now coadd ", lp3_dict['number_of_G130M_exposures'], " at LP = 3"
#        lp3_dict['G130M_final_lp3'] = create_coadd_format(final_wave_130)
#        lp3_dict = cos_counts_coadd(lp3_dict, 'G130M_final_lp3')
#        if (lp3_dict['number_of_G130M_exposures'] > 0): out_dict['G130M_final_lp3'] = lp3_dict['G130M_final_all']
#    if (lp3_dict['number_of_G160M_exposures'] > 0):
#        print "We will now coadd ", lp3_dict['number_of_G160M_exposures'], " at LP = 3"
#        lp3_dict['G160M_final_lp3'] = create_coadd_format(final_wave_160)
#        lp3_dict = cos_counts_coadd(lp3_dict, 'G160M_final_lp3')
#        if (lp3_dict['number_of_G160M_exposures'] > 0): out_dict['G160M_final_lp3'] = lp3_dict['G160M_final_all']
#    if (lp3_dict['number_of_G140L_exposures'] > 0):
#        print "We will now coadd ", lp3_dict['number_of_G140L_exposures'], " at LP = 3"
#        lp3_dict['G140L_final_lp3'] = create_coadd_format(final_wave_140)
#        lp3_dict = cos_counts_coadd(lp3_dict, 'G140L_final_lp3')
#        if (lp3_dict['number_of_G140L_exposures'] > 0): out_dict['G140L_final_lp3'] = lp3_dict['G140L_final_all']


#    lp12_dict = copy.deepcopy(out_dict)                             #### now we're going to do LP1+2 with cos_counts_coadd
#    count_exposures(lp12_dict)
#    print "Screening ", np.size(lp12_dict['exptable'].keys()), ' exposures to get LP = 1+2'
#    tt_lp12, h0_lp12, h1_lp12 = screen_dict_of_exposures(lp12_dict['exptable'], lp12_dict['h0table'], lp12_dict['h1table'], screen=screen, lifetime=12) 
#    get_interp_exptable(lp12_dict) 			     
#    exposure_interpolate(lp12_dict)			     
#    count_exposures(lp12_dict)
#    number_of_remaining_exposures = np.size(h0_lp12.keys())
#    print "There are ", number_of_remaining_exposures, ' exposures at LP = 3'
#    if (lp12_dict['number_of_G130M_exposures'] > 0):
#        print "We will now coadd ", lp12_dict['number_of_G130M_exposures'], " at LP = 1+2"
#        lp12_dict['G130M_final_lp12'] = create_coadd_format(final_wave_130)
#        lp12_dict = cos_counts_coadd(lp12_dict, 'G130M_final_lp12')
#        if (lp12_dict['number_of_G130M_exposures'] > 0): out_dict['G130M_final_lp12'] = lp12_dict['G130M_final_all']
#    if (lp12_dict['number_of_G160M_exposures'] > 0):
#        print "We will now coadd ", lp12_dict['number_of_G160M_exposures'], " at LP = 1+2"
#        lp12_dict['G160M_final_lp12'] = create_coadd_format(final_wave_160)
#        lp12_dict = cos_counts_coadd(lp12_dict, 'G160M_final_lp12')
#        if (lp12_dict['number_of_G160M_exposures'] > 0): out_dict['G160M_final_lp12'] = lp12_dict['G160M_final_all']
#    if (lp12_dict['number_of_G140L_exposures'] > 0):
#        print "We will now coadd ", lp12_dict['number_of_G140L_exposures'], " at LP = 1+2"
#        lp12_dict['G140L_final_lp12'] = create_coadd_format(final_wave_140)
#        lp12_dict = cos_counts_coadd(lp12_dict, 'G140L_final_lp12')
#        if (lp12_dict['number_of_G140L_exposures'] > 0): out_dict['G140L_final_lp12'] = lp12_dict['G140L_final_all']




    #######################################################################################################
    #               STEP 10:  obtain wavelength shifts for each exposure relatiive to the coadd, apply    #
    #######################################################################################################

    #if (out_dict['number_of_G130M_exposures'] > 0):
    #    shifts_to_apply = get_wavelength_shifts(out_dict) 
    #    print "WE DERIVED THE FOLLOWING WAVELENGTH SHIFTS:", shifts_to_apply 
    #    apply_wavelength_shifts(out_dict, shifts_to_apply) 

    #######################################################################################################
    #               STEP 11:  redo LP=ALL coaadd after shifts have been applied                           #
    #######################################################################################################

    #if (out_dict['number_of_G130M_exposures'] > 0):
    #    out_dict['G130M_final_shifted'] = create_coadd_format(final_wave_130)
    #    out_dict = cos_counts_coadd(out_dict, 'G130M_final_shifted')

    #######################################################################################################
    #               STEP 12:  get the chi2 comparison between inputs and outputs                          #
    #######################################################################################################

    if (out_dict['number_of_G130M_exposures'] > 0):
        if (out_dict.__contains__('G130M_final_all')): obtain_chi2(out_dict, 'G130M_final_all') 
    if (out_dict['number_of_G160M_exposures'] > 0):
        if (out_dict.__contains__('G160M_final_all')): obtain_chi2(out_dict, 'G160M_final_all') 

    if (out_dict.__contains__('G130M_final_lp1')): obtain_chi2(out_dict, 'G130M_final_lp1') 
    if (out_dict.__contains__('G160M_final_lp1')): obtain_chi2(out_dict, 'G160M_final_lp1') 

    if (out_dict.__contains__('G130M_final_lp2')): obtain_chi2(out_dict, 'G130M_final_lp2') 
    if (out_dict.__contains__('G160M_final_lp2')): obtain_chi2(out_dict, 'G160M_final_lp2') 

    if (out_dict.__contains__('G130M_final_lp3')): obtain_chi2(out_dict, 'G130M_final_lp3') 
    if (out_dict.__contains__('G160M_final_lp3')): obtain_chi2(out_dict, 'G160M_final_lp3') 

    if (out_dict.__contains__('G130M_final_lp12')): obtain_chi2(out_dict, 'G130M_final_lp12') 
    if (out_dict.__contains__('G160M_final_lp12')): obtain_chi2(out_dict, 'G160M_final_lp12') 

    #######################################################################################################
    #               STEP 13:  if you have both G130M and and G160M, splice them and store output          #
    #######################################################################################################

    if out_dict.__contains__('G130M_final_all') and out_dict.__contains__('G160M_final_all'):  
        splice_fuv(out_dict, 'all') 
    if out_dict.__contains__('G130M_final_lp1') and out_dict.__contains__('G160M_final_lp1'):  
        splice_fuv(out_dict, 'lp1') 
    if out_dict.__contains__('G130M_final_lp2') and out_dict.__contains__('G160M_final_lp2'):  
        splice_fuv(out_dict, 'lp2') 
    if out_dict.__contains__('G130M_final_lp3') and out_dict.__contains__('G160M_final_lp3'):  
        splice_fuv(out_dict, 'lp3') 
    if out_dict.__contains__('G130M_final_lp12') and out_dict.__contains__('G160M_final_lp12'):  
        splice_fuv(out_dict, 'lp12') 

    #######################################################################################################
    #               STEP 1X: write out the output file that the quick_look etc is going to use             #
    #######################################################################################################

    if (total_number_of_exposures > 0):
        write_output(out_dict, Linda)

    #######################################################################################################
    #               STEPS WHATEVER+1:                                                                     #
    #######################################################################################################

    return out_dict



######## 
######## 
######## 
######## END OF THE MAIN PROGRAM HERE
######## 
######## 
######## 





def splice_fuv(out_dict, key_to_combine): 

    g130m = out_dict['G130M_final_'+key_to_combine] 
    g160m = out_dict['G160M_final_'+key_to_combine] 

    where_160 = np.where((g160m['DQ'] == 8.) & (g160m['WAVE'] <  1450)) 
    if (np.size(where_160) > 1): 
       splitwave = np.max(g160m['WAVE'][np.where((g160m['DQ'] == 8.) & (g160m['WAVE'] <  1450))]) + 1. 
       i_130 = np.max(np.where(g130m['WAVE'] < splitwave))  
       i_160 = np.min(np.where(g160m['WAVE'] > splitwave))  
       out_dict['FUVM_final_'+key_to_combine] = vstack([ g130m[0:i_130], g160m[i_160:-1] ]) 
       print "SPLICE_FUV", splitwave 
    elif (np.size(where_160) < 1): 
       where_160 = np.where((g160m['DQ'] > 0.) & (g160m['WAVE'] <  1450)) 
       splitwave = np.max(g160m['WAVE'][np.where((g160m['DQ'] > 0.) & (g160m['WAVE'] <  1450))]) + 1. 
       i_130 = np.max(np.where(g130m['WAVE'] < splitwave))  
       i_160 = np.min(np.where(g160m['WAVE'] > splitwave))  
       out_dict['FUVM_final_'+key_to_combine] = vstack([ g130m[0:i_130], g160m[i_160:-1] ]) 
       print "SPLICE_FUV", splitwave 
    else: 
       print 'Problem with splitwave, nothing will be done' 
       






def apply_wavelength_shifts(out_dict, shifts_to_apply): 

    exp_interp_table = out_dict['exp_interp_table'] 
    exp_keys = exp_interp_table.keys()
    number_of_exposures = np.size(exp_interp_table.keys())

    columns_to_shift = ['FLUX', 'ERROR', 'GROSS', 'GCOUNTS','NET','BACKGROUND','DQ','DQ_WGT','EXP_PIX','DQ_WEIGHT','FLUXFACTOR'] 

    for nn in np.arange(number_of_exposures):
       for col in columns_to_shift: 
          exp_interp_table[exp_keys[nn]][1][col] = np.roll(exp_interp_table[exp_keys[nn]][1][col], int(-1*shifts_to_apply[nn])) 


def get_wavelength_shifts(out_dict): 
 
    exp_interp_table = out_dict['exp_interp_table'] 
    exp_keys = exp_interp_table.keys()
    number_of_exposures = np.size(exp_interp_table.keys())
 
    shifts_to_apply = np.zeros(number_of_exposures) 
 
    for nn in np.arange(number_of_exposures):

        exp_wave = exp_interp_table[exp_keys[nn]][1]['WAVELENGTH'] 				#### this is just the first exposure in the list 
        exp_flux = 1e14*exp_interp_table[exp_keys[nn]][1]['FLUX'] 					#### this is just the first exposure in the list 

        grating = exp_interp_table[exp_keys[nn]][1]['GRATING']
        
        if '130M' in grating: 
            wavelimits = [1255, 1265]  
            final_table = out_dict['G130M_final_all'] 
            coadd_wave = final_table['WAVE']							#### this is the final_table 
            coadd_flux = 1e14*final_table['FLUX']							#### this is the final_table 
        if '160M' in grating: 
            wavelimits = [1520, 1530]  
            final_table = out_dict['G160M_final_all'] 
            coadd_wave = final_table['WAVE']							#### this is the final_table 
            coadd_flux = 1e14*final_table['FLUX']							#### this is the final_table 

        if (exp_keys[nn] == 4): 
            exp_interp_table[exp_keys[nn]].write('exp.fits',overwrite=True) 
            final_table.write('final_table.fits',overwrite=True) 

        where1 = np.where((coadd_wave > wavelimits[0]) & (coadd_wave < wavelimits[1])) 
        where2 = np.where((exp_wave > wavelimits[0]) & (exp_wave < wavelimits[1])) 
 
        lags = np.arange(51)+1 
        shifts = np.concatenate((-1*lags[::-1], lags),axis=0) 
        coeffs = shifts * 0.0 

        i = 0 
        for shifter in shifts: 
          #  fig = plt.figure(figsize=(8,4),dpi=300)
          #  ax = fig.add_subplot(111) 
            coeff_here = np.correlate(np.roll(exp_flux[where2],shifter), coadd_flux[where1]) 
            #print "GET_WAVELENGTH SHIFTS", shifter, coeff_here 
            print 'EXP', exp_keys[nn], grating, wavelimits, shifter, coeff_here[0]
            coeffs[i] = coeff_here[0] 
    	    i = i + 1 
          #  ax.step(coadd_wave[where1], coadd_flux[where1]) 
          #  ax.step(exp_wave[where2], np.roll(exp_flux[where2],shifter), color='red') 
          #  plt.savefig('exp'+str(exp_keys[nn])+'_'+str(shifter)+'.png')
          #  plt.close() 
 
        best_shift = shifts[np.where(coeffs == np.max(coeffs))] 
        print 'BEST SHIFT', exp_keys[nn], grating, wavelimits, best_shift, coeffs[np.where(coeffs == np.max(coeffs))]
        shifts_to_apply[nn] = best_shift 

    return shifts_to_apply 


     
def obtain_chi2(out_dict, final_spectrum_key):					   ##### 

    out_dict[final_spectrum_key]['EXP_COUNTER'] = out_dict[final_spectrum_key]['FLUX'] * 0.0 

    print 'OBTAIN CHI2 for: ', final_spectrum_key 
    final_spectrum = out_dict[final_spectrum_key] 
   
    exp_interp_table = out_dict['exp_interp_table']       #### NOTE: exptable[i][0] refers to segment A, exptable[i][1] refers to segment B
    exp_keys = exp_interp_table.keys()

    temp_chi2_array = final_spectrum['WAVE'] * 0.0 
    temp_chi2_array_alt = final_spectrum['WAVE'] * 0.0 

    for pixel in np.arange(np.size(temp_chi2_array)): # for each pixel 
        exp_counter = 0. 
        pixsum = 0.0 
         
        for exp in exp_keys: 

            number_of_segments = np.size(exp_interp_table[exp]['SEGMENT'])

            segment = 0 
            if (exp_interp_table[exp]['GRATING'][segment] in final_spectrum_key) and (exp_interp_table[exp]['EXP_PIX'][segment][pixel] > 0.): 
               if (exp_interp_table[exp]['NET'][segment][pixel] >= final_spectrum['NETCOUNTS'][pixel]): 
                   NET_ERR = exp_interp_table[exp]['NET_ERR_DOWN'][segment][pixel]
               else: 
                   NET_ERR = exp_interp_table[exp]['NET_ERR_UP'][segment][pixel]
               if (exp_interp_table[exp]['NET'][segment][pixel] <= 0.): 
                   NET_ERR = exp_interp_table[exp]['NET_ERR_UP'][segment][pixel]
               pixsum = pixsum + (exp_interp_table[exp]['NET'][segment][pixel] - final_spectrum['NETCOUNTS'][pixel])**2 / NET_ERR**2 
               exp_counter = exp_counter + 1. 
  
            if (number_of_segments > 1): 
                segment = 1 
                if (exp_interp_table[exp]['GRATING'][segment] in final_spectrum_key) and (exp_interp_table[exp]['EXP_PIX'][segment][pixel] > 0.): 
                   NET_ERR_UP = exp_interp_table[exp]['NET_ERR_UP'][segment][pixel]
                   NET_ERR_DOWN = exp_interp_table[exp]['NET_ERR_DOWN'][segment][pixel]
                   if (exp_interp_table[exp]['NET'][segment][pixel] >= final_spectrum['NETCOUNTS'][pixel]): 
                       NET_ERR = NET_ERR_DOWN 
                   else: 
                       NET_ERR = NET_ERR_UP 
                   if (exp_interp_table[exp]['NET'][segment][pixel] <= 0.): 
                       NET_ERR = NET_ERR_UP 
                   pixsum = pixsum + (exp_interp_table[exp]['NET'][segment][pixel] - final_spectrum['NETCOUNTS'][pixel])**2 / NET_ERR**2 
                   exp_counter = exp_counter + 1. 

        if (exp_counter > 1): 
            temp_chi2_array_alt[pixel] = pixsum / (exp_counter - 1.) 
            out_dict[final_spectrum_key]['CHI2'][pixel] = temp_chi2_array_alt[pixel]
            out_dict[final_spectrum_key]['EXP_COUNTER'][pixel] = exp_counter 
#
    out_dict[final_spectrum_key]['CHI2'] = temp_chi2_array_alt 
  
 






def count_exposures(out_dict):				##### counts the number of exposures with certain properties and
                              				##### places these in the dictionary
    counter130 = 0
    counter160 = 0
    counter140 = 0
    counter185 = 0
    counter225 = 0

    for i in out_dict['h0table'].keys():
        print 'COUNT_EXPOSURES: OPT_ELEM: ', out_dict['h0table'][i]['ROOTNAME'], out_dict['h0table'][i]['OPT_ELEM']
        if(out_dict['h0table'][i]['OPT_ELEM'] == 'G130M'):
            counter130 = counter130 + 1
        if(out_dict['h0table'][i]['OPT_ELEM'] == 'G160M'):
            counter160 = counter160 + 1
        if(out_dict['h0table'][i]['OPT_ELEM'] == 'G140L'):
            counter140 = counter140 + 1
        if(out_dict['h0table'][i]['OPT_ELEM'] == 'G185M'):
            counter185 = counter185 + 1
        if(out_dict['h0table'][i]['OPT_ELEM'] == 'G225M'):
            counter225 = counter225 + 1

    out_dict['number_of_G130M_exposures'] = counter130
    out_dict['number_of_G160M_exposures'] = counter160
    out_dict['number_of_G140L_exposures'] = counter140
    out_dict['number_of_G185M_exposures'] = counter185
    out_dict['number_of_G225M_exposures'] = counter225 

    print 'COUNT_EXPOSURES: There are ', counter130, out_dict['number_of_G130M_exposures']
    print 'COUNT_EXPOSURES: There are ', counter160, out_dict['number_of_G160M_exposures']
    print 'COUNT_EXPOSURES: There are ', counter140, out_dict['number_of_G140L_exposures']
    print 'COUNT_EXPOSURES: There are ', counter185, out_dict['number_of_G185M_exposures']
    print 'COUNT_EXPOSURES: There are ', counter225, out_dict['number_of_G225M_exposures']


def create_coadd_format(wavegrid):                                ##### this creates a dictionary that holds the canonical outputs of the coadd
                                                ##### this has to be a dictionary because it contains a bunch of tables
                                                ##### with different levels of the coadd hierarchy in them
    flux = wavegrid * 0.0                                     ##### with the input wavelength grid
    error = wavegrid * 0.0

    print 'CREATE_OUTPUT: creating master output table'
    t = Table([wavegrid, flux, error], names=('WAVE','FLUX','ERROR'))                ##### these are the basics

    t['FLUXERR_UP'] = t['WAVE'] * 0.0
    t['FLUXERR_DOWN'] = t['WAVE'] * 0.0
    t['GROSSCOUNTS'] = t['WAVE'] * 0.0
    t['NETCOUNTS'] = t['WAVE'] * 0.0                                 ##### will contain the S/N per pixel
    t['NETCOUNTSERR_UP'] = t['WAVE'] * 0.0                                 ##### will contain the Gehrels (1986) error on the combined NETCOUNTS 
    t['NETCOUNTSERR_DOWN'] = t['WAVE'] * 0.0                                 ##### will contain the Gehrels (1986) error on the combined NETCOUNTS 
    t['BACK'] = t['WAVE'] * 0.0
    t['EXP_PIX'] = t['WAVE'] * 0.0
    t['DQ'] = t['WAVE'] * 0.0
    t['SN'] = t['WAVE'] * 0.0                                     ##### this is the factor multipled by gross count rate to obtain flux
    t['FLUXFACTOR'] = t['WAVE'] * 0.0
    t['FLUXWGT'] = t['WAVE'] * 0.0				  #### these will be the exposure time weighted fluxes, which are alternatives to the counts-based coadd 
    t['FLUXWGT_ERR'] = t['WAVE'] * 0.0
    t['CHI2'] = t['WAVE'] * 0.0                                    ##### will contain a measure of statistical fluctuation between coadd and inputs
    ##### add any other desired contents of the main output final data table here. . . .

    return t                                                ##### return the output dictionary to the calling routine


def get_pix_to_wave(out_dict, screen):

    exptable = out_dict['exptable']                                        #### we'll just use this variable inside this routine

    n_pix_per_segment = 16384 
    n_segments = 2 
    if ('NUV' in screen): 
       n_segments = 3 
       n_pix_per_segment = 1274  

    print 'Screen = ', screen, ' so n_segments = ', n_segments, 'and n_pix_per_segment = ', n_pix_per_segment

    number_of_exposures = np.size(exptable.keys())

    oldpix = np.empty([ n_segments, n_pix_per_segment, number_of_exposures ])
    newpix = np.empty([ n_segments, n_pix_per_segment, number_of_exposures ])

    if (out_dict.__contains__('G130M_final_all')):
        final_vector_130 = out_dict['G130M_final_all']['WAVE']                        #### this is the "Final" G130 wavelength vector to interpolate to
    if (out_dict.__contains__('G160M_final_all')):
        final_vector_160 = out_dict['G160M_final_all']['WAVE']                        #### this is the "Final" G160 wavelength vector to interpolate to
    if (out_dict.__contains__('G140L_final_all')):
        final_vector_140 = out_dict['G140L_final_all']['WAVE']                        #### this is the "Final" G160 wavelength vector to interpolate to

    exp_keys = exptable.keys()

    for nn in np.arange(number_of_exposures):
        exptable[exp_keys[nn]]['NEWPIX'] = exptable[exp_keys[nn]]['FLUX'] * 0.0 +  np.arange(n_pix_per_segment)  #### add NEWPIX vectors to the table for each exposure, no units
        exptable[exp_keys[nn]]['NEWPIX'].unit = ' '

        for i in np.arange(n_segments):                                    #### we may not actually need oldpix
            oldpix[i, :, nn] = np.arange(n_pix_per_segment)

        number_of_segments = np.size(exptable[exp_keys[nn]]['SEGMENT'])

        if (out_dict['h0table'][exp_keys[nn]]['OPT_ELEM'] == 'G130M'):
            for nseg in np.arange(number_of_segments):
                newpix[nseg,:,nn] = np.round(np.interp(exptable[exp_keys[nn]]['WAVELENGTH'][nseg], final_vector_130, np.arange(50000)))         #### G130M 
                for ipix in np.arange(np.size(newpix[nseg,:,nn])-3):
  		    if newpix[nseg,ipix+1,nn] - newpix[nseg,ipix,nn] > 1.: exptable[exp_keys[nn]]['DQ'][nseg][ipix] = 2.**14. 
                exptable[exp_keys[nn]]['NEWPIX'][nseg] = newpix[nseg,:,nn]
        elif (out_dict['h0table'][exp_keys[nn]]['OPT_ELEM'] == 'G160M'):
            for nseg in np.arange(number_of_segments):
                newpix[nseg,:,nn] = np.round(np.interp(exptable[exp_keys[nn]]['WAVELENGTH'][nseg], final_vector_160, np.arange(50000)))         #### G160M 
                for ipix in np.arange(np.size(newpix[nseg,:,nn])-3):
  		    if newpix[nseg,ipix+1,nn] - newpix[nseg,ipix,nn] > 1.: exptable[exp_keys[nn]]['DQ'][nseg][ipix] = 2.**14. 
                exptable[exp_keys[nn]]['NEWPIX'][nseg] = newpix[nseg,:,nn]
        elif (out_dict['h0table'][exp_keys[nn]]['OPT_ELEM'] == 'G140L'):
            for nseg in np.arange(number_of_segments):
                newpix[nseg,:,nn] = np.round(np.interp(exptable[exp_keys[nn]]['WAVELENGTH'][nseg], final_vector_140, np.arange(50000)))         #### G140M 
                for ipix in np.arange(np.size(newpix[nseg,:,nn])-3):
  		    if newpix[nseg,ipix+1,nn] - newpix[nseg,ipix,nn] > 1.: exptable[exp_keys[nn]]['DQ'][nseg][ipix] = 2.**14. 
                exptable[exp_keys[nn]]['NEWPIX'][nseg] = newpix[nseg,:,nn]
        elif (out_dict['h0table'][exp_keys[nn]]['OPT_ELEM'] == 'G185M'):
            for nseg in np.arange(number_of_segments):
                newpix[nseg,:,nn] = np.round(np.interp(exptable[exp_keys[nn]]['WAVELENGTH'][nseg], final_vector_140, np.arange(50000)))         #### G185M 
                for ipix in np.arange(np.size(newpix[nseg,:,nn])-3):
  		    if newpix[nseg,ipix+1,nn] - newpix[nseg,ipix,nn] > 1.: exptable[exp_keys[nn]]['DQ'][nseg][ipix] = 2.**14. 
                exptable[exp_keys[nn]]['NEWPIX'][nseg] = newpix[nseg,:,nn]
        elif (out_dict['h0table'][exp_keys[nn]]['OPT_ELEM'] == 'G225M'):
            for nseg in np.arange(number_of_segments):
                newpix[nseg,:,nn] = np.round(np.interp(exptable[exp_keys[nn]]['WAVELENGTH'][nseg], final_vector_140, np.arange(50000)))         #### G225M 
                for ipix in np.arange(np.size(newpix[nseg,:,nn])-3):
  		    if newpix[nseg,ipix+1,nn] - newpix[nseg,ipix,nn] > 1.: exptable[exp_keys[nn]]['DQ'][nseg][ipix] = 2.**14. 
                exptable[exp_keys[nn]]['NEWPIX'][nseg] = newpix[nseg,:,nn]
        else:
            print 'OOPS GET_PIX_TO_WAVE DID NOT FIND THE GRATING YOU WANT'

    print "STILL NEED TO IMPLEMENT DQ FLAGGING FOR HITCHES, HERE, HERE, AND HERE"




def get_wavelength_grid(grating):

    wavefile = "xaa18189l_disp.fits"
    print "looking for wavelength grid file",wavefile,"......should really be passing in a reference/ dir........"
    if os.path.isfile("/grp/hst/HST_spectro/datapile/"+wavefile):
        a = Table.read('/grp/hst/HST_spectro/datapile/'+wavefile)
    elif os.path.isfile("../../reference/"+wavefile):
        a = Table.read('../../reference/'+wavefile)
    elif os.path.isfile("../../msp-code/reference/"+wavefile):
        a = Table.read('../../msp-code/reference/'+wavefile)

    if (grating == 'G130M' or grating == 'G160M'):
        a = a[np.where(a['OPT_ELEM'] == grating)]
        a = a[np.where(a['CENWAVE'] > 1250)]
        a = a[np.where(a['APERTURE'] == 'PSA')]

        n_settings = np.size(a['COEFF'][:,1])

        w = np.float64(np.empty([16384, n_settings]))      #### create array for default wavelength solutions

        for i in np.arange(n_settings):
            w[:,i] = (np.arange(16384) - a['D'][i] + a['D_TV03'][i]) * np.float64(a['COEFF'][i,1]) + a['COEFF'][i,0]

        mean_dispersion = np.mean(a['COEFF'][:,1])     ##### this will be the "dispersion" of the final wavelength vector
        nominal_wavelength = np.mean(a['COEFF'][:,1] * (16384/2.) + a['COEFF'][:,0])            ##### this will be the "zeropoint" of the final wavelength vector
        print 'Will use mean dispersion: ', mean_dispersion, '  and nominal wavelengths : ', nominal_wavelength
        final_vector = np.float64(np.arange(50000)-50000/2) * mean_dispersion + nominal_wavelength    ##### create the final wavelength vector
    else:
        a = a[np.where(a['OPT_ELEM'] == grating)]
        a = a[np.where(a['APERTURE'] == 'PSA')]
        a = a[np.where(a['SEGMENT'] == 'FUVA')]

        n_settings = np.size(a['COEFF'][:,1])
        print 'there are N ', n_settings

        w = np.float64(np.empty([16384, n_settings]))   #### create array for default wavelength solutions

        for i in np.arange(n_settings):
            w[:,i] = (np.arange(16384) - a['D'][i] + a['D_TV03'][i]) * np.float64(a['COEFF'][i,1]) + a['COEFF'][i,0]
            print 'tell ', a['COEFF'][i,0], a['COEFF'][i,1], a['COEFF'][i,2]

        mean_dispersion = np.mean(a['COEFF'][:,1])    ##### this will be the "dispersion" of the final wavelength vector
        nominal_wavelength = np.mean(a['COEFF'][:,1] * (16384/2.) + a['COEFF'][:,0])            ##### this will be the "zeropoint" of the final wavelength vector
        print 'Will use mean dispersion: ', mean_dispersion, '  and nominal wavelengths : ', nominal_wavelength
        final_vector = np.float64(np.arange(50000)-50000/2) * mean_dispersion + nominal_wavelength    ##### create the final wavelength vector

    return final_vector




def get_interp_exptable(out_dict): 

    #### this takes in the exposure table and outputs one in the same structure but with the output wavelength grid

    exptable = out_dict['exptable']       		#### NOTE: exptable[i][0] refers to segment A, exptable[i][1] refers to segment B
    exp_keys = exptable.keys()         #### such that exptable[i][0]['WAVELENGTH'] is the segment A wavelength vector and so on
    number_of_exposures = np.size(exptable.keys())
    print "GET_INTERP_EXPTABLE: There are ", number_of_exposures, ' exposures to interpolate into exp_int_table' 

    exp_interp_table = copy.deepcopy(exptable) 

    for exp in exptable.keys(): 

       number_of_segments = np.size(exptable[exp]['SEGMENT'])

       if number_of_segments == 2: 
           blankarr = np.zeros((2,50000)) 
       elif number_of_segments == 1: 
           blankarr = np.zeros((1,50000)) 
       else: 
           print 'BAD NEWS, I DO NOT KNOW HOW TO HANDLE ', number_of_segments, ' SEGMENTS IN GET_INTERP_EXPTABLE' 
          

       print 'GET_INTERP_EXPTABLE: Deleting table entries:'  
       keys_to_delete = ['WAVELENGTH', 'FLUX', 'ERROR', 'GROSS', 'GCOUNTS', 'NET', 'BACKGROUND', 'DQ', 'DQ_WGT', \
                         'BACKGROUND_PER_ROW', 'NUM_EXTRACT_ROWS', 'ACTUAL_EE', 'Y_LOWER_OUTER', \
                         'Y_UPPER_OUTER', 'Y_LOWER_INNER', 'Y_UPPER_INNER', 'EXP_PIX', 'DQ_WEIGHT', 'FLUXFACTOR', 'NEWPIX']
 
       for key_here in keys_to_delete: 
          if key_here in exp_interp_table[exp].keys(): del exp_interp_table[exp][key_here] 

       exp_interp_table[exp]['WAVELENGTH'] = blankarr 
       exp_interp_table[exp]['WAVELENGTH'].unit = 'Angstrom' 
       exp_interp_table[exp]['FLUX'] = blankarr 
       exp_interp_table[exp]['FLUX'].unit = 'erg / (Angstrom cm2 s)' 
       exp_interp_table[exp]['ERROR'] = blankarr 
       exp_interp_table[exp]['ERROR'].unit = 'erg / (Angstrom cm2 s)' 
       exp_interp_table[exp]['GROSS'] = blankarr 
       exp_interp_table[exp]['GROSS'].unit = 'ct / s' 
       exp_interp_table[exp]['GCOUNTS'] = blankarr 
       exp_interp_table[exp]['GCOUNTS'].unit = 'ct'
       exp_interp_table[exp]['BACK'] = blankarr 
       exp_interp_table[exp]['BACK'].unit = 'ct / s' 
       exp_interp_table[exp]['NET'] = blankarr 
       exp_interp_table[exp]['NET'].unit = 'ct / s' 
       exp_interp_table[exp]['NET_ERR_UP'] = blankarr 
       exp_interp_table[exp]['NET_ERR_UP'].unit = 'ct / s' 
       exp_interp_table[exp]['NET_ERR_DOWN'] = blankarr 
       exp_interp_table[exp]['NET_ERR_DOWN'].unit = 'ct / s' 
       exp_interp_table[exp]['BACKGROUND'] = blankarr 
       exp_interp_table[exp]['BACKGROUND'].unit = 'ct / s' 
       exp_interp_table[exp]['DQ'] = blankarr 
       exp_interp_table[exp]['DQ_WGT'] = blankarr 
       exp_interp_table[exp]['EXP_PIX'] = blankarr 
       exp_interp_table[exp]['EXP_PIX'].unit = 'ct / s' 
       exp_interp_table[exp]['DQ_WEIGHT'] = blankarr 
       exp_interp_table[exp]['FLUXFACTOR'] = blankarr 

    out_dict['exp_interp_table'] = exp_interp_table 






def exposure_interpolate(out_dict): 

    #### here we are going to go through the exptable and "interpolate" the exposures onto the exp_interp_table 
    #### this should be done carefully so as to preserve DQ and such. 

    exptable = out_dict['exptable']       		   #### NOTE: exptable[i][0] refers to segment A, exptable[i][1] refers to segment B
    exp_interp_table = out_dict['exp_interp_table']        #### the interpolated arrays 

    number_of_exposures = np.size(exptable.keys())
    print "EXPOSURE_INTERPOLATE: There are ", number_of_exposures, ' exposures to interpolate'

    exp_keys = exptable.keys()                          

    exp_counter = 0 
    for i in np.arange(number_of_exposures):

        number_of_segments = np.size(exptable[exp_keys[i]]['SEGMENT'])

        for nseg in np.arange(number_of_segments):                #### loop over segments
            newpix = out_dict['exptable'][exp_keys[i]]['NEWPIX'][nseg,:]

            grating_here = exptable[exp_keys[i]][nseg]['GRATING']
            exp_interp_table[exp_keys[i]][nseg]['WAVELENGTH'] = out_dict[grating_here+'_final_all']['WAVE'] 

	    exp_interp_table[exp_keys[i]][nseg]['FLUX'][newpix.astype(int)] = exp_interp_table[exp_keys[i]][nseg]['FLUX'][newpix.astype(int)] + exptable[exp_keys[i]][nseg]['FLUX'] 
	    exp_interp_table[exp_keys[i]][nseg]['EXP_PIX'][newpix.astype(int)] = exp_interp_table[exp_keys[i]][nseg]['EXP_PIX'][newpix.astype(int)] + exptable[exp_keys[i]][nseg]['EXP_PIX'] 
            exp_interp_table[exp_keys[i]][nseg]['ERROR'][newpix.astype(int)] = exp_interp_table[exp_keys[i]][nseg]['ERROR'][newpix.astype(int)] + exptable[exp_keys[i]][nseg]['ERROR'] 
            exp_interp_table[exp_keys[i]][nseg]['GROSS'][newpix.astype(int)] = exp_interp_table[exp_keys[i]][nseg]['GROSS'][newpix.astype(int)] + exptable[exp_keys[i]][nseg]['GROSS'] 
            exp_interp_table[exp_keys[i]][nseg]['NET'][newpix.astype(int)] = exp_interp_table[exp_keys[i]][nseg]['NET'][newpix.astype(int)] + exptable[exp_keys[i]][nseg]['NET'] 
            exp_interp_table[exp_keys[i]][nseg]['BACK'][newpix.astype(int)] = exp_interp_table[exp_keys[i]][nseg]['BACK'][newpix.astype(int)] + exptable[exp_keys[i]][nseg]['BACKGROUND'] 

            counts_err_down, counts_err_up = counts_error.counts_error(exp_interp_table[exp_keys[i]][nseg]['NET'] * exp_interp_table[exp_keys[i]][nseg]['EXP_PIX']) / exp_interp_table[exp_keys[i]][nseg]['EXP_PIX'] 
            exp_interp_table[exp_keys[i]][nseg]['NET_ERR_UP'][newpix.astype(int)] = counts_err_up[newpix.astype(int)] 
            exp_interp_table[exp_keys[i]][nseg]['NET_ERR_DOWN'][newpix.astype(int)] = counts_err_down[newpix.astype(int)] 

            exp_interp_table[exp_keys[i]][nseg]['DQ'][newpix.astype(int)] = exptable[exp_keys[i]][nseg]['DQ'] # this may not be right unless it's done bitwise 
            exp_interp_table[exp_keys[i]][nseg]['FLUXFACTOR'][newpix.astype(int)] = exp_interp_table[exp_keys[i]][nseg]['FLUXFACTOR'][newpix.astype(int)] + exptable[exp_keys[i]][nseg]['FLUXFACTOR'] 
 
        exp_interp_table[exp_keys[i]]['ROOTNAME'] = out_dict['exptable'][exp_keys[i]]['ROOTNAME'] 

        exp_counter = exp_counter + 1                        #### increment the number of exposures we've added

    return out_dict






def cos_counts_coadd(out_dict, key_to_combine): 

    #### this preserves the functionality of cos_counts_coadd but operates on exp_interp_table since they have already been interpolated to the new grid 
    #### this is the new official cos_counts_coadd as of 111115 

    if(out_dict.__contains__(key_to_combine)):
        final_table = out_dict[key_to_combine]

    exp_interp_table = out_dict['exp_interp_table']       #### NOTE: exptable[i][0] refers to segment A, exptable[i][1] refers to segment B
    exp_keys = exp_interp_table.keys()         
    h0 = out_dict['h0table']
    number_of_exposures = np.size(exp_keys) 
    print "COS_COUNTS_COADD_REMIX: There are ", number_of_exposures, ' exposures to coadd for key_to_combine', key_to_combine
    print exp_keys 

    exp_counter = 0

    for i in np.arange(number_of_exposures):

        number_of_segments = np.size(exp_interp_table[exp_keys[i]]['SEGMENT'])

        if h0[exp_keys[i]]['OPT_ELEM'] in key_to_combine:

            print 'Stacking up (remix): ', h0[exp_keys[i]]['ROOTNAME'], h0[exp_keys[i]]['OPT_ELEM'], key_to_combine

            for nseg in np.arange(number_of_segments):                #### loop over segments

                final_table['EXP_PIX'] = final_table['EXP_PIX'] + exp_interp_table[exp_keys[i]][nseg]['EXP_PIX']

                final_table['GROSSCOUNTS'] = final_table['GROSSCOUNTS'] + exp_interp_table[exp_keys[i]][nseg]['GROSS'] * exp_interp_table[exp_keys[i]][nseg]['EXP_PIX']
                            
                final_table['NETCOUNTS'] = final_table['NETCOUNTS'] + exp_interp_table[exp_keys[i]][nseg]['NET'] * exp_interp_table[exp_keys[i]][nseg]['EXP_PIX']
                            
                final_table['BACK'] = final_table['BACK'] + exp_interp_table[exp_keys[i]][nseg]['BACK'] * exp_interp_table[exp_keys[i]][nseg]['EXP_PIX']

                final_table['FLUXWGT'] = final_table['FLUXWGT'] + exp_interp_table[exp_keys[i]][nseg]['FLUX'] * exp_interp_table[exp_keys[i]][nseg]['EXP_PIX']

                final_table['FLUXWGT_ERR'] = final_table['FLUXWGT_ERR'] + (exp_interp_table[exp_keys[i]][nseg]['EXP_PIX'] * exp_interp_table[exp_keys[i]][nseg]['ERROR'])**2 
                            
                final_table['DQ'] = (final_table['DQ']).astype(int) | (exp_interp_table[exp_keys[i]][nseg]['DQ']).astype(int) 
     
                final_table['FLUXFACTOR'] = final_table['FLUXFACTOR'] + exp_interp_table[exp_keys[i]][nseg]['FLUXFACTOR']*exp_interp_table[exp_keys[i]][nseg]['EXP_PIX']

            exp_counter = exp_counter + 1                        #### increment the number of exposures we've added
        else:
            print "THE GRATING ", h0[exp_keys[i]]['OPT_ELEM'], "DOES NOT MATCH WITH KEY_TO_COMBINE ", key_to_combine, " FOR EXPOSURE ", h0[exp_keys[i]]['ROOTNAME']

    i_good_pixels = np.where(final_table['EXP_PIX'] > 0.1)

    if (exp_counter > 0):

        final_table['GROSSCOUNTS'][i_good_pixels] = final_table['GROSSCOUNTS'][i_good_pixels] / final_table['EXP_PIX'][i_good_pixels]
        final_table['NETCOUNTS'][i_good_pixels] = final_table['NETCOUNTS'][i_good_pixels] / final_table['EXP_PIX'][i_good_pixels]
        final_table['BACK'][i_good_pixels] = final_table['BACK'][i_good_pixels] / final_table['EXP_PIX'][i_good_pixels]
        final_table['FLUXFACTOR'][i_good_pixels] = final_table['FLUXFACTOR'][i_good_pixels] / final_table['EXP_PIX'][i_good_pixels]

        final_table['FLUX'][i_good_pixels] = final_table['NETCOUNTS'][i_good_pixels] * final_table['FLUXFACTOR'][i_good_pixels]

        # obtain the error vector for the number of counts in each bin from the Gehrels (1986) tables, which are in counts_error.py 
        err_down, err_up = counts_error.counts_error(final_table['NETCOUNTS'][i_good_pixels] * final_table['EXP_PIX'][i_good_pixels]) 
        final_table['NETCOUNTSERR_UP'][i_good_pixels] = err_up / final_table['EXP_PIX'][i_good_pixels] 
        final_table['NETCOUNTSERR_DOWN'][i_good_pixels] = err_down / final_table['EXP_PIX'][i_good_pixels] 
        final_table['FLUXERR_UP'][i_good_pixels] = final_table['NETCOUNTSERR_UP'][i_good_pixels] * final_table['FLUXFACTOR'][i_good_pixels] 
        final_table['FLUXERR_DOWN'][i_good_pixels] = final_table['NETCOUNTSERR_DOWN'][i_good_pixels] * final_table['FLUXFACTOR'][i_good_pixels] 

##### WATCH OT THERE MAY BE A PROBLEM HERE!!!!! 
        print "WATCH OUT THERE MAY BE A PROBLEM HERE!!!!" 
        final_table['ERROR'][i_good_pixels] = final_table['FLUXERR_DOWN'][i_good_pixels] 
        final_table['ERROR'][i_good_pixels] = final_table['FLUXERR_UP'][i_good_pixels] 

        final_table['SN'][i_good_pixels] = final_table['FLUX'][i_good_pixels] / final_table['FLUXERR_UP'][i_good_pixels]        #### S/N per "pixel" or "bin", that is, not combining bins into resels 

        # process the exposure-weighted fluxes and errors from the original x1d "FLUX" and "ERROR" vectors 
        final_table['FLUXWGT'][i_good_pixels] = final_table['FLUXWGT'][i_good_pixels] / final_table['EXP_PIX'][i_good_pixels]     #### exposure-weighted mean
        final_table['FLUXWGT_ERR'][i_good_pixels] = np.sqrt(final_table['FLUXWGT_ERR'][i_good_pixels]) / final_table['EXP_PIX'][i_good_pixels]

    if(out_dict.__contains__(key_to_combine)):    out_dict[key_to_combine] = final_table

    return out_dict








def get_exposure_list(exposure_list):

    catalog = ascii.read(exposure_list)
    return catalog


def get_dict_of_exposures():

    tt = {}                            # this will hold the data
    h0 = {}                            # this will hold the headers [fits extention 0]
    h1 = {}                            # this will hold the headers [fits extention 1]

    if os.path.exists("all_exposures.txt"): 
        dataset_list =  get_exposure_list("all_exposures.txt")['Rootname','Flag']
    else: 
        print "Oh no, all_exposures.txt does not exist for this target, you might want to try something else" 
        return 0, 0, 0 

    if "Flag" in dataset_list.colnames:
        dataset_list.remove_rows(np.where(dataset_list['Flag'] == 0))    #### delete exposures flagged with 0 in the all_exposures.txt file

    if len(dataset_list) == 0: 
        print 'GET_DICT_OF_EXPOSURES: It appears that no exposures were flagged as 1, will try to exit gracefully!' 
        return 0, 0, 0 

    dataset_list = dataset_list['Rootname'] 

    for i in range(np.size(dataset_list)):

        #get the header
        filename = dataset_list[i]+'_x1d.fits.gz'
        h0[i] = fitsio.read_header(filename, 0)
        h1[i] = fitsio.read_header(filename, 1)

        t_new = Table.read(dataset_list[i]+'_x1d.fits.gz')
        #print 'GET_DICT_OF_EXPOSURES: We have read in ', dataset_list[i]+'_x1d.fits  for targname  ', h0[0]['TARGNAME']
        tt[i] = t_new
        tt[i]['ROOTNAME'] = dataset_list[i]

    print 'GET_DICT_OF_EXPOSURES: ', dataset_list
    print 'GET_DICT_OF_EXPOSURES: We have read in i = ', i+1, '  exposures '
    
    for key in h0.keys(): 
       h0[key]['OPT_ELEM'] = str(h0[key]['OPT_ELEM']).strip() 
       h0[key]['TARGNAME'] = str(h0[key]['TARGNAME']).strip() 
       h0[key]['ROOTNAME'] = str(h0[key]['ROOTNAME']).strip() 

    return tt, h0, h1


def screen_dict_of_exposures(tt, h0, h1, screen='FUVM', lifetime='ALL'):         # screen out the dictionary entries we don't want
    print 'SCREEN_DICT_OF_EXPOSURES: The requested screen is: ', screen
    print 'SCREEN_DICT_OF_EXPOSURES: The requested LP is: ', lifetime

    if (screen == 'FUVM' and lifetime == 'ALL'):                # for FUVM and all LPs
         print '1' 
         for i in h0.keys():
            if (h0[i]['DETECTOR'] == 'FUV' and (h0[i]['OPT_ELEM'] == 'G130M' or h0[i]['OPT_ELEM'] == 'G160M') and h0[i]['CENWAVE'] > 1250):
                print 'SCREEN_DICT_OF_EXPOSURES 1: Keeping dictionary entry for ', h0[i]['ROOTNAME'], h0[i]['DETECTOR'], h0[i]['OPT_ELEM'], h0[i]['CENWAVE'], h0[i]['LIFE_ADJ']
            else:
                print 'SCREEN_DICT_OF_EXPOSURES 1: Deleting dictionary entry for ', h0[i]['ROOTNAME'], h0[i]['DETECTOR'], h0[i]['OPT_ELEM'], h0[i]['CENWAVE'], h0[i]['LIFE_ADJ']
                del tt[i]
                del h0[i]
                del h1[i]
    elif (screen == 'FUVM' and lifetime < 10):                     # for FUVM and a specific LP
         print '2' 
         for i in h0.keys():
            if (h0[i]['DETECTOR'] == 'FUV' and (h0[i]['OPT_ELEM'] == 'G130M' or h0[i]['OPT_ELEM'] == 'G160M') and h0[i]['CENWAVE'] > 1250 and h0[i]['LIFE_ADJ'] == lifetime):
                print 'SCREEN_DICT_OF_EXPOSURES 2: Keeping dictionary entry for ', h0[i]['ROOTNAME'], h0[i]['DETECTOR'], h0[i]['OPT_ELEM'], h0[i]['CENWAVE'], h0[i]['LIFE_ADJ']
            else:
                print 'SCREEN_DICT_OF_EXPOSURES 2: Deleting dictionary entry for ', h0[i]['ROOTNAME'], h0[i]['DETECTOR'], h0[i]['OPT_ELEM'], h0[i]['CENWAVE'], h0[i]['LIFE_ADJ']
                del tt[i]
                del h0[i]
                del h1[i]
    elif (screen == 'FUVM' and lifetime == 12):				# for FUVM and 1+2 combo 
         print '3' 
         for i in h0.keys():
            if (h0[i]['DETECTOR'] == 'FUV' and (h0[i]['OPT_ELEM'] == 'G130M' or h0[i]['OPT_ELEM'] == 'G160M') and h0[i]['CENWAVE'] > 1250 and h0[i]['LIFE_ADJ'] < 3 and h0[i]['LIFE_ADJ'] > 0):
                print 'SCREEN_DICT_OF_EXPOSURES 2: Keeping dictionary entry for ', h0[i]['ROOTNAME'], h0[i]['DETECTOR'], h0[i]['OPT_ELEM'], h0[i]['CENWAVE'], h0[i]['LIFE_ADJ']
            else:
                print 'SCREEN_DICT_OF_EXPOSURES 2: Deleting dictionary entry for ', h0[i]['ROOTNAME'], h0[i]['DETECTOR'], h0[i]['OPT_ELEM'], h0[i]['CENWAVE'], h0[i]['LIFE_ADJ']
                del tt[i]
                del h0[i]
                del h1[i]
    elif (screen == 'FUVL' and lifetime == 'ALL'):                     # for FUVL and a specific LP
         print '4' 
         for i in h0.keys():
            if (h0[i]['DETECTOR'] == 'FUV' and h0[i]['OPT_ELEM'] == 'G140L'):
                print 'SCREEN_DICT_OF_EXPOSURES 2: Keeping dictionary entry for ', h0[i]['ROOTNAME'], h0[i]['DETECTOR'], h0[i]['OPT_ELEM'], h0[i]['CENWAVE'], h0[i]['LIFE_ADJ']
            else:
                print 'SCREEN_DICT_OF_EXPOSURES 2: Deleting dictionary entry for ', h0[i]['ROOTNAME'], h0[i]['DETECTOR'], h0[i]['OPT_ELEM'], h0[i]['CENWAVE'], h0[i]['LIFE_ADJ']
                del tt[i]
                del h0[i]
                del h1[i]
    elif (screen == 'FUVL'):                        # for FUVL 
         print '5' 
         for i in h0.keys():
            if (h0[i]['DETECTOR'] == 'FUV' and h0[i]['OPT_ELEM'] == 'G140L' and h0[i]['LIFE_ADJ'] == lifetime):
                print 'SCREEN_DICT_OF_EXPOSURES 3: Keeping dictionary entry for ', h0[i]['ROOTNAME'], h0[i]['DETECTOR'], h0[i]['OPT_ELEM'], h0[i]['CENWAVE'], h0[i]['LIFE_ADJ']
            else:
                print 'SCREEN_DICT_OF_EXPOSURES 3: Deleting dictionary entry for ', h0[i]['ROOTNAME'], h0[i]['DETECTOR'], h0[i]['OPT_ELEM'], h0[i]['CENWAVE'], h0[i]['LIFE_ADJ']
                del tt[i]
                del h0[i]
                del h1[i]
    elif (screen == 'FUVL'):                        # for FUVL and combo of 1+2 
         print '6' 
         for i in h0.keys():
            if (h0[i]['DETECTOR'] == 'FUV' and h0[i]['OPT_ELEM'] == 'G140L' and h0[i]['LIFE_ADJ'] < 3):
                print 'SCREEN_DICT_OF_EXPOSURES 3: Keeping dictionary entry for ', h0[i]['ROOTNAME'], h0[i]['DETECTOR'], h0[i]['OPT_ELEM'], h0[i]['CENWAVE'], h0[i]['LIFE_ADJ']
            else:
                print 'SCREEN_DICT_OF_EXPOSURES 3: Deleting dictionary entry for ', h0[i]['ROOTNAME'], h0[i]['DETECTOR'], h0[i]['OPT_ELEM'], h0[i]['CENWAVE'], h0[i]['LIFE_ADJ']
                del tt[i]
                del h0[i]
                del h1[i]
    elif (screen == 'NUVM' and lifetime == 'ALL'):                        # for NUVM 
         print '7' 
         for i in h0.keys():
            if (h0[i]['DETECTOR'] == 'NUV' and (h0[i]['OPT_ELEM'] == 'G185M' or h0[i]['OPT_ELEM'] == 'G225M')): 
                print 'SCREEN_DICT_OF_EXPOSURES 3: Keeping dictionary entry for ', h0[i]['ROOTNAME'], h0[i]['DETECTOR'], h0[i]['OPT_ELEM'], h0[i]['CENWAVE'], h0[i]['LIFE_ADJ']
            else:
                print 'SCREEN_DICT_OF_EXPOSURES 3: Deleting dictionary entry for ', h0[i]['ROOTNAME'], h0[i]['DETECTOR'], h0[i]['OPT_ELEM'], h0[i]['CENWAVE'], h0[i]['LIFE_ADJ']
                del tt[i]
                del h0[i]
    else:
        print "SORRY BUDDY, BUT I CAN'T DO THAT SCREEN:", screen, " YOUR CHOICE IS INVALID"

    return tt, h0, h1




def write_output(out_dict, Linda):

    ##### note this currently does not write out the fits header, prob requires astropy.io.fits
    if (out_dict.__contains__('G130M_final_all')): write_fits_output(out_dict,'G130M_final_all')
    if (out_dict.__contains__('G160M_final_all')): write_fits_output(out_dict,'G160M_final_all')
    if (out_dict.__contains__('G140L_final_all')): write_fits_output(out_dict,'G140L_final_all') 

    ##### note this currently does not write out the fits header, prob requires astropy.io.fits
    #if (out_dict.__contains__('G130M_final_shifted')): write_fits_output(out_dict,'G130M_final_shifted')

    #### LP1 
    if (out_dict.__contains__('G130M_final_lp1')): write_fits_output(out_dict,'G130M_final_lp1')
    if (out_dict.__contains__('G160M_final_lp1')): write_fits_output(out_dict,'G160M_final_lp1')
    if (out_dict.__contains__('G140L_final_lp1')): write_fits_output(out_dict,'G140L_final_lp1')

    #### LP2 
    if (out_dict.__contains__('G130M_final_lp2')): write_fits_output(out_dict,'G130M_final_lp2')
    if (out_dict.__contains__('G160M_final_lp2')): write_fits_output(out_dict,'G160M_final_lp2')
    if (out_dict.__contains__('G140L_final_lp2')): write_fits_output(out_dict,'G140L_final_lp2')

    #### LP3 
    if (out_dict.__contains__('G130M_final_lp3')): write_fits_output(out_dict,'G130M_final_lp3')
    if (out_dict.__contains__('G160M_final_lp3')): write_fits_output(out_dict,'G160M_final_lp3')
    if (out_dict.__contains__('G140L_final_lp3')): write_fits_output(out_dict,'G140L_final_lp3')

    #### LP1+2 
    if (out_dict.__contains__('G130M_final_lp12')): write_fits_output(out_dict,'G130M_final_lp12')
    if (out_dict.__contains__('G160M_final_lp12')): write_fits_output(out_dict,'G160M_final_lp12')
    if (out_dict.__contains__('G140L_final_lp12')): write_fits_output(out_dict,'G140L_final_lp12')

    #### FUV SPLICE 
    if (out_dict.__contains__('FUVM_final_all')): write_fits_output(out_dict,'FUVM_final_all')
    if (out_dict.__contains__('FUVM_final_lp1')): write_fits_output(out_dict,'FUVM_final_lp1')
    if (out_dict.__contains__('FUVM_final_lp2')): write_fits_output(out_dict,'FUVM_final_lp2')
    if (out_dict.__contains__('FUVM_final_lp3')): write_fits_output(out_dict,'FUVM_final_lp3')
    if (out_dict.__contains__('FUVM_final_lp12')): write_fits_output(out_dict,'FUVM_final_lp12')

    if (out_dict['lifetime'] == 1): print 'HELL YES ITS LP1', out_dict['lifetime'] 
    if (out_dict['lifetime'] == 2): print 'HELL YES ITS LP2', out_dict['lifetime'] 
    if (out_dict['lifetime'] == 3): print 'HELL YES ITS LP3', out_dict['lifetime'] 
    if (out_dict['lifetime'] == 12): print 'HELL YES ITS LP12', out_dict['lifetime'] 



    #### THIS IS A KLUDGE TO COPY THE SINGLE-GRATING FUV COADD INTO THE FUVM FILE IF ONLY ONE EXISTS 
    if (out_dict.__contains__('G130M_final_lp1') and not out_dict.__contains__('G160M_final_lp1')): 
        filename1 = out_dict['targname']+'_coadd_G130M_final_all.fits.gz'
        filename2 = out_dict['targname']+'_coadd_FUVM_final_all.fits.gz'
        command = 'cp -rp '+filename1+' '+filename2 
        os.system(command) 
        print 'DAMMIT I HAVE G130M BUT NOT G160M', command 
    if (out_dict.__contains__('G160M_final_lp1') and not out_dict.__contains__('G130M_final_lp1')):  
        filename1 = out_dict['targname']+'_coadd_G160M_final_all.fits.gz'
        filename2 = out_dict['targname']+'_coadd_FUVM_final_all.fits.gz'
        command = 'cp -rp '+filename1+' '+filename2 
        os.system(command) 
        print 'DAMMIT I HAVE G160M BUT NOT G130M', command 
    

    #### SPECIAL PROCESSING 
    if (Linda): 
        print ' this is for Linda' 
        if (out_dict.__contains__('FUVM_final_all')): 
            copyt = out_dict['FUVM_final_all']['WAVE','FLUX','ERROR']
            ascii.write(copyt, out_dict['targname'] + '_coadd_FUVM_final_all.txt') 
        if (out_dict.__contains__('G130M_final_all')): 
            copyt = out_dict['G130M_final_all']['WAVE','FLUX','ERROR']
            ascii.write(copyt, out_dict['targname'] + '_coadd_G130M_final_all.txt') 
        if (out_dict.__contains__('G160M_final_all')): 
            copyt = out_dict['G160M_final_all']['WAVE','FLUX','ERROR']
            ascii.write(copyt, out_dict['targname'] + '_coadd_G160M_final_all.txt') 
        if (out_dict.__contains__('G140L_final_all')): 
            copyt = out_dict['G140L_final_all']['WAVE','FLUX','ERROR']
            ascii.write(copyt, out_dict['targname'] + '_coadd_G140L_final_all.txt') 





def write_fits_output(out_dict, key_to_write_out):
    # can eventually expand this to include multiple extendsion for each LP instead of keeping those in separate files

    table_to_write_out = out_dict[key_to_write_out]
    prihdr = f.Header()
    prihdr['CREATOR'] = '50MP'                     # OK, here's how we can add keywords to the primary header
    prihdr['TIMESTAMP'] = str(datetime.now())
    prihdr['LEGACY'] = 'HST-50MP'
    prihdu = f.PrimaryHDU(header=prihdr)

    table_hdu = f.BinTableHDU.from_columns(np.array(table_to_write_out.filled()))
    table_hdu.header['LP'] = 'ALL'                    # this is how to add header keywords to the table extension!
    if 'lp1' in key_to_write_out: table_hdu.header['LP'] = 'LP1'
    if 'lp2' in key_to_write_out: table_hdu.header['LP'] = 'LP2'
    if 'lp3' in key_to_write_out: table_hdu.header['LP'] = 'LP3'

    thdulist = f.HDUList([prihdu, table_hdu])

    thdulist.writeto(out_dict['targname']+'_coadd_'+key_to_write_out+'.fits.gz',clobber=True)








