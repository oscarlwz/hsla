from astropy.io import fits as f
#JRD
#import fitsio
from astropy.io import fits
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
        #JRD python 3
        print("It appears there is no all_exposures.txt file in this directory, attemping to exit gracefully." )
        return 0  

    #JRD python3
    targname = h0[list(h0.keys())[0]]['TARGNAME']
    #JRD python 3
    print('ADD.MAIN: We will now perform a coadd Target ', targname)

    #######################################################################################################
    #               STEP 2: screen out the exposures that we don't want to have                           #
    #######################################################################################################
    #JRD
    print('We will perform coadds for lifetime position : ', lifetime)
    tt, h0, h1 = screen_dict_of_exposures(tt, h0, h1, screen=screen, lifetime=lifetime)
    number_of_remaining_exposures = np.size(list(h0.keys()))
    #JRD
    print('ADD.MAIN: Exposure cuts are: ', screen, ' and ', lifetime)
    print('ADD.MAIN: Number of exposures surviving the cuts: ', number_of_remaining_exposures)

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

    for nn in list(tt.keys()):                                    #### Add a few columns to the exposure tables to include per-pixel exptime and DQ weights
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

    # This is where the "final_all"  and "final_lp"-specific dictionary entries are created 
    if (out_dict['number_of_G130M_exposures'] > 0):
        final_wave_130 = get_wavelength_grid('G130M')
        #out_dict['G130M_final_all'] = create_coadd_format(final_wave_130)
        out_dict['G130M_final_lp'+lifetime] = create_coadd_format(final_wave_130)

    if (out_dict['number_of_G160M_exposures'] > 0):
        final_wave_160 = get_wavelength_grid('G160M')
        #out_dict['G160M_final_all'] = create_coadd_format(final_wave_160)
        out_dict['G160M_final_lp'+lifetime] = create_coadd_format(final_wave_160)

    if (out_dict['number_of_G140L_exposures'] > 0):
        final_wave_140 = get_wavelength_grid('G140L')
        #out_dict['G140L_final_all'] = create_coadd_format(final_wave_140)
        out_dict['G140L_final_lp'+lifetime] = create_coadd_format(final_wave_140)

    #JRD changing "all" to LP1 for NUV
    if (out_dict['number_of_G185M_exposures'] > 0):
        final_wave_185 = get_wavelength_grid('G185M')
        out_dict['G185M_final_lp1'] = create_coadd_format(final_wave_185)

    if (out_dict['number_of_G225M_exposures'] > 0):
        final_wave_225 = get_wavelength_grid('G225M')
        out_dict['G225M_final_lp1'] = create_coadd_format(final_wave_225)

    #JRD adding G285M
    if (out_dict['number_of_G285M_exposures'] > 0):
        final_wave_285 = get_wavelength_grid('G285M')
        out_dict['G285M_final_lp1'] = create_coadd_format(final_wave_285)
        
    #JRD adding G230L
    if (out_dict['number_of_G230L_exposures'] > 0):
        final_wave_230 = get_wavelength_grid('G230L')
        out_dict['G230L_final_lp1'] = create_coadd_format(final_wave_230)
        

    #######################################################################################################
    #               STEP 6: obtain the new pixel coordinates of the pixels in each x1d epxosure          #
    #######################################################################################################

    for nn in list(tt.keys()):                                    ##### Add a few column to the exposure table for the "newpix"

        if ('FUV' in screen): tt[nn]['NEWPIX'] = tt[nn]['FLUX'] * 0.0 +  np.arange(16384)
        #JRD uncommenting to include NUV
        if ('NUV' in screen): tt[nn]['NEWPIX'] = tt[nn]['FLUX'] * 0.0 +  np.arange(1274)

        tt[nn]['NEWPIX'].unit = ' '

    a = get_pix_to_wave(out_dict, screen, lifetime)


    #######################################################################################################
    #               STEP 7: replicate the structure of the exposure table and interpolate to it.          # 
    #######################################################################################################

    get_interp_exptable(out_dict, lifetime) 			     #### this creates an exposure table paralleling "exptable" 
							     #### but with the new wavelength grid instead of the old one 

    #############################################################################################################
    #               STEP 8: "interpolate" the input exposures to the output grid, populate exp_interp_table    # 
    #############################################################################################################

    exposure_interpolate(out_dict, lifetime)			     #### this interpolates the input exposures (exptable) 
							     #### onto the new wavelength grid (exp_interp_table) 

    #######################################################################################################
    #               STEPS 9: go ahead and do the final coadd                                              #
    #######################################################################################################

    #JRD PYTHON 3
    print("Passing ", np.size(list(tt.keys())), ' exposures to cos_counts_coadd at LP = ALL ')
    number_of_remaining_exposures = np.size(h0.keys())
    total_number_of_exposures = number_of_remaining_exposures
    print('Prior to calling cos_counts_coadd, the exposure counts are: ')
    print('       G130M: ', out_dict['number_of_G130M_exposures'])
    print('       G160M: ', out_dict['number_of_G160M_exposures'])
    print('       G140L: ', out_dict['number_of_G140L_exposures'])
    print('       G185M: ', out_dict['number_of_G185M_exposures'])
    print('       G225M: ', out_dict['number_of_G225M_exposures'])
    #JRD ADDING THE OTHER NUV SETTINGS
    print('       G285M: ', out_dict['number_of_G285M_exposures'])
    print('       G230L: ', out_dict['number_of_G230L_exposures'])


    if (out_dict['number_of_G130M_exposures'] > 0):
        #out_dict = cos_counts_coadd(out_dict, 'G130M_final_all')
        out_dict = cos_counts_coadd(out_dict, 'G130M_final_lp'+lifetime)
    if (out_dict['number_of_G160M_exposures'] > 0):
        #out_dict = cos_counts_coadd(out_dict, 'G160M_final_all')
        out_dict = cos_counts_coadd(out_dict, 'G160M_final_lp'+lifetime)
    if (out_dict['number_of_G140L_exposures'] > 0):
        #out_dict = cos_counts_coadd(out_dict, 'G140L_final_all')
        out_dict = cos_counts_coadd(out_dict, 'G140L_final_lp'+lifetime)
    if (out_dict['number_of_G185M_exposures'] > 0):
        out_dict = cos_counts_coadd(out_dict, 'G185M_final_lp1')
    if (out_dict['number_of_G225M_exposures'] > 0):
        out_dict = cos_counts_coadd(out_dict, 'G225M_final_lp1')
    #JRD ADDING OTHER NUV SETTINGS
    if (out_dict['number_of_G285M_exposures'] > 0):
        out_dict = cos_counts_coadd(out_dict, 'G285M_final_lp1')
    if (out_dict['number_of_G230L_exposures'] > 0):
        out_dict = cos_counts_coadd(out_dict, 'G230L_final_lp1')


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

    #JRD COMMENTING OUT BECAUSE THE FINAL_ALL ACTUALLY ONLY CONTAINS THE GIVEN LIFETIME AND IS IDENTICAL TO FINAL_LP
    #if (out_dict['number_of_G130M_exposures'] > 0):
    #    if (out_dict.__contains__('G130M_final_all')): obtain_chi2(out_dict, 'G130M_final_all') 
    #if (out_dict['number_of_G160M_exposures'] > 0):
    #    if (out_dict.__contains__('G160M_final_all')): obtain_chi2(out_dict, 'G160M_final_all') 

    #JRD ADD NUV and G140L
    #if (out_dict['number_of_G140L_exposures'] > 0):
    #    if (out_dict.__contains__('G140L_final_all')): obtain_chi2(out_dict, 'G140L_final_all') 

    if (out_dict['number_of_G185M_exposures'] > 0):
        if (out_dict.__contains__('G185M_final_lp1')): obtain_chi2(out_dict, 'G185M_final_lp1') 

    if (out_dict['number_of_G225M_exposures'] > 0):
        if (out_dict.__contains__('G225M_final_lp1')): obtain_chi2(out_dict, 'G225M_final_lp1')

    if (out_dict['number_of_G285M_exposures'] > 0):
        if (out_dict.__contains__('G285M_final_lp1')): obtain_chi2(out_dict, 'G285M_final_lp1')

    if (out_dict['number_of_G230L_exposures'] > 0):
        if (out_dict.__contains__('G230L_final_lp1')): obtain_chi2(out_dict, 'G230L_final_lp1') 

    #JRD REPLACING WITH LIFETIME
    if (out_dict.__contains__('G130M_final_lp'+lifetime)): obtain_chi2(out_dict, 'G130M_final_lp'+lifetime)
    if (out_dict.__contains__('G160M_final_lp'+lifetime)): obtain_chi2(out_dict, 'G160M_final_lp'+lifetime)
    if (out_dict.__contains__('G140L_final_lp'+lifetime)): obtain_chi2(out_dict, 'G140L_final_lp'+lifetime)
    
    #if (out_dict.__contains__('G130M_final_lp1')): obtain_chi2(out_dict, 'G130M_final_lp1') 
    #if (out_dict.__contains__('G160M_final_lp1')): obtain_chi2(out_dict, 'G160M_final_lp1') 
 
    #if (out_dict.__contains__('G130M_final_lp2')): obtain_chi2(out_dict, 'G130M_final_lp2') 
    #if (out_dict.__contains__('G160M_final_lp2')): obtain_chi2(out_dict, 'G160M_final_lp2') 
 
    #if (out_dict.__contains__('G130M_final_lp3')): obtain_chi2(out_dict, 'G130M_final_lp3') 
    #if (out_dict.__contains__('G160M_final_lp3')): obtain_chi2(out_dict, 'G160M_final_lp3')

    #if (out_dict.__contains__('G130M_final_lp4')): obtain_chi2(out_dict, 'G130M_final_lp4') 
    #if (out_dict.__contains__('G160M_final_lp4')): obtain_chi2(out_dict, 'G160M_final_lp4') 
 
    #if (out_dict.__contains__('G130M_final_lp12')): obtain_chi2(out_dict, 'G130M_final_lp12') 
    #if (out_dict.__contains__('G160M_final_lp12')): obtain_chi2(out_dict, 'G160M_final_lp12') 

    #######################################################################################################
    #               STEP 13:  if you have both G130M and and G160M, splice them and store output          #
    #######################################################################################################

    #JRD comment out since final_all is bogus. Need to call add.main("FUVM", "ALL")
    #if out_dict.__contains__('G130M_final_all') and out_dict.__contains__('G160M_final_all'):  
    #    splice_fuv(out_dict, 'all')
    if out_dict.__contains__('G130M_final_lp'+lifetime) and out_dict.__contains__('G160M_final_lp'+lifetime):  
        splice_fuv(out_dict, 'lp'+lifetime) 
    #if out_dict.__contains__('G130M_final_lp1') and out_dict.__contains__('G160M_final_lp1'):  
    #    splice_fuv(out_dict, 'lp1') 
    #if out_dict.__contains__('G130M_final_lp2') and out_dict.__contains__('G160M_final_lp2'):  
    #    splice_fuv(out_dict, 'lp2') 
    #if out_dict.__contains__('G130M_final_lp3') and out_dict.__contains__('G160M_final_lp3'):  
    #    splice_fuv(out_dict, 'lp3')
    #if out_dict.__contains__('G130M_final_lp4') and out_dict.__contains__('G160M_final_lp4'):  
    #    splice_fuv(out_dict, 'lp4')             
    #if out_dict.__contains__('G130M_final_lp12') and out_dict.__contains__('G160M_final_lp12'):  
    #    splice_fuv(out_dict, 'lp12')

    #JRD add NUV

    if (out_dict.__contains__('G185M_final_lp1') and out_dict.__contains__('G225M_final_lp1')) or (out_dict.__contains__('G185M_final_lp1') and out_dict.__contains__('G285M_final_lp1')) or (out_dict.__contains__('G225M_final_lp1') and out_dict.__contains__('G285M_final_lp1')):  
        #JRD NOT YET READY
        splice_nuv(out_dict)

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
       #JRD
       print("SPLICE_FUV ", key_to_combine, splitwave)
    elif (np.size(where_160) < 1): 
       where_160 = np.where((g160m['DQ'] > 0.) & (g160m['WAVE'] <  1450)) 
       splitwave = np.max(g160m['WAVE'][np.where((g160m['DQ'] > 0.) & (g160m['WAVE'] <  1450))]) + 1. 
       i_130 = np.max(np.where(g130m['WAVE'] < splitwave))  
       i_160 = np.min(np.where(g160m['WAVE'] > splitwave))  
       out_dict['FUVM_final_'+key_to_combine] = vstack([ g130m[0:i_130], g160m[i_160:-1] ])
       #JRD
       print("SPLICE_FUV ", key_to_combine, splitwave)
    else:
        #JRD
        print('Problem with splitwave, nothing will be done')
       

def splice_nuv(out_dict):

    f185 = False
    f225 = False
    f285 = False

    lambda_max = 0.
    spec = None
    
    if (out_dict.__contains__('G185M_final_lp1')):
        g185m = out_dict['G185M_final_lp1']
        f185 = True
        i185_max = np.max(np.where(g185m['FLUXWGT']>0.))
        spec = g185m[0:i185_max]
        lambda_max = g185m['WAVE'][i185_max]
    if (out_dict.__contains__('G225M_final_lp1')):
        g225m = out_dict['G225M_final_lp1']
        f225=True
        i225_min = np.min(np.where(g225m['FLUXWGT']>0.))
        lmin = g225m['WAVE'][i225_min]
        i225_max = np.max(np.where(g225m['FLUXWGT']> 0.))
        lmax = g225m['WAVE'][i225_max]
        if f185==True:
            if lmin > lambda_max:
                spec = vstack([spec, g225m[i225_min:i225_max]])
            else:
                good225 = np.where((g225m['WAVE'] > lambda_max) & (g225m['FLUXWGT'] > 0))
                i225_min = np.min(good225)
                spec = vstack([spec, g225m[i225_min:i225_max]])
        else:
            spec = g225m[0:i225_max]
        lambda_max = lmax

    if (out_dict.__contains__('G285M_final_lp1')):
        g285m = out_dict['G285M_final_lp1']
        f285=True
        i285_min = np.min(np.where(g285m['FLUXWGT']>0.))
        i285_max = np.max(np.where(g285m['FLUXWGT']> 0.))
        lmin = g285m['WAVE'][i285_max]
        lmax = g285m['WAVE'][i285_max]

        if spec != None:
            if lmin > lambda_max:
                spec = vstack([spec, g285m[i285_min:i285_max]])
            else:
                good285 = np.where((spec['WAVE'] > lambda_max) & (spec['FLUXWGT'] > 0.))
                i285_min = np.min(good285)
                spec = vstack([spec, g285m[i285_min:i285_max]])
        else:
            spec = g285m[0:i285_max]

        
    if spec !=None:
        out_dict['NUVM_final_lp1'] = spec
   






def apply_wavelength_shifts(out_dict, shifts_to_apply): 

    exp_interp_table = out_dict['exp_interp_table']
    #JRD python 3
    exp_keys = list(exp_interp_table.keys())
    number_of_exposures = np.size(list(exp_interp_table.keys()))

    columns_to_shift = ['FLUX', 'ERROR', 'GROSS', 'GCOUNTS','NET','BACKGROUND','DQ','DQ_WGT','EXP_PIX','DQ_WEIGHT','FLUXFACTOR'] 

    for nn in np.arange(number_of_exposures):
       for col in columns_to_shift: 
          exp_interp_table[exp_keys[nn]][1][col] = np.roll(exp_interp_table[exp_keys[nn]][1][col], int(-1*shifts_to_apply[nn])) 


def get_wavelength_shifts(out_dict): 
 
    exp_interp_table = out_dict['exp_interp_table']
    #JRD python 3
    exp_keys = list(exp_interp_table.keys())
    number_of_exposures = np.size(list(exp_interp_table.keys()))
 
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
            coeff_here = np.correlate(np.roll(exp_flux[where2],shifter), coadd_flux[where1])
            #JRD
            print('EXP', exp_keys[nn], grating, wavelimits, shifter, coeff_here[0])
            coeffs[i] = coeff_here[0]
            i = i + 1
 
        best_shift = shifts[np.where(coeffs == np.max(coeffs))]
        #JRD
        print('BEST SHIFT', exp_keys[nn], grating, wavelimits, best_shift, coeffs[np.where(coeffs == np.max(coeffs))])
        shifts_to_apply[nn] = best_shift 

    return shifts_to_apply 


     
def obtain_chi2(out_dict, final_spectrum_key):					   ##### 

    out_dict[final_spectrum_key]['EXP_COUNTER'] = out_dict[final_spectrum_key]['FLUX'] * 0.0 

    #JRD
    print('OBTAIN CHI2 for: ', final_spectrum_key)
    final_spectrum = out_dict[final_spectrum_key] 
   
    exp_interp_table = out_dict['exp_interp_table']       #### NOTE: exptable[i][0] refers to segment A, exptable[i][1] refers to segment B
    #JRD python 3
    exp_keys = list(exp_interp_table.keys())

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
                #JRD FOR NUV changing this (number_of_segments = 3)
                for segment in range(1,number_of_segments):
                #segment = 1 
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
    counter285 = 0
    counter230 = 0

    for i in out_dict['h0table'].keys():
        #JRD
        print('COUNT_EXPOSURES: OPT_ELEM: ', out_dict['h0table'][i]['ROOTNAME'], out_dict['h0table'][i]['OPT_ELEM'])
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
        #JRD adding NUV gratings
        if(out_dict['h0table'][i]['OPT_ELEM'] == 'G285M'):
            counter285 = counter285 + 1
        if(out_dict['h0table'][i]['OPT_ELEM'] == 'G230L'):
            counter230 = counter230 + 1
    out_dict['number_of_G130M_exposures'] = counter130
    out_dict['number_of_G160M_exposures'] = counter160
    out_dict['number_of_G140L_exposures'] = counter140
    out_dict['number_of_G185M_exposures'] = counter185
    out_dict['number_of_G225M_exposures'] = counter225
    out_dict['number_of_G285M_exposures'] = counter285
    out_dict['number_of_G230L_exposures'] = counter230

    #JRD
    print('COUNT_EXPOSURES: There are ', counter130, out_dict['number_of_G130M_exposures'])
    print('COUNT_EXPOSURES: There are ', counter160, out_dict['number_of_G160M_exposures'])
    print('COUNT_EXPOSURES: There are ', counter140, out_dict['number_of_G140L_exposures'])
    print('COUNT_EXPOSURES: There are ', counter185, out_dict['number_of_G185M_exposures'])
    print('COUNT_EXPOSURES: There are ', counter225, out_dict['number_of_G225M_exposures'])
    print('COUNT_EXPOSURES: There are ', counter285, out_dict['number_of_G285M_exposures'])
    print('COUNT_EXPOSURES: There are ', counter230, out_dict['number_of_G230L_exposures'])


def create_coadd_format(wavegrid):                                ##### this creates a dictionary that holds the canonical outputs of the coadd
                                                ##### this has to be a dictionary because it contains a bunch of tables
                                                ##### with different levels of the coadd hierarchy in them
    flux = wavegrid * 0.0                                     ##### with the input wavelength grid
    error = wavegrid * 0.0

    #JRD
    print('CREATE_OUTPUT: creating master output table')
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


def get_pix_to_wave(out_dict, screen, lifetime):

    exptable = out_dict['exptable']                                        #### we'll just use this variable inside this routine

    n_pix_per_segment = 16384 
    n_segments = 2 
    if ('NUV' in screen): 
       n_segments = 3 
       n_pix_per_segment = 1274  

    #JRD python 3
    print('Screen = ', screen, ' so n_segments = ', n_segments, 'and n_pix_per_segment = ', n_pix_per_segment)

    #JRD python 3
    number_of_exposures = np.size(list(exptable.keys()))

    oldpix = np.empty([ n_segments, n_pix_per_segment, number_of_exposures ])
    #JRD this is confusing so explaining here: newpix will contain the indices of the new wavelength grid corresponding to the old wavelengths. 
    newpix = np.empty([ n_segments, n_pix_per_segment, number_of_exposures ])

    #JRD PASSING THE LIFETIME ARGUMENT HERE, NO MORE FINAL_ALL
    if (out_dict.__contains__('G130M_final_lp' + lifetime)):
        final_vector_130 = out_dict['G130M_final_lp'+lifetime]['WAVE']                        #### this is the "Final" G130 wavelength vector to interpolate to
    if (out_dict.__contains__('G160M_final_lp'+lifetime)):
        final_vector_160 = out_dict['G160M_final_lp'+lifetime]['WAVE']                        #### this is the "Final" G160 wavelength vector to interpolate to
    if (out_dict.__contains__('G140L_final_lp'+lifetime)):
        final_vector_140 = out_dict['G140L_final_lp'+lifetime]['WAVE']                        #### this is the "Final" G160 wavelength vector to interpolate to
    #JRD ADDING NUV
    if (out_dict.__contains__('G185M_final_lp1')):
        final_vector_185 = out_dict['G185M_final_lp1']['WAVE']                        #### this is the "Final" G160 wavelength vector to interpolate to
    if (out_dict.__contains__('G225M_final_lp1')):
        final_vector_225 = out_dict['G225M_final_lp1']['WAVE']                        #### this is the "Final" G160 wavelength vector to interpolate to
    if (out_dict.__contains__('G285M_final_lp1')):
        final_vector_285 = out_dict['G285M_final_lp1']['WAVE']                        #### this is the "Final" G160 wavelength vector to interpolate to
    if (out_dict.__contains__('G230L_final_lp1')):
        final_vector_230 = out_dict['G230L_final_lp1']['WAVE']                        #### this is the "Final" G160 wavelength vector to interpolate to



    #JRD python 3
    exp_keys = list(exptable.keys())

    for nn in np.arange(number_of_exposures):
        exptable[exp_keys[nn]]['NEWPIX'] = exptable[exp_keys[nn]]['FLUX'] * 0.0 +  np.arange(n_pix_per_segment)  #### add NEWPIX vectors to the table for each exposure, no units
        exptable[exp_keys[nn]]['NEWPIX'].unit = ' '

        for i in np.arange(n_segments):                                    #### we may not actually need oldpix
            oldpix[i, :, nn] = np.arange(n_pix_per_segment)

        number_of_segments = np.size(exptable[exp_keys[nn]]['SEGMENT'])

        if (out_dict['h0table'][exp_keys[nn]]['OPT_ELEM'] == 'G130M'):
            #JRD adding numer of pixel of new grid
            final_npix = len(final_vector_130)
            for nseg in np.arange(number_of_segments):
                newpix[nseg,:,nn] = np.round(np.interp(exptable[exp_keys[nn]]['WAVELENGTH'][nseg], final_vector_130, np.arange(final_npix)))         #### G130M 
                for ipix in np.arange(np.size(newpix[nseg,:,nn])-3):
                    if newpix[nseg,ipix+1,nn] - newpix[nseg,ipix,nn] > 1.: exptable[exp_keys[nn]]['DQ'][nseg][ipix] = 2.**14. 
                exptable[exp_keys[nn]]['NEWPIX'][nseg] = newpix[nseg,:,nn]
        elif (out_dict['h0table'][exp_keys[nn]]['OPT_ELEM'] == 'G160M'):
            #JRD adding numer of pixel of new grid
            final_npix = len(final_vector_160)
            for nseg in np.arange(number_of_segments):
                newpix[nseg,:,nn] = np.round(np.interp(exptable[exp_keys[nn]]['WAVELENGTH'][nseg], final_vector_160, np.arange(final_npix)))         #### G160M 
                for ipix in np.arange(np.size(newpix[nseg,:,nn])-3):
                    if newpix[nseg,ipix+1,nn] - newpix[nseg,ipix,nn] > 1.: exptable[exp_keys[nn]]['DQ'][nseg][ipix] = 2.**14. 
                exptable[exp_keys[nn]]['NEWPIX'][nseg] = newpix[nseg,:,nn]
        elif (out_dict['h0table'][exp_keys[nn]]['OPT_ELEM'] == 'G140L'):
            #JRD adding numer of pixel of new grid
            final_npix = len(final_vector_140)
            for nseg in np.arange(number_of_segments):
                newpix[nseg,:,nn] = np.round(np.interp(exptable[exp_keys[nn]]['WAVELENGTH'][nseg], final_vector_140, np.arange(final_npix)))         #### G140M 
                for ipix in np.arange(np.size(newpix[nseg,:,nn])-3):
                    if newpix[nseg,ipix+1,nn] - newpix[nseg,ipix,nn] > 1.: exptable[exp_keys[nn]]['DQ'][nseg][ipix] = 2.**14. 
                exptable[exp_keys[nn]]['NEWPIX'][nseg] = newpix[nseg,:,nn]

        #JRD ADDING ALL NUV GRATINGS
        elif (out_dict['h0table'][exp_keys[nn]]['OPT_ELEM'] == 'G185M'):
            #JRD adding numer of pixel of new grid
            final_npix = len(final_vector_185)
            for nseg in np.arange(number_of_segments):
                newpix[nseg,:,nn] = np.round(np.interp(exptable[exp_keys[nn]]['WAVELENGTH'][nseg], final_vector_185, np.arange(final_npix)))         #### G185M 
                for ipix in np.arange(np.size(newpix[nseg,:,nn])-3):
                    if newpix[nseg,ipix+1,nn] - newpix[nseg,ipix,nn] > 1.: exptable[exp_keys[nn]]['DQ'][nseg][ipix] = 2.**14. 
                exptable[exp_keys[nn]]['NEWPIX'][nseg] = newpix[nseg,:,nn]
        elif (out_dict['h0table'][exp_keys[nn]]['OPT_ELEM'] == 'G225M'):
            #JRD adding numer of pixel of new grid
            final_npix = len(final_vector_225)
            for nseg in np.arange(number_of_segments):
                newpix[nseg,:,nn] = np.round(np.interp(exptable[exp_keys[nn]]['WAVELENGTH'][nseg], final_vector_225, np.arange(final_npix)))         #### G225M 
                for ipix in np.arange(np.size(newpix[nseg,:,nn])-3):
                    if newpix[nseg,ipix+1,nn] - newpix[nseg,ipix,nn] > 1.: exptable[exp_keys[nn]]['DQ'][nseg][ipix] = 2.**14. 
                exptable[exp_keys[nn]]['NEWPIX'][nseg] = newpix[nseg,:,nn]

        elif (out_dict['h0table'][exp_keys[nn]]['OPT_ELEM'] == 'G285M'):
            #JRD adding numer of pixel of new grid
            final_npix = len(final_vector_285)
            for nseg in np.arange(number_of_segments):
                newpix[nseg,:,nn] = np.round(np.interp(exptable[exp_keys[nn]]['WAVELENGTH'][nseg], final_vector_285, np.arange(final_npix)))         #### G225M 
                for ipix in np.arange(np.size(newpix[nseg,:,nn])-3):
                    if newpix[nseg,ipix+1,nn] - newpix[nseg,ipix,nn] > 1.: exptable[exp_keys[nn]]['DQ'][nseg][ipix] = 2.**14. 
                exptable[exp_keys[nn]]['NEWPIX'][nseg] = newpix[nseg,:,nn]

        elif (out_dict['h0table'][exp_keys[nn]]['OPT_ELEM'] == 'G230L'):
            #JRD adding numer of pixel of new grid
            final_npix = len(final_vector_230)
            for nseg in np.arange(number_of_segments):
                newpix[nseg,:,nn] = np.round(np.interp(exptable[exp_keys[nn]]['WAVELENGTH'][nseg], final_vector_230, np.arange(final_npix)))         #### G225M 
                for ipix in np.arange(np.size(newpix[nseg,:,nn])-3):
                    if newpix[nseg,ipix+1,nn] - newpix[nseg,ipix,nn] > 1.: exptable[exp_keys[nn]]['DQ'][nseg][ipix] = 2.**14. 
                exptable[exp_keys[nn]]['NEWPIX'][nseg] = newpix[nseg,:,nn]
                
        else:
            #JRD
            print('OOPS GET_PIX_TO_WAVE DID NOT FIND THE GRATING YOU WANT')

    print("STILL NEED TO IMPLEMENT DQ FLAGGING FOR HITCHES, HERE, HERE, AND HERE")




def get_wavelength_grid(grating):

    #JRD modified this code significantly to resolve some issues with FUV and enable the NUV
    
    #JRD point to NUV reference file
    if (grating=='G185M') or (grating == 'G225M') or (grating == 'G285M') or (grating == 'G230L'):
        wavefile = '/grp/crds/hst/references/hst/12e1345gl_disp.fits'
        npix_det = 1274
    else:
        #JRD updates path to disp reference file to be the new one for LP1. Should be LP dependent
        wavefile = "/grp/crds/hst/references/hst/05i1639ml_disp.fits"
        npix_det = 16384
    #JRD python 3
    print("looking for wavelength grid file",wavefile,"......should really be passing in a reference/ dir........", wavefile)
    #JRD path is now in wavefile, removing
    if os.path.isfile(wavefile):
        a = Table.read(wavefile)
    else:
        print("CRASHING: CRDS HAS LOST MY DISP REFERENCE FILE !!!")
    #JRD This is now useless, removing
    #elif os.path.isfile("../../reference/"+wavefile):
    #    a = Table.read('../../reference/'+wavefile)
    #elif os.path.isfile("../../msp-code/reference/"+wavefile):
    #    a = Table.read('../../msp-code/reference/'+wavefile)

    #JRD
    #The split in grating is meant to only do FUVA for G140L. For NUV, we need all stripes as well, so modifying statement to be if not G140L, then do this, else do that.
    #if (grating == 'G130M' or grating == 'G160M'):
    if grating != 'G140L':
        a = a[np.where(a['OPT_ELEM'] == grating)]
        #JRD: We want 1222, changing threshold in wavelength. Shoudl we include the blue modes? I think yes.
        #a = a[np.where(a['CENWAVE'] > 1200)]
        a = a[np.where(a['APERTURE'] == 'PSA')]

    #JRD moving the "else..." here since teh only difference betwen the two cases is whether or not we coadd FUVA only.
    else:
        a = a[np.where(a['OPT_ELEM'] == grating)]
        a = a[np.where(a['APERTURE'] == 'PSA')]
        a = a[np.where(a['SEGMENT'] == 'FUVA')]
        
    n_settings = np.size(a['COEFF'][:,1])

    w = np.float64(np.empty([npix_det, n_settings]))      #### create array for default wavelength solutions
    xfull = np.arange(npix_det)


    #JRD
    #NUV has quadratic dispersion, but we just wnat the final grid, which does not ahve to be quadratic as long as the resampling takes into account the quadratic terms in the original spectra . however, to know the min/max wavelength of each spectrum included in teh co-add, we need the quadratic solution.
    for i in np.arange(n_settings):
        w[:,i] = (xfull - a['D'][i] + a['D_TV03'][i])**2 * np.float64(a['COEFF'][i,2]) +  (xfull - a['D'][i] + a['D_TV03'][i]) * np.float64(a['COEFF'][i,1]) + a['COEFF'][i,0]

        mean_dispersion = np.mean(a['COEFF'][:,1])     ##### this will be the "dispersion" of the final wavelength vector
        #JRD: computing mean quadratic term, may not be used
        #mean_quad_disp = np.mean(a['COEFF'][:,2]) ##### this will be the "quadratic dispersion term" of the final wavelength vector
        
        #JRD THIS WAS WRONG (FORGOT the TV03/on orbit OFFSET). CORRECTING ... BUT NOW USELESS WITH NEW APPROACH, SO COMMENTING
        #nominal_wavelength_old = np.mean(a['COEFF'][:,1] * (16384/2.) + a['COEFF'][:,0])            ##### this will be the "zeropoint" of the final wavelength vector
        #nominal_wavelength = np.mean(w[16384/2, :])
        #print("JULIA CHECK WAVELENGTH ", nominal_wavelength_old, nominal_wavelength)
        
        #JRD python 3
        #print('Will use mean dispersion: ', mean_dispersion, '  and nominal wavelengths : ', nominal_wavelength)

        #JRD changing the approach here. 50000 pixels does not work for NUV. ALso, it does not work with FUV if blue modes are included. Need to instead compute min/max wavelenghs of all exposures and resample on grid with npix = (max-min)/mean_dispersion.

    min_w = np.min(w[0,:])
    max_w = np.max(w[npix_det-1,:])
        
    final_npix = int(np.round( (max_w - min_w)/mean_dispersion))
    
    final_vector = np.float64(np.arange(final_npix)) * mean_dispersion + min_w  #nominal_wavelength    ##### create the final wavelength vector
    nominal_wavelength = final_vector[final_npix//2]


    return final_vector




def get_interp_exptable(out_dict, lifetime): 

    #### this takes in the exposure table and outputs one in the same structure but with the output wavelength grid

    exptable = out_dict['exptable']       		#### NOTE: exptable[i][0] refers to segment A, exptable[i][1] refers to segment B
    #JRD python 3
    exp_keys = list(exptable.keys())        #### such that exptable[i][0]['WAVELENGTH'] is the segment A wavelength vector and so on
    number_of_exposures = np.size(list(exptable.keys()))
    #JRD
    print("GET_INTERP_EXPTABLE: There are ", number_of_exposures, ' exposures to interpolate into exp_int_table')

    exp_interp_table = copy.deepcopy(exptable)

    #JRD python 3
    for exp in list(exptable.keys()): 

        #number_of_segments = np.size(exptable[exp]['SEGMENT'])
        nseg = np.size(exptable[exp]['SEGMENT'])
        
        grating   = exptable[exp][0]['GRATING'] #I may beed exptable[nn][0]['gRATING'] for the segment

        if grating == 'G130M':
            final_npix = len(out_dict['G130M_final_lp'+lifetime]['WAVE'])
            #nseg = 2
        elif grating == 'G160M':
            final_npix = len(out_dict['G160M_final_lp'+lifetime]['WAVE'])
            #nseg = 2
        elif grating == 'G140L':
            final_npix = len(out_dict['G140L_final_lp'+lifetime]['WAVE'])
            #nseg = 1
        elif grating == 'G185M':
            final_npix = len(out_dict['G185M_final_lp1']['WAVE'])
            #nseg = 3
        elif grating == 'G225M':
            final_npix = len(out_dict['G225M_final_lp1']['WAVE'])
            #nseg = 3
        elif grating == 'G285M':
            final_npix = len(out_dict['G285M_final_lp1']['WAVE'])
            #nseg = 3
        elif grating == 'G230L':
            final_npix = len(out_dict['G230L_final_lp1']['WAVE'])
            #nseg = 3                        
        else:
            print('BAD NEWS, I DO NOT KNOW HOW TO HANDLE ', grating, 'IN GET_INTERP_EXPTABLE')

        
        #JRD this does nto work anymore
        #if number_of_segments == 2: 
        #    blankarr = np.zeros((2,50000)) 
        #elif number_of_segments == 1: 
        #    blankarr = np.zeros((1,50000)) 
        #else:
            #JRD
        #    print('BAD NEWS, I DO NOT KNOW HOW TO HANDLE ', number_of_segments, ' SEGMENTS IN GET_INTERP_EXPTABLE')

        #REPLACING WITH

        blankarr = np.zeros([nseg, final_npix])
          
        #JRD
        print('GET_INTERP_EXPTABLE: Deleting table entries:')
        keys_to_delete = ['WAVELENGTH', 'FLUX', 'ERROR', 'GROSS', 'GCOUNTS', 'NET', 'BACKGROUND', 'DQ', 'DQ_WGT', \
                         'BACKGROUND_PER_ROW', 'NUM_EXTRACT_ROWS', 'ACTUAL_EE', 'Y_LOWER_OUTER', \
                         'Y_UPPER_OUTER', 'Y_LOWER_INNER', 'Y_UPPER_INNER', 'EXP_PIX', 'DQ_WEIGHT', 'FLUXFACTOR', 'NEWPIX']

        for key_here in keys_to_delete:
            #JRD python 3

            if key_here in list(exp_interp_table[exp].keys()): del exp_interp_table[exp][key_here] 

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






def exposure_interpolate(out_dict, lifetime): 

    #### here we are going to go through the exptable and "interpolate" the exposures onto the exp_interp_table 
    #### this should be done carefully so as to preserve DQ and such. 

    exptable = out_dict['exptable']       		   #### NOTE: exptable[i][0] refers to segment A, exptable[i][1] refers to segment B
    exp_interp_table = out_dict['exp_interp_table']        #### the interpolated arrays 

    #JRD python 3
    number_of_exposures = np.size(list(exptable.keys()))
    #JRD
    print("EXPOSURE_INTERPOLATE: There are ", number_of_exposures, ' exposures to interpolate')

    #JRD python 3
    exp_keys = list(exptable.keys())                    

    exp_counter = 0 
    for i in np.arange(number_of_exposures):

        number_of_segments = np.size(exptable[exp_keys[i]]['SEGMENT'])

        for nseg in np.arange(number_of_segments):                #### loop over segments
            newpix = out_dict['exptable'][exp_keys[i]]['NEWPIX'][nseg,:]

            grating_here = exptable[exp_keys[i]][nseg]['GRATING']
            if (grating_here =='G185M') or (grating_here =='G225M') or (grating_here =='G285M') or (grating_here =='G230L'):
                this_lp= 'lp1'
            else:
                this_lp = 'lp'+lifetime                                                                                 
                            
            exp_interp_table[exp_keys[i]][nseg]['WAVELENGTH'] = out_dict[grating_here+'_final_' + this_lp]['WAVE']

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
    #JRD python 3
    exp_keys = list(exp_interp_table.keys())      
    h0 = out_dict['h0table']
    number_of_exposures = np.size(exp_keys)
    #JRD
    print("COS_COUNTS_COADD_REMIX: There are ", number_of_exposures, ' exposures to coadd for key_to_combine', key_to_combine)
    print(exp_keys)

    exp_counter = 0

    for i in np.arange(number_of_exposures):

        number_of_segments = np.size(exp_interp_table[exp_keys[i]]['SEGMENT'])

        if h0[exp_keys[i]]['OPT_ELEM'] in key_to_combine:

            #JRD
            print('Stacking up (remix): ', h0[exp_keys[i]]['ROOTNAME'], h0[exp_keys[i]]['OPT_ELEM'], key_to_combine)

            for nseg in np.arange(number_of_segments):                #### loop over segments

                final_table['EXP_PIX'] = final_table['EXP_PIX'] + exp_interp_table[exp_keys[i]][nseg]['EXP_PIX']

                final_table['GROSSCOUNTS'] = final_table['GROSSCOUNTS'] + exp_interp_table[exp_keys[i]][nseg]['GROSS'] * exp_interp_table[exp_keys[i]][nseg]['EXP_PIX']
                            
                final_table['NETCOUNTS'] = final_table['NETCOUNTS'] + exp_interp_table[exp_keys[i]][nseg]['NET'] * exp_interp_table[exp_keys[i]][nseg]['EXP_PIX']
                            
                final_table['BACK'] = final_table['BACK'] + exp_interp_table[exp_keys[i]][nseg]['BACK'] * exp_interp_table[exp_keys[i]][nseg]['EXP_PIX']

                final_table['FLUXWGT'] = final_table['FLUXWGT'] + exp_interp_table[exp_keys[i]][nseg]['FLUX'] * exp_interp_table[exp_keys[i]][nseg]['EXP_PIX']

                final_table['FLUXWGT_ERR'] = final_table['FLUXWGT_ERR'] + (exp_interp_table[exp_keys[i]][nseg]['EXP_PIX'] * exp_interp_table[exp_keys[i]][nseg]['ERROR'])**2 
                            
                final_table['DQ'] = (final_table['DQ']).astype(int) | (exp_interp_table[exp_keys[i]][nseg]['DQ']).astype(int) 
     
                final_table['FLUXFACTOR'] = final_table['FLUXFACTOR'] + exp_interp_table[exp_keys[i]][nseg]['FLUXFACTOR']*exp_interp_table[exp_keys[i]][nseg]['EXP_PIX']

                #print("CHECK EXP ", h0[exp_keys[i]]['CENWAVE'], np.max(exp_interp_table[exp_keys[i]][nseg]['EXP_PIX']), np.max(exp_interp_table[exp_keys[i]][nseg]['FLUX']), np.max(exp_interp_table[exp_keys[i]][nseg]['FLUXFACTOR']), np.max(exp_interp_table[exp_keys[i]][nseg]['NET']))

            exp_counter = exp_counter + 1                        #### increment the number of exposures we've added
        else:
            #JRD
            print("THE GRATING ", h0[exp_keys[i]]['OPT_ELEM'], "DOES NOT MATCH WITH KEY_TO_COMBINE ", key_to_combine, " FOR EXPOSURE ", h0[exp_keys[i]]['ROOTNAME'])

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
        #print "WATCH OUT THERE MAY BE A PROBLEM HERE!!!!" 
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
        #JRD python 3
        print("Oh no, all_exposures.txt does not exist for this target, you might want to try something else")
        return 0, 0, 0 

    if "Flag" in dataset_list.colnames:
        dataset_list.remove_rows(np.where(dataset_list['Flag'] == 0))    #### delete exposures flagged with 0 in the all_exposures.txt file

    if len(dataset_list) == 0:
        #JRD python 3
        print('GET_DICT_OF_EXPOSURES: It appears that no exposures were flagged as 1, will try to exit gracefully!')
        return 0, 0, 0 

    dataset_list = dataset_list['Rootname'] 

    for i in range(np.size(dataset_list)):

        #get the header
        filename = dataset_list[i]+'_x1d.fits.gz'
        #JRD
        h0[i] = fits.open(filename)[0].header #fitsio.read_header(filename, 0)
        h1[i] = fits.open(filename)[1].header #fitsio.read_header(filename, 1)

        t_new = Table.read(dataset_list[i]+'_x1d.fits.gz')
        #print 'GET_DICT_OF_EXPOSURES: We have read in ', dataset_list[i]+'_x1d.fits  for targname  ', h0[0]['TARGNAME']
        tt[i] = t_new
        tt[i]['ROOTNAME'] = dataset_list[i]

    #JRD python 3
    print('GET_DICT_OF_EXPOSURES: ', dataset_list)
    print('GET_DICT_OF_EXPOSURES: We have read in i = ', i+1, '  exposures ')

    #JRD python 3
    for key in list(h0.keys()): 
       h0[key]['OPT_ELEM'] = str(h0[key]['OPT_ELEM']).strip() 
       h0[key]['TARGNAME'] = str(h0[key]['TARGNAME']).strip() 
       h0[key]['ROOTNAME'] = str(h0[key]['ROOTNAME']).strip() 

    return tt, h0, h1


def screen_dict_of_exposures(tt, h0, h1, screen='FUVM', lifetime='ALL'):         # screen out the dictionary entries we don't want
    #JRD
    print("JRD update")
    print('SCREEN_DICT_OF_EXPOSURES: The requested screen is: ', screen)
    print('SCREEN_DICT_OF_EXPOSURES: The requested LP is: ', lifetime)

    if (screen == 'FUVM' and lifetime == 'ALL'):                # for FUVM and all LPs
        #JRD
         print('condition 1')
         #JRD python 3
         for i in list(h0.keys()):
             #JRD: We want 1222, removing and h0[i]['CENWAVE'] > 1250)
             if (h0[i]['DETECTOR'] == 'FUV' and (h0[i]['OPT_ELEM'] == 'G130M' or h0[i]['OPT_ELEM'] == 'G160M')) :
                #JRD
                 print('SCREEN_DICT_OF_EXPOSURES 1: Keeping dictionary entry for ', h0[i]['ROOTNAME'], h0[i]['DETECTOR'], h0[i]['OPT_ELEM'], h0[i]['CENWAVE'], h0[i]['LIFE_ADJ'])
             else:
                 #JRD
                 print('SCREEN_DICT_OF_EXPOSURES 1: Deleting dictionary entry for ', h0[i]['ROOTNAME'], h0[i]['DETECTOR'], h0[i]['OPT_ELEM'], h0[i]['CENWAVE'], h0[i]['LIFE_ADJ'])
                 del tt[i]
                 del h0[i]
                 del h1[i]
    elif (screen == 'FUVM' and ((lifetime == '1') or (lifetime == '2') or (lifetime == '3') or (lifetime=='4'))):                     # for FUVM and a specific LP
        #JRD
        print('condition 2')
        #JRD python 3
        for i in list(h0.keys()):
            #JRD: We want 1222, removing and h0[i]['CENWAVE'] > 1250)

            if (h0[i]['DETECTOR'] == 'FUV' and (h0[i]['OPT_ELEM'] == 'G130M' or h0[i]['OPT_ELEM'] == 'G160M') and h0[i]['LIFE_ADJ'] == int(lifetime)):
                #JRD
                print('SCREEN_DICT_OF_EXPOSURES 2: Keeping dictionary entry for ', h0[i]['ROOTNAME'], h0[i]['DETECTOR'], h0[i]['OPT_ELEM'], h0[i]['CENWAVE'], h0[i]['LIFE_ADJ'])
            else:
                #JRD
                print('SCREEN_DICT_OF_EXPOSURES 2: Deleting dictionary entry for ', h0[i]['ROOTNAME'], h0[i]['DETECTOR'], h0[i]['OPT_ELEM'], h0[i]['CENWAVE'], h0[i]['LIFE_ADJ'])
                del tt[i]
                del h0[i]
                del h1[i]
    elif (screen == 'FUVM' and lifetime == '12'):				# for FUVM and 1+2 combo
        #JRD
        print('condition 3')
        #JRD python 3
        for i in list(h0.keys()):
            #JRD: We want 1222, removing and h0[i]['CENWAVE'] > 1250)
            if (h0[i]['DETECTOR'] == 'FUV' and (h0[i]['OPT_ELEM'] == 'G130M' or h0[i]['OPT_ELEM'] == 'G160M') and h0[i]['LIFE_ADJ'] < 3 and h0[i]['LIFE_ADJ'] > 0):
                #JRD
                print('SCREEN_DICT_OF_EXPOSURES 2: Keeping dictionary entry for ', h0[i]['ROOTNAME'], h0[i]['DETECTOR'], h0[i]['OPT_ELEM'], h0[i]['CENWAVE'], h0[i]['LIFE_ADJ'])
            else:
                #JRD
                print('SCREEN_DICT_OF_EXPOSURES 2: Deleting dictionary entry for ', h0[i]['ROOTNAME'], h0[i]['DETECTOR'], h0[i]['OPT_ELEM'], h0[i]['CENWAVE'], h0[i]['LIFE_ADJ'])
                del tt[i]
                del h0[i]
                del h1[i]
    elif (screen == 'FUVL' and lifetime == 'ALL'):                     # for FUVL and a specific LP
        #JRD
        print('condition 4')
        #JRD python 3
        for i in list(h0.keys()):
            if (h0[i]['DETECTOR'] == 'FUV' and h0[i]['OPT_ELEM'] == 'G140L'):
                #JRD
                print('SCREEN_DICT_OF_EXPOSURES 2: Keeping dictionary entry for ', h0[i]['ROOTNAME'], h0[i]['DETECTOR'], h0[i]['OPT_ELEM'], h0[i]['CENWAVE'], h0[i]['LIFE_ADJ'])
            else:
                #JRD
                print('SCREEN_DICT_OF_EXPOSURES 2: Deleting dictionary entry for ', h0[i]['ROOTNAME'], h0[i]['DETECTOR'], h0[i]['OPT_ELEM'], h0[i]['CENWAVE'], h0[i]['LIFE_ADJ'])
                del tt[i]
                del h0[i]
                del h1[i]
    elif (screen == 'FUVL') and ((lifetime == '1') or (lifetime == '2') or (lifetime == '3') or (lifetime=='4')):                        # for FUVL
        #JRD
        print('condition 5')
        #JRD python 3
        for i in list(h0.keys()):
            if (h0[i]['DETECTOR'] == 'FUV' and h0[i]['OPT_ELEM'] == 'G140L' and h0[i]['LIFE_ADJ'] == int(lifetime)):
                #JRD
                print('SCREEN_DICT_OF_EXPOSURES 3: Keeping dictionary entry for ', h0[i]['ROOTNAME'], h0[i]['DETECTOR'], h0[i]['OPT_ELEM'], h0[i]['CENWAVE'], h0[i]['LIFE_ADJ'])
            else:
                #JRD
                print('SCREEN_DICT_OF_EXPOSURES 3: Deleting dictionary entry for ', h0[i]['ROOTNAME'], h0[i]['DETECTOR'], h0[i]['OPT_ELEM'], h0[i]['CENWAVE'], h0[i]['LIFE_ADJ'])
                del tt[i]
                del h0[i]
                del h1[i]
    elif (screen == 'FUVL'):                        # for FUVL and combo of 1+2
        #JRD
        print('condition 6')
        for i in list(h0.keys()):
            if (h0[i]['DETECTOR'] == 'FUV' and h0[i]['OPT_ELEM'] == 'G140L' and h0[i]['LIFE_ADJ'] < 5):
                #JRD
                print('SCREEN_DICT_OF_EXPOSURES 3: Keeping dictionary entry for ', h0[i]['ROOTNAME'], h0[i]['DETECTOR'], h0[i]['OPT_ELEM'], h0[i]['CENWAVE'], h0[i]['LIFE_ADJ'])
            else:
                #JRD
                print('SCREEN_DICT_OF_EXPOSURES 3: Deleting dictionary entry for ', h0[i]['ROOTNAME'], h0[i]['DETECTOR'], h0[i]['OPT_ELEM'], h0[i]['CENWAVE'], h0[i]['LIFE_ADJ'])
                del tt[i]
                del h0[i]
                del h1[i]

    #JRD lifetime = 1 for NUV
    elif (screen == 'NUVM' and (lifetime == 'ALL' or lifetime =='1')):                        # for NUVM
        #JRD
        print('condition 7')
        for i in list(h0.keys()):
            #JRD adding the G285M
            if (h0[i]['DETECTOR'] == 'NUV' ) and (h0[i]['OPT_ELEM'] == 'G185M' or h0[i]['OPT_ELEM'] == 'G225M' or h0[i]['OPT_ELEM'] == 'G285M'):
                #JRD
                print('SCREEN_DICT_OF_EXPOSURES 3: Keeping dictionary entry for ', h0[i]['ROOTNAME'], h0[i]['DETECTOR'], h0[i]['OPT_ELEM'], h0[i]['CENWAVE'], h0[i]['LIFE_ADJ'])
            else:
                #JRD
                print('SCREEN_DICT_OF_EXPOSURES 3: Deleting dictionary entry for ', h0[i]['ROOTNAME'], h0[i]['DETECTOR'], h0[i]['OPT_ELEM'], h0[i]['CENWAVE'], h0[i]['LIFE_ADJ'])
                del tt[i]
                del h0[i]

    #JRD ADDING NUVL
    elif (screen == 'NUVL' and (lifetime == 'ALL' or lifetime =='1')):                        # for NUVL
        #JRD
        print('condition 7')
        for i in list(h0.keys()):
            #JRD adding the G285M
            if (h0[i]['DETECTOR'] == 'NUV') and (h0[i]['OPT_ELEM'] == 'G230L'):
                #JRD
                print('SCREEN_DICT_OF_EXPOSURES 3: Keeping dictionary entry for ', h0[i]['ROOTNAME'], h0[i]['DETECTOR'], h0[i]['OPT_ELEM'], h0[i]['CENWAVE'], h0[i]['LIFE_ADJ'])
            else:
                #JRD
                print('SCREEN_DICT_OF_EXPOSURES 3: Deleting dictionary entry for ', h0[i]['ROOTNAME'], h0[i]['DETECTOR'], h0[i]['OPT_ELEM'], h0[i]['CENWAVE'], h0[i]['LIFE_ADJ'])
                del tt[i]
                del h0[i]
    else:
        #JRD
        print("SORRY BUDDY, BUT I CAN'T DO THAT SCREEN:", screen, " YOUR CHOICE IS INVALID")

    return tt, h0, h1




def write_output(out_dict, Linda):

    lifetime = out_dict['lifetime'] 

    if (out_dict.__contains__('G130M_final_lp'+lifetime)): write_fits_output(out_dict,'G130M_final_lp'+lifetime)
    if (out_dict.__contains__('G160M_final_lp'+lifetime)): write_fits_output(out_dict,'G160M_final_lp'+lifetime)
    if (out_dict.__contains__('G140L_final_lp'+lifetime)): write_fits_output(out_dict,'G140L_final_lp'+lifetime)
    if (out_dict.__contains__('FUVM_final_lp'+lifetime)): write_fits_output(out_dict,'FUVM_final_lp'+lifetime)
    if (out_dict.__contains__('FUVL_final_lp'+lifetime)): write_fits_output(out_dict,'FUVL_final_lp'+lifetime)

    #JRD adding NUV
    if (out_dict.__contains__('G185M_final_lp1')): write_fits_output(out_dict,'G185M_final_lp1')
    if (out_dict.__contains__('G225M_final_lp1')):
        print("JULIA, WRITING THE FILES!")
        write_fits_output(out_dict,'G225M_final_lp1')
    if (out_dict.__contains__('G285M_final_lp1')): write_fits_output(out_dict,'G285M_final_lp1')
    if (out_dict.__contains__('G230L_final_lp1')): write_fits_output(out_dict,'G230L_final_lp1')
    if (out_dict.__contains__('NUVM_final_lp1')): write_fits_output(out_dict,'NUVM_final_lp1')
    if (out_dict.__contains__('NUVL_final_lp1')): write_fits_output(out_dict,'NUVL_final_lp1')

    
    #### THIS IS A KLUDGE TO COPY THE SINGLE-GRATING FUV COADD INTO THE FUVM FILE IF ONLY ONE EXISTS 
    if (out_dict.__contains__('G130M_final_lp'+lifetime) and not out_dict.__contains__('G160M_final_lp'+lifetime)): 
        filename1 = out_dict['targname']+'_coadd_G130M_final_lp'+lifetime+'.fits.gz'
        filename2 = out_dict['targname']+'_coadd_FUVM_final_lp'+lifetime+'.fits.gz'
        command = 'cp -rp '+filename1+' '+filename2 
        os.system(command)
        #JRD
        print('I HAVE G130M BUT NOT G160M', command)
    if (out_dict.__contains__('G160M_final_lp'+lifetime) and not out_dict.__contains__('G130M_final_lp'+lifetime)):  
        filename1 = out_dict['targname']+'_coadd_G160M_final_lp'+lifetime+'.fits.gz'
        filename2 = out_dict['targname']+'_coadd_FUVM_final_lp'+lifetime+'.fits.gz'
        command = 'cp -rp '+filename1+' '+filename2 
        os.system(command)
        #JRD
        print('I HAVE G160M BUT NOT G130M', command)

    #if (out_dict.__contains__('FUVM_final_lp3')): # fourth best choice for "all"
    #    #JRD
    #    print("Using 4th best choice for all")
    #    filename1 = out_dict['targname']+'_coadd_FUVM_final_lp3.fits.gz'
    #    filename2 = out_dict['targname']+'_coadd_FUVM_final_all.fits.gz'
    #    command = 'cp -rp '+filename1+' '+filename2 
    #    os.system(command) 
    
    #if (out_dict.__contains__('FUVM_final_lp2')): # third best choice for "all"
    #    #JRD
    #    print("Using 3rd best choice for all")
    #    filename1 = out_dict['targname']+'_coadd_FUVM_final_lp2.fits.gz'
    #    filename2 = out_dict['targname']+'_coadd_FUVM_final_all.fits.gz'
    #    command = 'cp -rp '+filename1+' '+filename2 
    #    os.system(command) 
        
    #if (out_dict.__contains__('FUVM_final_lp1')): # second best choice for "all"
    #    #JRD
    #    print("Using 2nd best choice for all")
    #    filename1 = out_dict['targname']+'_coadd_FUVM_final_lp1.fits.gz'
    #    filename2 = out_dict['targname']+'_coadd_FUVM_final_all.fits.gz'
    #    command = 'cp -rp '+filename1+' '+filename2 
    #    os.system(command) 

    #if (out_dict.__contains__('FUVM_final_lp12')): # best choice for "all"
    #    #JRD
    #    print("Using 1st best choice for all")
    #    filename1 = out_dict['targname']+'_coadd_FUVM_final_lp12.fits.gz'
    #    filename2 = out_dict['targname']+'_coadd_FUVM_final_all.fits.gz'
    #    command = 'cp -rp '+filename1+' '+filename2 
    #    os.system(command) 


    #### SPECIAL PROCESSING 
    if (Linda):
        #JRD
        print(' this is for Linda')
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
    if 'lp4' in key_to_write_out: table_hdu.header['LP'] = 'LP4'

    thdulist = f.HDUList([prihdu, table_hdu])

    thdulist.writeto(out_dict['targname']+'_coadd_'+key_to_write_out+'.fits.gz',clobber=True)








