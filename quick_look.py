#! /usr/bin/env python

from astropy.io import fits as f
from astropy.io import ascii

import fitsio
from astropy.table import Table
from astropy.time import Time
from astropy.coordinates import SkyCoord
from datetime import datetime
import astropy.units as u
import collections
import scipy.io as scio
import numpy as np

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
import scipy.interpolate as sp
import glob
import os
import sys

def get_quick_look(): 

    mpl.rcParams['font.family'] = 'stixgeneral'

    window, wc, w_used = get_default_windows()
    
    #### Find the desired x1d's in the directory
    # first, check if all_exposures.txt exists. if it does, read it in
    if os.path.exists("all_exposures.txt"):
        exposure_cat = ascii.read("all_exposures.txt")
    else:
        # if all_exposures does not exist, import scrape_headers and make it.
        print "---> all_exposures.txt doesn't exist!!!!! so I'm going to make it now"
        from scrape_headers import make_exposure_catalog as make_exposure_catalog
        DATA_DIR = '.'
        dataset_list = glob.glob(os.path.join(DATA_DIR, '*x1d.fits.gz'))
        if len(dataset_list) == 0:
            print "did not find any x1d.fits.gz.gz files, trying for unzipped ones"
            print "I'm also going to gzip any x1d.fits.gz files I do find..."
            os.system("gzip *x1d.fits.gz")
            dataset_list = glob.glob(os.path.join(DATA_DIR, '*x1d.fits.gz'))
        if len(dataset_list) == 0:
            return "there's nothing here to make all_exposures with, trying to exit gracefully :-("
        exposure_cat = make_exposure_catalog(dataset_list)
    # second, cull for anything with Flag = 1
    mask = exposure_cat['Flag'] == 1
    exposure_cat = exposure_cat[mask]
    exposure_cat.sort('Cenwave')
    dataset_list = []
    for root in exposure_cat['Rootname']:
        dataset_list.append(root+"_x1d.fits.gz")
    if len(exposure_cat) == 0:
        return "the flags in all_exposures.txt told me to not do anything :-("
        
    # pathname = "file://"+os.getcwd()+"/"
    pathname = '' 
 
    ## does alias.txt exist?
    alias_file = "alias.txt"
    if (os.path.exists(alias_file)):
        af = open(alias_file, 'r')
        targname = str(af.readline()) 
        targname = targname.strip()
        af.close()
    else:
        print "---->>>> get_quick_look can't find "+alias_file+"!!!!! making one instead.....  <<<<----------"
        hdr0 = fitsio.read_header(dataset_list[0], 0) 
        targname = str(hdr0['TARGNAME']).strip()
        af = open(alias_file, 'w')
        af.write(targname)
        af.close()

    outfilename = targname + "_quicklook.html"

    exists = os.path.isfile(outfilename)
    if exists:
        command = "rm -f "+outfilename
        os.system(command)
        outfile = open(outfilename, "w")
    else:
        outfile = open(outfilename, "w")

    ra = exposure_cat['RA'][0]
    dec = exposure_cat['DEC'][0]
    mast_url = '"https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html?searchQuery='+str(ra)+','+str(dec)+'"'  
    tardescr = exposure_cat['Target Description'][0]
    coord = SkyCoord(ra=ra*u.degree, dec=dec*u.degree)
    targnamelink ="""<a href=".">"""+targname+"""</a>"""
    info = """<html>
    <head><title>"""+targname+"""</title>
    <h1 style="text-align:center;font-size:350%">"""+targnamelink+"""</h1></head>
    <body><p style="text-align:center;font-size=250%">"""+tardescr+"""<br>
    &alpha; = """+str(ra)+""", &delta; = """+str(dec)+""" (<a href="""+mast_url+""">"""+coord.to_string('hmsdms')+"""</a>)</font></p>
    <hr />
    """
    outfile.write(info)

    #### Now loop through all the files to figure out the demographics
    t, r = get_demographics(dataset_list)
    r_det = r.group_by(['detector'])
    mask = r_det.groups.keys['detector'] == 'FUV'
    r_fuv = r_det.groups[mask]
    mask = r_det.groups.keys['detector'] == 'NUV'
    r_nuv = r_det.groups[mask]

    if (exposure_cat['Detector'] == 'FUV').any():
        LAMBDA_MIN_FUV = np.min(r_fuv['minwave']-2)
        LAMBDA_MAX_FUV = np.max(r_fuv['maxwave']+2)
        MAX_FLUX_FUV = 1.05*np.max(r_fuv['maxflux'])
    if (exposure_cat['Detector'] == 'NUV').any():
        LAMBDA_MIN_NUV = np.min(r_nuv['minwave']-2)
        LAMBDA_MAX_NUV = np.max(r_nuv['maxwave']+2)
        MAX_FLUX_NUV = 1.05*np.max(r_nuv['maxflux'])    

    print 'MAX_FLUX', r['maxflux'] 

    t_pid = t.group_by(['PID'])
    info="""<p style="text-align:center;font-size=200%">Programs: """
    for k in t_pid.groups.keys:
        info = info+"""<a href="http://www.stsci.edu/cgi-bin/get-proposal-info?id="""+str(k[0])+"""&observatory=HST">"""+str(k[0])+"""</a> """
    info = info + """</p>"""
    outfile.write(info)

    width = 0.7
    colors = get_default_pid_colors(t_pid.groups.keys)
    print colors

    ### first do demographics for just the FUV exposures, then repeat first two for any NUV observations.
    
    if (exposure_cat['Detector'] == 'FUV').any():
        t_det = t.group_by(['detector'])
        mask = t_det.groups.keys['detector'] == 'FUV'
        t_fuv = t_det.groups[mask]
        t_fuv_pid = t_fuv.group_by(['PID'])
        print t_det
        print t_fuv
        
        figname = "fuv_exptime_histogram.png"
        title = "distribution of FUV exposure times by cenwave"
        plot_exptime_histogram(t_fuv, figname, colors=colors, t_pid=t_fuv_pid, width=width, title=title)
        addfig = r"""<br><img src='"""+pathname+figname+r"""' style="width:35%">"""
        outfile.write(addfig)

        figname = "fuv_exptime_fppos_histogram.png"
        title = "distribution of FUV exposure times by cenwave and FP-POS"
        plot_exptime_fppos_histogram(t_fuv, figname, colors=colors, t_pid=t_fuv_pid, width=width, title=title)
        addfig = r"""<img src='"""+pathname+figname+r"""' style="width:35%">"""
        outfile.write(addfig)
        
        figname = "fuv_lifetime_position_histogram.png"
        plot_lifetime_position_histogram(t_fuv, figname, colors=colors, t_pid=t_fuv_pid, width=width)
        addfig = r"""<img src='"""+pathname+figname+r"""' style="width:25%"><br>"""
        outfile.write(addfig)

    if (exposure_cat['Detector'] == 'NUV').any():
        t_det = t.group_by(['detector'])
        mask = t_det.groups.keys['detector'] == 'NUV'
        t_nuv = t_det.groups[mask]
        t_nuv_pid = t_nuv.group_by(['PID'])
        print t_det
        print t_nuv
        
        figname = "nuv_exptime_histogram.png"
        title = "distribution of NUV exposure times by cenwave"
        plot_exptime_histogram(t_nuv, figname, colors=colors, t_pid=t_nuv_pid, width=width, title=title)
        addfig = r"""<br><img src='"""+pathname+figname+r"""' style="width:35%">"""
        outfile.write(addfig)

        figname = "nuv_exptime_eppos_histogram.png"
        title = "distribution of NUV exposure times by cenwave and FP-POS"
        plot_exptime_fppos_histogram(t_nuv, figname, colors=colors, t_pid=t_nuv_pid, width=width, title=title)
        addfig = r"""<img src='"""+pathname+figname+r"""' style="width:35%">"""
        outfile.write(addfig)
        

    if (exposure_cat['Detector'] == 'FUV').any():
        print ">>>>>> I AM ASSUMING YOU HAVE COADDS ONLY IF YOU HAVE FUV DATA !!!!!!! <<<<<<<<<<<<<<<<"
        add_coadd = find_and_plot_coadds(targname, pathname, LAMBDA_MIN_FUV, LAMBDA_MAX_FUV, 0, MAX_FLUX_FUV, window=window, wc=wc)
        outfile.write(add_coadd)
        
    #### add legend
    info = """<p style="font-size=300%">Individual exposures<br>Legend: <b><font color="black">flux in black</font></b>, <b><font color="grey">errors in grey</font></b>, both smoothed over 6 pixels (~1 resel). S/N&equiv;median(flux/error), per ~1 resel, in shaded window.</p>"""
    outfile.write(info)

    #### Loop through files and make plots!
    time_flux = []
    for filename in dataset_list:
        with f.open(filename) as hdulist:   ##### want to change this to fitsio !!!!! ##### 
            hdr = hdulist[0].header
            data = hdulist[1].data
            if (hdulist[1].data == None):
                print "no data:",filename
                continue
            if (np.shape(data['flux'])[0] == 0):
                print np.shape(data['flux']),filename
                continue
            output_name = targname+"_"+str(hdr['proposid'])+"_"+str(hdr['cenwave'])+"_"+str(hdr['rootname'])+".png"
            print "FILENAME", output_name
            labeltext = str(hdr['rootname'])+ """ PID """+str(hdr['proposid'])+""" Visit """+str(hdr['linenum'].split('.')[0])+"""
            """+str(hdr['primesi'])+"""/"""+str(hdr['opt_elem']) + """/"""+str(hdr['cenwave']) + """/FPPOS"""+str(hdr['FPPOS'])+"""
            """+ hdulist[1].header['date-obs']+""", Exptime = """ +str(int(data['EXPTIME'][0]))+"""s"""
            time = np.average([hdulist[1].header['expstart'], hdulist[1].header['expend']])
            wavelength = data['wavelength']
            flux = data['flux']
            error = data['error']
            wgt = data['dq_wgt']
            if hdr['detector'] == 'FUV':
                LAMBDA_MIN, LAMBDA_MAX, MAX_FLUX = LAMBDA_MIN_FUV, LAMBDA_MAX_FUV, MAX_FLUX_FUV
            if hdr['detector'] == 'NUV':
                LAMBDA_MIN, LAMBDA_MAX, MAX_FLUX = LAMBDA_MIN_NUV, LAMBDA_MAX_NUV, MAX_FLUX_NUV
                if hdr['opt_elem'] == 'G230L':
                    ### want to get rid of stripe c
                    wavelength = data['wavelength'][0:2]
                    flux = data['flux'][0:2]
                    error = data['error'][0:2]
                    wgt = data['dq_wgt'][0:2]
            time_flux = plot_spectrum(output_name, wavelength, flux, LAMBDA_MIN, LAMBDA_MAX, 0, MAX_FLUX, \
                                      window=window, wc=wc, labeltext=labeltext, error=error, wgt=wgt, smooth=6, time=time, time_flux=time_flux)

            addfig = r"""<br><img src='"""+pathname+output_name+r"""' style="width:100%">"""
            print "writing out file ", addfig 
            outfile.write(addfig)


    if np.size(time_flux) > 0: 
        tf = Table(rows=time_flux, names=('mjd','flux','error','window'))
        plot_time_flux(tf, window=window, wc=wc)
        addfig = r"""<br><img src='"""+pathname+r"""time_flux.png' style="width:100%">"""
        outfile.write(addfig)
 
    outfile.write("</body></html>")
    outfile.close()
    message = """
    ~~~~~~~*~*~*~*~
    ~~~~~~~*~*~*~*~  all done!!!! spectra are fun!
    ~~~~~~~*~*~*~*~"""

#-----------------------------------------------------------------------------------------------------

def find_and_plot_coadds(targname, pathname, LAMBDA_MIN, LAMBDA_MAX, MIN_FLUX, MAX_FLUX, **kwargs):
    window = kwargs.get("window", get_default_windows()[0])
    wc = kwargs.get("wc", get_default_windows()[1])
    smooth = kwargs.get("smooth", 1)

    print "--->>>> assuming that coadds are named with ",targname,"!!! find_and_plot_coadds won't find them if not !!!!! <<<<------"
    print "--->>>> I'm gzipping your coadds because I won't find them if they aren't zipped !!!!! <<<<------"
    zipstring = "gzip "+targname+"_coadd*.fits.gz"
    os.system(zipstring)
    
    #### coadd legend
    # --->>>> should only be added if there is actually a coadd , fix!!!!! <<<<------ 
    info = """<p style="font-size=300%">Co-added spectra. Legend: <b><font color="black">flux in black</font></b>, <b><font color="grey">errors in grey</font></b>, both smoothed over 6 pixels (~1 resel). S/N&equiv;median(flux/error), per ~1 resel, in shaded window.</p>"""

    #### coadds????? ######
    addfig = ""
    coadd_exists = False
    if(os.path.exists(targname+'_coadd_G130M_final_all.fits.gz') or os.path.exists(targname+'_coadd_G160M_final_all.fits.gz')):
        coadd_exists = True
        output_name = targname+'_coadd_final_all.png'
        labeltext = """full coadd of """+targname+""" COS/FUV M"""
        if (os.path.exists(targname+'_coadd_G130M_final_all.fits.gz')): 
            # print '      ~~~ happy fuv data dance ~~~~' 
            # print '      YES!! I FOUND THE G130M coadd!'
            # print '      ~~~~ happy fuv data dance ~~~' 
            coadd = Table.read(targname+'_coadd_G130M_final_all.fits.gz') 
            print targname+':  quick_look opened  ' + targname+'_coadd_G130M_final_all.fits.gz for ', LAMBDA_MIN, LAMBDA_MAX 

    if(os.path.exists(targname + '_coadd_FUVM_final_all.fits.gz')):
        coadd_exists = True
        output_name = targname + '_coadd_final_all.png'
        labeltext = """full coadd of """+targname+""" COS/FUV M"""
        print '      ~~~ happy fuv data dance ~~~~' 
        print '      YES!! I FOUND THE FULL G130M+G160M COADD!'
        print '      ~~~~ happy fuv data dance ~~~' 
        coadd = Table.read(targname+'_coadd_FUVM_final_all.fits.gz')
        print targname+':  quick_look opened  ' + targname+'_coadd_FUVM_final_all.fits.gz for ', LAMBDA_MIN, LAMBDA_MAX
        plot_spectrum(output_name, coadd['WAVE'], coadd['FLUX'], 1100, 1900, 0, MAX_FLUX, window=window, wc=wc, labeltext=labeltext, error=coadd['ERROR'], smooth=6)
        addfig = addfig + r"""<br><img src='"""+pathname+output_name+r"""' style="width:100%">"""

    elif(os.path.exists(targname+'_coadd_G130M_final_all.fits.gz') or os.path.exists(targname+'_coadd_G160M_final_all.fits.gz')):
        coadd_exists = True
        output_name = targname+'_coadd_final_all.png'
        labeltext = """full coadd of """+str(targname)+""" COS/FUV M"""
        if (os.path.exists(targname+'_coadd_G130M_final_all.fits.gz')): 
            coadd = Table.read(targname+'_coadd_G130M_final_all.fits.gz') 
            print targname+':  quick_look opened  ' + targname+'_coadd_G130M_final_all.fits.gz for ', LAMBDA_MIN, LAMBDA_MAX 

            plot_spectrum(output_name, coadd['WAVE'], coadd['FLUX'], 1100, 1900, 0, MAX_FLUX, window=window, wc=wc, labeltext=labeltext, error=coadd['ERROR'], smooth=6)
            if (os.path.exists(targname+'_coadd_G160M_final_all.fits.gz')):
                sec_coadd = Table.read(targname+'_coadd_G160M_final_all.fits.gz') 
                print 'In the second if statement' 
                print targname+':  quick_look opened  ' + targname+'_coadd_G160M_final_all.fits.gz', LAMBDA_MIN, LAMBDA_MAX  
                plot_spectrum(output_name, coadd['WAVE'], coadd['FLUX'], 1100, 1900, 0, MAX_FLUX, \
                       window=window, wc=wc, labeltext=labeltext, error=coadd['ERROR'], smooth=6, \
                       overwave=sec_coadd['WAVE'], overflux=sec_coadd['FLUX'],overwgt=sec_coadd['DQ']+1.,overerror=sec_coadd['ERROR'],overcolor="black")
                print sec_coadd
            addfig = addfig + r"""<br><img src='"""+pathname+output_name+r"""' style="width:100%">"""
        else: 
            coadd = Table.read(targname+'_coadd_G160M_final_all.fits.gz') 
            print targname+':  quick_look opened  ' + targname+'_coadd_G160M_final_all.fits.gz', LAMBDA_MIN, LAMBDA_MAX  
            plot_spectrum(output_name, coadd['WAVE'], coadd['FLUX'], 1100, 1900, 0, MAX_FLUX, window=window, wc=wc, labeltext=labeltext, error=coadd['ERROR'], smooth=6)
            addfig = addfig + r"""<br><img src='"""+pathname+output_name+r"""' style="width:100%">"""

    if(os.path.exists(targname+'_coadd_G140L_final_all.fits.gz')): 
        coadd_exists = True
        output_name = targname+'_coadd_G140L_final_all.png'
        labeltext = """full coadd of """+str(targname)+""" COS/FUV L"""
        coadd = Table.read(targname+'_coadd_G140L_final_all.fits.gz') 
        print targname+':  quick_look opened  ' + targname+'_coadd_G140L_final_all.fits.gz', LAMBDA_MIN, LAMBDA_MAX  
        copy  = Table.read(targname+'_coadd_G140L_final_all.fits.gz') 
        copy['FLUX'] = 0.0 
        i_clip_short = np.where((copy['WAVE'] < 1000) & (copy['FLUX'] / copy['ERROR'] < 1.))  			#### screen out points at < 1100 with low S/N 
        i_clip_long  = np.where((copy['WAVE'] > 2000) & (copy['FLUX'] / copy['ERROR'] < 1.))  			#### screen out points at < 1100 with low S/N 
        print 'HELL', copy.keys() 
        copy['FLUX'][i_clip_short] = coadd['FLUX'][i_clip_short] 
        copy['FLUX'][i_clip_long] = coadd['FLUX'][i_clip_long] 
        plot_spectrum(output_name, coadd['WAVE'], coadd['FLUX'], 900, 2160, 0, MAX_FLUX, \
		window=window, wc=wc, labeltext=labeltext, error=coadd['ERROR'], smooth=6, color='red', \
		overwave=copy['WAVE'],overflux=copy['FLUX'],overwgt=copy['DQ']+1.,overerror=copy['ERROR'],overcolor='0.96' )
        addfig = addfig + r"""<br><img src='"""+pathname+output_name+r"""' style="width:100%">"""

    if(os.path.exists(targname+'_FUV_M_coadd.dat')):
        # if the FUV M coadd exists, use it 
        coadd_exists = True
        output_name = targname+"_FUV_M_coadd.png"
        coadd = scio.readsav(targname+'_FUV_M_coadd.dat')
        print targname+':  quick_look opened  ' + targname+'_FUV_M_coadd.dat'
        labeltext = """Colorado coadd of """+str(targname)+""" COS/FUV M"""    
        plot_spectrum(output_name, coadd['wave'], coadd['flux'], 1100, 1900, 0, MAX_FLUX, window=window, wc=wc, labeltext=labeltext, error=coadd['err'], smooth=6)
        addfig = addfig + r"""<br><img src='"""+pathname+output_name+r"""' style="width:100%">"""

    if(os.path.exists(targname+'_FUV_L_coadd.dat')):
        # if the FUV L coadd exists, use it
        coadd_exists = True
        output_name = targname+"_FUV_L_coadd.png"
        coadd = scio.readsav(targname+'_FUV_L_coadd.dat')
        print  'quick_look opened' + targname+'_FUV_L_coadd.dat'
        labeltext = """coadd of """+str(targname)+""" COS/FUV L"""
        plot_spectrum(output_name, coadd['wave'], coadd['flux'], 1100, 1900, 0, MAX_FLUX, window=window, wc=wc, labeltext=labeltext, error=coadd['err'], smooth=6)
        addfig = addfig + r"""<br><img src='"""+pathname+output_name+r"""' style="width:100%">"""
    
    if not coadd_exists:
        print "could not find a coadd, so sad"

    if coadd_exists:
        add_coadd = info + addfig
    else:
        add_coadd = ""
        
    return add_coadd

#-----------------------------------------------------------------------------------------------------

def smooth_spectrum(method, values, weights, scale):
    ## this theoretically allows users to choose if they want binning, Guassian smoothing, or moving average
    ## of course I only have moving average so far soooo.....:
    smooth = movingaverage(values, weights, scale)

    return smooth
      
#-----------------------------------------------------------------------------------------------------


def movingaverage(values, weights, window):
    indices = np.where(weights > 0)
    mvavg = np.convolve(values[indices], np.ones(window)/window,'same')                    
    return mvavg

#-----------------------------------------------------------------------------------------------------


def plot_spectrum(output_name, wavelength, flux, wavemin, wavemax, fluxmin, fluxmax, **kwargs):
    error = kwargs.get("error", None)
    labeltext = kwargs.get("labeltext", "")
    window = kwargs.get("window", get_default_windows()[0])
    wc = kwargs.get("wc", get_default_windows()[1])
    w_used = kwargs.get("w_used", get_default_windows()[2])
    smooth = kwargs.get("smooth", 1)
    wgt = kwargs.get("wgt", np.ones(np.shape(flux)))
    time = kwargs.get("time",-1)
    time_flux = kwargs.get("time_flux", [])

    if (fluxmax < 0): 
        fluxmax = 1.8 * np.mean(flux[np.where(flux > 0.)])
        print "fluxmax received by plot_spectrum is negative: adjusting it to: ", fluxmax 
    
    fig = plt.figure(figsize=(18, 6), dpi=300)
    ax = fig.add_subplot(111)
    if np.asarray(np.shape(np.shape(flux))) > 1:
        for i in range(np.shape(flux)[0]):
            indices = np.where(wgt[i] > 0)
            f = smooth_spectrum(1, flux[i], wgt[i], smooth)
            wave = wavelength[i][indices]
            ax.step(wave, f, lw=1, color='black')
    else:
        indices = np.where(wgt > 0)
        f = smooth_spectrum(1, flux, wgt, smooth)
        wave = wavelength[indices]
        ax.step(wave, f, lw=1, color='black')
    if "error" in kwargs:
        if np.asarray(np.shape(np.shape(flux))) > 1:
            for i in range(np.shape(flux)[0]):
                indices = np.where(wgt[i] > 0)
                e = smooth_spectrum(1, error[i], wgt[i], smooth)
                wave = wavelength[i][indices]
                ax.step(wave, e, lw=1, color='grey', alpha=0.6)
        else:
            e = smooth_spectrum(1, error, wgt, smooth)
            ax.step(wave, e, lw=1, color='grey', alpha=0.6)
        for w in range(np.shape(window)[0]):
            if np.asarray(np.shape(np.shape(flux))) > 1:
                for i in range(np.shape(flux)[0]):
                    sn_all = smooth_spectrum(1, np.ravel(flux[i]/error[i]), wgt[i], smooth)
                    f = smooth_spectrum(1, flux[i], wgt[i], smooth)
                    e = smooth_spectrum(1, error[i], wgt[i], smooth)
                    indices = np.where(wgt[i] > 0)
                    wave = wavelength[i][indices]
                    indices = np.where((f > 0) & (wave > window[w][0]) & (wave < window[w][1]))
                    if(np.shape(indices)[1] > 0):
                        #print "CALCULATING SN FOR ",w, wc[w]
                        medflux = np.median(f[indices])
                        err = np.sqrt(np.average(np.average((medflux-f)**2, weights=e)))/np.sqrt(np.size(indices))
                        sn = np.median(sn_all[indices]) * np.sqrt(smooth)  ## assuming smoothing by one resel
                        if (sn > 0):
                            plt.text(window[w][1], 0.75*fluxmax, "S/N="+"{:.1f}".format(sn), fontsize=10)
                            if w_used[w] == 0:
                                plt.axvspan(window[w][0], window[w][1], facecolor=wc[w], alpha=0.5)
                                w_used[w] = 1
                        print 'S to N 1:', time,medflux,err,w, sn 
                        if time > 0 and wc[w] != 'grey':
                            time_flux.append([time,medflux,err,w])
            else:
                f = smooth_spectrum(1, flux, wgt, smooth)
                e = smooth_spectrum(1, error, wgt, smooth)
                sn_all = smooth_spectrum(1, np.ravel(flux/error), wgt, smooth)
                indices = np.where((f > 0) & (wgt > 0) & (wavelength > window[w][0]) & (wavelength < window[w][1]))
            if(np.shape(indices)[1] > 0):
                print "CALCULATING SN FOR ",w, wc[w]
                medflux = np.median(f[indices])
                err = np.sqrt(np.average(np.average((medflux-f)**2, weights=e)))/np.sqrt(np.size(indices))
                sn = np.median(sn_all[indices]) * np.sqrt(smooth)
                print 'S to N 2a:', time, medflux, err, w, medflux / err, sn  
                if (sn > 0):
                    plt.text(window[w][1], 0.75*fluxmax, "S/N="+"{:.1f}".format(sn), fontsize=10)
                    if w_used[w] == 0:
                        plt.axvspan(window[w][0], window[w][1], facecolor=wc[w], alpha=0.5)
                        w_used[w] = 1
                print 'S to N 2b:', time,medflux,err,w, sn 
                if time > 0 and wc[w] != 'grey':
                    time_flux.append([time,medflux,err,w])
            # plt.axvspan(window[w][0], window[w][1], facecolor=wc[w], alpha=0.5)

    # do we want to overplot any other spectra?
    if "overwave" in kwargs and "overflux" in kwargs:
        print "overplotting something"
        overwavelength = kwargs.get("overwave", [])
        overflux = kwargs.get("overflux", [])
        overerror = kwargs.get("overerror", [])
        overwgt = kwargs.get("overwgt", np.ones(np.shape(overflux)))
        overcolor = kwargs.get("overcolor", "blue")
        if np.asarray(np.shape(np.shape(overflux))) > 1:
            for i in range(np.shape(overflux)[0]):
                indices = np.where(overwgt[i] > 0)
                f = smooth_spectrum(1, overflux[i], overwgt[i], smooth)
                overwave = overwavelength[i][indices]
                print "plot here", i
                ax.step(overwave, f, lw=1, color=overcolor, zorder=1)
            if "overerror" in kwargs:
                for i in range(np.shape(overflux)[0]):
                    e = smooth_spectrum(1, overerror[i], overwgt[i], smooth)
                    ax.step(overwave, e, lw=1, color=overcolor, zorder=10, alpha=0.2)
        else:
            print "or plot here"
            indices = np.where(overwgt > 0)
            f = smooth_spectrum(1, overflux, overwgt, smooth)
            overwave = overwavelength[indices]
            ax.step(overwave, f, lw=1, color=overcolor, zorder=2)
            if "overerror" in kwargs:
                e = smooth_spectrum(1, overerror, overwgt, smooth)
                ax.step(overwave, e, lw=1, color=overcolor, zorder=10, alpha=0.2)

    plt.xlim(wavemin, wavemax)
    plt.ylim(fluxmin, fluxmax)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)

    plt.ylabel(r'flux [erg / s / cm$^2/$ $\AA$]', fontsize=20)
    plt.xlabel(r'wavelength [$\AA$]',fontsize=20)
    plt.tight_layout()
    plt.text(0.99, 0.85, labeltext, fontsize=16, transform=ax.transAxes, horizontalalignment='right', color='black', bbox=dict(facecolor='white',linewidth=0.3,alpha=0.5), zorder=15)
    
    if os.path.isfile(output_name):   os.system("rm -f " + output_name) 
    plt.savefig(output_name)
    plt.close(fig)

    return time_flux

    
#-----------------------------------------------------------------------------------------------------
  

def get_default_windows():
    window = np.array([[1145,1155],
                       [1245,1255],
                       [1345,1355],
                       [1445,1455],
                       [1545,1555],
                       [1645,1655],
                       [1745,1755],
                       [1845,1855],
                       [1920,1930],
                       [2020,2030],
                       [2120,2130],
                       [2220,2230],
                       [2320,2330],
                       [2420,2430],
                       [2520,2530],
                       [2620,2630],
                       [2720,2730],
                       [2820,2830],
                       [2920,2930],
                       [3020,3030],
                       [3120,3130]])
    wc = ['indigo', 'violet', 'blue', 'green', 'gold', 'orange', 'red', 'darkred',
          'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey']

    w_used = np.zeros(21)
    
    return window, wc, w_used


def get_default_pid_colors(pids):
    color_list = ['purple','yellow','red','cyan','orange','green','blue','magenta','tan',
            'grey','pink','maroon','aqua','brown','greenyellow','olive','thistle']
    colors = {}
    for p in np.arange(len(pids)):
        colors[str(pids[p][0])] = color_list[p]
    return colors
    


#-----------------------------------------------------------------------------------------------------

def get_demographics(dataset_list):
    print "getting the demographics"

    # first, check if all_exposures.txt exists. if it does, read it in
    if os.path.exists("all_exposures.txt"):
        exposure_cat = ascii.read("all_exposures.txt")
    else:
        # if all_exposures does not exist, import scrape_headers and make it.
        print "---> all_exposures.txt doesn't exist!!!!! so I'm going to make it now"
        from scrape_headers import make_exposure_catalog as make_exposure_catalog
        exposure_cat = make_exposure_catalog(dataset_list)
    mask = exposure_cat['Flag'] == 1
    exposure_cat = exposure_cat[mask]

    # if scrape_headers can't be found, do it the hard way.
    # this can be done with a try except feature but things need re-arranging.
    if (False):
        LYA_MIN = 1206 ## should this depend on M vs L grating?
        LYA_MAX = 1226
        demographics = []
        ranges = []
        for filename in dataset_list:
            hdr = fitsio.read_header(filename, 0) 
            data, hdr1 = fitsio.read(filename, ext=1, header=True) 
            if (data == None):
                print "no data:", filename
                continue
            if (np.shape(data['FLUX'])[0] == 0):
                print np.shape(data['FLUX']), filename
                continue
            demographics.append((hdr['PROPOSID'], hdr['LINENUM'].split('.')[0], hdr['LINENUM'].split('.')[1], hdr['DETECTOR'], hdr['CENWAVE'], hdr['FPPOS'], hdr['LIFE_ADJ'], hdr['PR_INV_L'], hdr['ROOTNAME'], data['EXPTIME'][0]))
            indices = np.where((data["DQ_WGT"] > 0) & (data["DQ"] == 0) & ((data["WAVELENGTH"] > LYA_MAX) | (data["WAVELENGTH"] < LYA_MIN)))
            minwave = 1100 
            maxwave = 1900 
            maxflux = 1e-14 
            if(np.shape(indices)[1] > 0):
                if(hdr['OPT_ELEM'] == 'G140L'):
                    minwave = 900.
                    maxwave = 2160.
                    maxflux = 5.0 * np.mean(data["FLUX"][indices])   # changed from max of flux by JT 072015 
                if((hdr['OPT_ELEM'] == 'G130M') or (hdr['OPT_ELEM'] == 'G160M')): 
                    minwave = 1100
                    maxwave = 1900 
                    maxflux = 3.0 * np.median(data["FLUX"][indices])   # changed from max of flux by JT 072015 
                if hdr['DETECTOR'] == 'NUV': 
                    minwave = 1700
                    maxwave = 3000 
                    maxflux = 3.0 * np.mean(data["FLUX"][indices])   # changed from max of flux by JT 072015 
                ranges.append((minwave, maxwave, maxflux))
            print "Incorporating minwave, maxwave, maxflux from:    ", filename, hdr['DETECTOR'], hdr['OPT_ELEM'], hdr['CENWAVE'], minwave, maxwave, maxflux 
        t = Table(rows=demographics, names=('PID', 'Visit', 'expnum', 'detector', 'cenwave', 'FPPOS', 'lp', 'PI_NAME', 'rootname', 'exptime'))
        r = Table(rows=ranges, names=('minwave', 'maxwave', 'maxflux'), dtype=('f8', 'f8', 'f8'))


    # now from the exposure_cat table, find the stuff needed for t and r to be returned
    ranges = []
    demographics = []
    for i in range(len(exposure_cat)):
            demographics.append((exposure_cat['PropID'][i], exposure_cat['Detector'][i], exposure_cat['Cenwave'][i], exposure_cat['FPPOS'][i], exposure_cat['LP'][i], exposure_cat['PI Name'][i], exposure_cat['Rootname'][i], exposure_cat['Exptime'][i]))
            minwave = 1100 
            maxwave = 1900 
            maxflux = 1e-14 
            if(exposure_cat['Grating'][i] == 'G140L'):
                minwave = 900.
                maxwave = 2160.
                maxflux = 5.0 * exposure_cat['Mean Flux'][i]   # changed from max of flux by JT 072015 
            if((exposure_cat['Grating'][i] == 'G130M') or (exposure_cat['Grating'][i] == 'G160M')): 
                minwave = 1100
                maxwave = 1900 
                maxflux = 3.0 * exposure_cat['Median Flux'][i] # changed from max of flux by JT 072015 
            if exposure_cat['Detector'][i] == 'NUV': 
                minwave = 1680
                maxwave = 3250 
                maxflux = 5.0 * exposure_cat['Mean Flux'][i]   # changed from max of flux by JT 072015 
            ranges.append((minwave, maxwave, maxflux, exposure_cat['Detector'][i]))

    t = Table(rows=demographics, names=('PID', 'detector', 'cenwave', 'FPPOS', 'lp', 'PI_NAME', 'rootname', 'exptime'))
    r = Table(rows=ranges, names=('minwave', 'maxwave', 'maxflux', 'detector'), dtype=('f8', 'f8', 'f8','S3'))

    
    print "----------------------- THESE ARE THE DEMOGRAPHICS---------------------" 
    print t
    print r
    return t, r


#-----------------------------------------------------------------------------------------------------

def plot_exptime_histogram(t, figname, **kwargs):
    print "plotting exposure time by cenwave histogram"
    t_pid = kwargs.get("t_pid", t.group_by(['PID']))
    colors = kwargs.get("colors", get_default_pid_colors(t_pid.groups.keys))
    t_cen = kwargs.get("t_cen", t.group_by(['cenwave']))
    bins = kwargs.get("bins", t_cen.groups.keys['cenwave'])
    title = kwargs.get("title","distribution of exposure times by cenwave")
    width = kwargs.get("width", 0.7)

    fig = plt.figure(figsize=(6, 6), dpi=300)
    ax = fig.add_subplot(111)
    ind = np.arange(len(bins))
    bottom = np.zeros(len(ind))
    ## initialize cenwaves_dict
    cenwaves_dict = {}
    print range(len(t_pid.groups.keys)), len(colors)
    for k in range(len(t_pid.groups.keys)):
        for c in t_cen.groups.keys:
            cenwaves_dict[str(c[0])] = 0
        v = t_pid.groups[k].group_by(['cenwave'])
        for group in v.groups:
            cenwaves_dict[str(group['cenwave'][0])] += np.sum(group['exptime'])
        pid = str(v['PID'][0])
        p = plt.bar(ind, cenwaves_dict.values(), width, color=colors[pid], bottom=bottom, label=pid)
        bottom += cenwaves_dict.values()    
    plt.xlim(0-width/2., np.max(ind)+1.5*width)
    plt.xticks(ind+width/2., bins)
    plt.xlabel("cenwave", fontsize=16)
    plt.ylabel("exposure time [s]", fontsize=16)
    plt.title(title, fontsize=14)
    lg = ax.legend(loc='upper left', labelspacing=0.08, fontsize=16)
    lg.draw_frame(False)
    plt.tight_layout()
    if os.path.isfile(figname):
        command = "rm -f "+figname
        os.system(command)
    plt.savefig(figname)
    plt.close(fig)

    

#-----------------------------------------------------------------------------------------------------

def plot_exptime_fppos_histogram(t, figname, **kwargs): 
    print "plotting exposure time by cenwave+fppos histogram"
    t_pid = kwargs.get("t_pid", t.group_by(['PID']))
    colors = kwargs.get("colors", get_default_pid_colors(t_pid.groups.keys))
    t_cen = kwargs.get("t_cen", t.group_by(['cenwave']))
    bins = kwargs.get("bins", np.array(t_cen.groups.keys['cenwave']))
    title = kwargs.get("title","distribution of exposure times by cenwave and FP-POS")
    width = kwargs.get("width", 0.7)

    if len(bins) > 8:
        fig = plt.figure(figsize=(8, 6), dpi=300)
    else:
        fig = plt.figure(figsize=(6, 6), dpi=300)
    ax = fig.add_subplot(111)
    ind = np.arange(len(bins))
    bottom = {}
    for i in range(4):
        bottom[str(i)] = np.zeros(len(bins))
    ## initialize cenwaves_dict
    cenwaves_dict = {}
    for k in range(len(t_pid.groups.keys)):
        fp = t_pid.groups[k].group_by(['FPPOS'])
        for fg in fp.groups:
            for c in t_cen.groups.keys:
                cenwaves_dict[str(c[0])] = 0
            v = fg.group_by(['cenwave'])
            for group in v.groups:
                # print t_pid.groups.keys[k][0], group['cenwave'][0], np.sum(group['exptime'])
                cenwaves_dict[str(group['cenwave'][0])] += np.sum(group['exptime'])
                pid = str(v['PID'][0])
                p = plt.bar(ind+(fg['FPPOS'][0]-1)*width/4., np.array(cenwaves_dict.values()), width/4, color=colors[pid], bottom=bottom[str(fg['FPPOS'][0]-1)], label=pid)
            bottom[str(fg['FPPOS'][0]-1)] += cenwaves_dict.values()
    xticks = []
    xlbls = []
    xticks = np.concatenate((xticks, ind+width/2.))
    xlbls = np.concatenate((xlbls, np.round(np.array(bins))))
    va = np.zeros(len(ind))-0.03
    for i in range(len(ind)):
        for j in range(4):
            # print j, i+(2*j+1)*width/8.
            xticks = np.concatenate((xticks, [i+(2*j+1)*width/8.]))
            xlbls = np.concatenate((xlbls, [j+1]))
            va = np.concatenate((va, [0]))
    xticks = np.ravel(xticks)
    xlbls = np.ravel(xlbls)
    xlbls = xlbls.astype(int)
    va = np.ravel(va)
    # print xticks, xlbls
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlbls)
    for tick, y in zip( ax.get_xticklabels( ), va ):
        tick.set_y( y )
    plt.xlim(0-width/2., np.max(ind)+1.5*width)
    plt.xlabel("cenwave", fontsize=16)
    plt.ylabel("exposure time [s]", fontsize=16)
    plt.title(title, fontsize=14)
    plt.tight_layout()
    if os.path.isfile(figname):
        command = "rm -f "+figname
        os.system(command)
    plt.savefig(figname)
    plt.close(fig)

    
#-----------------------------------------------------------------------------------------------------

def plot_lifetime_position_histogram(t, figname, **kwargs): 
    print "plotting exposure time by lifetime position histogram"
    t_pid = kwargs.get("t_pid", t.group_by(['PID']))
    colors = kwargs.get("colors", get_default_pid_colors(t_pid.groups.keys))
    width = kwargs.get("width", 0.7)

      
    fig = plt.figure(figsize=(4, 6), dpi=300)
    ax = fig.add_subplot(111)
    t_lp = t.group_by(['lp'])
    if -1 in t_lp['lp']:
        bins = [-1, 1, 2, 3]
    else:
        bins = np.arange(3) + 1
    ind = np.arange(len(bins))
    bottom = np.zeros(len(bins))
    lp_dict = {}
    for k in range(len(t_pid.groups.keys)):
        for c in bins:
            lp_dict[str(c)] = 0
        v = t_pid.groups[k].group_by(['lp'])
        for group in v.groups:
            print t_pid.groups.keys[k][0], group['cenwave'][0], np.sum(group['exptime'])
            lp_dict[str(group['lp'][0])] += np.sum(group['exptime'])
        lp_sort = collections.OrderedDict(sorted(lp_dict.items()))
        pid = str(v['PID'][0])
        print k, ind, lp_sort.values(), colors[pid], bottom
        p = plt.bar(ind, lp_sort.values(), width, color=colors[pid], bottom=bottom, label=pid)
        bottom += lp_sort.values()
    plt.xlim(0-width/2., np.max(ind)+1.5*width)
    plt.xticks(ind+width/2., bins)
    plt.xlabel("lifetime position", fontsize=16)
    plt.ylabel("exposure time [s]", fontsize=16)
    plottitle = """distribution of FUV exposure
    times by lifetime position"""
    plt.title(plottitle, fontsize=14)
    plt.tight_layout()
    if os.path.isfile(figname):
        command = "rm -f "+figname
        os.system(command)
    plt.savefig(figname)
    plt.close(fig)

    
    
#-----------------------------------------------------------------------------------------------------

def plot_time_flux(tf, **kwargs):
    window = kwargs.get("window", get_default_windows()[0])
    wc = kwargs.get("wc", get_default_windows()[1])

    lpmoves = [56131.0, 57062.0]  ## COS FUV lifetime position moves, 2012-07-23 and 2015-02-09
    lpmoves_dec = Time(lpmoves, format='mjd').decimalyear
    
    tf['date'] = Time(tf['mjd'].data,format='mjd').datetime
    ## tf['date'] = Time(tf['mjd'].data,format='mjd').decimalyear
    
    tf_w = tf.group_by('window')
    print tf_w
        
    earliest = np.min(tf['date'])
    latest = np.max(tf['date'])

    fig = plt.figure(figsize=(18, 6), dpi=300)

    ax = fig.add_subplot(111)
    for w in range(np.shape(tf_w.groups.indices)[0]-1):
        ## print w, tf_w.groups[w]['window'],  np.shape(window)[0],  np.shape(wc),  np.shape(tf_w.groups.indices)[0]
        this_color_window = wc[tf_w.groups[w]['window'][0]]
        ## print this_color_window
        # ttt = Time(tf_w.groups[w]['mjd'].data,format='mjd').plot_date
        tf_w.groups[w].sort('date')
        ax.plot(list(Time(tf_w.groups[w]['mjd'].data,format='mjd').decimalyear), tf_w.groups[w]['flux'], color=wc[w])
        ##ax.plot_date(ttt, tf_w.groups[w]['flux'], color=wc[w])
        ax.errorbar(list(Time(tf_w.groups[w]['mjd'].data,format='mjd').decimalyear), tf_w.groups[w]['flux'], yerr=tf_w.groups[w]['error'], color=this_color_window)
        labeltext = str(window[w][0])+'$<\lambda<$'+str(window[w][1])+r'$\AA$'
        ax.scatter(list(Time(tf_w.groups[w]['mjd'].data,format='mjd').decimalyear), tf_w.groups[w]['flux'], s=50, color=this_color_window, alpha=0.5, label=labeltext)
    xr = np.array(ax.get_xlim())
    yr = np.array(ax.get_ylim())
    ax.plot([lpmoves_dec[0], lpmoves_dec[0]], [-1, 1], ls=':', color='k')
    ax.text(lpmoves_dec[0], 0.85*yr[1], r"move to LP2", fontsize=12)
    ax.plot([lpmoves_dec[1], lpmoves_dec[1]], [-1, 1], ls=':', color='k')
    ax.text(lpmoves_dec[1], 0.85*yr[1], r"move to LP3", fontsize=12)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.6f'))
    #lg = ax.legend(loc='upper left')
    #lg.draw_frame(False)
    plt.xlabel('year', fontsize=20)
    plt.ylabel(r'flux [erg/s/cm$^2$/$\AA$]', fontsize=20)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlim(xr[0], xr[1])
    plt.ylim(0, yr[1])
    plt.tight_layout()
    plt.savefig("time_flux.png")
    plt.close(fig)



#-----------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    message = get_quick_look()
    sys.exit(message)
