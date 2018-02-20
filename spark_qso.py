
"""Examples illustrating the use of plt.subplots().

This function creates a figure and a grid of subplots with a single call, while
providing reasonable control over how the individual plots are created.  For
very refined tuning of subplot creation, you can still use add_subplot()
directly on a new figure.
"""

import matplotlib.pyplot as plt
import numpy as np

from astropy.io import ascii
from astropy.table import Table 
import os 

x = np.linspace(0, 2 * np.pi, 400)
y = np.sin(x ** 2)

window = 16. 

fluxmax = 1e-14 

plt.close('all')
grid = np.random.rand(4, 4)
#fig, axes = plt.subplots(10, 10, figsize=(16, 9), # QSOs 
fig, axes = plt.subplots(2, 2, figsize=(24, 13.5), # QSOs 
                         subplot_kw={'xticks': [], 'yticks': [], 'axisbg':'black'})
fig.subplots_adjust(hspace=-0.1, wspace=0.00)
fig.patch.set_facecolor('black')

 
canonical = ascii.read('/grp/hst/HST_spectro/samples/QSOALS_spark.list')
canonical = ascii.read('/grp/hst/HST_spectro/samples/junk.list')
dirlist = canonical['targname'][np.where(canonical['flag'] == 1)]

for ax, dirname, nn in zip(axes.flat, dirlist[0:600], np.arange(576)): 

    os.chdir('/grp/hst/HST_spectro/datapile/'+dirname) 
    if(os.path.exists(dirname + '_coadd_FUVM_final_all.fits.gz')):
        print nn, '  found file', dirname + '_coadd_FUVM_final_all.fits.gz' 
        d = Table.read(dirname + '_coadd_FUVM_final_all.fits.gz') 
        wave = d['WAVE'][d['FLUX'] > 0.] 
        flux = d['FLUX'][d['FLUX'] > 0.] 
        flux[(wave > 1211) & (wave < 1218)] = 0. 
        flux = np.convolve(flux, np.ones(window,)/window,'same') 
        ax.plot(wave, flux, color='w') 
        ax.set_xlim([1150,1750])
        fluxmax = 2. * np.median(flux)  
        ax.set_ylim([0,fluxmax])
#        ax.text(1300, 0, dirname, fontsize=9, color='g')
    elif(os.path.exists(dirname + '_coadd_G140L_final_all.fits.gz')): 
        print nn, '  found file', dirname + '_coadd_G140L_final_all.fits.gz' 
        d = Table.read(dirname + '_coadd_G140L_final_all.fits.gz') 
        wave = d['WAVE'][d['FLUX'] > 0.] 
        flux = d['FLUX'][d['FLUX'] > 0.] 
        flux[(wave > 1211) & (wave < 1218)] = 0. 
        flux = np.convolve(flux, np.ones(window,)/window,'same') 
        ax.plot(wave, flux, color='w') 
        ax.set_xlim([1150,1750])
        fluxmax = 2. * np.median(flux)
        ax.set_ylim([0,fluxmax])
#        ax.text(1300, 0, dirname, fontsize=9, color='g')
    elif(os.path.exists(dirname + '_coadd_G130M_final_all.fits.gz')): 
        print nn, '  found file', dirname + '_coadd_G130M_final_all.fits.gz' 
        d = Table.read(dirname + '_coadd_G130M_final_all.fits.gz') 
        wave = d['WAVE'][d['FLUX'] > 0.] 
        flux = d['FLUX'][d['FLUX'] > 0.] 
        flux[(wave > 1211) & (wave < 1218)] = 0. 
        flux = np.convolve(flux, np.ones(window,)/window,'same') 
        ax.plot(wave, flux, color='w') 
        ax.set_xlim([1150,1750])
        fluxmax = 2. * np.median(flux) 
        ax.set_ylim([0,fluxmax])
#        ax.text(1300, 0, dirname, fontsize=9, color='r')
    else: 
        print 'OOPS, DID NOT FIND A FILE FOR ', dirname 
        ax.plot(x,y,color='r') 
#        ax.text(np.pi, 0, dirname, fontsize=9, color='g')


plt.tight_layout()
fig.tight_layout()
fig.patch.set_facecolor('black')
plt.show()
os.chdir('/grp/hst/HST_spectro/datapile') 
plt.savefig('all_qsos_sparks.png') 





