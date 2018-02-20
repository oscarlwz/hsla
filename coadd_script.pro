; IDL Version 8.2.3, Mac OS X (darwin x86_64 m64)
; Journal File for tumlinson@pancho.local
; Working directory: /Users/tumlinson/Dropbox/MESA/HSTSpecLegacy/datapile/white_dwarfs/WD1953-715
; Date: Sat Jul 25 12:08:20 2015
 
print, 'Running coadd_wd' 
pwd 
readcol, 'all_exposures.txt', flag, rootnames, targname, ra, dec, propid, piname, detector, segment, lifetime, grating, cenwave, fppos, exptime, nevents, extended, date, targdesc, $ 
		format='I,A,A,F,F,I,A,A,A,I,A,I,I,F,F,A,A,A', skip=1 

; identify and count the FUV M exposures 
	i_130 = where(grating eq 'G130M' and flag gt 0, n_130) 
	i_160 = where(grating eq 'G160M' and flag gt 0, n_160) 
	i_exposures = where((grating eq 'G130M' or grating  eq 'G160M') and cenwave ge 1100 and flag gt 0, n_exposures) 
	if (n_130 ge 1 and n_160 lt 1) then chan = 1 
	if (n_130 lt 1 and n_160 ge 1) then chan = 2 
	if (n_130 ge 1 and n_160 ge 1) then chan = 3 
	if (n_exposures ge 1) then coadd_x1d, wave, flux, err, chan=chan, files=rootnames[i_exposures] + '_x1d.fits', /no_scale, /no_align, method=1 
	if (n_exposures ge 1) then save, wave, flux, err, filename=targname[0]+'_FUV_M_coadd.dat' 
	if (n_exposures ge 1) then mwrfits, {wave:wave, flux:flux, err:err}, targname[0]+'_FUV_M_coadd.fits' 

; identify and count the FUV L exposures 
	i_140 = where(grating eq 'G140L' and flag gt 0, n_140) 
	if (n_140 ge 1) then chan = 4 
	if (n_140 ge 1) then coadd_x1d, wave, flux, err, chan=chan, /no_scale, files=rootnames[i_140] + '_x1d.fits', /no_scale, /no_align, method=1 
	if (n_140 ge 1) then save, wave, flux, err, filename=targname[0]+'_FUV_L_coadd.dat' 
	if (n_140 ge 1) then mwrfits, {wave:wave, flux:flux, err:err}, targname[0]+'_FUV_L_coadd.fits' 


