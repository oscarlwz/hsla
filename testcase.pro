
spawn, 'rm *eps *pdf' 
@../../fifty-megapictures/coadd_script.pro 

g130 = mrdfits(targname[0]+'_coadd_G130M_final_all.fits',1) 
g160 = mrdfits(targname[0]+'_coadd_G160M_final_all.fits',1) 

!p.multi=[0,1,4] 

psopen, targname[0]+'_testcase.eps', /color, /encap 


; Panel 1 
plot, [0], [0], xrange=[1200, 1220], yrange=[-1e-16, 2.*median(flux)], ysty=1, $ 
	title=targname[0]
oplot, wave, flux, psym=10 
oplot, wave, err, psym=10, color=30 
oplot, g130.wave, g130.flux, psym=10, color=80 
oplot, g130.wave, g130.error, psym=10, color=200 

; Panel 2 
plot, [0], [0], xrange=[1300, 1320], yrange=[-1e-16, 2.*median(flux)], ysty=1 
oplot, wave, flux, psym=10 
oplot, wave, err, psym=10, color=30 
oplot, g130.wave, g130.flux, psym=10, color=80 
oplot, g130.wave, g130.error, psym=10, color=200 

; Panel 3 
plot, [0], [0], xrange=[1390, 1410], yrange=[-1e-16, 2.*median(flux)], ysty=1 
oplot, wave, flux, psym=10 
oplot, wave, err, psym=10, color=30 
oplot, g130.wave, g130.flux, psym=10, color=80 
oplot, g130.wave, g130.error, psym=10, color=200 

; Panel 4 
plot, [0], [0], xrange=[1540, 1560], yrange=[-1e-16, 2.*median(flux)], ysty=1 
oplot, wave, flux, psym=10 
oplot, wave, err, psym=10, color=30 
oplot, g160.wave, g160.flux, psym=10, color=80 
oplot, g160.wave, g160.error, psym=10, color=200 


psclose 
spawn, '/Users/tumlinson/Dropbox/COS-Halos/idl/bin/idlepstopdf.sh '+targname[0]+'_testcase.eps' 
