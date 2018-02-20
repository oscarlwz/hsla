from astropy.io import ascii
import os 
import numpy as np 

def sample_download(canonical_filename): 

    os.chdir('/grp/hst/HST_spectro/datapile') 
	
    canonical = ascii.read('../samples/'+canonical_filename+'.list') 
    print canonical['flag']
    print canonical['targname']
    dirlist = canonical['targname'][np.where(canonical['flag'] == 1)]
    print "dirlist in driver"
    print dirlist 

    # create a directory 
    os.system('mkdir ../'+canonical_filename) 
 
    # move the sample files there 
    os.system('mv ../samples/'+canonical_filename+'.*  ../'+canonical_filename)

    for dirname in dirlist: 
        print "Moving target:  ", dirname
        if (os.path.isdir(dirname)): 
            print 'mv ' + dirname + ' ../' + canonical_filename  
            os.system('mv ' + dirname + ' ../' + canonical_filename) 

    os.chdir('..') 
    print 
    print 
    print 'tar --exclude="*_target.tar.gz" --exclude="DONE*" -cvf '+canonical_filename+'.tar ' + canonical_filename  
    os.system('tar --exclude="*_target.tar.gz" --exclude="DONE*" -cvf '+canonical_filename+'.tar ' + canonical_filename) 
    print 
    print 'gzip -f '+canonical_filename+'.tar '  
    os.system('gzip -f '+canonical_filename+'.tar ') 
  
    os.chdir(canonical_filename) 
    os.system('mv '+canonical_filename+'.* ../samples') 
    os.system('mv * ../datapile') 
    os.chdir('..') 
    os.system('rmdir '+canonical_filename) 
    os.chdir('/grp/hst/HST_spectro/distro/datapile') 

    #os.system("tar --exclude='*_target.tar.gz' --exclude='DONE*' -rvf ../samples/"+canonical_filename+".tar "+dirname) 

