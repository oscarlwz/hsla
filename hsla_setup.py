#! /usr/bin/env python

"""
Sets up a new HSLA version release in a new datapile directory that
includes older datasets from previous releases as well as datasets
that are new/reprocessed since this latest release. An updated alias 
file is also generated for the new targets. This should be manually
inspected, and the correct alias chosen before running hsla_finish.py.

The following outputs are created:
1. 

Authors
-------
    Ben Sunnquist, 2017

Use
---
    This script can be run via the command line as such:
        
        python hsla_setup.py --b <MM-DD-YYYY> --e <MM-DD-YYYY> --n <datapile_v#> --o <datapile_v#>

    --b [Required]: The date at which to start the new/reprocessed 
    file search from.
    --e [Required]: The date at which to end the new/reprocessed 
    file search from.
    --n [Required]: The name to use for the new datapile directory for
    this HSLA release.
    --o [Required]: The name of the datapile directory for the previous
    HSLA release.

    This script can also be run within python:

        >>> import hsla_setup as h
        >>> h.hsla_setup(begin, end, new_dir, old_dir)
    
    begin [Required]: str
        The date to start the search from (MM-DD-YYYY).
    end [Required]: str
        The date to end the search from (MM-DD-YYYY).
    new_dir [Required]: str
        The name to use for the new datapile directory for this HSLA release.
    old_dir [Required]: str
        The name of the datapile directory for the previous HSLA release.

Notes
-----
    Written in Python 2.7 to be consistent with the current HSLA software.
    It's assumed that the datapile directories, both old and new, are
    to be located in '/grp/hst/HST_spectro/hsla_releases/'.

    If things go wrong and you want to start over, you can quickly delete 
    a newly created datapile directory and all it contains with the 
    following commands (e.g. to delete the datapile_v4/ directory):

        >>> import shutil
        >>> shutil.rmtree('/grp/hst/HST_spectro/hsla_releases/datapile_v4/')
"""

import glob
import os
import shutil

import argparse
from astropy.io import ascii
from astropy.io import fits
from astropy.table import Table, vstack
import datetime as dt
import fitsio
import lxml
import numpy as np
import pandas as pd

import sort_targets as s
import target_alias as t

PUBLIC_DIR = '/ifs/archive/ops/hst/public/'
HSLA_DIR = '/grp/hst/HST_spectro/hsla_releases/'

# -----------------------------------------------------------------------------
def ban_visits(new_props, new_path, old_path):
    """
    Makes a table of failed COS/STIS visits according to the visit status page
    for each proposal. These visits will be excluded in the HSLA co-adds.
    For expediency, new failed visits are appended to the previous release's
    banned_visits.txt file.

    Parameters
    ----------
    new_props : array
        The proposal IDs of the new datasets.
    new_path : str
        The path to the new datapile directory used for this HSLA release.
    old_path : str
        The path to the old datapile directory used for the last HSLA release.

    Outputs
    -------
    banned_visits.txt
        A file containing all of the failed COS/STIS visits (with proposal 
        IDs).
    """

    # Read in the table of banned visits from the previous release
    banned_visits = ascii.read(os.path.join(old_path, 'banned_visits.txt'))

    for n,prop in enumerate(np.unique(new_props)):
        # Get the visit status info for this proposal
        url = 'http://www.stsci.edu/cgi-bin/get-visit-status?'
        url += 'id={}&markupFormat=html&observatory=HST'.format(prop)
        status_report = pd.read_html(url,header=0)[0]
        
        # Select only visits with COS/STIS data
        status_report = status_report[((status_report['Configs'].str.\
                                        contains('COS')) | 
                                       (status_report['Configs'].str.\
                                        contains('STIS')))]

        # Find any visits with a Status containing the word "Failed"
        failed_visits = status_report['Visit'][status_report['Status'].str.\
                                               contains('Failed')].values
        failed_visits = [str(v).zfill(2) for v in failed_visits]

        # Add any new failed visits to the old banned visits table
        for visit in failed_visits:
            select = ((banned_visits['Proposal']==prop) & 
                      (banned_visits['Visit']==visit))
            if len(banned_visits[select]) == 0:
                banned_visits.add_row([prop, visit])

        # progress tracker
        print('Done searching for failed visits in proposal {}.'.format(prop))
        print('{}/{} proposals checked.'.format(n, len(np.unique(new_props))))

    # Write out the updated banned visits table
    ascii.write(banned_visits, os.path.join(new_path, 'banned_visits.txt'))

# -----------------------------------------------------------------------------

def make_full_catalog(end, new_path):
    """
    Makes a full datapile catalog containing every COS/STIS spectrograph 
    observation from launch to the given end date.

    Parameters
    ----------
    end : str
        The date to end the search from (MM-DD-YYYY).
    new_path : str
        The path to write the CSV to.

    Outputs
    -------
    cos_exposures_{date}.txt
        A CSV file containing relevant information about each of the 
        observations.
    """

    # Query MAST for all COS spectrograph observations from launch to the 
    # given end date.
    url = 'https://archive.stsci.edu/hst/search.php?'
    url += '&sci_instrume=COS&sci_obs_type=spectrum&sci_aec=%' #for cal+sci
    url += '&sci_release_date=%3C{}'.format(end)
    url += '&max_records=50001&outputformat=CSV&coordformat=dec'
    url += '&selectedColumnsCsv=sci_data_set_name,sci_targname,sci_ra,sci_dec,'
    url += 'sci_expflag,sci_actual_duration,sci_instrume,sci_aper_1234,'
    url += 'sci_spec_1234,sci_central_wavelength,sci_pep_id,sci_status,'
    url += 'sci_target_descrip,sci_broad_category,sci_start_time&action=Search'
    print('Ignore the following 4 skipping line warnings. These will be added '
          'below. If there are more than 4, need to fix.')
    df_cos = pd.read_csv(url, error_bad_lines=False, skiprows=[1]) # skip dtype row

    # One person put an extra column in the target description for some 
    # asteroid obervations; this was skipped above with the 
    # 'error_bad_lines=False' option, so re-add it now.
    df_add = pd.read_csv(os.path.join(HSLA_DIR, 'code/additional_datasets.txt'),
                         names=df_cos.columns)
    df_cos = df_cos.append(df_add, ignore_index=True)

    # Query MAST for all STIS spectrograph observations from launch to the 
    # given end date.
    url = 'https://archive.stsci.edu/hst/search.php?'
    url += '&sci_instrume=STIS&sci_obs_type=spectrum&sci_aec=%' #for cal+sci
    url += '&sci_aper_1234=0.2X0.2,0.2x0.2FP*&sci_spec_1234=E230M,E140M'
    url += '&sci_release_date=%3C{}'.format(end)
    url += '&max_records=50001&outputformat=CSV&coordformat=dec'
    url += '&selectedColumnsCsv=sci_data_set_name,sci_targname,sci_ra,sci_dec,'
    url += 'sci_expflag,sci_actual_duration,sci_instrume,sci_aper_1234,'
    url += 'sci_spec_1234,sci_central_wavelength,sci_pep_id,sci_status,'
    url += 'sci_target_descrip,sci_broad_category,sci_start_time&action=Search'
    df_stis = pd.read_csv(url, error_bad_lines=False, skiprows=[1]) # skip dtype row

    # Combine the COS and STIS catalogs and write out the resulting catalog
    combined = pd.concat([df_cos, df_stis], ignore_index=True)
    t = dt.datetime.strptime(end, '%m-%d-%Y')
    t = t.strftime('%B%d_%Y')
    cat_name = os.path.join(new_path, 'all_exposures_{}.txt'.format(t))
    combined.to_csv(cat_name, index=False)

    return cat_name

# -----------------------------------------------------------------------------

def make_new_alias(new_targets, new_path, old_path):
    """
    Makes an alias file containing only new target names. This is the file
    that will be manually edited to select the correct alias name to use for
    the new targets.

    Parameters
    ----------
    new_targets : array
        Contains the target names of the new datasets (may have repeat names).
    new_path : str
        The path to the new datapile directory used for this HSLA release.
    old_path : str
        The path to the old datapile directory used for the last HSLA release.

    Outputs
    -------
    cos_exposures_{date}_new.alias
        An alias file only containing the alias entries for new target names
        that haven't been used in previous releases.
    """

    # Read in the reference alias from the previous release
    f = glob.glob(os.path.join(old_path, '*_final.alias'))
    reference_alias = ascii.read(f[0])
    
    # Read in the full alias file from this release
    f = glob.glob(os.path.join(new_path, '*_all.alias'))
    full_alias = ascii.read(f[0])

    # Create the new alias file
    cols = reference_alias.colnames
    dtypes = ['S' for c in range(len(cols))]
    new_alias = Table(names=cols, dtype=dtypes)

    # Append new target entries (those not contained in the alias file from
    # the previous release) to the new alias file
    for t in np.unique(new_targets):
        if t not in reference_alias['Target Name']:
            match = full_alias[full_alias['Target Name'] == t]
            new_alias = vstack([new_alias, match])

    # Write out the new alias file; add empty row if no entries
    if len(new_alias) == 0:
        new_alias.add_row(['' for c in range(len(cols))])
        print('None of the new datasets have new target names.')

    new_alias.write(f[0].replace('all.alias','new.alias'), 
                    format='ascii.fixed_width_two_line')

# -----------------------------------------------------------------------------

def make_new_datapile(begin, end, new, reprocessed, new_path, old_path):
    """
    Sets up a new datapile directory containing all relevant datasets,
    including new/reprocessed ones. The new/reprocessed lists from the new file
    query are checked to see if they are really new/reprocessed, and a text
    file is generated that summarizes what is new in this HSLA release.
   
    Parameters
    ----------
    begin : str
        The beginning date used in the date range for the new/reprocessed 
        file search (MM-DD-YYYY).
    end : str
        The end date used in the date range for the new/reprocessed file 
        search (MM-DD-YYYY).
    new : array
        Contains the new dataset names.
    reprocessed : array
        Contains the reprocessed dataset names. 
    new_path : str
        The path to the new datapile directory used for this HSLA release.
    old_path : str
        The path to the old datapile directory used for the last HSLA release.

    Outputs
    -------
    update_summary.txt
        A file that summarizes this new HSLA update. It includes a listing
        of all of the new/reprocessed datasets, their target names, and
        (for reprocessed files) the cal version and calibration files
        used for the old/new version of the dataset (i.e. for the version
        from the previous and new HSLA release). It also contains the
        version numbers of modules that are used in the HSLA creation.
    """

    # Open up the update_summary.txt file to summarize this HSLA update.
    outfile = open(os.path.join(new_path, 'update_summary.txt'), 'w')
    outfile.write('Update summary for HSLA version contained in %s.\n' 
                  % new_path)
    outfile.write('Date range used in the new/reprocessed file search for '
                  'this version: %s to %s.\n\n' % (begin, end))

    # Write the version numbers of modules that are used in the HSLA creation.
    import argparse
    outfile.write('argparse Version: %s\n' % argparse.__version__)
    import astropy
    outfile.write('astropy Version: %s\n' % astropy.__version__)
    import bokeh
    outfile.write('bokeh Version: %s\n' % bokeh.__version__)
    import fitsio
    outfile.write('fitsio Version: %s\n' % fitsio.__version__)
    import matplotlib
    outfile.write('matplotlib Version: %s\n' % matplotlib.__version__)
    import multiprocessing
    outfile.write('multiprocessing Version: %s\n' % multiprocessing.__version__)
    import numpy
    outfile.write('numpy Version: %s\n' % numpy.__version__)
    import pandas
    outfile.write('pandas Version: %s\n' % pandas.__version__)
    import scipy
    outfile.write('scipy Version: %s\n\n' % scipy.__version__)

    # Get the public archive paths to the x1d's of new datasets/asns. Copy
    # over any new files and record them in the update_summary.txt file.
    outfile.write('NEW FILES\n---------\n')
    n_new_files = 0
    for i,d in enumerate(new):
        d = d.lower()
        path = os.path.join(PUBLIC_DIR, d[0:4], d,  d + '_asn.fits')
        if os.path.exists(path):  # see if association exists
            data = fits.getdata(path)
            x1ds = data['MEMNAME'][data['MEMTYPE'] != 'PROD-FP']
            x1d_paths = [os.path.join(PUBLIC_DIR, d[0:4], dset.lower(), 
                         dset.lower()+'_x1d.fits') for dset in x1ds]
        else:  # for files without an association
            x1d_paths = [path.replace('_asn.fits', '_x1d.fits')]
        
        # If the x1d's didn't exist in the previous HSLA release, copy them
        # over to the new datapile directory and record their info.
        for f in x1d_paths:
            print('Working on {}.'.format(f))
            match = glob.glob(os.path.join(old_path, '*/', 
                                           os.path.basename(f)))
            if (len(match) == 0) & (os.path.exists(f)):
                shutil.copy(f, new_path)
                outfile.write('%s    %s\n' % (os.path.basename(f), 
                              fits.getheader(f, 0)['TARGNAME']))
                n_new_files += 1

        # progress tracker
        if i % 20 == 0:
            print('{} / {} new datasets checked/moved.'.format(i, len(new)))

    outfile.write('\n*** %s new x1ds in this release ***\n' % n_new_files)

    outfile.write('\n\nREPROCESSED FILES\n-----------------\n')

    # The keywords to search for why a COS/STIS file was reprocessed - these
    # include the cal version and calibration reference files.
    keywords_cos = ['FLATFILE','DEADTAB','BPIXTAB','SPOTTAB','GSAGTAB','HVTAB',
                    'BRFTAB','GEOFILE','DGEOFILE','TRACETAB','PROFTAB',
                    'TWOZXTAB','XWLKFILE','YWLKFILE','WALKTAB','PHATAB',
                    'PHAFILE','BADTTAB','XTRACTAB','LAMPTAB','DISPTAB',
                    'IMPHTTAB','FLUXTAB','WCPTAB','BRSTTAB','TDSTAB',
                    'SPWCSTAB','CAL_VER']

    keywords_stis = ['BPIXTAB','DARKFILE','PFLTFILE','LFLTFILE','PHOTTAB',
                     'IMPHTTAB','APERTAB','MLINTAB','WAVECAL','APDESTAB',
                     'SPTRCTAB','DISPTAB','INANGTAB','LAMPTAB','SDCTAB',
                     'XTRACTAB','PCTAB','MOFFTAB','WCPTAB','CDSTAB','ECHSCTAB',
                     'EXSTAB','RIPTAB','HALOTAB','TELTAB','SRWTAB','TDSTAB',
                     'TDCTAB','GACTAB','CAL_VER']

    # Get the public archive paths to the x1d's of reprocessed datasets/asns.
    # Copy over any reprocessed files and state the reason for the
    # reprocessing in the update_summary.txt file. Note: some files are both
    # new and reprocessed; these are accounted for here.
    n_rep_files = 0
    for i,d in enumerate(reprocessed):
        d = d.lower()
        path = os.path.join(PUBLIC_DIR, d[0:4], d,  d + '_asn.fits')
        if os.path.exists(path):  # see if association exists
            data = fits.getdata(path)
            x1ds = data['MEMNAME'][data['MEMTYPE'] != 'PROD-FP']
            x1d_paths = [os.path.join(PUBLIC_DIR, d[0:4], dset.lower(), 
                         dset.lower()+'_x1d.fits') for dset in x1ds]
        else:  # for files without an association
            x1d_paths = [path.replace('_asn.fits', '_x1d.fits')]
            
        # Copy over the reprocessed files and compare them to the old ones.
        for f in x1d_paths:
            print('Working on {}.'.format(f))
            if os.path.exists(f):
                # Some datasets that are both new/reprocessed were already
                # copied over in the 'new' step above.
                if not os.path.exists(os.path.join(new_path, 
                                                   os.path.basename(f))):
                    shutil.copy(f, new_path)
                
                # Find the matching dataset from the old HSLA release
                match = glob.glob(os.path.join(old_path, '*/', 
                                               os.path.basename(f)))
                if len(match) == 1:
                    old = fits.getheader(match[0], 0)
                    new = fits.getheader(f, 0)
                    diff = 0  # to see if anything changed

                    # Get relevant keywords based on instrument
                    if os.path.basename(f).startswith('l'):
                        keywords = keywords_cos
                    elif os.path.basename(f).startswith('o'):
                        keywords = keywords_stis
                    else:
                        print('Dataset doesnt appear to be COS or STIS.')

                    # Loop through the relevant keywords to see what changed
                    for k in keywords:
                        try:
                            new_kwd = new[k]
                            old_kwd = old[k]
                            if new_kwd != old_kwd:
                                if diff == 0:
                                    outfile.write('%s (%s)\n' % 
                                    (os.path.basename(f), new['TARGNAME']))
                                
                                outfile.write('%s: %s    %s\n' % (k, 
                                              old_kwd, new_kwd))
                                diff += 1
                        # not all keywords exist for every file
                        except KeyError:
                            pass
                    if diff != 0:
                        n_rep_files += 1  # count those with changes
                        outfile.write('\n')
                else:
                    print('No matching dataset (or, unlikely, multiple '
                          'datasets found) for {} in the previous HSLA '
                          'release directory. Its likely that this dataset '
                          'was both new and reprocessed in the specified '
                          'time range.'.format(os.path.basename(f)))

        # progress tracker
        if i % 100 == 0:
            print('{} / {} reprocessed datasets checked/moved.'.format(i, 
                  len(reprocessed)))

    outfile.write('*** %s x1ds with cal version/reference file '
                  'differences in this release ***\n' % n_rep_files)

    # Copy over all of the remaining files from the previous HSLA release.
    old_files = glob.glob(os.path.join(old_path, '*/*x1d.fits'))
    for i,f in enumerate(old_files):
        if not os.path.exists(os.path.join(new_path, os.path.basename(f))):
            shutil.copy(f, new_path)

        # progress tracker
        if i % 100 == 0:
            print('{} / {} old files checked/moved.'.format(i, len(old_files)))

    # Close the update_summary.txt file.
    outfile.close()

# -----------------------------------------------------------------------------

def new_file_query(begin, end):
    """
    Queries MAST for new and reprocessed COS and STIS files within a certain 
    date range.

    Parameters
    ----------
    begin : str
        The date to start the search from (MM-DD-YYYY).
    end : str
        The date to end the search from (MM-DD-YYYY).

    Returns
    -------
    new_datasets : array
        Contains the new dataset names.
    reprocessed_datasets : array
        Contains the reprocessed dataset names.
    new_targets : array
        Contains the target names of the new datasets.
    new_props : array
        Contains the proposal IDs of the new datasets.
    """

    # Query for all new COS files
    url = 'https://archive.stsci.edu/hst/search.php?'
    url += '&sci_instrume=COS&sci_obs_type=spectrum&sci_aec=%' #for cal+sci
    url += '&sci_release_date={}..{}'.format(begin, end)
    url += '&max_records=50001&outputformat=CSV'
    url += '&selectedColumnsCsv=sci_data_set_name,sci_targname,sci_pep_id'
    url += '&action=Search'
    new_cos = pd.read_csv(url, skiprows=[1]) # skip datatype row
    if len(new_cos)==0:
        new_cos = pd.DataFrame({'Dataset':[], 'Target Name':[], 'Proposal ID':[]})
        print('No new COS files in this time range.')

    # Query for all new STIS files
    url = 'https://archive.stsci.edu/hst/search.php?'
    url += '&sci_instrume=STIS&sci_obs_type=spectrum&sci_aec=%'
    url += '&sci_release_date={}..{}'.format(begin, end)
    url += '&sci_aper_1234=0.2X0.2,0.2x0.2FP*&sci_spec_1234=E230M,E140M'
    url += '&max_records=50001&outputformat=CSV'
    url += '&selectedColumnsCsv=sci_data_set_name,sci_targname,sci_pep_id'
    url += '&action=Search'
    new_stis = pd.read_csv(url, skiprows=[1]) # skip datatype row
    if len(new_stis)==0:
        new_stis = pd.DataFrame({'Dataset':[], 'Target Name':[], 'Proposal ID':[]})
        print('No new STIS files in this time range.')

    # Query for all COS reprocessed files
    url = 'https://archive.stsci.edu/hst/search.php?'
    url += '&sci_instrume=COS&sci_obs_type=spectrum&sci_aec=%'
    url += '&sci_archive_date={}..{}&sci_release_date=%3C{}'.format(begin, 
           end, end)
    url += '&max_records=50001&outputformat=CSV'
    url += '&selectedColumnsCsv=sci_data_set_name,sci_targname&action=Search'
    reprocessed_cos = pd.read_csv(url, skiprows=[1])
    if len(reprocessed_cos)==0:
        reprocessed_cos = pd.DataFrame({'Dataset':[], 'Target Name':[]})
        print('No reprocessed COS files in this time range.')

    # Query for all STIS reprocessed files
    url = 'https://archive.stsci.edu/hst/search.php?'
    url += '&sci_instrume=STIS&sci_obs_type=spectrum&sci_aec=%'
    url += '&sci_archive_date={}..{}&sci_release_date=%3C{}'.format(begin, 
           end, end)
    url += '&sci_aper_1234=0.2X0.2,0.2x0.2FP*&sci_spec_1234=E230M,E140M'
    url += '&max_records=50001&outputformat=CSV'
    url += '&selectedColumnsCsv=sci_data_set_name,sci_targname&action=Search'
    reprocessed_stis = pd.read_csv(url, skiprows=[1])
    if len(reprocessed_stis)==0:
        reprocessed_stis = pd.DataFrame({'Dataset':[], 'Target Name':[]})
        print('No reprocessed STIS files in this time range.')

    # Combine the COS and STIS results and get rid of bad targets
    new = pd.concat([new_cos, new_stis], ignore_index=True)
    if len(new) != 0:
        new = new[((new['Target Name'] != 'ANY') & 
                   (new['Target Name'] != 'WAVE') & 
                   (new['Target Name'] != 'CCDFLAT'))].reset_index(drop=True)
    reprocessed = pd.concat([reprocessed_cos, reprocessed_stis], ignore_index=True)
    if len(reprocessed) != 0:
        reprocessed = reprocessed[((reprocessed['Target Name'] != 'ANY') & 
                                   (reprocessed['Target Name'] != 'WAVE') & 
                                   (reprocessed['Target Name'] != 'CCDFLAT'))].reset_index(drop=True)

    # Return the new/reprocessed dataset names and new target names/prop IDs
    new_datasets = new['Dataset'].values
    reprocessed_datasets = reprocessed['Dataset'].values
    new_targets = new['Target Name'].values
    new_props = new['Proposal ID'].values

    return new_datasets, reprocessed_datasets, new_targets, new_props

# -----------------------------------------------------------------------------

def parse_args():
    """
    Parses command line arguments.
    
    Returns
    -------
    args : object
        Contains the begin/end new file search dates and the current/new HSLA 
        release directory names.
    """

    begin_help = 'The date to start searching for new COS files (MM-DD-YYYY)'
    end_help = 'The date to stop searching for new COS files (MM-DD-YYYY)'
    new_dir_help = ('The name to use for the new datapile directory for '
                    'this HSLA release.')
    old_dir_help = ('The name of the datapile directory for the previous '
                    'HSLA release.')

    parser = argparse.ArgumentParser()
    parser.add_argument('--b', dest='begin', action='store', type=str, 
                        required=True, help=begin_help)
    parser.add_argument('--e', dest='end', action='store', type=str, 
                        required=True, help=end_help)
    parser.add_argument('--n', dest='new_dir', action='store', type=str, 
                        required=True, help=new_dir_help)
    parser.add_argument('--o', dest='old_dir', action='store', type=str, 
                        required=True, help=old_dir_help)
    args = parser.parse_args()

    return args

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

def hsla_setup(begin, end, new_dir, old_dir):
    """
    The main function of the hsla_setup module. See module docstring for
    further details.

    Parameters
    ----------
    begin : str
        The date to start the search from (MM-DD-YYYY).
    end : str
        The date to end the search from (MM-DD-YYYY).
    new_dir : str
        The name to use for the new datapile directory for this HSLA release.
        Historically has been 'datapile_v#'.
    old_dir : str
        The name of the datapile directory for the previous HSLA release.
        Historically has been 'datapile'.
    """

    # Make a new datapile directory for this version and change working
    # directoy to it (other scripts depend on being in this directory).
    new_path = os.path.join(HSLA_DIR, new_dir)
    if not os.path.exists(new_path):
        os.mkdir(new_path)
    os.chdir(new_path)

    # The path to the last HSLA release
    old_path = os.path.join(HSLA_DIR, old_dir)

    # Get a list of the new and reprocessed COS spectrograph datasets
    # (and new target names/prop IDs) that were observed within a certain 
    # time range.
    new, reprocessed, new_targs, new_props = new_file_query(begin, end)
    print('{} new files.'.format(len(new)))
    print('{} reprocessed files.'.format(len(reprocessed)))

    # Make a table containing all failed COS visits. These visits will be
    # excluded in the co-adds.
    ban_visits(new_props, new_path, old_path)
    print('Banned visits table created.')

    # Set up the new datapile directory
    make_new_datapile(begin, end, new, reprocessed, new_path, old_path)
    print('New datapile created.')

    # Make a full exposure catalog of every COS spectrograph observation
    # from launch to the given end date.
    cat_name = make_full_catalog(end, new_path)
    print('Full catalog created.')

    # Sort the datasets into target directories using the full catalog
    s.sort_targets(cat_name)
    print('All datasets sorted into target directories.')

    # Make an alias file based on the full catalog
    t.target_alias(cat_name)
    print('Full alias file created.')

    # Make an alias file containing only the new target names
    make_new_alias(new_targs, new_path, old_path)
    print('New alias file created.')

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

if __name__ == '__main__':
    
    args = parse_args()

    hsla_setup(args.begin, args.end, args.new_dir, args.old_dir)
