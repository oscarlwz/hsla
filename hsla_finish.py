#! /usr/bin/env python

"""
Runs the necessary steps to complete the HSLA version that was
started with hsla_setup.py. It assumes that the new alias file
created in hsla_setup.py has already been manually inspected
and the correct aliases chosen.

The following outputs are created:
1. 

Authors
-------
    Ben Sunnquist, 2017

Use
---
    This script can be run via the command line as such:
        
        python hsla_finish.py --n <datapile_v#> --o <datapile_v#>

    --n [Required]: The name to use for the new datapile directory for
    this HSLA release.
    --o [Required]: The name of the datapile directory for the previous
    HSLA release.

    This script can also be run within python:

        >>> import hsla_finish as hf
        >>> hf.hsla_finish(new_dir, old_dir)
    
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
import numpy as np
import pandas as pd

import defeat_aliases as d
import scrape_headers as sc
import ban_programs as b
import drive as dr

HSLA_DIR = '/grp/hst/HST_spectro/hsla_releases/'

# -----------------------------------------------------------------------------

def combine_alias(new_path, old_path):
    """
    Combines the new target alias file with the reference alias file from
    the previous HSLA release.

    Parameters
    ----------
    new_path : str
        The path to the new datapile directory used for this HSLA release.
    old_path : str
        The path to the old datapile directory used for the last HSLA release.

    Outputs
    -------
    {COS/STIS}_exposures_{date}_final.alias
        The new reference alias file.
    """

    # Read in the reference alias from the previous release
    f = glob.glob(os.path.join(old_path, '*final.alias'))
    ref_alias = Table.read(f[0], format='ascii.fixed_width_two_line')

    # Read in the alias file containing the new target names for this release
    f = glob.glob(os.path.join(new_path, '*new.alias'))
    new_alias = Table.read(f[0], format='ascii.fixed_width_two_line')

    # Update the reference alias with the new target names
    for i,t in enumerate(new_alias['Target Name']):
        # Get the new target name and its chosen alias
        new_targ = new_alias['Target Name'][i]
        new_alias0 = new_alias['alias0'][i]
        
        # Find entries from reference alias that share this chosen alias
        indexes = np.where(ref_alias['alias0'] == new_alias0)[0]
        
        # For each of these entries, add the new target name if there is an 
        # empty alias column.
        for ix in indexes:
            # Keep those with full alias entries alone
            if not ref_alias[ix]['alias5'] == '. . .':
                print('Alias entries full for {}.'.format(
                      ref_alias[ix]['Target Name']))
            # Add new target name to first empty alias column
            else:
                aliases = np.array([ref_alias[ix]['alias{}'.format(n)] 
                                   for n in range(0,6)])
                first_empty = np.flatnonzero(aliases == '. . .')[0]
                ref_alias[ix]['alias{}'.format(first_empty)] = new_targ    

    # Combine this now updated reference alias with the alias file containing
    # the new target names.
    final_alias = vstack([ref_alias, new_alias])

    # Write out the final alias file for this release
    final_alias.write(f[0].replace('new.alias', 'final.alias'), 
                      format='ascii.fixed_width_two_line')

# -----------------------------------------------------------------------------

def make_target_list(new_path):
    """
    Makes a master target list containing all of the target names with x1ds.

    Parameters
    ----------
    new_path : str
        The path to the new datapile directory used for this HSLA release.
    ins : str
        The instrument to use for this datapile (COS or STIS).

    Outputs
    -------
    all_targets.list
        The master target list.
    """

    # Find all target directories with x1ds
    targs = []
    for f in glob.glob(new_path + '/*'):
        if (
            (os.path.isdir(f)) & 
            (len(glob.glob(os.path.join(f, '*x1d.fits'))) > 0)
           ):
            targs.append(os.path.basename(f))

    # Make a dataframe with these targets including flags (1=include)
    targets = pd.DataFrame({})
    flags = np.ones(len(targs), dtype=int)
    targets['flag'] = flags
    targets['targname'] = targs

    # Write out the master target list
    targets.to_csv(os.path.join(new_path, 'all_targets_{}.list'.format(ins)), 
                   index=0, sep=' ')

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

    new_dir_help = ('The name to use for the new datapile directory for '
                    'this HSLA release.')
    old_dir_help = ('The name of the datapile directory for the previous '
                    'HSLA release.')
    ins_help = 'The instrument to use for this datapile (COS or STIS).'

    parser = argparse.ArgumentParser()
    parser.add_argument('--n', dest='new_dir', action='store', type=str, 
                        required=True, help=new_dir_help)
    parser.add_argument('--o', dest='old_dir', action='store', type=str, 
                        required=True, help=old_dir_help)
    parser.add_argument('--i', dest='ins', action='store', type=str, 
                        required=True, help=ins_help)
    args = parser.parse_args()

    return args

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

def hsla_finish(new_dir, old_dir, ins):
    """
    The main function of the hsla_finish module. See module docstring for
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
    ins : str
        The instrument to use for this datapile (COS or STIS).
    """

    # CHANGE DATAPILE_V# TO DATAPILE

    # Change directories to that containing the new HSLA version (other
    # scripts depend on being in this directory)
    new_path = os.path.join(HSLA_DIR, '{}_{}'.format(new_dir, ins))
    os.chdir(new_path)

    # The path to the last HSLA release
    old_path = os.path.join(HSLA_DIR, '{}_{}'.format(old_dir, ins))

    # Combine the old reference alias file with the alias file for new target
    # names contained in this release.
    combine_alias(new_path, old_path)
    print('Final alias file created.')

    # Move x1ds to their chosen alias directory
    d.defeat_aliases(glob.glob(os.path.join(new_path, '*final.alias'))[0])

    # MAKE SAMPLES DIRECTORY AND MOVE OLD ONE

    # Make a master target list
    # make_target_list(new_path)
    # print('Master target list created.')

    # CHANGE CWD TO DATAPILE AGAIN?
    # Create list of exposures/useful info and tables for each target directory
    #sc.scrape_headers('../samples/all_targets', 'blank', False)

    # Flag x1ds from known bad programs and failed visits
    #b.ban_programs('../samples/all_targets')

    # Run the driver to generate the coadds/quicklooks
    #dr.drive_targets(('../samples/all_targets', True))

    # Re-run scrape headers to add proper links and regenerate the various
    # tables and quicklooks
    #sc.scrape_headers('../samples/all_targets', 'blank', False)

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

if __name__ == '__main__':
    
    args = parse_args()

    hsla_finish(args.new_dir, args.old_dir, args.ins)

