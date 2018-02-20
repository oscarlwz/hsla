#! /usr/bin/env python

"""
 return a list of targets in *.list form but only for a given PID
 takes in an exposure catalog
"""

import os
import sys 
from astropy.io import ascii 
from astropy.table import Table 
import fitsio
import glob
import argparse


def targs_by_program(catalog, pid_to_get): 

    #### first thing we do is open the MAST-provided data catalog
    targets = ascii.read(catalog) 
    targets['targname'] = targets['Target Name'] 
    targets['flag'] = 1 

    targets_for_this_pid = targets[targets['Proposal ID'] == pid_to_get] 

    unique_targets = targets_for_this_pid.group_by('Target Name') 
  
    bb = unique_targets[unique_targets.groups.indices[:-1]] 

    list_of_targs = Table([bb['flag'], bb['targname']]) 

    ascii.write(list_of_targs, '../samples/subsample_'+str(pid_to_get)+'.list')

    return unique_targets 


 
