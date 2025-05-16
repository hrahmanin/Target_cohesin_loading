import pickle
import os
import time
import numpy as np
import ast


import cooltools
import cooltools.lib.plotting
#import functions

import shutil

import pandas as pd
import warnings
import h5py 
import glob

#from scorefunctions import *
import matplotlib.pyplot as plt
import sys


path_dict = {}

directory='/scratch1/rahmanin/target_loading_cohesin/target/sims_chn/'

for fname  in glob.glob(directory+'folder*'):
    path_dict[fname.split('/sims_chn/')[1][:]]= fname
path_dict = dict(sorted(path_dict.items()))


def peak_positions(boundary_lst_t, window_sizes=[1]):
    peak_monomers = np.array([])
    for i in window_sizes:
        inds_to_add=[boundary_lst_t[j]+i for j in range(len(boundary_lst_t))]
        peak_monomers = np.hstack((peak_monomers,inds_to_add))
    return peak_monomers.astype(int)

def FRiP(num_sites_t, lef_positions, peak_positions ):
    
    hist,edges = np.histogram(  lef_positions  , np.arange(num_sites_t+1) )
    return np.sum(hist[peak_positions] )/len(lef_positions)



def peak_positions(boundary_list, window_sizes=[1]):
    """
    Calculate peak positions based on a boundary_list within window_sizes.

    Args:
        boundary_list (list): List of boundary values.
        window_sizes (list, optional): List of window sizes. Defaults to [1].

    Returns:
        np.ndarray: Array containing peak positions.
    """
    peak_monomers = np.array([])

    for i in window_sizes:
        inds_to_add = [boundary + i for boundary in boundary_list]
        peak_monomers = np.hstack((peak_monomers, inds_to_add))

    return peak_monomers.astype(int)

window_size = 1
numx,numy = 1,len(path_dict)
rep = 1 
mon = 1000
site = 10
min_time = 500
hist_dict = {}
c=1
frip_file = open('frip_target_respoints_w1_main_tau17_all_density_changing_chn.csv','w') #Aug 23 is for 100 sep, and Aug 24 is for 1000 sep
frip_file.write('birth,life,deltactcf,clife,cof,sep,face,frip\n')
for name in list(path_dict.keys())[:]:
    print(name)
    params=[ast.literal_eval(i) for i in name.split('folder_')[1].split('_')[1::2]]
    face, back, clife, cof, life, slife, birth, deltactcf, pause, sep, site, mon, rep, step, vel = params

    c+=1

    mapN=mon*site
    lefs = h5py.File(path_dict[name]+'/LEFPositions.h5','r')["positions"]
    print(len(lefs))

    ctcfrightlist = np.array(h5py.File(path_dict[name]+'/LEFPositions.h5','r')['CTCF_sites_right'])
    ctcfleftlist = np.array(h5py.File(path_dict[name]+'/LEFPositions.h5','r')['CTCF_sites_left'])
    lst = np.array(list(ctcfrightlist) + list(ctcfleftlist))
    ### list of boundary elements on all replications
    lst_t = []
    for i in range(rep):
        lst_t += list(np.array(lst))
    print(lst_t)
    
    
    lef_lefts = lefs[min_time:,:,0].flatten()
    lef_rights = lefs[min_time:,:,1].flatten()
    lef_positions = np.hstack((lef_lefts,lef_rights))

    peak_monomers = peak_positions(lst_t,window_sizes=np.arange(-window_size,(window_size)+1) )
    frip = FRiP(mapN * rep, lef_positions, peak_monomers)
    print(frip)
    frip_file.write('%s,%s,%s,%s,%s,%s,%s,%s\n'%(birth,life,deltactcf,clife,cof,sep,face,frip))

frip_file.close()


