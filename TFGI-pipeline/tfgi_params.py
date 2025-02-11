#!/usr/bin/env python

"""
10/Feb/2025   

@author: jalberto

Basic parameters for TFGI. 

HISTORY:
* 10/02/2025 - original version. 

"""

import numpy as np
import os


# Reads TFGI basic information from etc files
def read_tfgi_info_etc_file(jd, path='/net/calp-nas/proyectos/quijote2/tfgi/etc/'):

    lista = [f for f in os.listdir(path) if f.startswith('qt2_pixel_masterfile.') & f.endswith('.txt')]
    lista.sort()

    for i in np.arange(len(lista)):
        ff = path + lista[i]
        with open(ff, 'r') as file:
            jdref = float(file.readline().strip())
        tab = np.loadtxt(ff,skiprows=1)

        if (jdref <= jd):
            horn_fp_loc = tab[:,0] 
            pix_id      = tab[:,1] 
            fem1        = tab[:,2] 
            fem2        = tab[:,3] 
            bem         = tab[:,4] 
            das         = tab[:,5] 
            iuse        = i

    print(" * READ_TFGI_INFO_ETC_FILE: using ",lista[iuse])

    return das, horn_fp_loc, pix_id


# TFGI reference angles. 
def tfgi_ref_angle(jd_obs, nhorns=7, path='/net/calp-nas/proyectos/quijote2/tfgi/etc/',das=None):
    '''
    Returns the reference angle for each horn and channel. Values in DEGREES. 
    Keywords: 
        * jd_obs: julian date of the observation
    '''

    jd_change = 2459677.5
    if jd_obs < jd_change:
        ff = 'phic_required_values_all_pixels_all_channels_wo_considering_RJH_angles_mft_maps_nov24_beforeApr2022.txt'
    else:
        ff = 'phic_required_values_all_pixels_all_channels_wo_considering_RJH_angles_mft_maps_nov24_afterApr2022.txt'
    
    # Reference information
    if das is None:
        das, horn_fp_loc, pix_id = read_tfgi_info_etc_file(jd_obs)

    # Read data from Mateo's files. Remember that these values correspond to -REF_ANGLE.
    idas, ch1, ch2, ch3, ch4 =np.loadtxt(path+ff,skiprows=2,unpack=True)

    # Reference angle (DEG).
    ang_ref = np.zeros(nhorns*4).reshape((nhorns,4))
    fac     = -1.0 # values in Mateo's tables correspond to minus ref_angle
    for i in np.arange(len(idas)):
        if (das == idas[i]).any():
            j = np.where(das == idas[i])[0]
            if j<nhorns:
                ang_ref[j,0] = fac*ch1[i]
                ang_ref[j,1] = fac*ch2[i]
                ang_ref[j,2] = fac*ch3[i]
                ang_ref[j,3] = fac*ch4[i]

    return ang_ref


