#!/usr/bin/env python

"""
29/Jan/2026

@authors: jalberto, ubose

Routines for simulating CTOD files for TFGI

HISTORY:
* 29/01/2026 - original version. Based on previous version from Utsav

"""

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

from astropy.io import fits
from sancho_tfgi_io import read_tfgi_ctod2, read_tfgi_btod2, write_tfgi_ctod2
#from tfgi_params import *

import pysm3
import pysm3.units as u
from pysm3.models.cmb import CMBMap
from astropy.constants import k_B, c
import pysm3.models as models
import sys
import warnings


# Theoretical TFGI response
def sancho_tfgi_response(I, Q, U, parang, ichan, istgi=True):
    """
    Python translation of IDL function sancho_tfgi_response
    """

    # --- Argument check (IDL: n_params() lt 5)
    if any(v is None for v in (I, Q, U, parang, ichan)):
        print("")
        return -1

    # --- Phase states (degrees). Different order fo TGI and FGI
    phase = np.array([0.0, 90.0, 180.0, 270.0])
    if istgi==True:
        phase = np.array([0.0, 180.0, 90.0, 270.0])


    # --- Signs matrix: [Qsind, Qcosd, Usind, Ucosd]
    signs = np.zeros((4, 4), dtype=float)
    signs[0, :] = [-1,  0,  0,  1]
    signs[1, :] = [ 1,  0,  0, -1]
    signs[2, :] = [ 0,  1,  1,  0]
    signs[3, :] = [ 0, -1, -1,  0]

    # --- Response
    d2r = np.pi / 180.0
    ndata = I.size

    response = np.zeros((4, ndata), dtype=float)
    j = ichan

    for iph in range(4):
        phi = (phase[iph] + 2.0 * parang) * d2r
        response[iph, :] = 0.5 * (
            I
            + signs[j, 0] * Q * np.sin(phi)
            + signs[j, 1] * Q * np.cos(phi)
            + signs[j, 2] * U * np.sin(phi)
            + signs[j, 3] * U * np.cos(phi)
        )

    return response

# Generaring the PSYM model
def simulated_sky_pysm(nside=512,fwhm_arcmin=20*u.arcmin,freqGHz=30*u.GHz, model=['s1','d1','f1','a1'] ):
    '''
    -fwhm is in arcmin
    -freq is in Ghz'''
    warnings.filterwarnings("ignore")

    print(f'   Running PSYM for model {model} at nuGHz={freqGHz} and FWHM_arcmin={fwhm_arcmin}. ')
    sky = pysm3.Sky( nside=nside, preset_strings=model) #[ "s1",  "d1",  "f1",  "a1" ])
    # Output will be in uK_RJ by default
    maps = sky.get_emission(freqGHz)
    smoothed_map = pysm3.apply_smoothing_and_coord_transform(maps, fwhm=fwhm_arcmin)
    return smoothed_map

#Adding 1/f noise. Adaptaed from the routine of UBose
def simulate_1overf_noise_tod( tod, wn, fs, f_knee=0.1, alpha=1.0):

    nt = len(tod)
    dt = 1.0 / fs # time interval in seconds

    # white noise
    white_noise = wn*np.random.normal(0.0, 1.0, nt)

    # 1/f noise:
    if alpha==0:
         output = white_noise
    else:
        # fourier transofrm
        noise_fft = np.fft.rfft(white_noise)

        # Corresponding frequencies
        freqs = np.fft.rfftfreq(nt, d=dt)

        # making the white noise 1/f noise
        freqs_safe = np.maximum(freqs, freqs[1])  # avoid f=0

        filter_f = np.sqrt( 1.0 + (f_knee / freqs_safe) ** alpha)

        # Do not boost DC offset
        filter_f[0] = 0.0

        # Apply filter
        noise_fft *= filter_f

        #inverse transofrm
        output = np.fft.irfft(noise_fft, n=nt)

    return output


# Simulate one CTOD. Noise=0,1,2 for no noise, white and 1/f noise. Default is no noise
def simulate_one_ctod(ctod, map30, map40, nest=False, usenoise=0, verbose=True):

    # Info from the ctod file
    nhorns  = ctod['NHORNS']
    listpix = ctod['LISTPIX'] 
    data    = ctod['DATA']
    jd      = ctod['JD']
    gl      = ctod['GL']
    gb      = ctod['GB']
    parang  = ctod['PAR']
    msbin   = ctod['MSBIN']

    # Derived parameters for the CTOD
    istgi   = listpix<40
    dt      = msbin*1e-3 # time bin, s
    fs      = 1/dt       # sampling frequency, Hz

    # Display info
    if verbose:
        print( ' SIMULATING CTOD FILE. Parameters from input:')
        print(f'  nhorns = {nhorns}, msbin = {msbin} ms, f_samp = {fs} Hz. ')
        print(f'  noise case (0=no noise, 1=white, 2=correlated) = {usenoise}')

    # Info from the TFGI sky simulations
    npix   = len(map30[0])
    nside  = hp.npix2nside(npix)
    npix2  = len(map40[0])
    nside2 = hp.npix2nside(npix2)
    if not np.isclose( nside, nside2):
            raise ValueError("Inconsistent values of NSIDE in the simulations. ")

    # main part
    ctod_new = ctod.copy()
    ctod_new['DATA'] = ctod['DATA'].copy()
    for ihorn in range(nhorns):

        # IQU maps
        if istgi[ihorn]==True:
             II = map30[0]
             QQ = map30[1]
             UU = map30[2]
        else:
             II = map40[0]
             QQ = map40[1]
             UU = map40[2] 

        # Angles for that horn
        theta = np.deg2rad( 90.0 - gb[:,ihorn] )
        phi   = np.deg2rad( gl[:,ihorn] )
        ipix  = hp.ang2pix(nside, theta, phi, nest=nest)

        for ichan in np.arange(4):
                k = ihorn*4 + ichan
                
                # getting the noise level from the weights for Intensity data. 
                # In this case, sigma(I) coincides with sigma for a phase state
                noise_I = ctod['WEI_IQU'][:,0,k]
                sigma  = np.sqrt(1/noise_I) # Kelvin
                sigma_micro = sigma * 1e6   # microK

                # simulated sky signal
                sky = sancho_tfgi_response(II[ipix], QQ[ipix], UU[ipix], parang[:,ihorn], ichan, istgi=istgi[ihorn])

                # Simulated noise. Loop over phase states
                for iph in np.arange(4):

                    if usenoise == 0:  # No noise
                        noise = 0
                    elif usenoise == 1: # White noise
                        noise = simulate_1overf_noise_tod( II[ipix], sigma_micro, fs, f_knee=0.1, alpha=0.0)
                    elif usenoise == 2: # 1/f noise
                        noise = simulate_1overf_noise_tod( II[ipix], sigma_micro, fs, f_knee=10.0, alpha=1.75)
                    else:
                        print("Unexpected noise parameter. Assuming no noise")

                    # final output
                    data[:,iph,k] = sky[iph,:] + noise[:]
    ctod_new['DATA'] = data
    return ctod_new


# MAIN routine
def sancho_tfgi_simulate_ctod(ff, tail='ctod2', pathctod='/net/calp-nas/proyectos/quijote2/tfgi/ctod2/dec24/',
                              pathout='/net/calp-nas/proyectos/quijote2/tfgi/ctod2/sims/', usenoise=0):

    # Sky models
    nside = 2048 # 512 #2048
    freqs = [31.0, 41.0]*u.GHz
    fwhm  = [22.0, 18.0]*u.arcmin
    print(" (*) Preparing SKY models at TFGI frequencies. Units are uK_RJ by default. ")
    map30 = simulated_sky_pysm(nside=nside,fwhm_arcmin=fwhm[0],freqGHz=freqs[0] )
    map40 = simulated_sky_pysm(nside=nside,fwhm_arcmin=fwhm[1],freqGHz=freqs[1] )

    # Loop over files
    with open(ff, 'r') as index_file:
        for line in index_file:
            ffctod = pathctod + line.strip() + "." + tail
            ffout  = pathout + line.strip() + ".ctod2" # Simulated ctod2

            # Read ctod or btod
            if tail.startswith("ctod"):
                ctod = read_tfgi_ctod2(ffctod)
            elif tail.startswith("btod"):
                ctod = read_tfgi_btod2(ffctod)

            # Simulate new ctod
            simuctod = simulate_one_ctod(ctod, map30, map40, usenoise=usenoise)

            # Write ctod
            write_tfgi_ctod2(simuctod,ffout,overwrite=True)
    
    return


if __name__ == "__main__":
    
    # Test
    #sancho_tfgi_simulate_ctod('crab.txt', tail='ctod2', pathctod='/Users/jalberto/quijote/data/', pathout='./', usenoise=1)

    dothis = False
    if dothis:
        ff       = '/net/calp-nas/proyectos/quijote2/tfgi/list/crab/crab_2111_select_das24_pix23_allbtod.txt'
        pathctod = '/net/calp-nas/proyectos/quijote2/tfgi/ctod/dec24/' 
        pathout  = '/net/calp-nas/proyectos/quijote2/tfgi/ctod/sims/' 
        sancho_tfgi_simulate_ctod(ff, tail='ctod2', pathctod=pathctod, pathout=pathout, usenoise=1)

    dothis = True
    if dothis:
        ff       = '/net/calp-nas/proyectos/quijote2/tfgi/list/W44/W44_scans_for_destriper_05_17_times_RMS_allbtod2_das24.txt'
        pathctod = '/net/calp-nas/proyectos/quijote2/tfgi/ctod/dec24/' 
        pathout  = '/net/calp-nas/proyectos/quijote2/tfgi/ctod/sims/f_noise/' 
        sancho_tfgi_simulate_ctod(ff, tail='ctod2', pathctod=pathctod, pathout=pathout, usenoise=2)