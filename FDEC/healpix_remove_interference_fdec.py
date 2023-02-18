#!/usr/bin/env python3

"""
Created on Feb 18 2023

PURPOSE:  
Removes RFI interference by subtracting a constant pattern in declination, FDEC.
Used for the postprocessing stage of the QUIJOTE MFI wide survey maps (Rubino-Martin et al. 2023). 

Based on the original IDL code healpix_remove_mfi_interference_constdec_with_mask.pro used 
in the QUIJOTE MFI wide survey pipeline.

@author: jalberto, elenadelahoz
"""

import numpy as np
import healpy as hp
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt


def get_function_constdec_for_onemap(map_in, mask, dec):

    npix = len(map_in)

    # Copy map and remove bad pixels
    map_out = np.copy(map_in)
    mapa    = np.copy(map_in*mask) # masked
    mapa[mapa==hp.UNSEEN]=0
    
    # 1) Stack in declination
    cutno0  = np.argwhere(map_out != 0.0)[:,0]
    cutno0m = np.argwhere(mapa != 0.0)[:,0]
    
    decmin  = min(dec[cutno0m])
    decmax  = max(dec[cutno0m])
    deltad_arcmin  = 20.0 # in arcmin. Tested for 1deg beam (=1/3 of FWHM)
    deltad  = deltad_arcmin / 60.0
    ndecs   = int(np.rint( (decmax-decmin)/deltad ) )
    if ndecs == 0: ndecs = 300
        
    vdec    = np.arange(ndecs) * deltad + decmin
      
    stack  = np.zeros(ndecs)
    for j in range(ndecs): 
       cut = np.argwhere((dec >= (vdec[j]-deltad/2.0)) * (dec < (vdec[j]+deltad/2.0)))[:,0]
       if len(cut) >= 2:
          pixset = mapa[cut]
          ind = np.argwhere(pixset != 0)[:,0]
          if len(ind) >= 1: pixset = pixset[ind]
          if len(pixset) >= 2: stack[j] = np.median(pixset)

              
    # 2) Project template back onto sky
    template2 = np.zeros(npix)
    if (len(cutno0) >= 1):         
       f2 = interp1d(vdec, stack, kind='cubic', fill_value="extrapolate")
       cut = np.argwhere((dec >= decmin) * (dec <= decmax))[:,0]
       res = f2(dec[cut])
       if len(cut) >= 1: template2[cut] = res
    
    # 3) Substract template
    if len(cutno0) >= 1: map_out[cutno0] -= template2[cutno0]
    
    return map_out


# Remove FDEC from healpix map in intensity only.
def healpix_remove_interference_fdec_int(mapa, maskI):
    npix  = mapa.shape[-1]
    nside = hp.npix2nside(npix)
    
    # map with declination information. 
    ipix           = np.arange(npix)
    theta, phi     = hp.pix2ang(nside, ipix, nest=False)
    r              = hp.Rotator(coord=['G','C'])
    theta_c, phi_c = r(theta, phi)
    dec            = 90.0 - np.rad2deg(theta_c)
    
    # Correcting intensity maps.
    print(" --> Computing FDEC for intensity map: ")         
    map_out = get_function_constdec_for_onemap(mapa, maskI, dec)
    
    return map_out

# Rotate Stokes Q-U maps G-->C.
def rotate_stokes_healpixmap(Q_in, U_in, gal2eq=False, eq2gal=False):
    
    if gal2eq: signo= 1.
    if eq2gal: signo=-1.

    npix       = Q_in.shape[-1]
    nside      = hp.npix2nside(npix)
    ipix       = np.arange(npix)
    r          = hp.Rotator(coord=['G','C'])
    theta, phi = hp.pix2ang(nside, ipix, nest=False)
    ang        = r.angle_ref(theta,phi)

    Q_out =  Q_in*np.cos(2.*ang) - U_in*np.sin(2.*ang*signo)
    U_out =  Q_in*np.sin(2.*ang*signo) + U_in*np.cos(2.*ang)
    
    return Q_out, U_out
        

# Remove FDEC for healpix map in intensity and polarization
def healpix_remove_interference_fdec_pol(mapa, maskI, maskP):
    npix  = mapa.shape[-1]
    nside = hp.npix2nside(npix)
    
    # map with declination information. 
    ipix           = np.arange(npix)
    theta, phi     = hp.pix2ang(nside, ipix, nest=False)
    r              = hp.Rotator(coord=['G','C'])
    theta_c, phi_c = r(theta, phi)
    dec            = 90.0 - np.rad2deg(theta_c)
    
    # Correcting intensity maps.
    print(" --> Computing FDEC for the intensity map: ")         
    mapI_out = get_function_constdec_for_onemap(mapa[0], maskI, dec)

    # Correcting polarization maps.
    print(" --> Computing FDEC for the polarization maps: ")
    mapQ = np.copy(mapa[1])
    mapU = np.copy(mapa[2])
    
    mapQr, mapUr = rotate_stokes_healpixmap(mapQ, mapU, gal2eq=True)
    
    mapQ_out = get_function_constdec_for_onemap(mapQr, maskP, dec)
    mapU_out = get_function_constdec_for_onemap(mapUr, maskP, dec)

    mapQfinal, mapUfinal = rotate_stokes_healpixmap(mapQ_out, mapU_out, eq2gal=True)

    map_out = np.asarray([ mapI_out, mapQfinal, mapUfinal ])
    
    return map_out


# Run anafast for intensity map only
def run_anafast(m1, masc):
    npix  = m1.shape[-1]
    nside = hp.npix2nside(npix)
    lmax  = 2*nside
    m1m   = np.copy(m1)
    for j in np.arange(3): m1m[j,:] = m1[j,:]*masc[:]
    cl    = hp.anafast(m1m, lmax=lmax, pol=True)
    fsky  = np.sum(masc)/float(npix)
    ell   = np.arange(lmax + 1)
    cl   *= 1.0 / fsky  # approximately corrects for sky mask
    return ell, cl


# MAIN code
if __name__ == "__main__":

    # Paths
    path_wmap    = '/Users/jalberto/WMAP/WMAP9/'
    path_quijote = '../'

    # Global params
    nside = 512
    npix  = hp.nside2npix(nside)
    
    # FDEC in polarization
    # WMAP K-band map
    ff    = path_wmap + 'wmap_band_iqumap_r9_9yr_K_v5.fits'
    mapa  = hp.read_map(ff,field=[0,1,2],nest=False)
    nmaps = mapa.shape[0]
    
    # Masks.  
    maskI = hp.read_map(path_quijote+'fdec/quijote_mask_fdec_int_nside512_ring.fits',nest=False)
    maskP = hp.read_map(path_quijote+'fdec/quijote_mask_fdec_pol_nside512_ring.fits',nest=False)
    masc  = hp.read_map(path_quijote+'masks/mask_quijote_ncp_lowdec_satband_nside512.fits',nest=False)
    
    # Apply FDEC 
    map_out = healpix_remove_interference_fdec_pol(mapa, maskI, maskP)

    # Display maps
    qmax = np.max(mapa[1])
    qmin = np.min(mapa[1])
    hp.mollview(mapa[1], norm='hist',title="Q Map in (orig)",sub=(2,2,1),max=qmax,min=qmin)
    hp.mollview(map_out[1], norm='hist',title="Q FDEC corrected",sub=(2,2,2),max=qmax,min=qmin)
    hp.mollview(mapa[1]-map_out[1], norm='hist', title='Q FDEC',sub=(2,2,3),max=qmax,min=qmin)
    hp.mollview(masc,title='Standard QUIJOTE mask',sub=(2,2,4))
    plt.show()

    # Compute and plot power spectra.
    ell, cl_in   = run_anafast(mapa, masc)
    ell, cl_out  = run_anafast(map_out, masc)
    ell, cl_fdec = run_anafast(mapa-map_out, masc)
    
    plt.plot(ell,cl_in[0],label='TT input map')
    plt.plot(ell,cl_out[0],label='TT filtered map',linestyle='dashed')
    plt.plot(ell,cl_fdec[0],label='TT FDEC',linestyle='dotted')
    plt.plot(ell,cl_in[1],label='EE input map')
    plt.plot(ell,cl_out[1],label='EE filtered map',linestyle='dashed')
    plt.plot(ell,cl_fdec[1],label='EE FDEC',linestyle='dotted')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim([20,700])
    plt.xlabel(r'$\ell$', fontsize=16)
    plt.ylabel(r'$C_\ell$ ' , fontsize=16)
    plt.title(r'Spectra')
    plt.legend(loc='upper right', ncol=2, labelspacing=0.1)
    plt.show()
    
