#!/usr/bin/env python3

"""
Created on 15/Feb/2023
 
Script to illustrate the basic properties of QUIJOTE MFI wide maps in polarization:
A) Basic display of the maps 
B) Example of noise properties from half mission maps
C) Noise levels in power spectra.
D) Noise cross-correlation 11x13.

Maps can be downloaded from either the QUIJOTE web (https://research.iac.es/proyecto/quijote/), 
PLA (http://pla.esac.esa.int/pla/#maps) or LAMBDA (https://lambda.gsfc.nasa.gov/product/quijote/).

@author: jalberto
"""

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt


# Read half-mission maps and prepares noise map
def prepare_noise_map(path,txtfreq):
    comp = "IQU"
    ff1  = path+'quijote_mfi_skymap_'+txtfreq+'ghz_512_dr1_half1.fits'
    ff2  = path+'quijote_mfi_skymap_'+txtfreq+'ghz_512_dr1_half2.fits'
    h1   = hp.read_map(ff1,[c + "_STOKES" for c in comp],nest=False)
    h2   = hp.read_map(ff2,[c + "_STOKES" for c in comp],nest=False)

    w1  = hp.read_map(ff1,["WEI_"+c for c in comp],nest=False)
    w2  = hp.read_map(ff2,["WEI_"+c for c in comp],nest=False)
    w1[np.isnan(w1)]=0
    w2[np.isnan(w2)]=0
    w1[w1<0]=0  # Healpy bad values
    w2[w2<0]=0  
    
    w   = np.sqrt( (w1+w2)*(1./w1 + 1./w2) ) 
    n   = (h1-h2)/w
    n[w1*w2==0]=0
    return(n)

# Run anafast for intensity map only
def run_anafast(m1, m2, masc):
    npix  = np.size(m1,1)
    nside = hp.npix2nside(npix)
    lmax  = 2*nside
    fsky  = np.sum(masc)/float(npix)
    
    m1m   = np.copy(m1)
    m2m   = np.copy(m2)
    for i in np.arange(3):
        m1m[i,:] *= masc
        m2m[i,:] *= masc
    cl = hp.anafast(m1m, map2=m2m, lmax=lmax, pol=True)
    cl *= 1./fsky # approximately corrects for sky mask
    
    ell   = np.arange(lmax + 1)
    return ell, cl


# Main code
nside = 512
npix  = hp.nside2npix(nside)

# Path to QUIJOTE maps
path = '../'

# A) Display smoothed 1deg maps
mfi11s = hp.read_map(path+'quijote_mfi_smth_skymap_11ghz_512_dr1.fits',field=[0,1,2],nest=False)
mfi13s = hp.read_map(path+'quijote_mfi_smth_skymap_13ghz_512_dr1.fits',field=[0,1,2],nest=False)

hp.mollview(mfi11s[1],max=2,min=-2,norm='hist',title='Q 11GHz')
hp.mollview(mfi13s[1],max=2,min=-2,norm='hist',title='Q 13GHz')
plt.show()

# Analysis mask. No apodization applied. 
masc = hp.read_map(path+'masks/mask_quijote_ncp_lowdec_satband_nside512.fits',field=[0],nest=False)
hp.mollview(masc,title='Standard QUIJOTE mask')
plt.show()


# B) Noise maps for 11 and 13GHz, from half mission maps
n11 = prepare_noise_map(path,'11')
n13 = prepare_noise_map(path,'13')

hp.mollview(n11[1]*masc,max=5,min=-5,norm='hist',title='Noise Q 11GHz')
hp.mollview(n13[1]*masc,max=5,min=-5,norm='hist',title='Noise Q 13GHz')
plt.show()


# C) Noise levels. Compare with Fig. 15 and 16 in Rubino-Martin et al. (2023).
mfi11 = hp.read_map(path+'quijote_mfi_skymap_11ghz_512_dr1.fits',field=[0,1,2],nest=False)
ell, clsky_11 = run_anafast(mfi11, mfi11, masc) 
ell, cl_11    = run_anafast(n11, n11, masc)

# Simplified noise realization: anisotropic white noise.
wei   = hp.read_map(path+'quijote_mfi_skymap_11ghz_512_dr1.fits',['WEI_I','WEI_Q','WEI_U'],nest=False)
nsimu = np.random.randn(3,npix)*(1./np.sqrt(wei))
nsimu[wei<=0]=0 # zeros and bad values.
ell, cl_sim   = run_anafast(nsimu, nsimu, masc)


plt.plot(ell,clsky_11[1],label='signal EE 11GHz')
plt.plot(ell,cl_11[1],label='noise EE 11GHz half')
nbestfit = 6.13e-7 * (1 + (86.0/ell)**1.24) # Values from Table 11, Rubino-Martin et al. (2023).
plt.plot(ell,nbestfit,label='EE Fitted noise')
plt.plot(ell,cl_sim[1],label='EE anisotropic white noise sim')
plt.xscale('log')
plt.yscale('log')
plt.xlim([20,700])
plt.ylim([1e-7,3e-3])
plt.xlabel(r'$\ell$', fontsize=16)
plt.ylabel(r'$C_\ell^{EE}$ [mK$^2$]' , fontsize=16)
plt.title(r'EE spectra at 11 GHz')
plt.legend(loc='upper right', ncol=2, labelspacing=0.1)
plt.show()


# D) Cross-spectrum 11x13. Compare to Fig. 17 in Rubino-Martin et al. (2023).
ell, cl_11_13 = run_anafast(n11, n13, masc) # cross-spectrum
ell, cl_11 = run_anafast(n11, n11, masc)
ell, cl_13 = run_anafast(n13, n13, masc)

rho = cl_11_13/np.sqrt(cl_11*cl_13) * 100.0 # in percentage
print('Average correlation 11x13 EE (20<l<200) = ', np.mean( rho[1,np.argwhere((ell >=20) & (ell <=200))] ))
print('Average correlation 11x13 BB (20<l<200) = ', np.mean( rho[2,np.argwhere((ell >=20) & (ell <=200))] ))


plt.plot(ell,rho[1],label='EE')
plt.plot(ell,rho[2],label='BB')
plt.xscale('log')
plt.yscale('linear')
plt.xlim([20,1000])
plt.ylim([0,100])
plt.xlabel(r'$\ell$', fontsize=16)
plt.ylabel(r'$\rho$ (%)', fontsize=16)
plt.title(r'Polarization 11x13')
plt.legend(loc='upper right', ncol=2, labelspacing=0.1)
plt.show()
