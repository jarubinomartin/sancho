#!/usr/bin/env python3

"""
Created on 15/Feb/2023
 
Script to illustrate the basic properties of QUIJOTE MFI wide maps in intensity:
A) Basic display of the maps
B) Example of noise properties from half mission maps
C) Noise levels in power spectra
D) Noise cross-correlation 11x13

@author: jalberto
"""

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import pymaster as nmt


# Read half-mission maps and prepares noise map
def prepare_noise_map(path,txtfreq):
    h1 = hp.read_map(path+'quijote_mfi_skymap_'+txtfreq+'ghz_512_dr1_half1.fits',field=[0],nest=False)
    h2 = hp.read_map(path+'quijote_mfi_skymap_'+txtfreq+'ghz_512_dr1_half2.fits',field=[0],nest=False)

    w1  = hp.read_map(path+'quijote_mfi_skymap_'+txtfreq+'ghz_512_dr1_half1.fits',field=[5],nest=False)
    w2  = hp.read_map(path+'quijote_mfi_skymap_'+txtfreq+'ghz_512_dr1_half2.fits',field=[5],nest=False)
    w1[np.isnan(w1)]=0
    w2[np.isnan(w2)]=0
    w1[w1<0]=0 # Healpy bad values
    w2[w2<0]=0
    
    w   = np.sqrt( (w1+w2)*(1./w1 + 1./w2) )
    n   = (h1 - h2)/w
    n[w1*w2==0]=0
    return(n)

# Run anafast for intensity map only
def run_anafast_int(m1, m2, masc):
    npix  = np.size(m1,0)
    nside = hp.npix2nside(npix)
    lmax  = 2*nside
    cl    = hp.anafast(m1*masc, map2=m2*masc, lmax=lmax, pol=False)
    fsky  = np.sum(masc)/float(npix)
    ell   = np.arange(lmax + 1)
    cl   *= 1.0 / fsky  # approximately corrects for sky mask
    return ell, cl


# Main code
nside = 512

# Path to QUIJOTE maps
path = '../'

# A) Display smoothed 1deg maps
mfi11s = hp.read_map(path+'quijote_mfi_smth_skymap_11ghz_512_dr1.fits',field=[0],nest=False)
mfi13s = hp.read_map(path+'quijote_mfi_smth_skymap_13ghz_512_dr1.fits',field=[0],nest=False)
mfi17s = hp.read_map(path+'quijote_mfi_smth_skymap_17ghz_512_dr1.fits',field=[0],nest=False)
mfi19s = hp.read_map(path+'quijote_mfi_smth_skymap_19ghz_512_dr1.fits',field=[0],nest=False)

hp.mollview(mfi11s,max=100,min=-5,norm='hist',title='11GHz')
hp.mollview(mfi13s,max=100,min=-5,norm='hist',title='13GHz')
plt.show()

# Analysis mask, apodized
masc = hp.read_map(path+'masks/mask_quijote_ncp_lowdec_satband_nside512.fits',field=[0],nest=False)
#masc = nmt.mask_apodization(masc_raw, 5.0, apotype="C2")


# B) Noise maps for 11 and 13GHz, from half mission maps
n11 = prepare_noise_map(path,'11')
n13 = prepare_noise_map(path,'13')

hp.mollview(n11*masc,max=5,min=-5,norm='hist',title='Noise 11GHz')
hp.mollview(n13*masc,max=5,min=-5,norm='hist',title='Noise 13GHz')
plt.show()


# C) Noise levels. Compare with Fig. 15 and 16 in Rubino-Martin et al. (2023).
mfi11 = hp.read_map(path+'quijote_mfi_skymap_11ghz_512_dr1.fits',field=[0],nest=False)
ell, clsky_11 = run_anafast_int(mfi11, mfi11, masc) 
ell, cl_11 = run_anafast_int(n11, n11, masc)

plt.plot(ell,clsky_11,label='signal 11GHz')
plt.plot(ell,cl_11,label='noise 11GHz half')
nbestfit = 2.56e-6 * (1 + (221.4/ell)**1.27) # Values from Table 11, Rubino-Martin et al. (2023).
plt.plot(ell,nbestfit,label='Fitted noise')
plt.xscale('log')
plt.yscale('log')
plt.xlim([20,700])
plt.ylim([1e-7,1e0])
plt.xlabel(r'$\ell$', fontsize=16)
plt.ylabel(r'$C_\ell$ [mK$^2$]' , fontsize=16)
plt.title(r'Spectra 11 GHz')
plt.legend(loc='upper right', ncol=2, labelspacing=0.1)
plt.show()



# D) Cross-spectrum 11x13. Compare to Fig. 17 in Rubino-Martin et al. (2023).
ell, cl_11_13 = run_anafast_int(n11, n13, masc) # cross-spectrum
ell, cl_11 = run_anafast_int(n11, n11, masc)
ell, cl_13 = run_anafast_int(n13, n13, masc)

rho = cl_11_13/np.sqrt(cl_11*cl_13) * 100.0 # in percentage
print('Average correlation 11x13 TT (20<l<200) = ', np.mean( rho[np.argwhere((ell >=20) & (ell <=200))] ))

plt.plot(ell,rho)
plt.xscale('log')
plt.yscale('linear')
plt.xlim([20,1000])
plt.ylim([0,100])
plt.xlabel(r'$\ell$', fontsize=16)
plt.ylabel(r'$\rho$ (%)', fontsize=16)
plt.title(r'TT 11x13')
plt.show()
