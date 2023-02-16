#!/usr/bin/env python3

"""
Created on 16/Feb/2023
 
Script to read some basic properties of the QUIJOTE MFI DR1 RIMO file:
A) Basic information (frequencies, beam FWHM,...)
B) Bandpass
C) Beam transfer functions
D) Beam radial profile

RIMO file can be downloaded from either the QUIJOTE web (https://research.iac.es/proyecto/quijote/), 
or LAMBDA (https://lambda.gsfc.nasa.gov/product/quijote/). 

@author: jalberto
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

# Path to QUIJOTE RIMO files
path = '../rimo/'


# [A] Basic information. Values from Table 3 in Rubino-Martin et al. 2023.
ff      = path + 'rimo_quijote_mfi_params_dr1.fits'
hdulist = fits.open(ff)

cont = hdulist[1].data

freq = np.asarray(cont['FREQ_NOMINAL'])
fwhm = np.asarray(cont['FWHM_ARCMIN'])

print('Nominal central frequencies [GHz]',freq)
print('Beam FWHM [arcmin]',fwhm)



# [B] Bandpass
ff      = path + 'rimo_quijote_mfi_bandpass_dr1.fits'
hdulist = fits.open(ff)

cont = hdulist[1].data
freq_h3  = np.asarray(cont['freq_h3'])
bp_311 = np.asarray(cont['bp_311'])
bp_313 = np.asarray(cont['bp_313'])

freq_h2  = np.asarray(cont['freq_h2'])
bp_217 = np.asarray(cont['bp_217'])
bp_219 = np.asarray(cont['bp_219'])

freq_h4  = np.asarray(cont['freq_h4'])
bp_417 = np.asarray(cont['bp_417'])
bp_419 = np.asarray(cont['bp_419'])
hdulist.close()

# Plot
plt.plot(freq_h3[0,:],bp_311[0,:],linestyle='solid',label='311', color='black')
plt.plot(freq_h3[0,:],bp_313[0,:],linestyle='solid',label='313', color='red')
plt.plot(freq_h2[0,:],bp_217[0,:],linestyle='solid',label='217', color='cyan')
plt.plot(freq_h2[0,:],bp_219[0,:],linestyle='solid',label='219', color='orange')
plt.plot(freq_h4[0,:],bp_417[0,:],linestyle='solid',label='417', color='green')
plt.plot(freq_h4[0,:],bp_419[0,:],linestyle='solid',label='419', color='gray')

plt.xlim([9,21])
plt.ylim([0,4])
plt.xlabel(r'$\nu$ [GHz]', fontsize=10)
plt.ylabel('Transmission',fontsize=10)
plt.title('QUIJOTE MFI bandpasses')
plt.legend(loc='upper right', ncol=2, labelspacing=0.1)
#plt.savefig('rimo_quijote_mfi_bandpass_dr1.png',dpi=300, bbox_inches="tight")
plt.show()


# [C] Beam radial profiles
ff      = path + 'rimo_quijote_mfi_beamrp_dr1.fits'
hdulist = fits.open(ff)

cont = hdulist[1].data
theta  = np.asarray(cont['THETA'])
rp_311 = np.asarray(cont['rp_311'])
rp_313 = np.asarray(cont['rp_313'])
rp_217 = np.asarray(cont['rp_217'])
rp_219 = np.asarray(cont['rp_219'])
rp_417 = np.asarray(cont['rp_417'])
rp_419 = np.asarray(cont['rp_419'])
hdulist.close()

# Plot
plt.plot(theta[0,:],rp_311[0,:],linestyle='solid',label='311', color='black')
plt.plot(theta[0,:],rp_313[0,:],linestyle='dashed',label='313', color='red')
plt.plot(theta[0,:],rp_217[0,:],linestyle='solid',label='217', color='cyan')
plt.plot(theta[0,:],rp_219[0,:],linestyle='dashed',label='219', color='orange')
plt.plot(theta[0,:],rp_417[0,:],linestyle='solid',label='417', color='green')
plt.plot(theta[0,:],rp_419[0,:],linestyle='dashed',label='419', color='gray')

plt.yscale('log')
plt.xlim([0,7])
plt.ylim([1e-8,1])
plt.xlabel(r'$\theta$ [deg]', fontsize=10)
plt.ylabel(r'$r_p(\theta)$',fontsize=10)
plt.title('QUIJOTE MFI radial profiles')
plt.legend(loc='upper right', ncol=2, labelspacing=0.1)
#plt.savefig('rimo_quijote_mfi_beamrp_dr1.png',dpi=300, bbox_inches="tight")
plt.show()



# [D] Beam transfer functions
ff      = path + 'rimo_quijote_mfi_beamtf_dr1.fits'
hdulist = fits.open(ff)

cont = hdulist[1].data
ell  = np.asarray(cont['ell'])
bl_311 = np.asarray(cont['bl_311'])
bl_313 = np.asarray(cont['bl_313'])
bl_217 = np.asarray(cont['bl_217'])
bl_219 = np.asarray(cont['bl_219'])
bl_417 = np.asarray(cont['bl_417'])
bl_419 = np.asarray(cont['bl_419'])
hdulist.close()

# Plot
plt.plot(ell[0,:],bl_311[0,:]**2,linestyle='solid',label='311', color='black')
plt.plot(ell[0,:],bl_313[0,:]**2,linestyle='dashed',label='313', color='red')
plt.plot(ell[0,:],bl_217[0,:]**2,linestyle='solid',label='217', color='cyan')
plt.plot(ell[0,:],bl_219[0,:]**2,linestyle='dashed',label='219', color='orange')
plt.plot(ell[0,:],bl_417[0,:]**2,linestyle='solid',label='417', color='green')
plt.plot(ell[0,:],bl_419[0,:]**2,linestyle='dashed',label='419', color='gray')

plt.yscale('log')
plt.xlim([2,800])
plt.ylim([1e-12,1])
plt.xlabel(r'$\ell$', fontsize=10)
plt.ylabel(r'$W_\ell = B_\ell^2$',fontsize=10)
plt.title('QUIJOTE MFI window functions')
plt.legend(loc='lower left', ncol=2, labelspacing=0.1)
#plt.savefig('rimo_quijote_mfi_beamwf_dr1.png',dpi=300, bbox_inches="tight")
plt.show()

