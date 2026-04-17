'''
@author: jalberto

This program provides the unit conversions between usual quantities in CMB analyses:
       K_CMB (thermodynamic), K_RJ (rayleigh jeans), and Jy/sr (intensity).

Also provides the conversion between Planck brightness temperature (adopted by some codes at AM) and the 
usual RJ brightness temperarure (K_RJ) used in Radioastronomy.

HISTORY:
 25/Mar/2025 - original version, based on cmb_unit_conversion.pro
 04/Dec/2025 - added planck_uc_hfi(), which returns the HFI unit conversions.
 17/Apr/2026 - added routines to convert between different brightness temperature conventions (Planck brightness or RJ brightness)

'''

import numpy as np
from scipy.constants import c,h,k
Tcmb = 2.72548 # TCMB from Fixsen et al (2009).

def cmb_unit_conversion(nuGHz=None,option='KCMB2KRJ',help=False):

    casos = ['KCMB2KRJ', 'KRJ2KCMB', 'KCMB2Jysr', 'Jysr2KCMB', 'KRJ2Jysr', 'Jysr2KRJ']
    if nuGHz is None or help:
       print('  Syntax -- cmb_unit_conversion(nuGHz,option=)')
       print('  Possible options are',casos)
       print('  Default option is KCMB2KRJ.')
       return None

    # Basic computation
    nu  = nuGHz*1e9
    x   = h * nu/ (k*Tcmb)
    thermo = x**2 * np.exp(x)/(np.exp(x)-1.)**2
    rj     = ( 2.0 * k * nu**2 / c**2 ) * 1e26 

    # Identify case 
    if option == 'KCMB2KRJ':
       fac = thermo
    elif option == 'KRJ2KCMB':
       fac = 1/thermo
    elif option == 'KCMB2Jysr':
       fac = thermo * rj
    elif option == 'Jysr2KCMB':
       fac = 1 / (thermo*rj)
    elif option == 'KRJ2Jysr':
       fac = rj
    elif option == 'Jysr2KRJ':
       fac = 1/rj
    else:
        print("Units not identified. Returning -1")
        fac = -1

    return fac


# Unit corrections for PLANCK HFI. Values extracted from PLA, explanatory supplement.
#  https://wiki.cosmos.esa.int/planck-legacy-archive/index.php/UC_CC_Tables
# RGS values sent via email, 26-27/Feb/2025 
def planck_uc_hfi(use_bps=True):
   bands_hfi = np.array([100, 143, 217, 353, 545, 857], dtype=float) # GHz
  
   # Table 1. Coefficient MJy/sr/KCMB. Values correspond to "avg" entry.
   UC_HFI_KCMB2MJysr_PLA = np.array([244.0960, 371.7327, 483.6874, 287.4517, 58.0356, 2.2681])

   # Computed by RGS, including bandpass shift
   UC_HFI_KCMB2MJysr_rgs = np.array([242.09786, 370.53512, 481.93046, 287.22432, 56.659334, 2.1156277])

   # Table 2. KRJ/(MJy/sr).
   # Note: Coincides with Table 5 of Planck 2013, IX. HFI spectral response
   UC_HFI_MJysr2KRJ_PLA = np.array([0.0032548074, 0.0015916707, 0.00069120334, 0.00026120163, 0.00010958025, 4.4316316e-05 ])

   # Computed by RGS, evaluated at the center
   UC_HFI_MJysr2KRJ_rgs = 1./np.array([307.09143, 627.97125, 1446.0629, 3826.6356, 9121.3834, 22554.299])

   # Derived quantities (final outputs)
   UC_HFI_KCMB2KRJ = UC_HFI_KCMB2MJysr_rgs * UC_HFI_MJysr2KRJ_rgs # includes Bandpass shifts
   uc_hfi_no_bps   = UC_HFI_KCMB2MJysr_PLA * UC_HFI_MJysr2KRJ_PLA

   # Select output. Default is bandpass shift corrected value.
   output = UC_HFI_KCMB2KRJ 
   if use_bps==False: output = uc_hfi_no_bps 

   return output

# Routine to convert brightness temperature defintions, from Planck brightness temperature (as in Rybicki & Lightman)
# corresponding to I_nu = B_nu(T)), to the usual RJ brightness temperature used in Radioastronomy
def brightness_temperature_conversion(nuGHz=None, T=None, option='Planck2RJ',help=False): 
   casos = ['Planck2RJ', 'RJ2Planck']

   if nuGHz is None or T is None or help:
      print("""
Usage:
    brightness_temperature_conversion(nuGHz, T, option='Planck2RJ')

Parameters:
    nuGHz : frequency in GHz
    T     : temperature in K (either Planck brightness or K_RJ)
    option: 'Planck2RJ' o 'RJ2Planck'

Example:
    brightness_temperature_conversion(30.0, 2.7)
            """)
      return None

   nu  = nuGHz*1e9 # in Hz, SI units
   x   = h * nu / (k*T)

   if option == 'Planck2RJ':
      fac = x /(np.exp(x)-1.) * T
   elif option == 'RJ2Planck':
      fac = (h*nu/k) / np.log( x +1 ) 
   else:
      print("Temperature brightness units not identified. Returning -1")
      fac = -1

   return fac