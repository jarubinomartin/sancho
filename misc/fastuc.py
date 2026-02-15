'''
@author: jalberto

This program provides the unit conversions to be used with FASTCC.

HISTORY:
  9/Mar/2025 - original version.
 15/Feb/2026 - now using cmb_unit_conversion.py as external code. For the moment, just 
                using QUIJOTE-MFI, WMAP and PLANCK.

'''

import numpy as np
from cmb_unit_conversion import *

# Main code
def fastuc(freq, returnfreq=False):
    # Define dictionaries containing the coefficients for different detectors and frequencies.
    # Using same notation as in fastcc.
    nu = {
        'Q11': 11.1,
        'Q13': 12.9, 
        'Q17': 16.8, 
        'Q19': 18.8,
        'WK':  22.8,
        'WKa': 33.0,
        'P30': 28.4,
        'P44': 44.1, 
        'P70': 70.4,
        'P100': 100.0,
        'P143': 143.0,
        'P217': 217.0,
        'P353': 353.0,
        'P545': 545.0,
        'P857': 857.0
    }

    hfi = planck_uc_hfi(use_bps=True)

    uc_v1 = {
    # QUIJOTE. CC defined for central frequency
        'Q11': cmb_unit_conversion(nu['Q11'], option='KCMB2KRJ'),
        'Q13': cmb_unit_conversion(nu['Q13'], option='KCMB2KRJ'),
        'Q17': cmb_unit_conversion(nu['Q17'], option='KCMB2KRJ'),
        'Q19': cmb_unit_conversion(nu['Q19'], option='KCMB2KRJ'),
    # WMAP9. CC defined for central frequency.
        'WK': cmb_unit_conversion(nu['WK'], option='KCMB2KRJ'),
        'WKa': cmb_unit_conversion(nu['WKa'], option='KCMB2KRJ'),
    # LFI: CC defined for central frequency.
        'P30': cmb_unit_conversion(nu['P30'], option='KCMB2KRJ'),
        'P44': cmb_unit_conversion(nu['P44'], option='KCMB2KRJ'),
        'P70': cmb_unit_conversion(nu['P70'], option='KCMB2KRJ'),
    # HFI: using HFI convention. Values include BPS. Used for PR3 and NPIPE
        'P100': hfi[0],
        'P143': hfi[1],
        'P217': hfi[2],
        'P353': hfi[3],    
        'P545': hfi[4],
        'P857': hfi[5]
    }

    uc = uc_v1.get(freq, 0)

    if uc==0:
        print(' FASTUC error')
        print(f" > freq={freq} is not defined. Returning uc=0. ")
        print(f" > Available keys: {list(nu.keys())}")

    if returnfreq:
        return [nu[freq],uc]
    else:
        return uc