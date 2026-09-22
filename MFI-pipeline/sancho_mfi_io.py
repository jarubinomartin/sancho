#!/usr/bin/env python

"""
14/Sep/2026

@authors: jarm (Fortran original), adapted to Python following
          sancho_mfi2_io.py and sancho_tfgi_io.py

Basic I/O routines for the data processing of QUIJOTE MFI data.
Mirrors the FITS layout of btod_quijote / tod_quijote / cal_quijote
defined in sancho_fitsio_new.f90 (subroutines read_mfi_tod,
read_mfi_btod, read_mfi_cal, write_mfi_tod) from the PICASSO code.

* read_mfi_tod
* read_mfi_btod
* read_mfi_cal

* write_mfi_tod

Note on CTOD: in the current pipeline .ctod files share exactly the
same FITS layout as .btod files (gains already applied upstream) and
are read with read_mfi_btod -- there is no separate CTOD reader in
sancho_fitsio_new.f90, so none is added here either.

HISTORY:
* 14/09/2026 - original version. Mirrors sancho_fitsio_new.f90 (MFI part).

"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys

# Same parameters as in sancho_fitsio_new.f90
NCHAN_MFI = 32
NHORN_MFI = 4
JD_REF    = 2456244.5


# READ MFI data in TOD format, and returns a structure (dictionary).
# Mirrors read_mfi_tod() / tod_quijote in sancho_fitsio_new.f90:
#   JD          DOUBLE   Array[ndata]   (JD_REF is added, as in Fortran)
#   AZ          FLOAT    Array[ndata]
#   EL          FLOAT    Array[ndata]
#   DATA        FLOAT    Array[32, ndata]
#   MOD_ANGLE   FLOAT    Array[4, ndata]
#   CAL         LONG     Array[ndata]
def read_mfi_tod(filename):
    print(' READ_MFI_TOD: reading '+filename)
    hdulist = fits.open(filename)
    cont    = hdulist[1].data

    jd       = cont['JD'][0, :] + JD_REF   # Fortran adds jd_ref for TOD (not for BTOD)
    az       = cont['AZ'][0, :]
    el       = cont['EL'][0, :]
    data     = cont['DATA'][0, :, :]
    mod_ang  = cont['MOD_ANGLE'][0, :, :]
    cal      = cont['CAL'][0, :]

    # Output dictionary
    tod = {'JD': jd, 'AZ': az, 'EL': el, 'DATA': data,
           'MOD_ANGLE': mod_ang, 'CAL': cal}
    hdulist.close()

    return (tod)


# READ MFI data in BTOD format, and returns a structure (dictionary).
# Mirrors read_mfi_btod() / btod_quijote in sancho_fitsio_new.f90:
#   JD          DOUBLE   Array[ndata]
#   AZ          FLOAT    Array[ndata]
#   EL          FLOAT    Array[ndata]
#   GL          FLOAT    Array[4, ndata]
#   GB          FLOAT    Array[4, ndata]
#   PAR         FLOAT    Array[4, ndata]
#   DATA        FLOAT    Array[32, ndata]
#   MOD_ANG     FLOAT    Array[4, ndata]
#   WEI         FLOAT    Array[32, ndata]
#   FLAG        INT      Array[32, ndata]
#   MSBIN       FLOAT    scalar
#   AZHORN      FLOAT    Array[4, ndata]
#   ELHORN      FLOAT    Array[4, ndata]
#   SIGMACOV    FLOAT    Array[4, 2, 2, ndata]
def read_mfi_btod(filename):
    print(' READ_MFI_BTOD: reading '+filename)
    hdulist = fits.open(filename)
    cont    = hdulist[1].data

    jd       = cont['JD'][0, :]
    az       = cont['AZ'][0, :]
    el       = cont['EL'][0, :]
    gl       = cont['GL'][0, :, :]
    gb       = cont['GB'][0, :, :]
    par      = cont['PAR'][0, :, :]

    data     = cont['DATA'][0, :, :]
    mod_ang  = cont['MOD_ANG'][0, :, :]
    wei      = cont['WEI'][0, :, :]
    flag     = cont['FLAG'][0, :, :]

    msbin    = cont['MSBIN'][0]

    azhorn   = cont['AZHORN'][0, :, :]
    elhorn   = cont['ELHORN'][0, :, :]
    sigmacov = cont['SIGMACOV'][0, :, :, :, :]

    # Output dictionary
    btod = {'JD': jd, 'AZ': az, 'EL': el, 'GL': gl, 'GB': gb, 'PAR': par,
            'DATA': data, 'MOD_ANG': mod_ang, 'WEI': wei, 'FLAG': flag,
            'MSBIN': msbin, 'AZHORN': azhorn, 'ELHORN': elhorn,
            'SIGMACOV': sigmacov}
    hdulist.close()

    return (btod)


# READ MFI calibration data (CAL format), and returns a structure (dictionary).
# Mirrors read_mfi_cal() / cal_quijote in sancho_fitsio_new.f90:
#   JD        DOUBLE   Array[ndata]
#   GAIN      FLOAT    Array[32, ndata]
#   RAIN      FLOAT    Array[32, ndata]
#   BASE      FLOAT    Array[32, ndata]
#   ACTI      FLOAT    Array[32, ndata]
#   FLAG      INT      Array[32, ndata]
#   MOD_ANG   FLOAT    Array[4, ndata]
#   RMS       FLOAT    Array[32, ndata]
def read_mfi_cal(filename):
    print(' READ_MFI_CAL: reading '+filename)
    hdulist = fits.open(filename)
    cont    = hdulist[1].data

    jd      = cont['JD'][0, :]
    gain    = cont['GAIN'][0, :, :]
    rain    = cont['RAIN'][0, :, :]
    base    = cont['BASE'][0, :, :]
    acti    = cont['ACTI'][0, :, :]
    flag    = cont['FLAG'][0, :, :]
    mod_ang = cont['MOD_ANG'][0, :, :]
    rms     = cont['RMS'][0, :, :]

    # Output dictionary
    cal = {'JD': jd, 'GAIN': gain, 'RAIN': rain, 'BASE': base, 'ACTI': acti,
           'FLAG': flag, 'MOD_ANG': mod_ang, 'RMS': rms}
    hdulist.close()

    return (cal)


# Write TOD file. Mirrors write_mfi_tod() in sancho_fitsio_new.f90
# (the only MFI writer implemented in Fortran):
# columns JD, AZ, EL, MOD_ANGLE, CAL, DATA.
def write_mfi_tod(tod, ffout, overwrite=False):
    print(' WRITE_MFI_TOD: writing '+ffout)

    col1 = fits.Column(name='JD', format=str(len(tod['JD']))+'D', array=[tod['JD']])
    col2 = fits.Column(name='AZ', format=str(len(tod['AZ']))+'E', array=[tod['AZ']])
    col3 = fits.Column(name='EL', format=str(len(tod['EL']))+'E', array=[tod['EL']])

    col4 = fits.Column(name='MOD_ANGLE', format=str(tod['MOD_ANGLE'].size)+'E',
                        array=[tod['MOD_ANGLE']], dim=str(tod['MOD_ANGLE'].shape[::-1]))
    col5 = fits.Column(name='CAL', format=str(len(tod['CAL']))+'I', array=[tod['CAL']])
    col6 = fits.Column(name='DATA', format=str(tod['DATA'].size)+'E',
                        array=[tod['DATA']], dim=str(tod['DATA'].shape[::-1]))

    # Bin table
    hdu = fits.BinTableHDU.from_columns([col1, col2, col3, col4, col5, col6])

    # write file
    hdu.writeto(ffout, overwrite=overwrite)

    return
