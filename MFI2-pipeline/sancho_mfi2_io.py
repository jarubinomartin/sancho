#!/usr/bin/env python

"""
5/Jun/2023   

@authors: jalberto

Basic I/O routines for the data processing of MFI2 data.
* read_mfi2_tod3 
* read_mfi2_btod3

* write_mfi2_btod2


HISTORY:
* 05/06/2023 - original version, based on sancho_tfgi_io. JARM

"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys


# READ MFI2 data in TOD3 format, and returns a structure (dictionary). Example:
#** Structure <aa404638>, 4 tags, length=103680000, data length=103680000, refs=1:
#   JD              DOUBLE    Array[720000]
#   AZ              FLOAT     Array[720000]
#   EL              FLOAT     Array[720000]
#   DATA            FLOAT     Array[32, 720000]
#
def read_mfi2_tod3(filename):
    print(' READ_MFI2_TOD3: reading '+filename)
    hdulist = fits.open(filename)
    cont    = hdulist[1].data

    jd      = cont['JD'][0,:]
    az      = cont['AZ'][0,:]
    el      = cont['EL'][0,:]
    data    = cont['DATA'][0,:,:]

    # Output dictionary
    tod  = {'JD':jd, 'AZ':az, 'EL':el,'DATA':data}
    hdulist.close()
    
    return(tod)
    
# READ MFI2 data in BTOD3 format, and returns a structure (dictionary). Example:
#** Structure <af815008>, 15 tags, length=237816032, data length=237816022, refs=1:
#   NHORNS          INT              4
#   JD              DOUBLE    Array[495450]
#   AZ              FLOAT     Array[495450]
#   EL              FLOAT     Array[495450]
#   GL              FLOAT     Array[4, 495450]
#   GB              FLOAT     Array[4, 495450]
#   PAR             FLOAT     Array[4, 495450]
#   DATA            FLOAT     Array[32, 495450]
#   WEI             FLOAT     Array[32, 495450]
#   FLAG            INT       Array[32, 495450]
#   MSBIN           FLOAT           20.0000
#   POINTINGMODEL   STRING    'Version1'
#   AZHORN          FLOAT     Array[4, 495450]
#   ELHORN          FLOAT     Array[4, 495450]
#   SIGMACOV        FLOAT     Array[4, 2, 2, 495450]
def read_mfi2_btod3(filename):
    print(' READ_MFI2_BTOD3: reading '+filename)
    hdulist = fits.open(filename)
    cont    = hdulist[1].data
    names   = np.asarray(cont.names)

    nhorns   = cont['NHORNS'][0]
    jd       = cont['JD'][0,:]
    az       = cont['AZ'][0,:]
    el       = cont['EL'][0,:]
    gl       = cont['GL'][0,:,:]
    gb       = cont['GB'][0,:,:]
    par      = cont['PAR'][0,:,:]

    data     = cont['DATA'][0,:,:]
    wei      = cont['WEI'][0,:,:]
    flag     = cont['FLAG'][0,:,:]

    msbin    = cont['MSBIN'][0]
    pmodel   = cont['POINTINGMODEL'][0]
    
    azhorn   = cont['AZHORN'][0,:,:]
    elhorn   = cont['ELHORN'][0,:,:]
    sigmacov = cont['SIGMACOV'][0,:,:,:,:]
    
    # Output dictionary
    btod  = {'NHORNS':nhorns, 'JD':jd, 'AZ':az, 'EL':el,'GL':gl, 'GB':gb, 'PAR':par,
                 'DATA':data, 'FLAG':flag,'WEI':wei, 'MSBIN':msbin, 'POINTINGMODEL':pmodel,
                 'AZHORN':azhorn, 'ELHORN':elhorn, 'SIGMACOV':sigmacov }
    hdulist.close()
    
    return(btod)


# Write BTOD3 file. 
#
def write_mfi2_btod3(btod, ffout, overwrite=False):
    print(' WRITE_MFI2_BTOD3: writing '+ffout)

    # Definition of columns in BTOD3 files
    col1 = fits.Column(name='NHORNS', format='I', array = [btod['NHORNS']])
    
    col2 = fits.Column(name='JD', format=str(len(btod['JD']))+'D', array = [btod['JD']] )
    col3 = fits.Column(name='AZ', format=str(len(btod['AZ']))+'E', array = [btod['AZ']] )
    col4 = fits.Column(name='EL', format=str(len(btod['EL']))+'E', array = [btod['EL']] )

    col5 = fits.Column(name='GL', format=str(btod['GL'].size)+'E', array = [btod['GL']], dim=str(btod['GL'].shape[::-1]) )
    col6 = fits.Column(name='GB', format=str(btod['GB'].size)+'E', array = [btod['GB']], dim=str(btod['GB'].shape[::-1]) )
    col7 = fits.Column(name='PAR', format=str(btod['PAR'].size)+'E', array = [btod['PAR']], dim=str(btod['PAR'].shape[::-1]) )

    col8 = fits.Column(name='DATA', format=str(btod['DATA'].size)+'E', array = [btod['DATA']], dim=str(btod['DATA'].shape[::-1]) )
    col9 = fits.Column(name='WEI', format=str(btod['WEI'].size)+'E', array = [btod['WEI']], dim=str(btod['WEI'].shape[::-1]) )

    col10 = fits.Column(name='FLAG', format=str(btod['FLAG'].size)+'I', array = [btod['FLAG']], dim=str(btod['FLAG'].shape[::-1]) )

    col11 = fits.Column(name='MSBIN', format='E', array = [btod['MSBIN']])
    col12 = fits.Column(name='POINTINGMODEL', format=str(len(btod['POINTINGMODEL']))+'A', array = [btod['POINTINGMODEL']])
    
    col13 = fits.Column(name='AZHORN', format=str(btod['AZHORN'].size)+'E', array = [btod['AZHORN']], dim=str(btod['AZHORN'].shape[::-1]) )
    col14 = fits.Column(name='ELHORN', format=str(btod['ELHORN'].size)+'E', array = [btod['ELHORN']], dim=str(btod['ELHORN'].shape[::-1]) )
    col15 = fits.Column(name='SIGMACOV', format=str(btod['SIGMACOV'].size)+'E', array = [btod['SIGMACOV']], dim=str(btod['SIGMACOV'].shape[::-1]) )


    # Bin table
    hdu = fits.BinTableHDU.from_columns([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15])
    
    # write file
    hdu.writeto(ffout,overwrite=overwrite)
    
    return
