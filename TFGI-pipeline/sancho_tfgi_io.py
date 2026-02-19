#!/usr/bin/env python

"""
10/Jul/2023   

@authors: jalberto, mpeel, rjh

Basic I/O routines for the data processing of TFGI data.
* read_tfgi_pixel_masterfile
* read_tfgi_tod2 
* read_tfgi_btod2
* read_tfgi_ctod2
* read_tfgi_sci2

* write_tfgi_btod2
* write_tfgi_ctod2


HISTORY:
* 10/07/2023 - original version. JARM
* 25/10/2023 - updated read_tfgi_btod2, to include WEI_IQU. Added read_tfgi_ctod2
* 30/01/2026 - adding write_tfgi_ctod2

"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys


# Based on read_sci_FTGI_2018 from RJH.
def read_tfgi_sci2(filename):
    print(' Reading in SCI2 mode: ',filename)

    data=np.empty(shape=(124,60*4000),dtype=float)
    dat=np.empty(4000,dtype='i8')
    with open(filename, "rb") as f:
        for k in range(60):
            print(k)
            for i in range(124):
                campo1= np.fromfile(f,count=1, dtype='>a2')
                campo2= np.fromfile(f,count=1, dtype='>u2')            
                campo3= np.fromfile(f,count=1, dtype='>u2')
                campo4= np.fromfile(f,count=1, dtype='>u1')
                campo5= np.fromfile(f,count=1, dtype='>u1')
                campo6= np.fromfile(f,count=1, dtype='>u1')
                campo7= np.fromfile(f,count=1, dtype='>u1')
                Tstamp= np.fromfile(f,count=1, dtype='>u4')
                Time=   np.fromfile(f,count=1, dtype='>u4')
                SID=    np.fromfile(f,count=1, dtype='>u2')
                AqErr=  np.fromfile(f,count=1, dtype='>u1')
                channel_ID=  np.fromfile(f,count=1, dtype='>u1')
                cal_flag=  np.fromfile(f,count=1, dtype='>u1')
                cal_sw=  np.fromfile(f,count=1, dtype='>u1')
                ph_sw=  np.fromfile(f,count=1, dtype='>u1') 
                repuesto1=  np.fromfile(f,count=1, dtype='>u1') 
                Nphstates=  np.fromfile(f,count=1, dtype='>u4')
                PHseq= np.fromfile(f,count=16, dtype='>u1')
                Processtype= np.fromfile(f,count=1, dtype='>u1')
                repuesto2= np.fromfile(f,count=1, dtype='>u1')
                Pack_count=  np.fromfile(f,count=1, dtype='>u2')
                Sc_factor=np.fromfile(f,count=1, dtype='>f8')
                Samprate=  np.fromfile(f,count=1, dtype='>u4')
                NSample=  np.fromfile(f,count=1, dtype='u4')
                dat= np.fromfile(f,count=4000, dtype='>i4')
                dat=dat*Sc_factor
                data[i,(k*4000):(k*4000)+4000]=dat
           
                print(k,i)      
                print('campo1 = ' + str(campo1))
                print('campo2 = ' + str(campo2))
                print('campo3 = ' + str(campo3))
                print('campo4 = ' + str(campo4))
                print('campo5 = ' + str(campo5))
                print('campo6 = ' + str(campo6))
                print('campo7 = ' + str(campo7))
                print('Tstamp = ' + str(Tstamp))
                print('Time = ' + str(Time))
                print('SID = ' + str(SID))
                print('AqErr = ' + str(AqErr)) 
                print('channel_ID = ' + str(channel_ID))
                print('cal_flag = ' + str(cal_flag))
                print('cal__sw = ' + str(cal_sw))
                print('ph__sw = ' + str(ph_sw) + ' (1- 16KHz, 2- 8KHz)')
                print('repuesto1 = ' + str(repuesto1))
                print('Nphstates = ' + str(Nphstates))
                print('Phase sequence = ' + str(PHseq))
                print('Process type = ' + str(Processtype))
                print('repuesto2 = ' + str(repuesto2))
                print('packet counter = ' + str(Pack_count))
                print('Scale factor = ' + str(Sc_factor))
                print('Sampling rate = ' + str(Samprate))
                print('No of samples = ' + str(NSample))
                print(dat[0:10])

        return data


# Reads TFGI pixel masterfile in ascii format. Original version from M. Peel.
def read_tfgi_pixel_masterfile(filename,usefloat=False):
        jd    = 0
        array = []
        with open(filename) as f:
                for line in f:
                        if jd == 0:
                                jd = float(line.strip())
                        else:
                                val = line.strip().split()
                                # Convert to integers
                                if usefloat:
                                        val = [float(x) for x in val]
                                else:
                                        val = [int(x) for x in val]
                                array.append(val)
        return jd, array


# READ TFGI data in TOD2 format, and returns a structure (dictionary). Example:
#** Structure <ca809208>, 7 tags, length=31680224, data length=31680224, refs=1:
#   DAS             LONG      Array[28]
#   OUTPUT          LONG      Array[28]
#   JD              DOUBLE    Array[240000]
#   AZ              FLOAT     Array[240000]
#   EL              FLOAT     Array[240000]
#   DATA            FLOAT     Array[28, 240000]
#   PS              LONG      Array[240000]
def read_tfgi_tod2(filename):
    hdulist = fits.open(filename)
    cont    = hdulist[1].data

    das     = cont['DAS'][0,:]
    output  = cont['OUTPUT'][0,:] 
    
    jd      = cont['JD'][0,:]
    az      = cont['AZ'][0,:]
    el      = cont['EL'][0,:]

    data    = cont['DATA'][0,:,:]
    ps      = cont['PS'][0,:]

    # Output dictionary
    tod  = {'DAS':das, 'OUTPUT':output, 'JD':jd, 'AZ':az, 'EL':el,'DATA':data, 'PS':ps}
    hdulist.close()
    
    return(tod)
    
# READ TFGI data in BTOD2 format, and returns a structure (dictionary). Example: 
#** Structure <c6851208>, 17 tags, length=158077368, data length=158077366, refs=1:
#   NHORNS          INT              7
#   LISTCHAN        INT       Array[28]
#   LISTPIX         INT       Array[7]
#   LISTDAS         INT       Array[7]
#   JD              DOUBLE    Array[123885]
#   AZ              FLOAT     Array[123885]
#   EL              FLOAT     Array[123885]
#   GL              FLOAT     Array[7, 123885]
#   GB              FLOAT     Array[7, 123885]
#   PAR             FLOAT     Array[7, 123885]
#   DATA            FLOAT     Array[28, 4, 123885]
#   WEI             FLOAT     Array[28, 4, 123885]
#   FLAG            INT       Array[28, 4, 123885]
#   AZHORN          FLOAT     Array[7, 123885]
#   ELHORN          FLOAT     Array[7, 123885]
#   MSBIN           FLOAT           4.00000
#   POINTINGMODEL   STRING    'MFT-Feb2022'
#   WEI_IQU         FLOAT     Array[28, 3, 123885]  # This field is only present in binned data.
def read_tfgi_btod2(filename):
    print(' READ_TFGI_BTOD2: reading '+filename)
    hdulist = fits.open(filename)
    cont    = hdulist[1].data
    names   = np.asarray(cont.names)

    nhorns   = cont['NHORNS'][0]
    listchan = cont['LISTCHAN'][0,:] 
    listpix  = cont['LISTPIX'][0,:] 
    listdas  = cont['LISTDAS'][0,:] 
    
    jd       = cont['JD'][0,:]
    az       = cont['AZ'][0,:]
    el       = cont['EL'][0,:]
    gl       = cont['GL'][0,:,:]
    gb       = cont['GB'][0,:,:]
    par      = cont['PAR'][0,:,:]

    data     = cont['DATA'][0,:,:,:]
    flag     = cont['FLAG'][0,:,:,:]
    wei      = cont['WEI'][0,:,:,:]

    azhorn   = cont['AZHORN'][0,:,:]
    elhorn   = cont['ELHORN'][0,:,:]

    msbin    = cont['MSBIN'][0]
    pmodel   = cont['POINTINGMODEL'][0]

    # WEI_IQU. This field only appears in binned BTOD files.
    di       = np.argwhere(names == 'WEI_IQU')
    if len(di)==1:
        wei_iqu  = cont['WEI_IQU'][0,:,:,:]
    else:
        nsamp    = len(jd)
        wei_iqu  = np.zeros((nsamp,3,nhorns*4))
        print(' READ_TFGI_BTOD2: warning -- WEI_IQU field is not present. Setting values to zero. ')

    # Output dictionary
    btod  = {'NHORNS':nhorns, 'LISTCHAN':listchan, 'LISTPIX':listpix, 'LISTDAS':listdas,
                 'JD':jd, 'AZ':az, 'EL':el,'GL':gl, 'GB':gb, 'PAR':par, 'DATA':data, 'FLAG':flag,
                 'WEI':wei, 'AZHORN':azhorn, 'ELHORN':elhorn, 'MSBIN':msbin, 'POINTINGMODEL':pmodel, 'WEI_IQU':wei_iqu }
    hdulist.close()
    
    return(btod)


###################
# READ CTOD2 FILE. Same format at BTOD2, but with few extra fields related to gains.
def read_tfgi_ctod2(filename):
    print(' READ_TFGI_CTOD2: reading '+filename)
    hdulist = fits.open(filename)
    cont    = hdulist[1].data

    nhorns   = cont['NHORNS'][0]
    listchan = cont['LISTCHAN'][0,:] 
    listpix  = cont['LISTPIX'][0,:] 
    listdas  = cont['LISTDAS'][0,:] 
    
    jd       = cont['JD'][0,:]
    az       = cont['AZ'][0,:]
    el       = cont['EL'][0,:]
    gl       = cont['GL'][0,:,:]
    gb       = cont['GB'][0,:,:]
    par      = cont['PAR'][0,:,:]

    data     = cont['DATA'][0,:,:,:]
    flag     = cont['FLAG'][0,:,:,:]
    wei      = cont['WEI'][0,:,:,:]

    azhorn   = cont['AZHORN'][0,:,:]
    elhorn   = cont['ELHORN'][0,:,:]

    msbin    = cont['MSBIN'][0]
    pmodel   = cont['POINTINGMODEL'][0]

    wei_iqu  = cont['WEI_IQU'][0,:,:,:]
    gainmodel = cont['GAINMODEL'][0]
    gains    = cont['GAINS'][0,:,:]

    # Output dictionary
    ctod  = {'NHORNS':nhorns, 'LISTCHAN':listchan, 'LISTPIX':listpix, 'LISTDAS':listdas,
                 'JD':jd, 'AZ':az, 'EL':el,'GL':gl, 'GB':gb, 'PAR':par, 'DATA':data, 'FLAG':flag,
                 'WEI':wei, 'AZHORN':azhorn, 'ELHORN':elhorn, 'MSBIN':msbin, 'POINTINGMODEL':pmodel,
                 'WEI_IQU':wei_iqu, 'GAINMODEL':gainmodel, 'GAINS':gains }
    hdulist.close()
    
    return(ctod)



# Write BTOD2 file. An example of the format in header is:
#TFORM1  = 'I       '           /                                                
#TFORM2  = '28I     '           /                                                
#TFORM3  = '7I      '           /                                                
#TFORM4  = '7I      '           /                                                
#TFORM5  = '123885D '           /                                                
#TFORM6  = '123885E '           /                                                
#TFORM7  = '123885E '           /                                                
#TFORM8  = '867195E '           /                                                
#TFORM9  = '867195E '           /                                                
#TFORM10 = '867195E '           /                                                
#TFORM11 = '13875120E'          /                                                
#TFORM12 = '13875120E'          /                                                
#TFORM13 = '13875120I'          /                                                
#TFORM14 = '867195E '           /                                                
#TFORM15 = '867195E '           /                                                
#TFORM16 = 'E       '           /                                                
#TFORM17 = '11A     '           /                                                
#COMMENT                                                                         
#COMMENT  *** Column dimensions (2 D or greater) ***                             
#COMMENT                                                                         
#TDIM8   = '( 7, 123885)'       /                                                
#TDIM9   = '( 7, 123885)'       /                                                
#TDIM10  = '( 7, 123885)'       /                                                
#TDIM11  = '( 28, 4, 123885)'   /                                                
#TDIM12  = '( 28, 4, 123885)'   /                                                
#TDIM13  = '( 28, 4, 123885)'   /                                                
#TDIM14  = '( 7, 123885)'       /                                                
#TDIM15  = '( 7, 123885)'       /                                                
#COMMENT                                                                         
#COMMENT  *** Column names ***                                                   
#COMMENT                                                                         
#TTYPE1  = 'NHORNS  '           /                                                
#TTYPE2  = 'LISTCHAN'           /                                                
#TTYPE3  = 'LISTPIX '           /                                                
#TTYPE4  = 'LISTDAS '           /                                                
#TTYPE5  = 'JD      '           /                                                
#TTYPE6  = 'AZ      '           /                                                
#TTYPE7  = 'EL      '           /                                                
#TTYPE8  = 'GL      '           /                                                
#TTYPE9  = 'GB      '           /                                                
#TTYPE10 = 'PAR     '           /                                                
#TTYPE11 = 'DATA    '           /                                                
#TTYPE12 = 'WEI     '           /                                                
#TTYPE13 = 'FLAG    '           /                                                
#TTYPE14 = 'AZHORN  '           /                                                
#TTYPE15 = 'ELHORN  '           /                                                
#TTYPE16 = 'MSBIN   '           /                                                
#TTYPE17 = 'POINTINGMODEL'      / 
#
def write_tfgi_btod2(btod, ffout, overwrite=False):

    # Definition of columns in BTOD2 files
    col1 = fits.Column(name='NHORNS', format='I', array = [btod['NHORNS']])
    
    col2 = fits.Column(name='LISTCHAN', format=str(len(btod['LISTCHAN']))+'I', array = [btod['LISTCHAN']] )
    col3 = fits.Column(name='LISTPIX', format=str(len(btod['LISTPIX']))+'I', array = [btod['LISTPIX']] )
    col4 = fits.Column(name='LISTDAS', format=str(len(btod['LISTDAS']))+'I', array = [btod['LISTDAS']] )

    col5 = fits.Column(name='JD', format=str(len(btod['JD']))+'D', array = [btod['JD']] )
    col6 = fits.Column(name='AZ', format=str(len(btod['AZ']))+'E', array = [btod['AZ']] )
    col7 = fits.Column(name='EL', format=str(len(btod['EL']))+'E', array = [btod['EL']] )

    col8 = fits.Column(name='GL', format=str(btod['GL'].size)+'E', array = [btod['GL']], dim=str(btod['GL'].shape[::-1]) )
    col9 = fits.Column(name='GB', format=str(btod['GB'].size)+'E', array = [btod['GB']], dim=str(btod['GB'].shape[::-1]) )
    col10 = fits.Column(name='PAR', format=str(btod['PAR'].size)+'E', array = [btod['PAR']], dim=str(btod['PAR'].shape[::-1]) )

    col11 = fits.Column(name='DATA', format=str(btod['DATA'].size)+'E', array = [btod['DATA']], dim=str(btod['DATA'].shape[::-1]) )
    col12 = fits.Column(name='WEI', format=str(btod['WEI'].size)+'E', array = [btod['WEI']], dim=str(btod['WEI'].shape[::-1]) )

    col13 = fits.Column(name='FLAG', format=str(btod['FLAG'].size)+'I', array = [btod['FLAG']], dim=str(btod['FLAG'].shape[::-1]) )

    col14 = fits.Column(name='AZHORN', format=str(btod['AZHORN'].size)+'E', array = [btod['AZHORN']], dim=str(btod['AZHORN'].shape[::-1]) )
    col15 = fits.Column(name='ELHORN', format=str(btod['ELHORN'].size)+'E', array = [btod['ELHORN']], dim=str(btod['ELHORN'].shape[::-1]) )

    col16 = fits.Column(name='MSBIN', format='E', array = [btod['MSBIN']])
    col17 = fits.Column(name='POINTINGMODEL', format=str(len(btod['POINTINGMODEL']))+'A', array = [btod['POINTINGMODEL']])

    # Bin table
    hdu = fits.BinTableHDU.from_columns([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15,col16,col17])
    
    # write file
    hdu.writeto(ffout,overwrite=overwrite)
    
    return

# Write CTOD file
def write_tfgi_ctod2(ctod, ffout, overwrite=False):

    print(' WRITE_TFGI_CTOD2: writing '+ffout)
    # Definition of columns in CTOD2 files
    col1 = fits.Column(name='NHORNS', format='I', array = [ctod['NHORNS']])

    col2 = fits.Column(name='LISTCHAN', format=str(len(ctod['LISTCHAN']))+'I', array = [ctod['LISTCHAN']] )
    col3 = fits.Column(name='LISTPIX', format=str(len(ctod['LISTPIX']))+'I', array = [ctod['LISTPIX']] )
    col4 = fits.Column(name='LISTDAS', format=str(len(ctod['LISTDAS']))+'I', array = [ctod['LISTDAS']] )

    col5 = fits.Column(name='JD', format=str(len(ctod['JD']))+'D', array = [ctod['JD']] )
    col6 = fits.Column(name='AZ', format=str(len(ctod['AZ']))+'E', array = [ctod['AZ']] )
    col7 = fits.Column(name='EL', format=str(len(ctod['EL']))+'E', array = [ctod['EL']] )

    col8 = fits.Column(name='GL', format=str(ctod['GL'].size)+'E', array = [ctod['GL']], dim=str(ctod['GL'].shape[::-1]) )
    col9 = fits.Column(name='GB', format=str(ctod['GB'].size)+'E', array = [ctod['GB']], dim=str(ctod['GB'].shape[::-1]) )
    col10 = fits.Column(name='PAR', format=str(ctod['PAR'].size)+'E', array = [ctod['PAR']], dim=str(ctod['PAR'].shape[::-1]) )

    col11 = fits.Column(name='DATA', format=str(ctod['DATA'].size)+'E', array = [ctod['DATA']], dim=str(ctod['DATA'].shape[::-1]) )
    col12 = fits.Column(name='WEI', format=str(ctod['WEI'].size)+'E', array = [ctod['WEI']], dim=str(ctod['WEI'].shape[::-1]) )

    col13 = fits.Column(name='FLAG', format=str(ctod['FLAG'].size)+'I', array = [ctod['FLAG']], dim=str(ctod['FLAG'].shape[::-1]) )

    col14 = fits.Column(name='AZHORN', format=str(ctod['AZHORN'].size)+'E', array = [ctod['AZHORN']], dim=str(ctod['AZHORN'].shape[::-1]) )
    col15 = fits.Column(name='ELHORN', format=str(ctod['ELHORN'].size)+'E', array = [ctod['ELHORN']], dim=str(ctod['ELHORN'].shape[::-1]) )

    col16 = fits.Column(name='MSBIN', format='E', array = [ctod['MSBIN']])
    col17 = fits.Column(name='POINTINGMODEL', format=str(len(ctod['POINTINGMODEL']))+'A', array = [ctod['POINTINGMODEL']])

    # === CTOD-only columns ===
    col18 = fits.Column( name='WEI_IQU', format=str(ctod['WEI_IQU'].size)+'E', array=[ctod['WEI_IQU']], dim=str(ctod['WEI_IQU'].shape[::-1]))
    col19 = fits.Column( name='GAINMODEL', format=str(len(ctod['GAINMODEL']))+'A', array=[ctod['GAINMODEL']])
    col20 = fits.Column( name='GAINS', format=str(ctod['GAINS'].size)+'E', array=[ctod['GAINS']], dim=str(ctod['GAINS'].shape[::-1]))

    # Build binary table
    hdu = fits.BinTableHDU.from_columns([ col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, 
                                         col11, col12, col13, col14, col15, col16, col17, col18, col19, col20 ])

    # Write file
    hdu.writeto(ffout, overwrite=overwrite)

    return
