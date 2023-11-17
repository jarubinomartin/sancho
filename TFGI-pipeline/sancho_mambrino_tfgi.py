#!/usr/bin/env python

"""
26/Oct/2023   

@author: jalberto

Implementation of the naive map-making for TFGI btod2 or ctod2 lists of files

HISTORY:
* 26/10/2023 - original version. JARM.
*  3/11/2023 - added horns as input option. Limits the computation of the maps to a sublist.

"""

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from astropy.io import fits
from sancho_tfgi_io import read_tfgi_btod2
import sys


def mambrino_tfgi(root, path='/net/calp-nas/proyectos/quijote2/ctod/',tail='.ctod2', nside=512, nest=False, nhorns=7, 
                      horns=np.arange(7), dobaserm=True, t_base=6.0, usewei=True):
    '''
    Keywords:
       * dobaserm : if True, removes baseline computing the median of the data in scans of length t_base.
       * t_base   : baseline size in seconds.
       * usewei   : if True, uses the WEI_IQU information.
    '''

    nf   = len(root)

    # Display basic info
    print('*** MAMBRINO_TFGI code ***')
    print('')
    print(' Settings:')
    print('    Path                            = ',path)
    print('    Output NSIDE                    = ',nside)
    print('    Remove baseline (median filter) = ',dobaserm)
    if dobaserm: print('    Baseline length (secs)          = ',t_base)
    print('    Use WEI_IQU                     = ',usewei)
    print('    Making maps for horns           = ',horns)
    print('')

    # Main loop
    npix   = hp.nside2npix(nside)
    mapa   = np.zeros((npix,nhorns*4))
    deno   = np.zeros((npix,nhorns*4))
    nhits  = np.zeros((npix,nhorns*4))
    for i in np.arange(nf):
        ff   = path + root[i] + tail
        btod = read_tfgi_btod2(ff)

        data  = btod['DATA']
        flag  = btod['FLAG']
        wei_iqu=btod['WEI_IQU']
        gl    = btod['GL']
        gb    = btod['GB']
        msbin = btod['MSBIN'] 

        nsamp = len(gl[:,0])
        if nsamp != np.max(gl.shape):
            sys.exit('SYS.EXIT() -- incorrect nsamp. ')
        navg  = np.int32(t_base*1000/msbin)
        nbase = np.max([1,np.int32(nsamp/navg)])
    
        for ihorn in horns:
            theta = np.deg2rad( 90.0 - gb[:,ihorn] )
            phi   = np.deg2rad( gl[:,ihorn] )
            ipix  = hp.ang2pix(nside, theta, phi, nest=nest)

            for ichan in np.arange(4):
                k = ihorn*4 + ichan
                
                # Preparing Intensity maps only
                Vd    = ( data[:,0,k] + data[:,1,k] + data[:,2,k] + data[:,3,k] ) / 2.
                flagd = flag[:,0,k] + flag[:,1,k] + flag[:,2,k] + flag[:,3,k]
                weid  = wei_iqu[:,0,k]  # Intensity wei.

                # Use weights? 
                if usewei==False:
                    weid = np.ones_like(wei_iqu[:,0,k])
            
                # Removing median (if requested).
                if dobaserm==True:
                    Vd_smth = np.zeros_like(Vd)
                    for j in np.arange(nbase):
                        i1 = np.int64(j * navg)
                        i2 = i1 + navg #- 1 in python

                        datos = Vd[i1:i2]
                        flags = flagd[i1:i2]
                        usa   = np.where(flags == 0)
                        if len(usa[-1] > 2):
                            mediana = np.median(datos[usa])
                            Vd_smth[i1:i2] = datos[:] - mediana
                            
                else:
                    Vd_smth = Vd
            
                # Project onto a Healpy map. Exclude flags!=0 and Vd=0
                cut = np.where( (flagd == 0) & (Vd_smth != 0 ) )
                if (len(cut[-1]) != 0):
                    np.add.at(mapa[:,k], ipix[cut], Vd_smth[cut]*weid[cut])
                    np.add.at(deno[:,k], ipix[cut], weid[cut])
                    np.add.at(nhits[:,k], ipix[cut], np.ones(len(cut)))


    # Normalising and UNSEEN pixels
    print(' Normalising map')
    cut       = np.where(deno != 0)
    mapa[cut] = mapa[cut]/deno[cut]

    print(' Setting UNSEEN pixels')
    ceros        = np.where(deno == 0)
    mapa[ceros]  = hp.UNSEEN
    nhits[ceros] = hp.UNSEEN
    
    print('*** MAMBRINO_TFGI code ends ***')

    return mapa, nhits


# Write maps to FITS file
def write_mambrino_tfgi_maps(mapa,nhits,ffout='test_mambrino.fits'):
    col1 = fits.Column(name='map', format=str(mapa.size)+'E', array = [mapa], dim=str(mapa.shape[::-1]) )
    col2 = fits.Column(name='nhits', format=str(nhits.size)+'E', array = [nhits], dim=str(nhits.shape[::-1]) )
    hdu = fits.BinTableHDU.from_columns([col1,col2])
    hdu.writeto(ffout,overwrite=True)
    return


##############
# MAIN code
if __name__ == "__main__":

    # Example path and list of files
    path  = '/Users/jalberto/quijote/data/' 
    root  =  np.asarray(['CRAB-211130-0230-000','CRAB-211130-0230-001','CRAB-211130-0230-002','CRAB-211130-0230-003'])
    tail  = '.btod2'
    nside = 512
    nest  = False

    #mapa, nhits = mambrino_tfgi(root, path=path,tail=tail, nside=nside, nest=nest, nhorns=7, horns=np.asarray([2]), usewei=False)
    mapa, nhits = mambrino_tfgi(root, path=path,tail=tail, nside=nside, nest=nest, nhorns=7, usewei=False)

    # Display
    hp.gnomview(mapa[:,0],rot=[185,-5],norm='hist',reso=3,xsize=500)
    hp.gnomview(nhits[:,0],rot=[185,-5],norm='hist',reso=3,xsize=500)
    plt.show()

    hp.write_map('testmam.fits',mapa[:,0],nest=nest,coord='G',overwrite=True)

    
    # Write all maps to FITS file
    write_mambrino_tfgi_maps(mapa,nhits,ffout='test_mambrino.fits')

