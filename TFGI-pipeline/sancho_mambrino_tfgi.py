#!/usr/bin/env python

"""
26/Oct/2023   

@author: jalberto

Implementation of the naive map-making for TFGI btod2 or ctod2 lists of files

HISTORY:
* 26/10/2023 - original version. JARM. Only intensity.
*  3/11/2023 - added horns as input option. Limits the computation of the maps to a sublist.
* 12/02/2024 - Now including polarization map-making. Tested with TauA simulations.

"""

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from astropy.io import fits
from sancho_tfgi_io import read_tfgi_btod2
import sys

def mambrino_tfgi(root, path='/net/calp-nas/proyectos/quijote2/ctod/',tail='.ctod2', nside=512, nest=False, nhorns=7, horns=np.arange(7), dobaserm=True, t_base=6.0, usewei=True, istgi=np.ones(7,dtype=bool)):
    '''
    Keywords:
       * dobaserm : if True, removes baseline computing the median of the data in scans of length t_base.
       * t_base   : baseline size in seconds.
       * usewei   : if True, uses the WEI_IQU information.
    '''

    nf   = len(root)

    # Display basic info
    print('*** MAMBRINO_TFGI intensity and polarization code ***')
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
    aa     = np.zeros((npix,nhorns*4))
    bb     = np.zeros((npix,nhorns*4))
    cc     = np.zeros((npix,nhorns*4))
    f      = np.zeros((npix,nhorns*4))
    g      = np.zeros((npix,nhorns*4))
    nhitspol = np.zeros((npix,nhorns*4))

    for i in np.arange(nf):
        ff   = path + root[i] + tail
        btod = read_tfgi_btod2(ff)

        data   = btod['DATA']
        flag   = btod['FLAG']
        wei_iqu= btod['WEI_IQU']
        gl     = btod['GL']
        gb     = btod['GB']
        parang = btod['PAR']
        msbin  = btod['MSBIN'] 

        nsamp = len(gl[:,0])
        if nsamp != np.max(gl.shape):
            sys.exit('SYS.EXIT() -- incorrect nsamp.')
        navg  = np.int32(t_base*1000/msbin)
        nbase = np.max([1,np.int32(nsamp/navg)])
    
        for ihorn in horns:
            theta = np.deg2rad( 90.0 - gb[:,ihorn] )
            phi   = np.deg2rad( gl[:,ihorn] )
            ipix  = hp.ang2pix(nside, theta, phi, nest=nest)

            ang_ref = 0.0 ##-50.22
            psi     = (2*parang[:,ihorn] - 2*ang_ref)*np.pi/180. # in radians. No reference angle applied at the moment.

            for ichan in np.arange(4):
                k = ihorn*4 + ichan
                
                # Intensity timeline
                Vd    = ( data[:,0,k] + data[:,1,k] + data[:,2,k] + data[:,3,k] ) / 2.
                flagd = flag[:,0,k] + flag[:,1,k] + flag[:,2,k] + flag[:,3,k]
                weid  = wei_iqu[:,0,k]  # Intensity wei.

                # Polarization timelines. Assuming TGI ordering
                iph0    = 0
                iph90   = 1
                iph180  = 2
                iph270  = 3
                if istgi[ihorn]==True:
                    iph90  = 2
                    iph180 = 1
                    
                Qd      = data[:,iph0,k] - data[:,iph180,k]
                flagQd  = flag[:,iph0,k] + flag[:,iph180,k] 
                wei_qd  = wei_iqu[:,1,k] # Q wei.
                
                Ud      = data[:,iph90,k] - data[:,iph270,k]
                flagUd  = flag[:,iph90,k] + flag[:,iph270,k] 
                wei_ud  = wei_iqu[:,2,k] # U wei.
                
                # Use weights? 
                if usewei==False:
                    weid   = np.ones_like(wei_iqu[:,0,k])
                    wei_qd = np.ones_like(wei_iqu[:,1,k])
                    wei_ud = np.ones_like(wei_iqu[:,2,k])
            
                # Removing median (if requested), for the three time streams.
                if dobaserm==True:
                    Vd_smth = np.zeros_like(Vd)
                    Qd_smth = np.zeros_like(Vd)
                    Ud_smth = np.zeros_like(Vd)
                    for j in np.arange(nbase):
                        i1 = np.int64(j * navg)
                        i2 = i1 + navg #- 1 in python not needed

                        # Intensity
                        datos = Vd[i1:i2]
                        flags = flagd[i1:i2]
                        usa   = np.where(flags == 0)
                        if len(usa[-1] > 2):
                            mediana = np.median(datos[usa])
                            Vd_smth[i1:i2] = datos[:] - mediana

                        # Q in local coords
                        datos = Qd[i1:i2]
                        flags = flagQd[i1:i2]
                        usa   = np.where(flags == 0)
                        if len(usa[-1] > 2):
                            mediana = np.median(datos[usa])
                            Qd_smth[i1:i2] = datos[:] - mediana

                        # U in local coords.
                        datos = Ud[i1:i2]
                        flags = flagUd[i1:i2]
                        usa   = np.where(flags == 0)
                        if len(usa[-1] > 2):
                            mediana = np.median(datos[usa])
                            Ud_smth[i1:i2] = datos[:] - mediana
                            
                else:
                    Vd_smth = Vd
                    Qd_smth = Qd
                    Ud_smth = Ud

                # Project onto a Healpy map. Exclude flags!=0 and Vd=0    
                # a) Intensity. 
                cut = np.where( (flagd == 0) & (Vd_smth != 0 ) & (weid != 0) )
                if (len(cut[-1]) != 0):
                    np.add.at(mapa[:,k], ipix[cut], Vd_smth[cut]*weid[cut])
                    np.add.at(deno[:,k], ipix[cut], weid[cut])
                    np.add.at(nhits[:,k], ipix[cut], np.ones(len(cut)))

                # b) Polarization. Equations assume a response of the type r=1/2(I + Qcosd + Usind).
                cut = np.where( (flagQd == 0) & (Qd_smth != 0 ) & (wei_qd != 0) & (flagUd == 0) & (Ud_smth != 0 ) & (wei_ud != 0) )
                if (len(cut[-1]) != 0):
                    np.add.at(aa[:,k], ipix[cut], np.cos(psi[cut])**2 *wei_qd[cut] + np.sin(psi[cut])**2 *wei_ud[cut])
                    np.add.at(bb[:,k], ipix[cut], np.sin(psi[cut])**2 *wei_qd[cut] + np.cos(psi[cut])**2 *wei_ud[cut])
                    np.add.at(cc[:,k], ipix[cut], np.sin(psi[cut])*np.cos(psi[cut])* (wei_qd[cut] - wei_ud[cut]) )
                    np.add.at( f[:,k], ipix[cut], np.cos(psi[cut])*Qd_smth[cut]*wei_qd[cut] - np.sin(psi[cut])*Ud_smth[cut]*wei_ud[cut] )
                    np.add.at( g[:,k], ipix[cut], np.sin(psi[cut])*Qd_smth[cut]*wei_qd[cut] + np.cos(psi[cut])*Ud_smth[cut]*wei_ud[cut] )
                    np.add.at(nhitspol[:,k], ipix[cut], np.ones(len(cut)))


    # Normalising and UNSEEN pixels
    print(' Normalising maps')
    cut       = np.where(deno != 0)
    mapa[cut] = mapa[cut]/deno[cut]

    determ     = aa*bb -cc**2.
    cut        = np.where(determ !=0)
    mapaQ      = np.zeros_like(mapa)
    mapaU      = np.zeros_like(mapa)
    mapaQ[cut] = ( bb[cut]*f[cut] - cc[cut]*g[cut]) / determ[cut]
    mapaU[cut] = (-cc[cut]*f[cut] + aa[cut]*g[cut]) / determ[cut]

    print(' Setting UNSEEN pixels')
    ceros        = np.where(deno == 0)
    mapa[ceros]  = hp.UNSEEN
    nhits[ceros] = hp.UNSEEN
    
    cerosp       = np.where(determ == 0)
    mapaQ[cerosp]= hp.UNSEEN
    mapaU[cerosp]= hp.UNSEEN
    nhitspol[cerosp] = hp.UNSEEN
    
    
    print('*** MAMBRINO_TFGI code ends ***')
    #outmap = np.zeros((npix,nhorns*4,3))
    #outmap[:,:,0] = mapa
    #outmap[:,:,1] = mapaQ
    #outmap[:,:,2] = mapaU
    return mapa, nhits, mapaQ, mapaU, nhitspol


# Write maps to FITS file
def write_mambrino_tfgi_pol_maps(mapa,nhits,ffout='test_mambrino.fits'):
    col1 = fits.Column(name='map', format=str(mapa.size)+'E', array = [mapa], dim=str(mapa.shape[::-1]) )
    col2 = fits.Column(name='nhits', format=str(nhits.size)+'E', array = [nhits], dim=str(nhits.shape[::-1]) )
    hdu = fits.BinTableHDU.from_columns([col1,col2])
    hdu.writeto(ffout,overwrite=True)
    return


##############
# MAIN code
if __name__ == "__main__":

# Binned btod maps
    path  = '/net/calp-nas/proyectos/quijote2/ctod/sims/nonoise/' #'btod/test/'
    path  = '/net/calp-nas/proyectos/quijote2/ctod/' #'btod/test/'
    #path  = '/net/calp-nas/proyectos/quijote2/btod/' #'btod/test/'

    source   = 'crab'
    tail     = '.ctod2'
    nside    = 512
    nest     = False
    usewei   = True
    dobaserm = True

    ff     = '/net/calp-nas/proyectos/quijote2/list/'+source+'/'+source+'_good.txt'
    #ff     = '/net/calp-nas/proyectos/quijote2/list/crab/crab_2201.txt'
    root   = np.loadtxt(ff,unpack=True,dtype='U')
    mapa, nhits, mapaQ, mapaU, nhitspol = mambrino_tfgi(root, path=path,tail=tail, nside=nside, nest=nest, nhorns=7, usewei=usewei,
                                            horns=np.array([0]),istgi=np.array([True]), dobaserm=dobaserm)

    print('Writing maps to files')
    print(mapaQ.shape)
    txttail='_2201_ns'+str(nside)+'_ctod2'
    txttail='_good_ns'+str(nside)+'_ctod2_wei_cal'
    txttail='_good_ns'+str(nside)+'_ctod2_wei'
    hp.write_map(source+'I'+txttail+'_V1.fits', mapa[:,0], nest=False, coord='G', overwrite=True)
    hp.write_map(source+'Q'+txttail+'_V1.fits', mapaQ[:,0], nest=False, coord='G', overwrite=True)
    hp.write_map(source+'U'+txttail+'_V1.fits', mapaU[:,0], nest=False, coord='G', overwrite=True)

    print(' and writing all maps')
    write_mambrino_tfgi_pol_maps(mapa[:,0:4],nhits[:,0:4],ffout='mambrino_I_'+source+txttail+'.fits')
    write_mambrino_tfgi_pol_maps(mapaQ[:,0:4],nhitspol[:,0:4],ffout='mambrino_Q_'+source+txttail+'.fits')
    write_mambrino_tfgi_pol_maps(mapaU[:,0:4],nhitspol[:,0:4],ffout='mambrino_U_'+source+txttail+'.fits')
    
    sys.exit('Code ends.')

# Old code. Running perseus maps
  
    listpath    = '/net/calp-nas/proyectos/quijote2/list/perseus/Perseus_selected_observations_2/'
    listpath    = '/net/calp-nas/proyectos/quijote2/list/perseus/Perseus_half_maps_2/'
    txthalf     = 'half1' #'selection_2'
    root        = np.loadtxt(listpath+'Perseus_flagged_pix24_ch1_'+txthalf+'_clean.txt',unpack=True,dtype='U')
    mapa, nhits = mambrino_tfgi(root, path=path,tail=tail, nside=nside, nest=nest, nhorns=7, usewei=usewei, horns=np.array([0]))

    roots = ['Perseus_flagged_pix21_ch1_'+txthalf+'_clean.txt','Perseus_flagged_pix25_ch1_'+txthalf+'_clean.txt',
                 'Perseus_flagged_pix4_ch1_'+txthalf+'_clean.txt', 'Perseus_flagged_pix6_ch1_'+txthalf+'_clean.txt']
    horns = np.asarray([2,6,3,5])

    for j in np.arange(4):
        root = np.loadtxt(listpath+roots[j], unpack=True,dtype='U')
        mapai, nhitsi = mambrino_tfgi(root, path=path,tail=tail, nside=nside, nest=nest, nhorns=7, usewei=usewei, horns=np.array([horns[j]]))
        i1 = horns[j]*4
        i2 = i1+4
        mapa[:,i1:i2]  = np.copy(mapai[:,i1:i2])
        nhits[:,i1:i2] = np.copy(nhitsi[:,i1:i2])
    
    print('Writing...')
    write_mambrino_tfgi_maps(mapa,nhits,ffout='mambrino_perseus_ctod_flagged_select17nov23_wei_half1.fits')

    



