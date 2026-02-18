"""
26/Jan/2025   

@author: jalberto

Implementation of the naive map-making for MFI2 btod3 or ctod3 lists of files. Based on sancho_mambrino2_tfgi.py.
The code accepts a different list for each pixel. Thus, the input is now the lists of files

HISTORY:
* 26/01/2025 - original version. JARM. Only intensity.
* 27/01/2025 - now it creates a merged list of input files, which is read only once. This avoids multiple 
               readings of the same file.
* 07/02/2025 - including polarization. 
               
"""

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from astropy.io import fits
from sancho_mfi2_io import read_mfi2_btod3
from info_quijote import mfi2
import sys

# Given one list per horn, it creates the union list, and the mask to identify which horn is using what.
def get_merged_list(fflist=np.array(['test.txt']*4)):

    nff = len(fflist)
    # 1) Merged list
    for i in np.arange(nff):
        sublist = np.loadtxt(fflist[i],unpack=True,dtype='U')
        if i == 0:
            lista = sublist
        else:
            unionset = set(lista).union(set(sublist))
            lista = np.asarray( list(unionset) )

    lista.sort() # Reorder

    # 2) Mask for the different horns
    n    = len(lista)
    masc = np.zeros((n,nff),dtype='bool')
    for i in np.arange(nff):
        sublist = np.loadtxt(fflist[i],unpack=True,dtype='U')
        for j in np.arange(n):
            masc[j,i] = np.any(sublist == lista[j])

    return lista, masc

# Main routine for intensity only
def mambrino2_mfi2_intensity(fflist=np.array(['/net/calp-nas/proyectos/quijote2/mfi2/list/crab/crab_2201.txt']*4), 
                   path='/net/calp-nas/proyectos/quijote2/mfi2/ctod/', tail='.ctod3',nside=512, nest=False, 
                   nhorns=4, horns=np.arange(4)+1, dobaserm=True, t_base=6.0, usewei=True ):
    '''
    Keywords:
       * dobaserm : if True, removes baseline computing the median of the data in scans of length t_base.
       * t_base   : baseline size in seconds. Default is 6.0 seconds.
       * usewei   : if True, uses the WEI information.
    '''


    # Display basic info
    print('*** MAMBRINO2_MFI2 code for naive map-making in intensity ***')
    print('')
    print(' Settings:')
    print('    Path                            = ',path)
    print('    Output NSIDE                    = ',nside)
    print('    Remove baseline (median filter) = ',dobaserm)
    if dobaserm: print('    Baseline length (secs)          = ',t_base)
    print('    Use WEI                         = ',usewei)
    print('    Making maps for horns           = ',horns)
    print('')

    # Basic checks
    if len(horns) != len(fflist):
        sys.exit('SYS.EXIT() -- incorrect dimensions of fflist.')
 
    # Main loop over horns
    npix   = hp.nside2npix(nside)
    mapa   = np.zeros((npix,nhorns,2,2))
    deno   = np.zeros((npix,nhorns,2,2))
    nhits  = np.zeros((npix,nhorns,2,2))
 

    # First, we prepare the merged list, and masc for horns
    root, masc = get_merged_list(fflist)
    nf = len(root)
    print('    File lists                      = ',fflist)
    print('    Number of files in merged list  = ',nf)
    print('')

    # Main loop for reading CTOD3 files
    print('*** MAMBRINO2_MFI2: starting main loop, reading CTOD3/BTOD3 files')
    for ifile in np.arange(nf):
        ff   = path + root[ifile] + tail 
        btod = read_mfi2_btod3(ff)
        print(f' > File {ifile+1} out of {nf}.')

        data   = btod['DATA']
        flag   = btod['FLAG']
        wei    = btod['WEI']
        gl     = btod['GL']
        gb     = btod['GB']
        parang = btod['PAR']
        msbin  = btod['MSBIN'] 

        nsamp = len(gl[:,0])
        if nsamp != np.max(gl.shape):
            sys.exit('SYS.EXIT() -- incorrect nsamp.')
        navg  = np.int32(t_base*1000/msbin)
        nbase = np.max([1,np.int32(nsamp/navg)])
        #print(nsamp,navg,nbase)
        #print(data.shape)

    # Loop over horns
        for idxh in np.arange(len(horns)):
            ihorn   = horns[idxh]-1
            usehorn = masc[ifile,idxh]

            if usehorn==True:
                print(f' > File used in HORN={ihorn+1}, ihorn={ihorn}, and list {fflist[idxh]}')
                
                # Angles for that horn. Only needed for polarization. Not used.
                theta = np.deg2rad( 90.0 - gb[:,ihorn] )
                phi   = np.deg2rad( gl[:,ihorn] )
                ipix  = hp.ang2pix(nside, theta, phi, nest=nest)

                ang_ref = 0.0
                psi     = (2*parang[:,ihorn] - 2*ang_ref)*np.pi/180. # in radians. No reference angle applied at the moment.

                # Loop over frequencies. IFREQ
                ic = mfi2['ic']
                for ifreq in np.arange(2):
                    # Correlated or uncorrelated channels. ICORR
                    for icorr in np.arange(2):

                        ichan_s = ic[ihorn,ifreq,2*icorr]   # X+Y or X
                        ichan_d = ic[ihorn,ifreq,2*icorr+1] # X-Y or Y
                
                        # Intensity timeline.
                        Vd    = data[:,ichan_s] + data[:,ichan_d]
                        flagd = flag[:,ichan_s] + flag[:,ichan_d]
                        weid  = np.zeros_like(Vd)
                        cutd  = np.where( (wei[:,ichan_s] != 0) & (wei[:,ichan_d] != 0))
                        weid[cutd] = 1.0 / (1/wei[cutd,ichan_s] + 1/wei[cutd,ichan_d])

                        # Use weights? 
                        if usewei==False:
                            weid   = np.ones_like(wei[:,ichan_s])
            
                        # Removing median (if requested).
                        if dobaserm==True:
                            Vd_smth = np.zeros_like(Vd)
                            for j in np.arange(nbase):
                                i1 = np.int64(j * navg)
                                i2 = i1 + navg #- 1 in python not needed

                                # Intensity
                                datos = Vd[i1:i2]
                                flags = flagd[i1:i2]
                                usa   = np.where( (flags == 0) & (datos != 0))
                                if len(usa[-1] > 2):
                                    mediana = np.median(datos[usa])
                                    Vd_smth[i1:i2] = datos[:] - mediana                            
                        else:
                            Vd_smth = Vd
                

                        # Project onto a Healpy map. Exclude flags!=0 and Vd=0    
                        # a) Intensity. 
                        cut = np.where( (flagd == 0) & (Vd_smth != 0 ) & (weid != 0) )
                        if (len(cut[-1]) != 0):
                            np.add.at(mapa[:,ihorn,ifreq,icorr], ipix[cut], Vd_smth[cut]*weid[cut])
                            np.add.at(deno[:,ihorn,ifreq,icorr], ipix[cut], weid[cut])
                            np.add.at(nhits[:,ihorn,ifreq,icorr], ipix[cut], np.ones(len(cut)))

            

    # Normalising and UNSEEN pixels
    print(' Normalising maps')
    cut       = np.where(deno != 0)
    mapa[cut] = mapa[cut]/deno[cut]

    print(' Computing listpix (non zero pixels)')
    coaddmap = np.sum(mapa,axis=1)
    listpix  = np.where( coaddmap != 0 )[0]

    print(' Setting UNSEEN pixels in map ')
    ceros        = np.where(deno == 0)
    mapa[ceros]  = hp.UNSEEN
    nhits[ceros] = 0
    
    print('*** MAMBRINO2_MFI2: map-making ends ***')
    return mapa, nhits, listpix


# Main routine for intensity and polarization
def mambrino2_mfi2_intpol(fflist=np.array(['/net/calp-nas/proyectos/quijote2/mfi2/list/crab/crab_2201.txt']*4), 
                   path='/net/calp-nas/proyectos/quijote2/mfi2/ctod/', tail='.ctod3',nside=512, nest=False, 
                   nhorns=4, horns=np.arange(4)+1, dobaserm=True, t_base=6.0, usewei=True ):
    '''
    Keywords:
       * dobaserm : if True, removes baseline computing the median of the data in scans of length t_base.
       * t_base   : baseline size in seconds. Default is 6.0 seconds.
       * usewei   : if True, uses the WEI information.
    '''


    # Display basic info
    print('*** MAMBRINO2_MFI2 code for naive map-making in intensity and polarization ***')
    print('')
    print(' Settings:')
    print('    Path                            = ',path)
    print('    Output NSIDE                    = ',nside)
    print('    Remove baseline (median filter) = ',dobaserm)
    if dobaserm: print('    Baseline length (secs)          = ',t_base)
    print('    Use WEI                         = ',usewei)
    print('    Making maps for horns           = ',horns)
    print('')

    # Basic checks
    if len(horns) != len(fflist):
        sys.exit('SYS.EXIT() -- incorrect dimensions of fflist.')
 
    # Main loop over horns
    npix   = hp.nside2npix(nside)
    mapa   = np.zeros((npix,nhorns,2,2))
    deno   = np.zeros((npix,nhorns,2,2))
    nhits  = np.zeros((npix,nhorns,2,2))
    aa     = np.zeros((npix,nhorns,2,2))
    bb     = np.zeros((npix,nhorns,2,2))
    cc     = np.zeros((npix,nhorns,2,2))
    f      = np.zeros((npix,nhorns,2,2))
    g      = np.zeros((npix,nhorns,2,2))
    nhitspol = np.zeros((npix,nhorns,2,2))
 
    # First, we prepare the merged list, and masc for horns
    root, masc = get_merged_list(fflist)
    nf = len(root)
    print('    File lists                      = ',fflist)
    print('    Number of files in merged list  = ',nf)
    print('')

    # Main loop for reading CTOD3 files
    print('*** MAMBRINO2_MFI2: starting main loop, reading CTOD3/BTOD3 files')
    for ifile in np.arange(nf):
        ff   = path + root[ifile] + tail 
        btod = read_mfi2_btod3(ff)
        print(f' > File {ifile+1} out of {nf}.')

        data   = btod['DATA']
        flag   = btod['FLAG']
        wei    = btod['WEI']
        gl     = btod['GL']
        gb     = btod['GB']
        parang = btod['PAR']
        msbin  = btod['MSBIN'] 

        nsamp = len(gl[:,0])
        if nsamp != np.max(gl.shape):
            sys.exit('SYS.EXIT() -- incorrect nsamp.')
        navg  = np.int32(t_base*1000/msbin)
        nbase = np.max([1,np.int32(nsamp/navg)])
        #print(nsamp,navg,nbase)
        #print(data.shape)

    # Loop over horns
        for idxh in np.arange(len(horns)):
            ihorn   = horns[idxh]-1
            usehorn = masc[ifile,idxh]

            if usehorn==True:
                print(f' > File used in HORN={ihorn+1}, ihorn={ihorn}, and list {fflist[idxh]}')
                
                # Angles for that horn. Only needed for polarization. 
                theta = np.deg2rad( 90.0 - gb[:,ihorn] )
                phi   = np.deg2rad( gl[:,ihorn] )
                ipix  = hp.ang2pix(nside, theta, phi, nest=nest)

                ang_ref = 0.0
                psi     = (2*parang[:,ihorn] - 2*ang_ref)*np.pi/180. # radians. No reference angle applied.

                # Loop over frequencies. IFREQ
                ic = mfi2['ic']
                for ifreq in np.arange(2):
                    # Correlated or uncorrelated channels. ICORR
                    for icorr in np.arange(2):

                        ichan_s = ic[ihorn,ifreq,2*icorr]   # X+Y or X
                        ichan_d = ic[ihorn,ifreq,2*icorr+1] # X-Y or Y
                
                        # Intensity timeline.
                        Vd    = data[:,ichan_s] + data[:,ichan_d]
                        flagd = flag[:,ichan_s] + flag[:,ichan_d]
                        weid  = np.zeros_like(Vd)
                        cutd  = np.where( (wei[:,ichan_s] != 0) & (wei[:,ichan_d] != 0))
                        weid[cutd] = 1.0 / (1/wei[cutd,ichan_s] + 1/wei[cutd,ichan_d])

                        # Polarization timelines.
                        Qd      = data[:,ichan_s] - data[:,ichan_d]
                        flagQd  = flag[:,ichan_s] + flag[:,ichan_d] 
                        wei_qd  = np.zeros_like(Vd)
                        cut_qd  = np.where( (wei[:,ichan_s] != 0) & (wei[:,ichan_d] != 0))
                        wei_qd[cutd] = 1.0 / (1/wei[cut_qd,ichan_s] + 1/wei[cut_qd,ichan_d])

                        # Use weights? 
                        if usewei==False:
                            weid   = np.ones_like(wei[:,ichan_s])
                            wei_qd = np.ones_like(wei[:,ichan_s])
            
                        # Removing median (if requested).
                        if dobaserm==True:
                            Vd_smth = np.zeros_like(Vd)
                            Qd_smth = np.zeros_like(Qd)
                            for j in np.arange(nbase):
                                i1 = np.int64(j * navg)
                                i2 = i1 + navg #- 1 in python not needed

                                # Intensity
                                datos = Vd[i1:i2]
                                flags = flagd[i1:i2]
                                usa   = np.where( (flags == 0) & (datos != 0))
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
                        else:
                            Vd_smth = Vd
                            Qd_smth = Qd
                

                        # Project onto a Healpy map. Exclude flags!=0 and Vd=0    
                        # a) Intensity. 
                        cut = np.where( (flagd == 0) & (Vd_smth != 0 ) & (weid != 0) )
                        if (len(cut[-1]) != 0):
                            np.add.at(mapa[:,ihorn,ifreq,icorr], ipix[cut], Vd_smth[cut]*weid[cut])
                            np.add.at(deno[:,ihorn,ifreq,icorr], ipix[cut], weid[cut])
                            np.add.at(nhits[:,ihorn,ifreq,icorr], ipix[cut], np.ones(len(cut)))
                        # b) Polarization.
                        cut = np.where( (flagQd == 0) & (Qd_smth != 0 ) & (wei_qd != 0) )
                        if (len(cut[-1]) != 0):
                            np.add.at(aa[:,ihorn,ifreq,icorr], ipix[cut], np.cos(psi[cut])**2 *wei_qd[cut] )
                            np.add.at(bb[:,ihorn,ifreq,icorr], ipix[cut], np.sin(psi[cut])**2 *wei_qd[cut] )
                            np.add.at(cc[:,ihorn,ifreq,icorr], ipix[cut], np.sin(psi[cut])*np.cos(psi[cut])* wei_qd[cut] )
                            np.add.at( f[:,ihorn,ifreq,icorr], ipix[cut], np.cos(psi[cut])*Qd_smth[cut]*wei_qd[cut] )
                            np.add.at( g[:,ihorn,ifreq,icorr], ipix[cut], np.sin(psi[cut])*Qd_smth[cut]*wei_qd[cut] )
                            np.add.at(nhitspol[:,ihorn,ifreq,icorr], ipix[cut], np.ones(len(cut)))


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

    print(' Computing listpix (non zero pixels)')
    coaddmap = np.sum(mapa,axis=1)
    listpix  = np.where( coaddmap != 0 )[0]

    print(' Setting UNSEEN pixels in maps (int and pol) ')
    ceros        = np.where(deno == 0)
    mapa[ceros]  = hp.UNSEEN
    nhits[ceros] = 0

    cerosp       = np.where(determ == 0)
    mapaQ[cerosp]= hp.UNSEEN
    mapaU[cerosp]= hp.UNSEEN
    nhitspol[cerosp] = hp.UNSEEN
    
    print('*** MAMBRINO2_MFI2: map-making for intensity and polarization ends ***')
    return mapa, nhits, listpix, mapaQ, mapaU, nhitspol


# Write maps to FITS files
def write_mambrino2_mfi2_maps(mapa, nhits, listpix, ffout='test_mambrino_mfi2.fits'):

    print('*** MAMBRINO2_MFI2: Writing intensity maps to files *** ')
    
    nlist     = len(listpix)
    nchan     = len(mapa[0,:])
    submapa   = np.zeros((nlist,4,2,2))
    submapaQ  = np.zeros((nlist,4,2,2))
    submapaU  = np.zeros((nlist,4,2,2))
    subnhits  = np.zeros((nlist,4,2,2))
    subnhitspol = np.zeros((nlist,4,2,2))

    for ihorn in np.arange(4):
        for ifreq in np.arange(2):
            for icorr in np.arange(2):
                submapa[:,ihorn,ifreq,icorr]  = mapa[listpix[:],ihorn,ifreq,icorr]
                subnhits[:,ihorn,ifreq,icorr] = nhits[listpix[:],ihorn,ifreq,icorr]


    col1 = fits.Column(name='nside',format='K',array=[nside])
    col2 = fits.Column(name='listpix', format=str(listpix.size)+'E', array = [listpix]) 
    col3 = fits.Column(name='map', format=str(submapa.size)+'E', array = [submapa], dim=str(submapa.shape[::-1]) )
    col4 = fits.Column(name='nhits', format=str(subnhits.size)+'E', array = [subnhits], dim=str(subnhits.shape[::-1]) )
    
    hdu = fits.BinTableHDU.from_columns([col1,col2,col3,col4])
    hdu.writeto(ffout,overwrite=True)
    return

def write_mambrino2_mfi2_maps_intpol(mapa, nhits, listpix, mapaQ, mapaU, nhitspol, ffout='test_mambrino_mfi2_intpol.fits'):

    print('*** MAMBRINO2_MFI2: Writing int and pol maps to files *** ')
    
    nlist     = len(listpix)
    nchan     = len(mapa[0,:])
    submapa   = np.zeros((nlist,4,2,2))
    submapaQ  = np.zeros((nlist,4,2,2))
    submapaU  = np.zeros((nlist,4,2,2))
    subnhits  = np.zeros((nlist,4,2,2))
    subnhitspol = np.zeros((nlist,4,2,2))

    for ihorn in np.arange(4):
        for ifreq in np.arange(2):
            for icorr in np.arange(2):
                submapa[:,ihorn,ifreq,icorr]  = mapa[listpix[:],ihorn,ifreq,icorr]
                subnhits[:,ihorn,ifreq,icorr] = nhits[listpix[:],ihorn,ifreq,icorr]
                submapaQ[:,ihorn,ifreq,icorr]  = mapaQ[listpix[:],ihorn,ifreq,icorr]
                submapaU[:,ihorn,ifreq,icorr]  = mapaU[listpix[:],ihorn,ifreq,icorr]
                subnhitspol[:,ihorn,ifreq,icorr] = nhitspol[listpix[:],ihorn,ifreq,icorr]

    col1 = fits.Column(name='nside',format='K',array=[nside])
    col2 = fits.Column(name='listpix', format=str(listpix.size)+'E', array = [listpix]) 
    col3 = fits.Column(name='mapI', format=str(submapa.size)+'E', array = [submapa], dim=str(submapa.shape[::-1]) )
    col4 = fits.Column(name='nhits', format=str(subnhits.size)+'E', array = [subnhits], dim=str(subnhits.shape[::-1]) )
    col5 = fits.Column(name='mapQ', format=str(submapaQ.size)+'E', array = [submapaQ], dim=str(submapaQ.shape[::-1]) )
    col6 = fits.Column(name='mapU', format=str(submapaU.size)+'E', array = [submapaU], dim=str(submapaU.shape[::-1]) )
    col7 = fits.Column(name='nhitspol', format=str(subnhitspol.size)+'E', array = [subnhitspol], dim=str(subnhitspol.shape[::-1]) )


    hdu = fits.BinTableHDU.from_columns([col1,col2,col3,col4,col5,col6,col7])
    hdu.writeto(ffout,overwrite=True)
    return



##############
# MAIN code
if __name__ == "__main__":

# Binned btod maps
    path  = '/net/calp-nas/proyectos/quijote2/mfi2/ctod/'
    pathout = '/net/calp-nas/proyectos/quijote2/mfi2/maps/'
    
    tail     = '.ctod3'
    nside    = 512
    nest     = False
    usewei   = True
    dobaserm = True
    nhorns   = 4

    # Testing intensity
    dothis = 0
    if dothis !=0:
        ff = 'test2.txt'
        mapa, nhits, listpix = mambrino2_mfi2_intensity(np.array([ff]*nhorns), path=path,tail=tail, nside=nside, nest=nest, 
                                              nhorns=nhorns, usewei=usewei,horns=np.arange(nhorns)+1, dobaserm=dobaserm)
        #hp.write_map('test2.fits', mapa[:,0], nest=False, coord='G', overwrite=True)

        write_mambrino2_mfi2_maps(mapa, nhits, listpix, ffout='test2_mambrino_mfi2.fits')
        
    # Testing intensity and polarization
    dothis = 0
    if dothis !=0:
        ff = 'test2.txt'
        mapa, nhits, listpix, mapaQ, mapaU, nhitspol = mambrino2_mfi2_intpol(np.array([ff]*nhorns), path=path,tail=tail, nside=nside, nest=nest, 
                                              nhorns=nhorns, usewei=usewei,horns=np.arange(nhorns)+1, dobaserm=dobaserm)
        #hp.write_map('test2.fits', mapa[:,0], nest=False, coord='G', overwrite=True)

        write_mambrino2_mfi2_maps(mapa, nhits, listpix, ffout='test3_mambrino_mfi2_int.fits')
        write_mambrino2_mfi2_maps(mapaQ, nhitspol, listpix, ffout='test3_mambrino_mfi2_q.fits')
        write_mambrino2_mfi2_maps(mapaU, nhitspol, listpix, ffout='test3_mambrino_mfi2_u.fits')
        write_mambrino2_mfi2_maps_intpol(mapa, nhits, listpix, mapaQ, mapaU, nhitspol, ffout='test3_mambrino_mfi2_iqu.fits')


    # Tau A observations Nov-Dec2024
    dothis = 0
    if dothis !=0:
        ff    = '/net/calp-nas/proyectos/quijote2/mfi2/list/crab/crab_nov_dic2024.txt'
        ffout = pathout+'mambrino2_mfi2_crab_nov_dic2024.fits'
        nhorn  = 4
        mapa, nhits, listpix = mambrino2_mfi2_intensity(np.array([ff]*nhorns), path=path,tail=tail, nside=nside, nest=nest, 
                                              nhorns=nhorns, usewei=usewei,horns=np.arange(nhorns)+1, dobaserm=dobaserm)
        write_mambrino2_mfi2_maps(mapa, nhits, listpix, ffout=ffout)


    # NOMINAL60 observations until Jan2025
    dothis = 0
    if dothis !=0:
        ff    = '/net/calp-nas/proyectos/quijote2/mfi2/list/nominal/nominal60_allbtod.txt'
        ffout = pathout+'mambrino2_mfi2_nominal60_until_20250113_jflag.fits'
        path  = '/net/calp-nas/proyectos/quijote2/mfi2/ctod/jflag/'
        mapa, nhits, listpix = mambrino2_mfi2_intensity(np.array([ff]*nhorns), path=path,tail=tail, nside=nside, nest=nest, 
                                              nhorns=nhorns, usewei=usewei,horns=np.arange(nhorns)+1, dobaserm=dobaserm)
        write_mambrino2_mfi2_maps(mapa, nhits, listpix, ffout=ffout)

    dothis = 0
    if dothis !=0:
        path  = '/net/calp-nas/proyectos/quijote2/mfi2/ctod/jflag/'

        ff    = '/net/calp-nas/proyectos/quijote2/mfi2/list/nominal/nominal60_allbtod_half1.txt'
        ffout = pathout+'mambrino2_mfi2_nominal60_until_20250113_jflag_half1.fits'
        mapa, nhits, listpix = mambrino2_mfi2_intensity(np.array([ff]*nhorns), path=path,tail=tail, nside=nside, nest=nest, 
                                              nhorns=nhorns, usewei=usewei,horns=np.arange(nhorns)+1, dobaserm=dobaserm)
        write_mambrino2_mfi2_maps(mapa, nhits, listpix, ffout=ffout)

        ff    = '/net/calp-nas/proyectos/quijote2/mfi2/list/nominal/nominal60_allbtod_half2.txt'
        ffout = pathout+'mambrino2_mfi2_nominal60_until_20250113_jflag_half2.fits'
        mapa, nhits, listpix = mambrino2_mfi2_intensity(np.array([ff]*nhorns), path=path,tail=tail, nside=nside, nest=nest, 
                                              nhorns=nhorns, usewei=usewei,horns=np.arange(nhorns)+1, dobaserm=dobaserm)
        write_mambrino2_mfi2_maps(mapa, nhits, listpix, ffout=ffout)


    # NOMINAL60 observations until Jan2025, intensity and polarization
    dothis = 1
    if dothis !=0:
        ff    = '/net/calp-nas/proyectos/quijote2/mfi2/list/nominal/nominal60_allbtod.txt'
        path  = '/net/calp-nas/proyectos/quijote2/mfi2/ctod/jflag/'
        mapa, nhits, listpix, mapaQ, mapaU, nhitspol = mambrino2_mfi2_intpol(np.array([ff]*nhorns), path=path,tail=tail, nside=nside, nest=nest, 
                                              nhorns=nhorns, usewei=usewei,horns=np.arange(nhorns)+1, dobaserm=dobaserm)
        #ffout1 = pathout+'mambrino2_mfi2_nominal60_until_20250113_jflag_int.fits'
        #ffout2 = pathout+'mambrino2_mfi2_nominal60_until_20250113_jflag_q.fits'
        #ffout3 = pathout+'mambrino2_mfi2_nominal60_until_20250113_jflag_u.fits'
        #write_mambrino2_mfi2_maps(mapa, nhits, listpix, ffout=ffout1)
        #write_mambrino2_mfi2_maps(mapaQ, nhitspol, listpix, ffout=ffout2)
        #write_mambrino2_mfi2_maps(mapaU, nhitspol, listpix, ffout=ffout3)
        ffout4 = pathout+'mambrino2_mfi2_nominal60_until_20250113_jflag_iqu.fits'
        write_mambrino2_mfi2_maps_intpol(mapa, nhits, listpix, mapaQ, mapaU, nhitspol, ffout=ffout4)


    sys.exit('*** MAMBRINO2_MFI2: Code ends. ***')

