#!/usr/bin/env python

"""
26/Oct/2023   

@author: jalberto

Implementation of the naive map-making for TFGI btod2 or ctod2 lists of files

HISTORY:
* 26/10/2023 - original version. JARM. Only intensity.
*  3/11/2023 - added horns as input option. Limits the computation of the maps to a sublist.
* 12/02/2024 - Now including polarization map-making. Tested with TauA simulations.
* 17/12/2024 - sancho_mambrino2_tfgi. Now accepts a different list for each pixel/DAS number, as it is the usual case. 
               Thus, the input is now the lists of files
* 10/02/2025 - using tfgi_params.py, and reading reference angles.
* 12/02/2025 - new keyword usezeroangle. If True, the pol maps are done assuming ref_angle=0. To be used in calibration.
* 19/02/2025 - Using poleff from TauA. New keyword usepoleff added, to replicate old values.
*  6/07/2026 - Updated to compute the error map (errmap[npix,3,nmaps]) and the covariance covQU.

"""

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from astropy.io import fits
from sancho_tfgi_io import read_tfgi_btod2
from tfgi_params import *
import sys


def mambrino2_tfgi(fflist=np.array(['/net/calp-nas/proyectos/quijote2/tfgi/list/crab/crab_2201.txt']*7), path='/net/calp-nas/proyectos/quijote2/tfgi/ctod/',tail='.ctod2',
                   nside=512, nest=False, nhorns=7, horns=np.arange(7), dobaserm=True, t_base=6.0, usewei=True, istgi=np.ones(7,dtype=bool), 
                   usezeroangle=False, useflag=True, usepoleff=False, usephasestates=0):
    '''
    MAMBRINO2_TFGI map-making code. Makes maps for TFGI data.
    Based on a naive map-making approach, using baseline removal with median filter.

    Inputs:
        * fflist. Input list of files to be read. Array of length equal to nhorns. Default is 7.
        * nhorns. Number of horns in the TFGI configuration.

    Keywords:
        * dobaserm : if True, removes baseline computing the median of the data in scans of length t_base.
        * t_base   : baseline size in seconds. Default is 6.0 seconds.
        * usewei   : if True, uses the WEI_IQU information.
        * usephasestates: default is 0 (uses all phase states). If 1 or 2, uses only one of the combinations.
    '''


    # Display basic info
    print('******************************************************')
    print('*** MAMBRINO2_TFGI intensity and polarization code ***')
    print('*** Version: July 6th, 2026')
    print('')
    print(' Settings:')
    print('    Path                             = ',path)
    print('    Output NSIDE                     = ',nside)
    print('    Remove baseline (median filter)  = ',dobaserm)
    if dobaserm: print('    Baseline length (secs)           = ',t_base)
    print('    usezeroangle (if true, ref_ang=0) = ',usezeroangle)
    print('    usepoleff                        = ',usepoleff)
    print('    Use WEI_IQU                      = ',usewei)
    print('    Use FLAGS (default is True)      = ',useflag)
    print('    Making maps for horns            = ',horns)
    print('    Use_phase_states (0=both,1,2)    = ',usephasestates)
    print('')

    # Basic checks
    if len(horns) != len(fflist):
        sys.exit('SYS.EXIT() -- incorrect dimensions of fflist.')
    if len(horns) != len(istgi):
        sys.exit('SYS.EXIT() -- incorrect dimensions of istgi.')

        
    # Main loop over horns
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
    mapaQ  = np.zeros_like(mapa)
    mapaU  = np.zeros_like(mapa)
    errI   = np.zeros((npix,nhorns*4))
    errQ   = np.zeros((npix,nhorns*4))
    errU   = np.zeros((npix,nhorns*4))
    covQU  = np.zeros((npix,nhorns*4))

    for idxh in np.arange(len(horns)):
        ihorn = horns[idxh]
        
        # Read list of CTOD files for that horn.
        print(f'(*) HORN: ihorn={ihorn}, and file {fflist[idxh]}')
        root = np.loadtxt(fflist[idxh],unpack=True,dtype='U')
        nf   = len(root)

        # Loop over tod files
        use_das = True
        for i in np.arange(nf):
            ff   = path + root[i] + tail 
            btod = read_tfgi_btod2(ff)

            jd     = btod['JD']
            data   = btod['DATA']
            flag   = btod['FLAG']
            wei_iqu= btod['WEI_IQU']
            gl     = btod['GL']
            gb     = btod['GB']
            parang = btod['PAR']
            msbin  = btod['MSBIN'] 

            # Use flags? For planets, we might want to deactivate this.
            if useflag==False:
                flag[:,:,:] = 0

            nsamp = len(gl[:,0])
            if nsamp != np.max(gl.shape):
                sys.exit('SYS.EXIT() -- incorrect nsamp.')
            navg  = np.int32(t_base*1000/msbin)
            nbase = np.max([1,np.int32(nsamp/navg)])
            #print(nsamp,navg,nbase)

            # Initialize das values
            if use_das:
                das, horn_fp_loc, pix_id = read_tfgi_info_etc_file(jd[0])
                use_das = False
                
            # Angles for that horn
            theta = np.deg2rad( 90.0 - gb[:,ihorn] )
            phi   = np.deg2rad( gl[:,ihorn] )
            ipix  = hp.ang2pix(nside, theta, phi, nest=nest)

            #ang_ref_matrix1, ang_ref_matrix2 = tfgi_ref_angle_2ps(jd[0],das=das)
            #ang_ref_matrix = tfgi_ref_angle(jd[0],das=das,jarm=True)
            ang_ref_matrix = tfgi_ref_angle_monthly(jd[0],das=das,jarm=False) # use MFT values
            #print(ang_ref_matrix)

            # Polarization efficiencies
            #poleff = tfgi_poleff(jd[0],das=das)
            #poleff1, poleff2 = tfgi_poleff_2ps(jd[0],das=das)
            poleff1 = tfgi_poleff(jd[0],das=das)
            poleff2 = poleff1
            #print(poleff1)
            

            # Loop over ichan   
            for ichan in np.arange(4):
                k = ihorn*4 + ichan
                ang_ref = ang_ref_matrix[ihorn,ichan]
                if usezeroangle: 
                    ang_ref = 0.0 
                psi     = 2*(parang[:,ihorn] - ang_ref)*np.pi/180. # in radians. 

                rho1 = 1.0
                rho2 = 1.0
                if usepoleff:
                    rho1 = 1/poleff1[ihorn,ichan]
                    rho2 = 1/poleff2[ihorn,ichan]

                # Intensity timeline
                Vd    = ( data[:,0,k] + data[:,1,k] + data[:,2,k] + data[:,3,k] ) / 2.
                flagd = flag[:,0,k] + flag[:,1,k] + flag[:,2,k] + flag[:,3,k]
                weid  = wei_iqu[:,0,k]  # Intensity wei.

                # Polarization timelines. Assuming TGI ordering
                iph0    = 0
                iph90   = 1
                iph180  = 2
                iph270  = 3
                if istgi[idxh]==True:
                    iph90  = 2
                    iph180 = 1
                    
                Qd      = rho1*(data[:,iph0,k] - data[:,iph180,k])
                flagQd  = flag[:,iph0,k] + flag[:,iph180,k] 
                wei_qd  = wei_iqu[:,1,k] # Q wei.
                
                Ud      = rho2*(data[:,iph90,k] - data[:,iph270,k])
                flagUd  = flag[:,iph90,k] + flag[:,iph270,k] 
                wei_ud  = wei_iqu[:,2,k] # U wei.
                
                # Use weights? 
                if usewei==False:
                    weid   = np.ones_like(wei_iqu[:,0,k])
                    wei_qd = np.ones_like(wei_iqu[:,1,k])
                    wei_ud = np.ones_like(wei_iqu[:,2,k])


                # Use both Qd and Ud, or only one?
                if usephasestates==1: 
                    flagUd += 1
                    wei_ud *= 0
                if usephasestates==2:
                    flagQd += 1
                    wei_qd *= 0      
            
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
                if usephasestates==1:
                    cut = np.where( (flagQd == 0) & (Qd_smth != 0 ) & (wei_qd != 0) )
                elif usephasestates==2:
                    cut = np.where( (flagUd == 0) & (Ud_smth != 0 ) & (wei_ud != 0) )
                else:
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
    # Computing error map I
    errI[cut] = 1./np.sqrt(deno[cut])

    # Polarization
    determ     = aa*bb -cc*cc
    cut        = np.where(determ !=0)
    mapaQ[cut] = ( bb[cut]*f[cut] - cc[cut]*g[cut]) / determ[cut]
    mapaU[cut] = (-cc[cut]*f[cut] + aa[cut]*g[cut]) / determ[cut]

    # Computing error maps Q and U, and covariance QU
    errQ[cut]  = np.sqrt(bb[cut]/determ[cut]) 
    errU[cut]  = np.sqrt(aa[cut]/determ[cut]) 
    covQU[cut] = -cc[cut]/determ[cut] 

    print(' Computing listpix (non zero pixels)')
    coaddmap = np.sum(mapa+mapaQ+mapaU,axis=1)
    listpix  = np.where( coaddmap != 0 )[0]


    print(' Setting UNSEEN pixels')
    ceros        = np.where(deno == 0)
    mapa[ceros]  = hp.UNSEEN
    nhits[ceros] = hp.UNSEEN
    errI[ceros]  = hp.UNSEEN
    
    cerosp       = np.where(determ == 0)
    mapaQ[cerosp]= hp.UNSEEN
    mapaU[cerosp]= hp.UNSEEN
    nhitspol[cerosp] = hp.UNSEEN
    errQ[cerosp] = hp.UNSEEN
    errU[cerosp] = hp.UNSEEN
    covQU[cerosp] = hp.UNSEEN
    
    print('*** MAMBRINO_TFGI code ends ***')

    return mapa, nhits, mapaQ, mapaU, nhitspol, listpix, errI, errQ, errU, covQU


# Write maps to FITS file
def write_mambrino2_tfgi_maps(mapa, nhits, mapaQ, mapaU, nhitspol, listpix, errI, errQ, errU, covQU, ffout='test_mambrino.fits'):
#    mapa2    = np.copy(mapa)
#    mapa2[np.where(mapa == hp.UNSEEN)]=0.0
    
#    coaddmap = np.sum(mapa2,axis=1)
#    listpix  = np.where( coaddmap != 0 )[0]

    print('*** MAMBRINO_TFGI: Writing maps to files *** ')
    
    nlist     = len(listpix)
    nchan     = len(mapa[0,:])
    submapa   = np.zeros((nlist,nchan))
    submapaQ  = np.zeros((nlist,nchan))
    submapaU  = np.zeros((nlist,nchan))
    subnhits  = np.zeros((nlist,nchan))
    subnhitspol = np.zeros((nlist,nchan))
    suberrI   = np.zeros((nlist,nchan))
    suberrQ   = np.zeros((nlist,nchan))
    suberrU   = np.zeros((nlist,nchan))
    subcovQU  = np.zeros((nlist,nchan))


    for k in np.arange(nchan):
        submapa[:,k]  = mapa[listpix[:],k]
        submapaQ[:,k] = mapaQ[listpix[:],k]
        submapaU[:,k] = mapaU[listpix[:],k]
        subnhits[:,k] = nhits[listpix[:],k]
        subnhitspol[:,k] = nhitspol[listpix[:],k]
        suberrI[:,k]  = errI[listpix[:],k]
        suberrQ[:,k]  = errQ[listpix[:],k]
        suberrU[:,k]  = errU[listpix[:],k]
        subcovQU[:,k] = covQU[listpix[:],k]


    col1 = fits.Column(name='nside',format='K',array=[nside])
    col2 = fits.Column(name='listpix', format=str(listpix.size)+'E', array = [listpix]) 
    col3 = fits.Column(name='map', format=str(submapa.size)+'E', array = [submapa], dim=str(submapa.shape[::-1]) )
    col4 = fits.Column(name='nhits', format=str(subnhits.size)+'E', array = [subnhits], dim=str(subnhits.shape[::-1]) )
    col5 = fits.Column(name='mapQ', format=str(submapaQ.size)+'E', array = [submapaQ], dim=str(submapaQ.shape[::-1]) )
    col6 = fits.Column(name='mapU', format=str(submapaU.size)+'E', array = [submapaU], dim=str(submapaU.shape[::-1]) )
    col7 = fits.Column(name='nhitspol', format=str(subnhitspol.size)+'E', array = [subnhitspol], dim=str(subnhitspol.shape[::-1]) )

    col8  = fits.Column(name='errI', format=str(suberrI.size)+'E', array = [suberrI], dim=str(suberrI.shape[::-1]) )
    col9  = fits.Column(name='errQ', format=str(suberrQ.size)+'E', array = [suberrQ], dim=str(suberrQ.shape[::-1]) )
    col10 = fits.Column(name='errU', format=str(suberrU.size)+'E', array = [suberrU], dim=str(suberrU.shape[::-1]) )
    col11 = fits.Column(name='covQU', format=str(subcovQU.size)+'E', array = [subcovQU], dim=str(subcovQU.shape[::-1]) )

    hdu = fits.BinTableHDU.from_columns([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11])
    hdu.writeto(ffout,overwrite=True)
    return


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
    path  = '/net/calp-nas/proyectos/quijote2/tfgi/ctod/dec24/'
    pathout = '/net/calp-nas/proyectos/quijote2/tfgi/maps/'
    
    tail     = '.ctod2'
    nside    = 512
    nest     = False
    usewei   = True
    dobaserm = True


    # Testing
    dothis = 0
    if dothis !=0:
        ff     = '/net/calp-nas/proyectos/quijote2/tfgi/list/crab/crab_2111_select_das24_pix23_allbtod.txt'
        ff     = 'test.txt'
        mapa, nhits, mapaQ, mapaU, nhitspol, listpix = mambrino2_tfgi(np.array([ff]*7), path=path,tail=tail, nside=nside, nest=nest, nhorns=7, usewei=usewei,
                                                                      horns=np.arange(7),istgi=np.array([True]*7), dobaserm=dobaserm)
        #    hp.write_map('test1.fits', mapa[:,0], nest=False, coord='G', overwrite=True)

        
    # 18/Dic72024. CRAB maps. 
    dothis = 0
    if dothis !=0:
        fflist = np.array(['/net/calp-nas/proyectos/quijote2/tfgi/list/crab/crab_select_allbtod_das24.txt',
                           '/net/calp-nas/proyectos/quijote2/tfgi/list/crab/crab_select_allbtod_das21.txt',
                           '/net/calp-nas/proyectos/quijote2/tfgi/list/crab/crab_select_allbtod_das5.txt'] )
        horns = np.array([0,2,4])
        istgi = np.array([True,True,False])
        mapa, nhits, mapaQ, mapaU, nhitspol, listpix = mambrino2_tfgi(fflist, path=path,tail=tail, nside=nside, nest=nest, nhorns=7, usewei=usewei,
                                                                      horns=horns,istgi=istgi, dobaserm=dobaserm)
        write_mambrino2_tfgi_maps(mapa, nhits, mapaQ, mapaU, nhitspol, listpix, ffout=pathout+'mambrino_crab_select_dec2024_wei_baserm.fits')


    # 18/Dic72024. W63 maps.
    dothis = 0
    if dothis !=0:
        fflist = np.array(['/net/calp-nas/proyectos/quijote2/tfgi/list/W63_cygnus/W63_scans_for_destriper_05_17_times_RMS_allbtod2_das24_W63.lst',
                           '/net/calp-nas/proyectos/quijote2/tfgi/list/W63_cygnus/W63_scans_for_destriper_05_17_times_RMS_allbtod2_das21_W63.lst',
                           '/net/calp-nas/proyectos/quijote2/tfgi/list/W63_cygnus/W63_scans_for_destriper_05_17_times_RMS_allbtod2_das5_W63.lst'])
        horns = np.array([0,2,4])
        istgi = np.array([True,True,False])
        mapa, nhits, mapaQ, mapaU, nhitspol, listpix = mambrino2_tfgi(fflist, path=path,tail=tail, nside=nside, nest=nest, nhorns=7, usewei=usewei,
                                                                      horns=horns,istgi=istgi, dobaserm=dobaserm)
        write_mambrino2_tfgi_maps(mapa, nhits, mapaQ, mapaU, nhitspol, listpix, ffout=pathout+'mambrino_w63_select_dec2024_wei_baserm.fits')

    
    # 18/Dic/2024. W44 maps.
    dothis = 0
    if dothis !=0:
        fflist = np.array(['/net/calp-nas/proyectos/quijote2/tfgi/list/W44/W44_scans_for_destriper_05_17_times_RMS_allbtod2_das24.lst',
                           '/net/calp-nas/proyectos/quijote2/tfgi/list/W44/W44_scans_for_destriper_05_17_times_RMS_allbtod2_das21.lst',
                           '/net/calp-nas/proyectos/quijote2/tfgi/list/W44/W44_scans_for_destriper_05_17_times_RMS_allbtod2_das5.lst' ] )
        horns = np.array([0,2,4])
        istgi = np.array([True,True,False])
        mapa, nhits, mapaQ, mapaU, nhitspol, listpix = mambrino2_tfgi(fflist, path=path,tail=tail, nside=nside, nest=nest, nhorns=7, usewei=usewei,
                                                                      horns=horns,istgi=istgi, dobaserm=dobaserm)
        write_mambrino2_tfgi_maps(mapa, nhits, mapaQ, mapaU, nhitspol, listpix, ffout=pathout+'mambrino_w44_select_dec2024_wei_baserm.fits')

    # 10/Feb/2025. W44 maps but with calibrated angle
    dothis = 0
    if dothis !=0:
        fflist = np.array(['/net/calp-nas/proyectos/quijote2/tfgi/list/W44/W44_scans_for_destriper_05_17_times_RMS_allbtod2_das24.lst',
                           '/net/calp-nas/proyectos/quijote2/tfgi/list/W44/W44_scans_for_destriper_05_17_times_RMS_allbtod2_das21.lst',
                           '/net/calp-nas/proyectos/quijote2/tfgi/list/W44/W44_scans_for_destriper_05_17_times_RMS_allbtod2_das5.lst' ] )
        horns = np.array([0,2,4])
        istgi = np.array([True,True,False])
        mapa, nhits, mapaQ, mapaU, nhitspol, listpix = mambrino2_tfgi(fflist, path=path,tail=tail, nside=nside, nest=nest, nhorns=7, usewei=usewei,
                                                                      horns=horns,istgi=istgi, dobaserm=dobaserm)
        write_mambrino2_tfgi_maps(mapa, nhits, mapaQ, mapaU, nhitspol, listpix, ffout=pathout+'mambrino_w44_select_feb2025_wei_baserm_angle_final.fits')

    # 11/Feb/2025. W63 maps with calibrated angle.
    dothis = 0
    if dothis !=0:
        fflist = np.array(['/net/calp-nas/proyectos/quijote2/tfgi/list/W63_cygnus/W63_scans_for_destriper_05_17_times_RMS_allbtod2_das24_W63.lst',
                           '/net/calp-nas/proyectos/quijote2/tfgi/list/W63_cygnus/W63_scans_for_destriper_05_17_times_RMS_allbtod2_das21_W63.lst',
                           '/net/calp-nas/proyectos/quijote2/tfgi/list/W63_cygnus/W63_scans_for_destriper_05_17_times_RMS_allbtod2_das5_W63.lst'])
        horns = np.array([0,2,4])
        istgi = np.array([True,True,False])
        mapa, nhits, mapaQ, mapaU, nhitspol, listpix = mambrino2_tfgi(fflist, path=path,tail=tail, nside=nside, nest=nest, nhorns=7, usewei=usewei,
                                                                      horns=horns,istgi=istgi, dobaserm=dobaserm)
        write_mambrino2_tfgi_maps(mapa, nhits, mapaQ, mapaU, nhitspol, listpix, ffout=pathout+'mambrino_w63_select_feb20245_wei_baserm_angle_final.fits')

    # 12/Feb/2025. CRAB maps, before and after Apr2022 
    dothis = 0
    if dothis !=0:
        fflist = np.array(['/net/calp-nas/proyectos/quijote2/tfgi/list/crab/crab_select_allbtod_das24_beforeApr2022.txt',
                           '/net/calp-nas/proyectos/quijote2/tfgi/list/crab/crab_select_allbtod_das21_beforeApr2022.txt',
                           '/net/calp-nas/proyectos/quijote2/tfgi/list/crab/crab_select_allbtod_das5_beforeApr2022.txt'] )
        horns = np.array([0,2,4])
        istgi = np.array([True,True,False])
        mapa, nhits, mapaQ, mapaU, nhitspol, listpix = mambrino2_tfgi(fflist, path=path,tail=tail, nside=nside, nest=nest, nhorns=7, usewei=usewei,
                                                                      horns=horns,istgi=istgi, dobaserm=dobaserm, usezeroangle=True)
        write_mambrino2_tfgi_maps(mapa, nhits, mapaQ, mapaU, nhitspol, listpix, ffout=pathout+'mambrino_crab_select_feb2025_wei_baserm_beforeApr2022.fits')

    if dothis !=0:
        fflist = np.array(['/net/calp-nas/proyectos/quijote2/tfgi/list/crab/crab_select_allbtod_das24_afterApr2022.txt',
                           '/net/calp-nas/proyectos/quijote2/tfgi/list/crab/crab_select_allbtod_das21_afterApr2022.txt',
                           '/net/calp-nas/proyectos/quijote2/tfgi/list/crab/crab_select_allbtod_das5_afterApr2022.txt'] )
        horns = np.array([0,2,4])
        istgi = np.array([True,True,False])
        mapa, nhits, mapaQ, mapaU, nhitspol, listpix = mambrino2_tfgi(fflist, path=path,tail=tail, nside=nside, nest=nest, nhorns=7, usewei=usewei,
                                                                      horns=horns,istgi=istgi, dobaserm=dobaserm, usezeroangle=True)
        write_mambrino2_tfgi_maps(mapa, nhits, mapaQ, mapaU, nhitspol, listpix, ffout=pathout+'mambrino_crab_select_feb2025_wei_baserm_afterApr2022.fits')


    # 12/Feb/2025. W44 half1/2 maps with calibrated angle
    dothis = 0
    if dothis !=0:
        txtnull = ['half1','half2']
        for txt in txtnull:
            fflist = np.array(['/net/calp-nas/proyectos/quijote2/tfgi/list/W44/W44_scans_for_destriper_05_17_times_RMS_allbtod2_das24_'+txt+'.txt', 
                               '/net/calp-nas/proyectos/quijote2/tfgi/list/W44/W44_scans_for_destriper_05_17_times_RMS_allbtod2_das21_'+txt+'.txt',
                               '/net/calp-nas/proyectos/quijote2/tfgi/list/W44/W44_scans_for_destriper_05_17_times_RMS_allbtod2_das5_'+txt+'.txt' ] )
            horns = np.array([0,2,4])
            istgi = np.array([True,True,False])
            mapa, nhits, mapaQ, mapaU, nhitspol, listpix = mambrino2_tfgi(fflist, path=path,tail=tail, nside=nside, nest=nest, nhorns=7, usewei=usewei,
                                                                      horns=horns,istgi=istgi, dobaserm=dobaserm)
            write_mambrino2_tfgi_maps(mapa, nhits, mapaQ, mapaU, nhitspol, listpix, ffout=pathout+'mambrino_w44_'+txt+'_select_feb2025_wei_baserm_angle_final.fits')
        
    # 12/Feb/2025. W63 half1/2 maps with calibrated angle
    dothis = 0
    if dothis !=0:
        txtnull = ['half1','half2']
        for txt in txtnull:
            fflist = np.array(['/net/calp-nas/proyectos/quijote2/tfgi/list/W63_cygnus/W63_scans_for_destriper_05_17_times_RMS_allbtod2_das24_'+txt+'.txt', 
                               '/net/calp-nas/proyectos/quijote2/tfgi/list/W63_cygnus/W63_scans_for_destriper_05_17_times_RMS_allbtod2_das21_'+txt+'.txt',
                               '/net/calp-nas/proyectos/quijote2/tfgi/list/W63_cygnus/W63_scans_for_destriper_05_17_times_RMS_allbtod2_das5_'+txt+'.txt' ] )
            horns = np.array([0,2,4])
            istgi = np.array([True,True,False])
            mapa, nhits, mapaQ, mapaU, nhitspol, listpix = mambrino2_tfgi(fflist, path=path,tail=tail, nside=nside, nest=nest, nhorns=7, usewei=usewei,
                                                                      horns=horns,istgi=istgi, dobaserm=dobaserm)
            write_mambrino2_tfgi_maps(mapa, nhits, mapaQ, mapaU, nhitspol, listpix, ffout=pathout+'mambrino_w63_'+txt+'_select_feb2025_wei_baserm_angle_final.fits')
        

    # 19/Feb/2025. CRAB maps, with corrected angle before and after Apr2022 
    dothis = 0
    if dothis !=0:
        fflist = np.array(['/net/calp-nas/proyectos/quijote2/tfgi/list/crab/crab_select_allbtod_das24.txt',
                           '/net/calp-nas/proyectos/quijote2/tfgi/list/crab/crab_select_allbtod_das21.txt',
                           '/net/calp-nas/proyectos/quijote2/tfgi/list/crab/crab_select_allbtod_das5.txt'] )
        horns = np.array([0,2,4])
        istgi = np.array([True,True,False])
        mapa, nhits, mapaQ, mapaU, nhitspol, listpix = mambrino2_tfgi(fflist, path=path,tail=tail, nside=nside, nest=nest, nhorns=7, usewei=usewei,
                                                                      horns=horns,istgi=istgi, dobaserm=dobaserm, usezeroangle=False)
        write_mambrino2_tfgi_maps(mapa, nhits, mapaQ, mapaU, nhitspol, listpix, ffout=pathout+'mambrino_crab_select_feb2025_wei_baserm_angle_final.fits')

    dothis = 0
    if dothis !=0:
        fflist = np.array(['/net/calp-nas/proyectos/quijote2/tfgi/list/crab/crab_select_allbtod_das24.txt',
                           '/net/calp-nas/proyectos/quijote2/tfgi/list/crab/crab_select_allbtod_das21.txt',
                           '/net/calp-nas/proyectos/quijote2/tfgi/list/crab/crab_select_allbtod_das5.txt'] )
        horns = np.array([0,2,4])
        istgi = np.array([True,True,False])
        mapa, nhits, mapaQ, mapaU, nhitspol, listpix = mambrino2_tfgi(fflist, path=path,tail=tail, nside=nside, nest=nest, nhorns=7, usewei=usewei,
                                                                      horns=horns,istgi=istgi, dobaserm=dobaserm, usezeroangle=False, usepoleff=True)
        write_mambrino2_tfgi_maps(mapa, nhits, mapaQ, mapaU, nhitspol, listpix, ffout=pathout+'mambrino_crab_select_feb2025_wei_baserm_poleff_angle_final.fits')


    # 19/Feb/2025. W44 maps with calibrated angle and poleff applied
    dothis = 0
    if dothis !=0:
        fflist = np.array(['/net/calp-nas/proyectos/quijote2/tfgi/list/W44/W44_scans_for_destriper_05_17_times_RMS_allbtod2_das24.lst',
                           '/net/calp-nas/proyectos/quijote2/tfgi/list/W44/W44_scans_for_destriper_05_17_times_RMS_allbtod2_das21.lst',
                           '/net/calp-nas/proyectos/quijote2/tfgi/list/W44/W44_scans_for_destriper_05_17_times_RMS_allbtod2_das5.lst' ] )
        horns = np.array([0,2,4])
        istgi = np.array([True,True,False])
        mapa, nhits, mapaQ, mapaU, nhitspol, listpix = mambrino2_tfgi(fflist, path=path,tail=tail, nside=nside, nest=nest, nhorns=7, usewei=usewei,
                                                                      horns=horns,istgi=istgi, dobaserm=dobaserm, usezeroangle=True, usepoleff=False)
        write_mambrino2_tfgi_maps(mapa, nhits, mapaQ, mapaU, nhitspol, listpix, ffout=pathout+'mambrino_w44_select_feb2025_wei_baserm_poleff_angle_final.fits')

    # 24/Mar/2025. Single TauA observation (Fasano et al. paper). CRAB-211130-0230
    dothis = 0
    if dothis !=0:
        path  = '/net/calp-nas/proyectos/quijote2/tfgi/btod/test/'
        tail     = '.btod2'
        nside    = 512
        nest     = False
        usewei   = True
        dobaserm = True

        ff = 'onecrab.txt'
        #np.savetxt(ff,['CRAB-211130-0230-000','CRAB-211130-0230-001','CRAB-211130-0230-002','CRAB-211130-0230-003','CRAB-211130-0230-004'],fmt='%s')
        fflist = np.array([ff]*6)
        horns = np.array([0,2,3,4,5,6])
        istgi = np.array([True,True,False,False,False,True])
        mapa, nhits, mapaQ, mapaU, nhitspol, listpix = mambrino2_tfgi(fflist, path=path,tail=tail, nside=nside, nest=nest, nhorns=7, usewei=usewei,
                                                                      horns=horns,istgi=istgi, dobaserm=dobaserm, usezeroangle=True, usepoleff=False)
        write_mambrino2_tfgi_maps(mapa, nhits, mapaQ, mapaU, nhitspol, listpix, ffout=pathout+'mambrino_onecrab_nside512_btod_wei_baserm_zeroangle_nopoleff.fits')


# 8/Apr/2025. W44 maps with calibrated angle and poleff applied, per month
    dothis = 0
    if dothis !=0:
        tails = ['','_2112_2202','_2204','_2205','_2204_2205']
        horns = np.array([0,2,4])
        istgi = np.array([True,True,False])
        pathlist = '/net/calp-nas/proyectos/quijote2/tfgi/list/W44/'
        for xx in tails:
            fflist = np.array([pathlist+'W44_scans_for_destriper_05_17_times_RMS_allbtod2_das24'+xx+'.lst',
                               pathlist+'W44_scans_for_destriper_05_17_times_RMS_allbtod2_das21'+xx+'.lst',
                               pathlist+'W44_scans_for_destriper_05_17_times_RMS_allbtod2_das5'+xx+'.lst' ] )

            mapa, nhits, mapaQ, mapaU, nhitspol, listpix = mambrino2_tfgi(fflist, path=path,tail=tail, nside=nside, nest=nest, nhorns=7, usewei=usewei,
                                                                      horns=horns,istgi=istgi, dobaserm=dobaserm, usezeroangle=False, usepoleff=True)
            write_mambrino2_tfgi_maps(mapa, nhits, mapaQ, mapaU, nhitspol, listpix, ffout=pathout+'mambrino_w44_select_apr2025_wei_baserm_poleff_angle_final'+xx+'.fits')


# 08/Apr/2025. CRAB maps, before and after Apr2022. But selecting phase states
    dothis = 0
    if dothis !=0:
        fflist = np.array(['/net/calp-nas/proyectos/quijote2/tfgi/list/crab/crab_select_allbtod_das24_beforeApr2022.txt',
                           '/net/calp-nas/proyectos/quijote2/tfgi/list/crab/crab_select_allbtod_das21_beforeApr2022.txt',
                           '/net/calp-nas/proyectos/quijote2/tfgi/list/crab/crab_select_allbtod_das5_beforeApr2022.txt'] )
        horns = np.array([0,2,4])
        istgi = np.array([True,True,False])
        mapa, nhits, mapaQ, mapaU, nhitspol, listpix = mambrino2_tfgi(fflist, path=path,tail=tail, nside=nside, nest=nest, nhorns=7, usewei=usewei,
                                                                      horns=horns,istgi=istgi, dobaserm=dobaserm, usezeroangle=True, usepoleff=False,usephasestates=1)
        write_mambrino2_tfgi_maps(mapa, nhits, mapaQ, mapaU, nhitspol, listpix, ffout=pathout+'mambrino_crab_select_apr2025_wei_baserm_beforeApr2022_phasest1.fits')

        mapa, nhits, mapaQ, mapaU, nhitspol, listpix = mambrino2_tfgi(fflist, path=path,tail=tail, nside=nside, nest=nest, nhorns=7, usewei=usewei,
                                                                      horns=horns,istgi=istgi, dobaserm=dobaserm, usezeroangle=True, usepoleff=False,usephasestates=2)
        write_mambrino2_tfgi_maps(mapa, nhits, mapaQ, mapaU, nhitspol, listpix, ffout=pathout+'mambrino_crab_select_apr2025_wei_baserm_beforeApr2022_phasest2.fits')

    if dothis !=0:
        fflist = np.array(['/net/calp-nas/proyectos/quijote2/tfgi/list/crab/crab_select_allbtod_das24_afterApr2022.txt',
                           '/net/calp-nas/proyectos/quijote2/tfgi/list/crab/crab_select_allbtod_das21_afterApr2022.txt',
                           '/net/calp-nas/proyectos/quijote2/tfgi/list/crab/crab_select_allbtod_das5_afterApr2022.txt'] )
        horns = np.array([0,2,4])
        istgi = np.array([True,True,False])
        mapa, nhits, mapaQ, mapaU, nhitspol, listpix = mambrino2_tfgi(fflist, path=path,tail=tail, nside=nside, nest=nest, nhorns=7, usewei=usewei,
                                                                      horns=horns,istgi=istgi, dobaserm=dobaserm, usezeroangle=True, usepoleff=False, usephasestates=1)
        write_mambrino2_tfgi_maps(mapa, nhits, mapaQ, mapaU, nhitspol, listpix, ffout=pathout+'mambrino_crab_select_apr2025_wei_baserm_afterApr2022_phasest1.fits')

        mapa, nhits, mapaQ, mapaU, nhitspol, listpix = mambrino2_tfgi(fflist, path=path,tail=tail, nside=nside, nest=nest, nhorns=7, usewei=usewei,
                                                                      horns=horns,istgi=istgi, dobaserm=dobaserm, usezeroangle=True, usepoleff=False, usephasestates=2)
        write_mambrino2_tfgi_maps(mapa, nhits, mapaQ, mapaU, nhitspol, listpix, ffout=pathout+'mambrino_crab_select_apr2025_wei_baserm_afterApr2022_phasest2.fits')

    # 11/Apr/2025. W44 maps with calibrated angle and poleff applied, per month
    dothis = 0
    if dothis !=0:
        tails = ['','_2112_2202','_2204','_2205','_2204_2205']
        horns = np.array([0,2,4])
        istgi = np.array([True,True,False])
        pathlist = '/net/calp-nas/proyectos/quijote2/tfgi/list/W44/'
        for xx in tails:
            fflist = np.array([pathlist+'W44_scans_for_destriper_05_17_times_RMS_allbtod2_das24'+xx+'.lst',
                               pathlist+'W44_scans_for_destriper_05_17_times_RMS_allbtod2_das21'+xx+'.lst',
                               pathlist+'W44_scans_for_destriper_05_17_times_RMS_allbtod2_das5'+xx+'.lst' ] )

            mapa, nhits, mapaQ, mapaU, nhitspol, listpix = mambrino2_tfgi(fflist, path=path,tail=tail, nside=nside, nest=nest, nhorns=7, usewei=usewei,
                                                                      horns=horns,istgi=istgi, dobaserm=dobaserm, usezeroangle=False, usepoleff=True)
            write_mambrino2_tfgi_maps(mapa, nhits, mapaQ, mapaU, nhitspol, listpix, ffout=pathout+'mambrino_w44_select_apr2025_wei_baserm_poleff12_angle_final'+xx+'.fits')


    # 14/Apr/2025. CRAB maps, before and after Apr2022, but with poleff per ph state included
    dothis = 0
    if dothis !=0:
        fflist = np.array(['/net/calp-nas/proyectos/quijote2/tfgi/list/crab/crab_select_allbtod_das24_beforeApr2022.txt',
                           '/net/calp-nas/proyectos/quijote2/tfgi/list/crab/crab_select_allbtod_das21_beforeApr2022.txt',
                           '/net/calp-nas/proyectos/quijote2/tfgi/list/crab/crab_select_allbtod_das5_beforeApr2022.txt'] )
        horns = np.array([0,2,4])
        istgi = np.array([True,True,False])
        mapa, nhits, mapaQ, mapaU, nhitspol, listpix = mambrino2_tfgi(fflist, path=path,tail=tail, nside=nside, nest=nest, nhorns=7, usewei=usewei,
                                                                      horns=horns,istgi=istgi, dobaserm=dobaserm, usezeroangle=True, usepoleff=True)
        write_mambrino2_tfgi_maps(mapa, nhits, mapaQ, mapaU, nhitspol, listpix, ffout=pathout+'mambrino_crab_select_apr2025_zeroang_poleff12_wei_baserm_beforeApr2022.fits')

    if dothis !=0:
        fflist = np.array(['/net/calp-nas/proyectos/quijote2/tfgi/list/crab/crab_select_allbtod_das24_afterApr2022.txt',
                           '/net/calp-nas/proyectos/quijote2/tfgi/list/crab/crab_select_allbtod_das21_afterApr2022.txt',
                           '/net/calp-nas/proyectos/quijote2/tfgi/list/crab/crab_select_allbtod_das5_afterApr2022.txt'] )
        horns = np.array([0,2,4])
        istgi = np.array([True,True,False])
        mapa, nhits, mapaQ, mapaU, nhitspol, listpix = mambrino2_tfgi(fflist, path=path,tail=tail, nside=nside, nest=nest, nhorns=7, usewei=usewei,
                                                                      horns=horns,istgi=istgi, dobaserm=dobaserm, usezeroangle=True, usepoleff=True)
        write_mambrino2_tfgi_maps(mapa, nhits, mapaQ, mapaU, nhitspol, listpix, ffout=pathout+'mambrino_crab_select_apr2025_zeroang_poleff12_wei_baserm_afterApr2022.fits')


    # 14/Apr/2025. W44 maps with calibrated angle and poleff applied, per month
    dothis = 0
    if dothis !=0:
        tails = ['','_2112_2202','_2204','_2205','_2204_2205']
        horns = np.array([0,2,4])
        istgi = np.array([True,True,False])
        pathlist = '/net/calp-nas/proyectos/quijote2/tfgi/list/W44/'
        for xx in tails:
            fflist = np.array([pathlist+'W44_scans_for_destriper_05_17_times_RMS_allbtod2_das24'+xx+'.lst',
                               pathlist+'W44_scans_for_destriper_05_17_times_RMS_allbtod2_das21'+xx+'.lst',
                               pathlist+'W44_scans_for_destriper_05_17_times_RMS_allbtod2_das5'+xx+'.lst' ] )

            mapa, nhits, mapaQ, mapaU, nhitspol, listpix = mambrino2_tfgi(fflist, path=path,tail=tail, nside=nside, nest=nest, nhorns=7, usewei=usewei,
                                                                      horns=horns,istgi=istgi, dobaserm=dobaserm, usezeroangle=False, usepoleff=True)
            write_mambrino2_tfgi_maps(mapa, nhits, mapaQ, mapaU, nhitspol, listpix, ffout=pathout+'mambrino_w44_select_14apr2025_wei_baserm_poleff12_angle_final'+xx+'.fits')



   # 3/Jun/2025. Moon maps with no calibrated angle and no poleff applied
    dothis = 0
    if dothis !=0:
        path  = '/net/calp-nas/proyectos/quijote2/tfgi/btod/'
        tail  = '.btod2'
        tails = ['211222','211223','211224','211227','211228','220123','220917']
        mfiles = ['moon_'+s+'.txt' for s in tails]
        pathlist = '/net/calp-nas/proyectos/quijote2/tfgi/list/moon/'
        horns = np.array([0,2,3,4,5,6])
        istgi = np.array([True,True,False,False,False,True])
        usewei = False
        
        for xx in tails:
            fflist = np.array([pathlist+'moon_'+xx+'.txt']*6)

            mapa, nhits, mapaQ, mapaU, nhitspol, listpix = mambrino2_tfgi(fflist, path=path,tail=tail, nside=nside, nest=nest, nhorns=7, usewei=usewei,
                                                                      horns=horns,istgi=istgi, dobaserm=dobaserm, usezeroangle=True, usepoleff=False, useflag=False)
            write_mambrino2_tfgi_maps(mapa, nhits, mapaQ, mapaU, nhitspol, listpix, ffout=pathout+'mambrino_moon_3jun2025_noflag_baserm_nopoleff_noangle_'+xx+'.fits')



    # 4/June/2025. W44 maps with calibrated angle PER MONTH and poleff applied.
    dothis = 0
    if dothis !=0:
        tails = ['','_half1','_half2'] #,'_2112_2202','_2204','_2205','_2204_2205']
        horns = np.array([0,2,4])
        istgi = np.array([True,True,False])
        pathlist = '/net/calp-nas/proyectos/quijote2/tfgi/list/W44/'
        for xx in tails:
            fflist = np.array([pathlist+'W44_scans_for_destriper_05_17_times_RMS_allbtod2_das24'+xx+'.txt',
                               pathlist+'W44_scans_for_destriper_05_17_times_RMS_allbtod2_das21'+xx+'.txt',
                               pathlist+'W44_scans_for_destriper_05_17_times_RMS_allbtod2_das5'+xx+'.txt' ] )

            mapa, nhits, mapaQ, mapaU, nhitspol, listpix = mambrino2_tfgi(fflist, path=path,tail=tail, nside=nside, nest=nest, nhorns=7, usewei=usewei,
                                                                      horns=horns,istgi=istgi, dobaserm=dobaserm, usezeroangle=False, usepoleff=True)
            write_mambrino2_tfgi_maps(mapa, nhits, mapaQ, mapaU, nhitspol, listpix, ffout=pathout+'mambrino_w44_select_4june2025_wei_baserm_poleff_monthly_angle'+xx+'.fits')


   # 4/July/2025. W44 maps with calibrated angle PER MONTH and poleff applied. Using orig_reso.
    dothis = 0
    if dothis !=0:
        tails = ['','_half1','_half2'] 
        horns = np.array([0,2,4])
        istgi = np.array([True,True,False])
        pathlist = '/net/calp-nas/proyectos/quijote2/tfgi/list/W44/'
        for xx in tails:
            fflist = np.array([pathlist+'W44_scans_for_destriper_05_17_times_RMS_allbtod2_das24'+xx+'.txt',
                               pathlist+'W44_scans_for_destriper_05_17_times_RMS_allbtod2_das21'+xx+'.txt',
                               pathlist+'W44_scans_for_destriper_05_17_times_RMS_allbtod2_das5'+xx+'.txt' ] )

            mapa, nhits, mapaQ, mapaU, nhitspol, listpix = mambrino2_tfgi(fflist, path=path,tail=tail, nside=nside, nest=nest, nhorns=7, usewei=usewei,
                                                                      horns=horns,istgi=istgi, dobaserm=dobaserm, usezeroangle=False, usepoleff=True)
            write_mambrino2_tfgi_maps(mapa, nhits, mapaQ, mapaU, nhitspol, listpix, ffout=pathout+'mambrino_w44_select_4july2025_wei_baserm_poleff_monthly_angle_orig_reso'+xx+'.fits')

   # 6/July/2026. repeat W44 maps with calibrated angle PER MONTH and poleff applied. Using orig_reso. But now generates err maps
    dothis = 1
    if dothis !=0:
        tails = ['','_half1','_half2'] 
        horns = np.array([0,2,4])
        istgi = np.array([True,True,False])
        pathlist = '/net/calp-nas/proyectos/quijote2/tfgi/list/W44/'
        for xx in tails:
            fflist = np.array([pathlist+'W44_scans_for_destriper_05_17_times_RMS_allbtod2_das24'+xx+'.txt',
                               pathlist+'W44_scans_for_destriper_05_17_times_RMS_allbtod2_das21'+xx+'.txt',
                               pathlist+'W44_scans_for_destriper_05_17_times_RMS_allbtod2_das5'+xx+'.txt' ] )

            mapa, nhits, mapaQ, mapaU, nhitspol, listpix, errI, errQ, errU, covQU = mambrino2_tfgi(fflist, path=path,tail=tail, nside=nside, nest=nest, nhorns=7, usewei=usewei,
                                                                      horns=horns,istgi=istgi, dobaserm=dobaserm, usezeroangle=False, usepoleff=True)
            write_mambrino2_tfgi_maps(mapa, nhits, mapaQ, mapaU, nhitspol, listpix, errI, errQ, errU, covQU, ffout=pathout+'mambrino_w44_select_6july2026_wei_baserm_poleff_monthly_angle_withcov_orig_reso'+xx+'.fits')




    sys.exit('Code ends.')

