#!/usr/bin/env python

"""
2/Sep/2025

@author: jalberto

Averaging raw (PCAP) files and generating .RAW files for MFI2 FPGA.
Also includes the analysis of the CAL data.

HISTORY:
*  2/09/2025 - original version. JARM
* 22/09/2025 - treating CAL data. 

"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from sancho_mfi2_fpga_io import *
import sys

# Global variables
fpga_nu_samp_MHz = 4900.0
n_fpga           = 512
n_avg_fpga       = 8192
dnu_MHz          = fpga_nu_samp_MHz/n_fpga                   # 9.57 MHz
dt_sec           = n_fpga*n_avg_fpga / fpga_nu_samp_MHz*1e-6 # 856 us
nstokes          = 8
nsbin            = 35               # choice: bin every 35 samples, equivalent to ~30ms
msbin            = nsbin*dt_sec*1e3 # bin size in ms

# Averaging code
def sancho_bin_data_mfi2_fpga(alldata, verbose=False):

    print(f" SANCHO_BIN_DATA_MFI2_FPGA: averaging in bins of {nsbin} samples, equivalent to {msbin} ms.")

    nfreq     = alldata[0].shape[1]
    nsamp_all = alldata[0].shape[0]
    n   = (nsamp_all // nsbin) * nsbin
    nb  = n//nsbin 

    if verbose:
        print(f'   > Averaging in bins of {nsbin} samples, equivalent to {msbin} ms.')
        print(f'   > Number of Stokes  = {len(alldata)}')
        print(f'   > Number of samples = {nsamp_all}')
        print(f'   > Number of freqs   = {nfreq}')
        print(f'   > Size of bins      = {nsbin}')
        print(f'   > Number of binned samples = {nb}')
        print(f'   > Trimmed size of samples {nsbin} x {nb} = {n}')
        print(f'   > Extra samples = {nsamp_all-n}')

    data = np.zeros(nstokes*nb*nfreq).reshape((nstokes,nb,nfreq))
    wei  = np.copy(data)
    jd   = np.zeros(nb)
    nu   = np.zeros(nfreq)

    # Full time and freq vectors
    t_all = np.arange(0, nsamp_all, 1, dtype=np.float64) * dt_sec # in seconds
    nuMHz = np.arange(0, nfreq, 1, dtype=np.float64) * dnu_MHz    # in MHz

    # Averaged time vectors
    arr_trimmed = t_all[:n].astype(np.float64)
    arr_grouped = arr_trimmed.reshape(nb, nsbin)
    jd[:]       = arr_grouped.mean(axis=1)

    # Averaged data
    for j in range(nstokes):
        x = alldata[j]
        arr_trimmed = x[:n].astype(float) #make sure it is a float
        arr_grouped = arr_trimmed.reshape(nb, nsbin, nfreq)
        
        arr_avg = arr_grouped.mean(axis=1)
        data[j,:,:] = arr_avg[:,:]

        arr_std = arr_grouped.std(axis=1)
        arr_wei = np.zeros_like(arr_std, dtype=float)
        np.divide(1.0, arr_std**2, out=arr_wei, where=arr_std != 0)
        wei[j,:,:] = arr_wei[:,:]
        #wei[j,:,:]  = 1.0/arr_std[:,:]**2
        #wei[j,:,:]  = np.where(arr_std == 0, 0, 1.0 / arr_std**2)

    print(" SANCHO_BIN_DATA_MFI2_FPGA: end") 
    return jd, data, wei, nuMHz, t_all

# CAL diode lockin code
def sancho_cal_lockin_mfi2_fpga(cal_mask, cal_sgn, alldata, verbose=False):

    print(" SANCHO_CAL_LOCKIN_MFI2_FPGA.")
    nsamp = 1280 # Number of samples in a CAL cycle
    f_cal = 20   # patterns (repetitions) per cal cycle 
    nfreq = alldata[0].shape[1]
    for j in range(nstokes):
        cal2d = cal_mask[j]
        cal   = cal2d[:,0] # will be the same for all frequencies.
        copia_cal = cal.astype(np.int16).copy()
        copia_cal[-1] = 0
        steps = cal - np.roll(copia_cal, 1)
        on     = np.where(steps == 1)[0]
        ncals1 = len(on)
        off = np.where(steps == -1)[0]
        ncals2 = len(off)
        ncals = min(ncals1, ncals2)

        goahead = 1 
        if ncals==0:
            goahead = 0
            ncals   = 1
        
        if j==0:
            base  = np.zeros(ncals*nfreq*nstokes).reshape((nstokes,nfreq,ncals))
            acti  = np.zeros(ncals*nfreq*nstokes).reshape((nstokes,nfreq,ncals))
            gain  = np.zeros(ncals*nfreq*nstokes).reshape((nstokes,nfreq,ncals))
            gain3 = np.zeros(ncals*nfreq*nstokes).reshape((nstokes,nfreq,ncals))

        # Lockin pattern
        # 1) Default scheme ON-OFF
        nlock_one  = int(np.rint(nsamp / f_cal))
        if not np.isclose(nlock_one * f_cal, nsamp):
            raise ValueError("nsamp / f_cal is not an integer")
        lockin_one = np.zeros(nlock_one)
        nh         = nlock_one//2
        lockin_one[0:nh] = -1
        lockin_one[nh:nlock_one] = 1

        nlocks  = nsamp//nlock_one
        lockin  = np.tile(lockin_one, nlocks)
        nominal = np.where(lockin > 0, lockin, 0) 

        # 2) Symmetric scheme: OFF-ON-OFF used in MFI.
        lockin3 = np.tile(lockin_one, nlocks-2)

        lockin_one_first        = np.zeros(nlock_one)
        lockin_one_first[0:nh]  =  -0.5
        lockin_one_first[nh:nlock_one] = 1.0

        lockin_one_last         = np.zeros(nlock_one)
        lockin_one_last[0:nh]   =  -0.5
        lockin_one_last[nh:nlock_one]  = 0.0

        lockin3 = np.concatenate([lockin_one_first, lockin3, lockin_one_last])

        if (j==0) & verbose:
            print(f' Number of samples in a CAL cycle = {nsamp}')
            print(f' Number of repetitions of the cal diode = {f_cal}')
            print(f' CAL cycle has {f_cal} repetitions of {nlock_one} samples.')
            print(' Number of CALs in this file = ',ncals)
            print('   > ON vector  = ',on)
            print('   > OFF vector = ',off)
            print('i  i1  i2  ni error')
        for i in range(ncals):
            i1 = on[i]
            i2 = off[i]
            ni = i2-i1
            im = (i1+i2)//2
            delta_error = np.sum(cal_sgn[j][i1:i2,0]-nominal) 
            if not np.isclose(delta_error, 0):
                raise ValueError("The CAL signal template is not matching the expectation.") 
            if (j==0) & verbose: 
                print(i,i1,i2,ni, delta_error)

            altura = np.max(lockin) - np.min(lockin)
            norm   = np.sum(lockin**2) / altura
            norm3  = np.sum(lockin3[lockin3 > 0.0]**2)  

            data_cal = alldata[j][i1:i2,:]
            #data_cal = cal_sgn[j][i1:i2,:]
            # Gain = amplitude of the cal signal
            for k in range(nfreq):
                gain[j,k,i]  = np.sum(data_cal[:,k]*lockin)  / norm
                gain3[j,k,i] = np.sum(data_cal[:,k]*lockin3) / norm3
            # Base and activation levels
                yy    = data_cal[:,k]
                mask1 = (lockin < 0.0) & (yy != 0.0)
                ncut1 = np.count_nonzero(mask1)
                if ncut1>1:
                    base[j,k,i]  = np.mean( yy[mask1] )
                mask2 = (lockin > 0.0) & (yy != 0.0)
                ncut2 = np.count_nonzero(mask2)
                if ncut2>1:
                    acti[j,k,i] = np.mean( yy[mask2] )
        if (j==0) & verbose:
            print('gain  = ',gain[0,0,:])
            print('gain3 = ',gain3[0,0,:])
            print('base  = ',base[0,0,:])
            print('acti  = ',acti[0,0,:])

    return gain3, lockin


# Example codes
def test_sept2025():

    # Example data generated on Sept 1st, 2025.
    file_directory = '/Users/jalberto/quijote/RIM/example_data/'
    file_name = 'test.pcap'

    alldata, allhdr = read_mfi2_fpga_pcap_sept2025(file_name, dir=file_directory, verbose=True)

    # Display basic info
    for i, arr in enumerate(alldata):
        print(f"Values in Array {i}: shape={arr.shape}, dtype={arr.dtype}")

    for i, arr in enumerate(allhdr):
        print(f"Headers in Array {i}: shape={arr.shape}, dtype={arr.dtype}")

    # Bin data
    jd, data, wei, nuMHz, t_all = sancho_bin_data_mfi2_fpga(alldata)

    # Diode analysis. CAL lockin.
    gain, lockin = sancho_cal_lockin_mfi2_fpga(allhdr, allhdr, alldata) 

    # Output dictionary
    raw  = {'JD':jd, 'DATA':data, 'WEI':wei, 'NSBIN':nsbin, 'NUMHZ':nuMHz }

    # Write file:
    write_mfi2_fpga_raw_data(raw, '/Users/jalberto/quijote/RIM/example_data/test.fits', overwrite=True)


def test_dec2025():

    # Example data generated on December 2025.
    file_directory = '/Users/jalberto/quijote/RIM/example_data/'
    file_name = 'test1'

    alldata, cal_sgn, cal_mask, dindex = read_mfi2_fpga_pcap(file_name, dir=file_directory, verbose=True)
    check_dindex(dindex,nstokes) # checking dindex values

    # Display basic info
    for i, arr in enumerate(alldata):
        print(f"Values in Array {i}: shape={arr.shape}, dtype={arr.dtype}")

    # Binning
    jd, data, wei, nuMHz, t_all = sancho_bin_data_mfi2_fpga(alldata, verbose=True)
    
    # Write binned file
    raw  = {'JD':jd, 'DATA':data, 'WEI':wei, 'NSBIN':nsbin, 'NUMHZ':nuMHz }
    write_mfi2_fpga_raw_data(raw, '/Users/jalberto/quijote/RIM/example_data/test1.fits', overwrite=True)

    # CAL signal
    gain0, lockin0 = sancho_cal_lockin_mfi2_fpga(cal_mask, cal_sgn, cal_sgn, verbose=True) # Returns gain=1
    gain, lockin = sancho_cal_lockin_mfi2_fpga(cal_mask, cal_sgn, alldata, verbose=False)
    print(cal_mask[0][0:1280,0])
    print(cal_sgn[0][0:1280,0])
    print(gain[0,0,:])
    plt.plot(cal_sgn[0][0:1280,0])
    plt.plot(lockin,'r.')
    plt.show()

    # Plots
    doplot = True
    if doplot:
    # Plot1
        stokes  = 5
        channel = 0
        a = alldata[stokes][:,channel]
        b = cal_sgn[stokes][:,channel]
        c = cal_mask[stokes][:,channel]
        d = dindex[stokes][:,channel]

        plt.plot(t_all, a,'b')
        plt.plot(jd, data[stokes][:,channel],'r')
        sigma = 1/np.sqrt(wei[stokes][:,channel])
        plt.plot(jd, data[stokes][:,channel]-sigma,'g')
        plt.plot(jd, data[stokes][:,channel]+sigma,'g')
        #plt.errorbar(jd,data[stokes][:,channel],yerr=1/np.sqrt(wei[stokes][:,channel]),fmt='r*')
        plt.title("Data")
        plt.show()

        plt.plot(b, 'b')
        plt.plot(c, 'r')
        plt.title("Calibration Signal (blue) and Mask (red)", pad=20, loc='center', fontsize=24)
        plt.ylabel("Digital values [a.u.]", labelpad=20)
        plt.xlabel("Sample index [a.u.] @ sample time = 856us", labelpad=14)
        plt.show()

        plt.title("Samples counter ", pad=20, loc='center', fontsize=24)
        plt.ylabel("Counts [a.u.]", labelpad=20)
        plt.xlabel("Sample index [a.u.] @ sample time = 856us", labelpad=14)
        plt.plot(d, 'b')
        plt.show()
    
        #Plot 2. Espectrograma
        stokes = {
        "i_lowband":  0,
        "q_lowband":  1,
        "u_lowband":  2,
        "v_lowband":  3,
        "i_highband": 4,
        "q_highband": 5,
        "u_highband": 6,
        "v_highband": 7 }
        stokes_param_label = "i_highband"
        stokes_param_to_plot = stokes[stokes_param_label]

        a = []
        for i in range(0, len(alldata[stokes_param_to_plot])):
            a.append(alldata[stokes_param_to_plot][i])

        A = np.vstack([np.asarray(row, dtype=float) for row in a])
        A = A.T
        A_db = 10 * np.log10(A + 1e-12)
        plt.figure(figsize=(10, 5))
        plt.imshow(A_db, aspect='auto', origin='lower', cmap='viridis')
        plt.colorbar(label="Amplitude [dB]")
        plt.xlabel("spectrum index [a.u.] @ spectrum time = 856us", labelpad=14)
        plt.ylabel("Channel index [a.u.] @ channel freq = index * 9.57MHz", labelpad=14)
        plt.title("Averaged-" + stokes_param_label + " Spectrogram", pad=20, loc='center', fontsize=24)
        plt.show() 

    return

def test_jan2026():
    # Example data generated on Jan 2026, during testing phase of the calibration
    file_directory = '/Users/jalberto/quijote/RIM/example_data/'
    file_name = 'sweep.pcap'

    alldata, cal_sgn, cal_mask, dindex = read_mfi2_fpga_pcap(file_name, dir=file_directory, verbose=True)
    print('alldata = ',len(alldata))
    print('dindex = ',len(dindex))
    print(dindex[0].shape)
    print(dindex)
    #check_dindex(dindex,nstokes) # checking dindex values

    # Plot1
    stokes  = 0
    channel = 0
    a = alldata[stokes][:,channel]
    b = cal_sgn[stokes][:,channel]
    c = cal_mask[stokes][:,channel]
    d = dindex[stokes][:,channel]

    plt.plot(a,'b')
    plt.title("Data")
    plt.show()


    return


##############
# MAIN code
if __name__ == "__main__":

    # Datos de Sept 2025. Primer formato
    dothis = False
    if dothis:
        test_sept2025()

    # Datos Dic 2025. Nuevo formato
    dothis = False
    if dothis:
        test_dec2025()

    # Datos Enero 2026. Nuevo formato
    dothis = True
    if dothis:
        test_jan2026()