#!/usr/bin/env python

"""
2/Sep/2025

@author: jalberto

Averaging raw (PCAP) files and generating .RAW files for MFI2 FPGA.
Also includes the analysis of the CAL data.

HISTORY:
*  2/09/2025 - original version. JARM
* 22/09/2025 - treating CAL data. 
* 17/09/2026 - new version of the code, dealing with leftovers. Join binning and CAL extraction in a single routine.

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

# New versions of binning, accounting for leftovers
def sancho_bin_data_mfi2_fpga_leftovers(alldata, leftover=None, t_offset_sec=0.0, verbose=False):
    """
    Average data into bins of nsbin samples, propagating leftover samples between files.

    Parameters
    ----------
    alldata      : list of nstokes arrays of shape (nsamp, nfreq)
    leftover     : leftover from previous file, as {'data': [...], 't0': float}, or None
    t_offset_sec : absolute time (in seconds) of the first sample in alldata
    verbose      : if True, print diagnostic information

    Returns
    -------
    jd           : array of shape (nb,), mean absolute time of each bin (seconds)
    data         : array of shape (nstokes, nb, nfreq), binned averages
    wei          : array of shape (nstokes, nb, nfreq), weights (1/sigma^2)
    nuMHz        : array of shape (nfreq,), frequency axis in MHz
    leftover_new : dict with leftover samples and their absolute start time
    """

    # 1. Prepend leftover samples from the previous file
    if leftover is not None:
        merged = [np.concatenate([leftover['data'][j], alldata[j]], axis=0)
                  for j in range(nstokes)]
        t0 = leftover['t0']       # leftover starts at this absolute time
    else:
        merged = [alldata[j].copy() for j in range(nstokes)]
        t0 = t_offset_sec

    nfreq     = merged[0].shape[1]
    nsamp_all = merged[0].shape[0]
    n         = (nsamp_all // nsbin) * nsbin  # largest multiple of nsbin that fits
    nb        = n // nsbin

    if verbose:
        print(f'   > Averaging in bins of {nsbin} samples, equivalent to {msbin} ms.')
        print(f'   > Number of Stokes  = {len(alldata)}')
        print(f'   > Number of samples = {nsamp_all}')
        print(f'   > Number of freqs   = {nfreq}')
        print(f'   > Size of bins      = {nsbin}')
        print(f'   > Number of binned samples = {nb}')
        print(f'   > Trimmed size of samples {nsbin} x {nb} = {n}')
        print(f'   > Extra samples = {nsamp_all-n}')
        print(f'   > nsamp_all={nsamp_all}, n={n}, nb={nb}, leftover={nsamp_all - n}')

    # 2. Absolute time axis for all samples in this block
    t_all = t0 + np.arange(nsamp_all, dtype=np.float64) * dt_sec

    # 3. Mean time of each bin
    jd = t_all[:n].reshape(nb, nsbin).mean(axis=1)

    # 4. Frequency axis
    nuMHz = np.arange(nfreq, dtype=np.float64) * dnu_MHz

    # 5. Average data and compute weights
    data = np.zeros((nstokes, nb, nfreq))
    wei  = np.zeros_like(data)
    for j in range(nstokes):
        x   = merged[j][:n].astype(float)
        arr = x.reshape(nb, nsbin, nfreq)
        data[j] = arr.mean(axis=1)
        std = arr.std(axis=1)
        np.divide(1.0, std**2, out=wei[j], where=(std != 0))

    # 6. Store leftover samples with their absolute start time
    leftover_new = {
        'data': [merged[j][n:, :] for j in range(nstokes)],
        't0'  : t0 + n * dt_sec    # absolute time of the first leftover sample
    }

    return jd, data, wei, nuMHz, leftover_new


def sancho_bin_data_mfi2_fpga_new_v1(alldata, leftover=None, verbose=False):
    '''
    Docstring for sancho_bin_data_mfi2_fpga_new_v1. Not tested
    
    :param alldata: input data. 
    :param leftover: leftover from previous file
    :param verbose: if active, prints info on screen
    '''

    print(f" SANCHO_BIN_DATA_MFI2_FPGA_NEW: averaging in bins of {nsbin} samples, equivalent to {msbin} ms.")

    # Concatenar resto del fichero anterior
    if leftover is not None:
        for j in range(nstokes):
            alldata[j] = np.concatenate(
                [leftover[j], alldata[j]],
                axis=0
            )

    nfreq     = alldata[0].shape[1]
    nsamp_all = alldata[0].shape[0]

    # número máximo múltiplo de nsbin
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

    # save the leftover samples
    leftover_new = []
    for j in range(nstokes):
        leftover_new.append( alldata[j][n:, :])

    # Binned arrays
    data = np.zeros((nstokes,nb,nfreq))
    wei  = np.zeros_like(data)
    jd   = np.zeros(nb)
    nu   = np.zeros(nfreq)

    # Full time and freq vectors
    t_all = np.arange(0, nsamp_all, 1, dtype=np.float64) * dt_sec # in seconds
    nuMHz = np.arange(0, nfreq, 1, dtype=np.float64) * dnu_MHz    # in MHz

    # Averaged time vectors
    jd = t_all[:n].reshape(nb,nsbin).mean(axis=1)

    # Averaged data
    for j in range(nstokes):
        x   = alldata[j][:n].astype(float)
        arr = x.reshape(nb,nsbin,nfreq)
        data[j] = arr.mean(axis=1)

        arr_std=np.std(arr,axis=1)
        np.divide( 1.0, arr_std**2, out=wei[j], where=(arr_std!=0))

    print(" SANCHO_BIN_DATA_MFI2_FPGA_NEW: end") 
    return jd, data, wei, nuMHz, t_all, leftover_new

def read_and_bin_mfi2_fpga(file_list, dir='./', verbose=False):
    """
    Read and bin a time-ordered list of PCAP files, preserving leftover
    samples between files so no data is lost at file boundaries.

    Parameters
    ----------
    file_list : list of str, file names in time order (without directory)
    dir       : str, directory where the files are located
    verbose   : if True, print diagnostic information per file

    Returns
    -------
    jd    : array of shape (nb_total,), absolute time of each bin (seconds)
    data  : array of shape (nstokes, nb_total, nfreq), binned averages
    wei   : array of shape (nstokes, nb_total, nfreq), weights (1/sigma^2)
    nuMHz : array of shape (nfreq,), frequency axis in MHz

    Example of use
    --------------
    file_dir  = '/Users/jalberto/quijote/RIM/example_data/'
    file_list = ['obs_001.pcap', 'obs_002.pcap', 'obs_003.pcap']
    jd, data, wei, nuMHz = read_and_bin_mfi2_fpga(file_list, dir=file_dir, verbose=False)

    """

    leftover = None
    t_abs    = 0.0

    jd_list   = []
    data_list = []
    wei_list  = []

    for i, file_name in enumerate(file_list):
        print(f' [{i+1}/{len(file_list)}] Reading {file_name} ...')
        alldata, cal_sgn, cal_mask, dindex = read_mfi2_fpga_pcap(file_name, dir=dir, verbose=verbose)

        jd, data, wei, nuMHz, leftover = sancho_bin_data_mfi2_fpga_leftovers(
            alldata, leftover=leftover, t_offset_sec=t_abs, verbose=verbose
        )

        # Advance absolute time offset for the next file
        t_abs += alldata[0].shape[0] * dt_sec

        jd_list.append(jd)
        data_list.append(data)
        wei_list.append(wei)

    # Concatenate results from all files along the time axis
    jd   = np.concatenate(jd_list)
    data = np.concatenate(data_list, axis=1)
    wei  = np.concatenate(wei_list,  axis=1)

    print(f' Done. Total bins: {len(jd)}, total time: {jd[-1]:.1f} s')
    return jd, data, wei, nuMHz

# New version of lockin, accounting for leftovers
def sancho_cal_lockin_mfi2_fpga_leftovers(cal_mask, cal_sgn, alldata, leftover=None, verbose=False):
    """
    CAL diode lock-in with leftover support for file boundaries.

    Prepends the tail of the previous file (stored in leftover) so that CAL cycles
    split across file boundaries are correctly recovered.

    Parameters
    ----------
    cal_mask : list of nstokes arrays (nsamp, nfreq), calibration mask
    cal_sgn  : list of nstokes arrays (nsamp, nfreq), calibration signal
    alldata  : list of nstokes arrays (nsamp, nfreq), science data
    leftover : dict {'cal_mask', 'cal_sgn', 'alldata'} from previous file, or None
    verbose  : if True, print diagnostic information

    Returns
    -------
    gain3        : array (nstokes, nfreq, ncals), symmetric lock-in gain per cycle
    lockin       : lock-in template used
    leftover_new : dict with tail arrays for the next file call
    """

    nsamp = 1280  # samples in one CAL cycle
    f_cal = 20    # repetitions per cycle
    nfreq = alldata[0].shape[1]

    # 1. Prepend leftover from the previous file
    if leftover is not None:
        merged_mask = [np.concatenate([leftover['cal_mask'][j], cal_mask[j]], axis=0) for j in range(nstokes)]
        merged_sgn  = [np.concatenate([leftover['cal_sgn'][j],  cal_sgn[j]],  axis=0) for j in range(nstokes)]
        merged_data = [np.concatenate([leftover['alldata'][j],  alldata[j]],  axis=0) for j in range(nstokes)]
    else:
        merged_mask = list(cal_mask)
        merged_sgn  = list(cal_sgn)
        merged_data = list(alldata)

    nsamp_merged = merged_mask[0].shape[0]

    # 2. Build lock-in templates (unchanged from original)
    nlock_one = int(np.rint(nsamp / f_cal))
    if not np.isclose(nlock_one * f_cal, nsamp):
        raise ValueError("nsamp / f_cal is not an integer")
    lockin_one       = np.zeros(nlock_one)
    nh               = nlock_one // 2
    lockin_one[:nh]  = -1
    lockin_one[nh:]  =  1
    nlocks           = nsamp // nlock_one
    lockin           = np.tile(lockin_one, nlocks)
    altura           = np.max(lockin) - np.min(lockin)
    norm             = np.sum(lockin**2) / altura

    # Symmetric template (OFF-ON-OFF)
    lockin_one_first       = np.zeros(nlock_one)
    lockin_one_first[:nh]  = -0.5
    lockin_one_first[nh:]  =  1.0
    lockin_one_last        = np.zeros(nlock_one)
    lockin_one_last[:nh]   = -0.5
    lockin_one_last[nh:]   =  0.0
    lockin3 = np.concatenate([lockin_one_first, np.tile(lockin_one, nlocks - 2), lockin_one_last])
    norm3   = np.sum(lockin3[lockin3 > 0.0]**2)

    # 3. Find complete CAL cycles in merged data (using Stokes 0, same for all)
    cal   = merged_mask[0][:, 0].astype(np.int16)
    copia = cal.copy()
    copia[-1] = 0
    steps = cal - np.roll(copia, 1)
    on  = np.where(steps ==  1)[0]
    off = np.where(steps == -1)[0]

    # Only keep ON/OFF pairs where the full cycle fits within the merged data
    pairs = [(i_on, i_on + nsamp) for i_on in on if i_on + nsamp <= nsamp_merged]
    ncals = len(pairs)

    if verbose:
        print(f'   > Merged samples: {nsamp_merged}, complete CAL cycles found: {ncals}')

    # Handle the case with no complete cycles (carry everything as leftover)
    if ncals == 0:
        leftover_new = {'cal_mask': [merged_mask[j] for j in range(nstokes)],
                        'cal_sgn':  [merged_sgn[j]  for j in range(nstokes)],
                        'alldata':  [merged_data[j]  for j in range(nstokes)]}
        return None, lockin, leftover_new

    # 4. Compute lock-in gain (vectorised over nfreq)
    gain3 = np.zeros((nstokes, nfreq, ncals))
    base  = np.zeros((nstokes, nfreq, ncals))
    acti  = np.zeros((nstokes, nfreq, ncals))

    for j in range(nstokes):
        for i, (i1, i2) in enumerate(pairs):
            d = merged_data[j][i1:i2, :].astype(float)  # (nsamp, nfreq)
            gain3[j, :, i] = (d.T @ lockin3) / norm3    # dot product over time axis

            # Base and activation levels
            mask_off = lockin < 0.0
            mask_on  = lockin > 0.0
            d_off = d[mask_off, :]
            d_on  = d[mask_on,  :]
            valid_off = np.any(d_off != 0, axis=0)
            valid_on  = np.any(d_on  != 0, axis=0)
            base[j, valid_off, i] = d_off[:, valid_off].mean(axis=0)
            acti[j, valid_on,  i] = d_on[ :, valid_on ].mean(axis=0)

    # 5. Leftover: everything after the last complete cycle
    i_last = pairs[-1][1]
    leftover_new = {
        'cal_mask': [merged_mask[j][i_last:, :] for j in range(nstokes)],
        'cal_sgn':  [merged_sgn[j][i_last:, :]  for j in range(nstokes)],
        'alldata':  [merged_data[j][i_last:, :]  for j in range(nstokes)],
    }

    if verbose:
        print(f'   > Leftover size: {nsamp_merged - i_last} samples')

    return gain3, lockin, leftover_new

def read_and_cal_mfi2_fpga(file_list, dir='./', verbose=False):
    """
    Read a time-ordered list of PCAP files and compute the CAL lock-in gain
    for all files, correctly handling CAL cycles at file boundaries.

    Parameters
    ----------
    file_list : list of str, file names in time order (without directory)
    dir       : str, directory where the files are located
    verbose   : if True, print diagnostic information per file

    Returns
    -------
    gain3 : array of shape (nstokes, nfreq, ncals_total), lock-in gain per cycle
    lockin : lock-in template used
    """

    leftover_cal = None
    gain3_list   = []

    for i, file_name in enumerate(file_list):
        print(f' [{i+1}/{len(file_list)}] Reading {file_name} ...')
        alldata, cal_sgn, cal_mask, dindex = read_mfi2_fpga_pcap(file_name, dir=dir, verbose=verbose)

        gain3, lockin, leftover_cal = sancho_cal_lockin_mfi2_fpga_leftovers(
            cal_mask, cal_sgn, alldata, leftover=leftover_cal, verbose=verbose
        )

        if gain3 is not None:
            gain3_list.append(gain3)
        else:
            print(f'   > Warning: no complete CAL cycle found in {file_name}')

    if not gain3_list:
        raise RuntimeError('No complete CAL cycles found in any file.')

    # Concatenate along the CAL cycle axis
    gain3_all = np.concatenate(gain3_list, axis=2)  # (nstokes, nfreq, ncals_total)
    print(f' Done. Total CAL cycles: {gain3_all.shape[2]}')

    return gain3_all, lockin

# Main routine that carries out all tasks: read data, binning and CAL evaluation
def sancho_read_bin_cal_mfi2_fpga(file_list, dir='./', verbose=False):
    """
    Read a time-ordered list of PCAP files once per file, computing both
    the binned data and the CAL lock-in gain in a single pass.

    The two leftovers (binning and CAL) are tracked independently, as they
    cover different amounts of data and involve different arrays.

    Parameters
    ----------
    file_list : list of str, file names in time order (without directory)
    dir       : str, directory where the files are located
    verbose   : if True, print diagnostic information per file

    Returns
    -------
    jd    : array (nb_total,), absolute time of each bin (seconds)
    data  : array (nstokes, nb_total, nfreq), binned averages
    wei   : array (nstokes, nb_total, nfreq), weights (1/sigma^2)
    nuMHz : array (nfreq,), frequency axis in MHz
    gain3 : array (nstokes, nfreq, ncals_total), CAL lock-in gain per cycle
    lockin : lock-in template used
    """

    leftover_bin = None   # carries nsbin-1 samples of alldata at most
    leftover_cal = None   # carries nsamp-1 samples of alldata+cal_mask+cal_sgn at most
    t_abs        = 0.0

    jd_list    = []
    data_list  = []
    wei_list   = []
    gain3_list = []
    lockin     = None

    print(f" SANCHO_READ_BIN_CAL_MFI2_FPGA: reads data, averages over {nsbin}-sample bins, and extracts the CAL signal.")

    for i, file_name in enumerate(file_list):
        print(f' [{i+1}/{len(file_list)}] Reading {file_name} ...')

        # Single read per file
        alldata, cal_sgn, cal_mask, dindex = read_mfi2_fpga_pcap(file_name, dir=dir, verbose=verbose)

        # --- Binning (uses alldata, leftover_bin) ---
        jd, data, wei, nuMHz, leftover_bin = sancho_bin_data_mfi2_fpga_leftovers(
            alldata, leftover=leftover_bin, t_offset_sec=t_abs, verbose=verbose
        )
        jd_list.append(jd)
        data_list.append(data)
        wei_list.append(wei)

        # --- CAL lock-in (uses alldata + cal_mask + cal_sgn, leftover_cal) ---
        gain3, lockin, leftover_cal = sancho_cal_lockin_mfi2_fpga_leftovers(
            cal_mask, cal_sgn, alldata, leftover=leftover_cal, verbose=verbose
        )
        if gain3 is not None:
            gain3_list.append(gain3)
        else:
            print(f'   > Warning: no complete CAL cycle in {file_name}')

        # Advance absolute time offset for the next file
        t_abs += alldata[0].shape[0] * dt_sec

    # Concatenate results along time / CAL-cycle axes
    jd    = np.concatenate(jd_list)
    data  = np.concatenate(data_list,  axis=1)
    wei   = np.concatenate(wei_list,   axis=1)
    gain3 = np.concatenate(gain3_list, axis=2)

    print(f' Done. Total bins: {len(jd)}, total CAL cycles: {gain3.shape[2]}')
    return jd, data, wei, nuMHz, gain3, lockin


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


##################
# Examples of use
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

def test_sept2026():
    import time

    # Example data generated on December 2025.
    file_directory = '/Users/jalberto/quijote/RIM/example_data/'
    file_list = ['test1']

    # New code (allows leftovers)
    t0 = time.perf_counter()
    print('New code')
    jd2, data2, wei2, nuMHz2, gain2, lockin2 = sancho_read_bin_cal_mfi2_fpga(file_list, dir=file_directory, verbose=True)
    t_new = time.perf_counter() - t0
    print(f' New code elapsed time: {t_new:.3f} s')

    # OLD style
    t0 = time.perf_counter()
    print('Old code')
    alldata, cal_sgn, cal_mask, dindex = read_mfi2_fpga_pcap(file_list[0], dir=file_directory, verbose=True)
    jd1, data1, wei1, nuMHz1, t_all1 = sancho_bin_data_mfi2_fpga(alldata, verbose=True)
    gain1, lockin1 = sancho_cal_lockin_mfi2_fpga(cal_mask, cal_sgn, alldata, verbose=False)
    t_old = time.perf_counter() - t0
    print(f' Old code elapsed time: {t_old:.3f} s')

    # --- Speedup ---
    print(f'\n Speedup: {t_old/t_new:.2f}x  ({t_old:.3f} s -> {t_new:.3f} s)')

    # --- Numerical comparison ---
    print(jd1.shape,jd2.shape, np.std(jd1-jd2))
    print(data1.shape,data2.shape, np.std(data1-data2))
    print(wei1.shape, wei2.shape, np.std(wei1-wei2))
    print(gain1.shape, gain2.shape, np.std(gain1-gain2))
    print(lockin1.shape, lockin2.shape, np.std(lockin1-lockin2))

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
    dothis = False
    if dothis:
        test_jan2026()

    # Tests Septiembre 2026. Nuevos codigos para concatenar. Testing.
    dothis = True 
    if dothis:
        test_sept2026() 