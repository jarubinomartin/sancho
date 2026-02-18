#!/usr/bin/env python

"""
2/Sep/2025   

@authors: jalberto, david.diaz

Basic I/O routines for the data processing of raw MFI2 files with the digital BEM based on FPGAs.
* read_mfi2_fpga_pcap
* write_mfi2_fpga_raw_data


HISTORY:
* 02/09/2025 - original version. Using the codes from David, and some routines from sancho_mfi2_io.py.
* 19/01/2026 - New version of the reading codes, based on the latest format from David (email 5/Dic/2025)
* 21/01/2026 - New version of the format_pcap_data() code, using sign_extend_52(). Requires Python 1.13.

"""

import numpy as np
import matplotlib.pyplot as plt
from scapy.all import *
from astropy.io import fits
import sys
from datetime import datetime


# Reading PCAP files.
# Codes from Sept 2025
def read_mfi2_fpga_pcap_sept2025(file_name, dir='./', verbose=False):
    if verbose:
        print(' READ_MFI2_FPGA_PCAP_SEPT2025: ',file_name)
    data_raw = read_pcap_files_sept2025(dir, file_name)
    [values, headers] = format_pcap_data_sept2025(data_raw)
    return values, headers

def read_pcap_files_sept2025(directory, file_name):
    pcap_file = directory + "/" + file_name
    total_packets = []
    reader = RawPcapReader(pcap_file)
    for packet_data, _ in reader:
        total_packets.append(packet_data)
    reader.close()
    return total_packets

def format_pcap_data_sept2025(data_raw):
    dout = {
        0xA0: [], 0xA1: [], 0xA2: [], 0xA3: [], 0xA4: [], 0XA5: [], 0xA6: [], 0xA7: [],
        0xA8: [], 0xA9: [], 0xAA: [], 0xAB: [], 0xAC: [], 0XAD: [], 0xAE: [], 0xAF: []
    }
    hout = {
        0xA0: [], 0xA1: [], 0xA2: [], 0xA3: [], 0xA4: [], 0XA5: [], 0xA6: [], 0xA7: [],
        0xA8: [], 0xA9: [], 0xAA: [], 0xAB: [], 0xAC: [], 0XAD: [], 0xAE: [], 0xAF: []
    }
    for packet in data_raw:
        try:
            packet_id = int.from_bytes(packet[42:48], 'big')
        except Exception:
            continue
        if packet_id not in dout:
            continue
        buffer = np.frombuffer(packet, dtype=np.uint8).copy()
        samples_raw_value = buffer[128:]
        samples_raw_header_4msb = (samples_raw_value[::8] & 0xF0) >> 7
        samples_raw_header_4lsb = samples_raw_value[::8] & 0x0F
        samples_raw_value[::8] = ((samples_raw_header_4lsb >> 3) * 0xF0) | samples_raw_header_4lsb
        if len(samples_raw_value) % 4 != 0:
            continue
        try:
            samples_int_values = np.frombuffer(samples_raw_value, dtype='>i8').ravel()
            samples_int_header = np.frombuffer(samples_raw_header_4msb, dtype=np.uint8).ravel()
        except Exception:
            continue
        dout[packet_id].append(samples_int_values)
        hout[packet_id].append(samples_int_header)
    return [[np.array(dout[h]) for h in
             (0xA0, 0xA1, 0xA2, 0xA3, 0xA4, 0xA5, 0xA6, 0xA7, 0xA8, 0xA9, 0xAA, 0xAB, 0xAC, 0xAD, 0xAE, 0xAF)],
            [np.array(hout[h]) for h in
             (0xA0, 0xA1, 0xA2, 0xA3, 0xA4, 0xA5, 0xA6, 0xA7, 0xA8, 0xA9, 0xAA, 0xAB, 0xAC, 0xAD, 0xAE, 0xAF)]]

# Codes from Jan2026. Final versions
# data     -> valor de la muestra
# cal_sgn  -> valor de la señal para el diodo de calibracion
# cal_mask -> valor de la máscara de activación de la calibración
# dindex   -> indice de la muestra
def read_mfi2_fpga_pcap(file_name, dir='./', verbose=False):
    if verbose:
        print(' READ_MFI2_FPGA_PCAP_JAN2026: ',file_name)
    [data_raw, time_points] = read_pcap_files(dir, file_name)
    print(time_points[0])
    [data, cal_sgn, cal_mask, dindex] = format_pcap_data(data_raw)
    return data, cal_sgn, cal_mask, dindex

def check_dindex(dindex,nstokes):
    nfreq = dindex[0].shape[1]
    nsamp = dindex[0].shape[0]
    for i in range(nstokes):
        idx = dindex[0]
        for j in range(nfreq):
            teo = (idx[0,j] + np.arange(nsamp,dtype=np.int16)) % 256 
            derror = np.sum(idx[:,j] - teo)
            if not np.isclose(derror, 0):
                raise ValueError("There is a missing index in dindex.") 
    return

def read_pcap_files(directory, file_name):
    pcap_file = directory + "/" + file_name
    total_packets = []
    time_points = []
    reader = RawPcapReader(pcap_file)
    for packet_data, pkt_metadata in reader:
        # Leer segundos y microsegundos
        ts_sec = pkt_metadata.sec
        ts_usec = pkt_metadata.usec
        # Combinar en timestamp con precisión de microsegundos
        ts = ts_sec + ts_usec / 1e6
        time_points.append(datetime.fromtimestamp(ts))
        # Guardar tupla (data cruda, timestamp legible)
        total_packets.append(packet_data)
    reader.close()
    return [total_packets, time_points]

def sign_extend_52(u64: np.ndarray) -> np.ndarray:
    # Sign-extend desde 52 bits (bit 51 es el signo) a int64
    sign_bit = np.uint64(1) << np.uint64(51)
    mask = (np.uint64(1) << np.uint64(52)) - 1  # 0x000FFFFFFFFFFFFF
    u = u64 & mask
    return ((u ^ sign_bit) - sign_bit).astype(np.int64)

def format_pcap_data(data_raw):
    ids = (0xA0, 0xA1, 0xA2, 0xA3, 0xA4, 0xA5, 0xA6, 0xA7,
           0xA8, 0xA9, 0xAA, 0xAB, 0xAC, 0xAD, 0xAE, 0xAF)

    data = {k: [] for k in ids}
    cal_sgn = {k: [] for k in ids}
    cal_msk = {k: [] for k in ids}
    dindex = {k: [] for k in ids}

    for packet in data_raw:
        try:
            packet_id = int.from_bytes(packet[42:48], 'big')
        except Exception:
            continue

        if packet_id not in data:
            continue

        nbytes = len(packet) - 128
        if nbytes <= 0 or (nbytes % 8) != 0:
            continue

        nsamples = nbytes // 8
        words = np.frombuffer(packet, dtype='>u8', count=nsamples, offset=128)

        # --- Extraer cabecera ---
        cal_sgn_bit = ((words >> 63) & 0x1).astype(np.uint8)
        cal_msk_bit = ((words >> 62) & 0x1).astype(np.uint8)

        # bits 59..52 (8 bits)
        dindex_word = ((words >> 52) & 0xFF).astype(np.uint8)

        # --- Extraer data (bits 51..0) con signo ---
        samples_int_values = sign_extend_52(words)

        data[packet_id].append(samples_int_values)
        cal_sgn[packet_id].append(cal_sgn_bit)
        cal_msk[packet_id].append(cal_msk_bit)
        dindex[packet_id].append(dindex_word)

    return [
        [np.array(data[h], dtype=object) for h in ids],
        [np.array(cal_sgn[h], dtype=object) for h in ids],
        [np.array(cal_msk[h], dtype=object) for h in ids],
        [np.array(dindex[h], dtype=object) for h in ids],
    ]

def format_pcap_data_dec2025(data_raw):
    ids = (0xA0, 0xA1, 0xA2, 0xA3, 0xA4, 0xA5, 0xA6, 0xA7,
           0xA8, 0xA9, 0xAA, 0xAB, 0xAC, 0xAD, 0xAE, 0xAF)
    data = {k: [] for k in ids}
    cal_sgn = {k: [] for k in ids}
    cal_msk = {k: [] for k in ids}
    dindex = {k: [] for k in ids}
    for packet in data_raw:
        try:
            packet_id = int.from_bytes(packet[42:48], 'big')
        except Exception:
            continue
        if packet_id not in data:
            continue
        nbytes = len(packet) - 128
        if nbytes <= 0 or (nbytes % 8) != 0:
            continue
        nsamples = nbytes // 8
        words = np.frombuffer(packet, dtype='>u8', count=nsamples, offset=128)
        byte1_msb = (words >> 56) & 0xFF
        byte0_msb = (words >> 48) & 0xFF
        cal_sgn_bit = ((byte1_msb >> 7) & 0x1).astype(np.uint8)
        cal_msk_bit = ((byte1_msb >> 6) & 0x1).astype(np.uint8)
        dindex_word = (((byte1_msb & 0x0F) << 4) | (byte0_msb >> 4)).astype(np.uint8)
        v = ((byte0_msb & 0x0F) * 0x11).astype(words.dtype)
        low48_mask = (1 << 48) - 1
        new_words = (words & low48_mask) | (v << 56) | (v << 48)
        samples_int_values = new_words.astype('>i8', copy=False)
        data[packet_id].append(samples_int_values)
        cal_sgn[packet_id].append(cal_sgn_bit)
        cal_msk[packet_id].append(cal_msk_bit)
        dindex[packet_id].append(dindex_word)
    return [
        [np.array(data[h], dtype=object) for h in ids],
        [np.array(cal_sgn[h], dtype=object) for h in ids],
        [np.array(cal_msk[h], dtype=object) for h in ids],
        [np.array(dindex[h], dtype=object) for h in ids]]



# Write RAW files for MFI2 with FPGAs. 
#
def write_mfi2_fpga_raw_data(raw, ffout, overwrite=False):
    print(' WRITE_MFI2_FPGA_RAW_DATA: writing '+ffout)

    # Definition of columns in RAW files    
    col1 = fits.Column(name='JD', format=str(len(raw['JD']))+'D', array = [raw['JD']] )
    col2 = fits.Column(name='DATA', format=str(raw['DATA'].size)+'E', array = [raw['DATA']], dim=str(raw['DATA'].shape[::-1]) )
    col3 = fits.Column(name='WEI', format=str(raw['WEI'].size)+'E', array = [raw['WEI']], dim=str(raw['WEI'].shape[::-1]) )
    col4 = fits.Column(name='NSBIN', format='E', array = [raw['NSBIN']])
    col5 = fits.Column(name='NUMHZ', format=str(len(raw['NUMHZ']))+'D', array = [raw['NUMHZ']] )

    # Bin table
    hdu = fits.BinTableHDU.from_columns([col1,col2,col3,col4,col5])
    
    # write file
    hdu.writeto(ffout,overwrite=overwrite)
    
    return
