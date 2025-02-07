'''
@author: jalberto

Basic information for QUIJOTE. It defines one dictionary per instrument.
Available dictionaries: mfi, tfgi, mfi2.

HISTORY:
 26/Jan/2025 - original version, based on info_quijote.pro


'''

import numpy as np

# JD reference. Taken as the first day of the Commissioning (13/Nov/2012, 0.00h)
jd_ref = 2456244.5

#=============================================
# MFI information
# Nominal values
freq_mfi  = [ 11.0, 13.0, 17.0, 19.0] # GHz
beam_mfi  = [ 0.92, 0.92, 0.60, 0.60] # fwhm, deg
nchan_mfi = [    8,    8,    8,    8] # number of channels

# Focal plane: nominal values
horn_mfi  = ["Horn10-14A", "Horn16-20A","Horn10-14B", "Horn16-20B"]
xpos_mfi  = [90.0,  -155.88, -90.0, 155.88] # mm
ypos_mfi  = [155.88, 90.0, -155.88, -90.0]  # mm
focal_mfi = 3700.0                          # mm New fitted value

# MFI channel identification (for software)
nhorn_mfi          = 4  
nfreq_per_horn_mfi = 2  
nchan_per_horn_mfi = 4  
nchan_mfi          = nhorn_mfi * nfreq_per_horn_mfi * nchan_per_horn_mfi
ic        = np.zeros(nchan_mfi,dtype='I').reshape(nhorn_mfi,nfreq_per_horn_mfi,nchan_per_horn_mfi)
ic[0,0,:] = np.asarray([8,2,4,6])-1
ic[0,1,:] = np.asarray([7,1,3,5])-1
ic[1,0,:] = np.asarray([16,10,12,14])-1
ic[1,1,:] = np.asarray([15,9,11,13])-1
ic[2,0,:] = np.asarray([24,18,20,22])-1
ic[2,1,:] = np.asarray([23,17,19,21])-1
ic[3,0,:] = np.asarray([32,26,28,30])-1
ic[3,1,:] = np.asarray([31,25,27,29])-1

chan_name  = np.asarray( [" 113d"," 111d", " 113x"," 111x"," 113y"," 111y"," 113s"," 111s",
                          " 219d"," 217d", " 219x"," 217x"," 219y"," 217y"," 219s"," 217s",
                          " 313d"," 311d", " 313x"," 311x"," 313y"," 311y"," 313s"," 311s", 
                          " 419d"," 417d", " 419x"," 417x"," 419y"," 417y"," 419s"," 417s" ] )


polarizer_id = np.asarray( [1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 
                            3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4 ] )

# Dictionary
mfi = {"ic":ic, "nhorn":nhorn_mfi, "nchan":nchan_mfi, "xpos":xpos_mfi, "ypos":ypos_mfi, "focal":focal_mfi,
       "chan_name":chan_name, "polarizer_id":polarizer_id, "jd_ref":jd_ref}


#=============================================
# TFGI. 
nhorn_tfgi = 31

x=np.zeros(nhorn_tfgi)
y=np.zeros(nhorn_tfgi)

x[1] = 86.639
y[1] = 0.0
  
x[2] = 43.319
y[2] = 75.027
  
x[3] = -43.319
y[3] = 75.027
  
x[4] = -86.639
y[4] = 0.0
  
x[5] = -43.319
y[5] = -75.027

x[6] = 43.319
y[6] = -75.027

x[7] = 129.958
y[7] = 75.027
  
x[8] = 0.0
y[8] = 150.065
  
x[9] = -129.958
y[9] = 75.027
  
x[10] = -129.958
y[10] = -75.027
  
x[11] = 0.0
y[11] = -150.065
  
x[12] = 129.958
y[12] = -75.027

x[13] = 173.278
y[13] = 0.0
  
x[14] = 86.639
y[14] = 150.065
  
x[15] = -86.639
y[15] = 150.065
  
x[16] = -173.278
y[16] = 0.0
  
x[17] = -86.639
y[17] = -150.065
  
x[18] = 86.639
y[18] = -150.065

x[19] = 216.597
y[19] = 75.027
  
x[20] = 173.278
y[20] = 150.065
  
x[21] = 43.319
y[21] = 225.092
  
x[22] = -43.319
y[22] = 225.092
  
x[23] = -173.278
y[23] = 150.065
  
x[24] = -216.597
y[24] = 75.027
  
x[25] = -216.597
y[25] = -75.027
  
x[26] = -173.278
y[26] = -150.065

x[27] = -43.319
y[27] = -225.092
  
x[28] = 43.319
y[28] = -225.092
  
x[29] = 173.278
y[29] = -150.065
  
x[30] = 216.597
y[30] = -75.027

focal_tgi = 3700.0  # mm

tfgi = {"nhorn":nhorn_tfgi, "xpos":x, "ypos":y, "focal":focal_tgi, "jd_ref":jd_ref}

#=============================================
# MFI2. 
# Nominal values
freq_mfi2  = [ 11.0, 13.0, 17.0, 19.0]  # GHz
beam_mfi2  = [ 0.92, 0.92, 0.60, 0.60]  # fwhm, deg
nchan_mfi2 = [   12,   12,    8,    8]  # number of channels

# Focal plane: nominal values. These are the numbers provided by
# Haroldo, email 8/May/2024.
horn_mfi2_all  = np.asarray(["Horn10-14A", "Horn16-20A","Horn10-14B", "Horn16-20B", "Horn10-14C"])
xpos_mfi2_all  = np.asarray([0.0,  97.50, -168.87,  -97.50, 168.87]) # mm
ypos_mfi2_all  = np.asarray([0.0, 168.87,   97.50, -168.87, -97.50]) # mm
# Using the ordering of connections to old DAS. Pixel 3 is the central
# one. Using X-coordinate as seen from sky, as in MFI and TFGI. Thus x->-x
idxmfi2    = [2,1,0,3,4]
horn_mfi2  = horn_mfi2_all[idxmfi2]
xpos_mfi2  = -1.0*xpos_mfi2_all[idxmfi2] # mm
ypos_mfi2  = ypos_mfi2_all[idxmfi2]      # mm
focal_mfi2 = 3700.0                      # mm

mfi2 = {"nhorn":5, "xpos":xpos_mfi2, "ypos":ypos_mfi2, "focal":focal_mfi2, "jd_ref":jd_ref, "ic":ic}
