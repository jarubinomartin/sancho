;   8/Nov/2012   J.A.Rubino, D. Tramonte
;    
;   HISTORY:
;    8/Nov/2012 - original version
;   21/Oct/2014 - added coordinates for TGI horns and TGI structure
;                 (D. Tramonte)
;   22/Jul/2015 - reads the modification dates  
;
;-
;


pro tgi_horn_positions, x, y

  x = dblarr(31)
  y = dblarr(31) 

;horn 1  -- in the center

  x[0] = 0.d
  y[0] = 0.d
  
;horns 2 to 7 -- distance from center r~86.64 mm

  x[1] = 86.639d 
  y[1] = 0.d
  
  x[2] = 43.319d
  y[2] = 75.027d
  
  x[3] = -43.319d
  y[3] = 75.027d
  
  x[4] = -86.639d
  y[4] = 0.d
  
  x[5] = -43.319d
  y[5] = -75.027d

  x[6] = 43.319d
  y[6] = -75.027d

;horns 8 to 13 -- distance from center r~150.07
  
  x[7] = 129.958d
  y[7] = 75.027d
  
  x[8] = 0.d
  y[8] = 150.065d
  
  x[9] = -129.958d
  y[9] = 75.027d
  
  x[10] = -129.958d
  y[10] = -75.027d
  
  x[11] = 0.d
  y[11] = -150.065d
  
  x[12] = 129.958d
  y[12] = -75.027d 

;horns 14 to 19 -- distance from center r~173.28

  x[13] = 173.278d
  y[13] = 0.d
  
  x[14] = 86.639d
  y[14] = 150.065d
  
  x[15] = -86.639d
  y[15] = 150.065d
  
  x[16] = -173.278d
  y[16] = 0.d
  
  x[17] = -86.639d
  y[17] = -150.065d
  
  x[18] = 86.639d
  y[18] = -150.065d
	
;horns 20 to 31 -- distance from center r~229.22

  x[19] = 216.597d
  y[19] = 75.027d
  
  x[20] = 173.278d
  y[20] = 150.065d
  
  x[21] = 43.319d
  y[21] = 225.092d
  
  x[22] = -43.319d
  y[22] = 225.092d
  
  x[23] = -173.278d
  y[23] = 150.065d
  
  x[24] = -216.597d
  y[24] = 75.027d
  
  x[25] = -216.597d
  y[25] = -75.027d
  
  x[26] = -173.278d
  y[26] = -150.065d
  
  x[27] = -43.319d
  y[27] = -225.092d
  
  x[28] = 43.319d
  y[28] = -225.092d
  
  x[29] = 173.278d
  y[29] = -150.065d
  
  x[30] = 216.597d
  y[30] = -75.027d

return
end

;******************************************************************************************


PRO info_quijote,mfi=mfi,tgi=tgi


;=============================================
; MFI
;
; MFI. Nominal values
freq_mfi  = [ 11.0, 13.0, 17.0, 19.0] ; GHz
beam_mfi  = [ 0.92, 0.92, 0.60, 0.60]  ; fwhm, deg
nchan_mfi = [    8,    8,    8,    8]  ; number of channels

; Focal plane: nominal values
horn_mfi  = ["Horn10-14A", "Horn16-20A","Horn10-14B", "Horn16-20B"]
xpos_mfi  = [90.0,  -155.88, -90.0, 155.88] ; mm
ypos_mfi  = [155.88, 90.0, -155.88, -90.0]  ; mm
;focal_mfi = 3637.0                          ; mm
focal_mfi = 3700.0                          ; mm New fitted value


; Channels of the MFI
nhorn          = 4  
nfreq_per_horn = 2  
nchan_per_horn = 4  
nchan          = nhorn * nfreq_per_horn * nchan_per_horn

chan_name  = strarr(nchan)
chan_name  = [" 113d"," 111d", " 113x"," 111x"," 113y"," 111y"," 113s"," 111s",$
              " 219d"," 217d", " 219x"," 217x"," 219y"," 217y"," 219s"," 217s",$
              " 313d"," 311d", " 313x"," 311x"," 313y"," 311y"," 313s"," 311s",$
              " 419d"," 417d", " 419x"," 417x"," 419y"," 417y"," 419s"," 417s" ]
chan_name = strtrim(chan_name,2)


polarizer_id = [1, 1, 1, 1, 1, 1, 1, 1, $ 
                2, 2, 2, 2, 2, 2, 2, 2, $
                3, 3, 3, 3, 3, 3, 3, 3, $
                4, 4, 4, 4, 4, 4, 4, 4 ]

chan_id      = ["d", "d", "x", "x", "y", "y", "s", "s", $
                "d", "d", "x", "x", "y", "y", "s", "s", $
                "d", "d", "x", "x", "y", "y", "s", "s", $
                "d", "d", "x", "x", "y", "y", "s", "s" ]



chan_freq_id = [" 12-14GHZ A"," 10-12GHZ A", " 12-14GHZ A"," 10-12GHZ A", " 12-14GHZ A"," 10-12GHZ A", " 12-14GHZ A"," 10-12GHZ A", $
                " 18-20GHz A"," 16-18GHz A", " 18-20GHz A"," 16-18GHz A", " 18-20GHz A"," 16-18GHz A", " 18-20GHz A"," 16-18GHz A",$ 
                " 12-14GHZ B"," 10-12GHZ B", " 12-14GHZ B"," 10-12GHZ B", " 12-14GHZ B"," 10-12GHZ B", " 12-14GHZ B"," 10-12GHZ B", $
                " 18-20GHz B"," 16-18GHz B", " 18-20GHz B"," 16-18GHz B", " 18-20GHz B"," 16-18GHz B", " 18-20GHz B"," 16-18GHz B"]
chan_freq_id = strtrim(chan_freq_id,2)


central_freq_id = [ 13,11, 13,11, 13,11, 13,11, $ 
                    19,17, 19,17, 19,17, 19,17, $
                    13,11, 13,11, 13,11, 13,11, $ 
                    19,17, 19,17, 19,17, 19,17 ] * 1.0



; Channel names (for plots)
chan_id_txt   = strarr(nchan)
chan_id_txt[where(chan_id eq "d")] = "(x-y)"
chan_id_txt[where(chan_id eq "x")] = "(x)"
chan_id_txt[where(chan_id eq "y")] = "(y)"
chan_id_txt[where(chan_id eq "s")] = "(x+y)"

chan_full_name = "Horn "+strtrim(chan_freq_id)+" -"+chan_id_txt


; Channel identification (for software)
ic        = intarr(nhorn,nfreq_per_horn,nchan_per_horn)
ic[0,0,*] = [8,2,4,6]-1
ic[0,1,*] = [7,1,3,5]-1
ic[1,0,*] = [16,10,12,14]-1
ic[1,1,*] = [15,9,11,13]-1
ic[2,0,*] = [24,18,20,22]-1
ic[2,1,*] = [23,17,19,21]-1
ic[3,0,*] = [32,26,28,30]-1
ic[3,1,*] = [31,25,27,29]-1


; JD reference. Taken as the first day of the Commissioning (13/Nov/2012, 0.00h)
jd_ref = 2456244.5d

; MFI modification dates. Taken from
; /net/nas4/quijote/etc/mfi_modifications_jd.txt
jd_mfi_modifications = [0.0, 300.0, 331.0, 514.0, 626.0, 735.0 ] + jd_ref




; Output structure
mfi = { freq:freq_mfi, beam:beam_mfi, fwhm_arcmin:beam_mfi*60.0, nhorn:nhorn, $
	nfreq_per_horn:nfreq_per_horn, nchan_per_horn:nchan_per_horn, $
        nchan_mfi:nchan_mfi, horn:horn_mfi, xpos:xpos_mfi, ypos:ypos_mfi, $
        focal:focal_mfi, ic:ic, nchan:nchan, chan_name:chan_name, $
        chan_full_name:chan_full_name, polarizer_id:polarizer_id, jd_ref:jd_ref,$
	central_freq_id:central_freq_id, jd_mfi_mod:jd_mfi_modifications, chan_id:chan_id }


;-----------------------------------------------------------------------------------------

; TGI. Nominal values

nhorn = 31
freq_tgi = dblarr(nhorn) + 30.d
beam_tgi  =  dblarr(nhorn) + 0.37d  ; fwhm, deg
nchan_tgi = intarr(nhorn) + 4  ; number of channels

; Focal plane: nominal values
horn_tgi = "Horn_" + strtrim(string(indgen(nhorn)+1),2)

tgi_horn_positions, xpos_tgi, ypos_tgi
focal_tgi = 3700.0  ; mm


; Channels of the TGI
nfreq_per_horn = 1
nchan_per_horn = 4  
nchan          = nhorn * nfreq_per_horn * nchan_per_horn

chan_name=["c01_30","c01_30","c01_30","c01_30","c02_30","c02_30","c02_30","c02_30","c03_30","c03_30","c03_30","c03_30",$ 
	  "c04_30","c04_30","c04_30","c04_30","c05_30","c05_30","c05_30","c05_30","c06_30","c06_30","c06_30","c06_30",$ 
	  "c07_30","c07_30","c07_30","c07_30","c08_30","c08_30","c08_30","c08_30","c09_30","c09_30","c09_30","c09_30",$ 
	  "c10_30","c10_30","c10_30","c10_30","c11_30","c11_30","c11_30","c11_30","c12_30","c12_30","c12_30","c12_30",$ 
	  "c13_30","c13_31","c13_30","c13_30","c14_30","c14_30","c14_30","c14_30","c15_30","c15_30","c15_30","c15_30",$
	  "c16_30","c16_30","c16_30","c16_30","c17_30","c17_30","c17_30","c17_30","c18_30","c18_30","c18_30","c18_30",$
	  "c19_30","c19_30","c19_30","c19_30","c20_30","c20_30","c20_30","c20_30","c21_30","c21_30","c21_30","c21_30",$	
	  "c20_30","c20_30","c20_30","c20_30","c21_30","c21_30","c21_30","c21_30","c22_30","c22_30","c22_30","c22_30",$
	  "c23_30","c23_30","c23_30","c23_30","c24_30","c24_30","c24_30","c24_30","c25_30","c25_30","c25_30","c25_30",$
	  "c26_30","c26_30","c26_30","c26_30","c27_30","c27_30","c27_30","c27_30","c28_30","c28_30","c28_30","c28_30",$
	  "c29_30","c29_30","c29_30","c29_30","c30_30","c30_30","c30_30","c30_30","c31_30","c31_30","c31_30","c31_30"]


chan_name = strtrim(chan_name,2)


polarizer_id = intarr(nchan)
polarizer_id = [1, 1, 1, 1, 2, 2, 2, 2, $
                3, 3, 3, 3, 4, 4, 4, 4, $
                5, 5, 5, 5, 6, 6, 6, 6, $
                7, 7, 7, 7, 8, 8, 8, 8, $
		9, 9, 9, 9, 10, 10, 10, 10, $
		11, 11, 11, 11, 12, 12, 12, 12, $
                13, 13, 13, 13, 14, 14, 14, 14, $
                15, 15, 15, 15, 16, 16, 16, 16, $
                17, 17, 17, 17, 18, 18, 18, 18, $
		19, 19, 19, 19, 20, 20, 20, 20, $
		21, 21, 21, 21, 22, 22, 22, 22, $
                23, 23, 23, 23, 24, 24, 24, 24, $
                25, 25, 25, 25, 26, 26, 26, 26, $
                27, 27, 27, 27, 28, 28, 28, 28, $
		29, 29, 29, 29, 30, 30, 30, 30, $
		31, 31, 31, 31]

chan_id      = intarr(nchan)
;TO BE INITIALIZED


chan_freq_id =indgen(nchan)/4 +1
chan_freq_id = strtrim(chan_freq_id,2)


central_freq_id = dblarr(nhorn) + 30.d


; Channel names (for plots)
chan_id_txt   = strarr(nchan)  
; TO BE INITIALIZED

chan_full_name = "Horn "+strtrim(chan_freq_id)+" -"+chan_id_txt


; Channel identification (for software)
ic        = intarr(nhorn,nfreq_per_horn,nchan_per_horn)
;TO BE INITIALIZED

; JD reference. Taken as the first day of the Commissioning (13/Nov/2012, 0.00h)
jd_ref_tgi = 2456244.5d   
 

; Output structure
tgi = { freq:freq_tgi, beam:beam_tgi, fwhm_arcmin:beam_tgi*60.0, nhorn:nhorn, $
	nfreq_per_horn:nfreq_per_horn, nchan_per_horn:nchan_per_horn, $
        nchan_tgi:nchan_tgi, horn:horn_tgi, xpos:xpos_tgi, ypos:ypos_tgi, $
        focal:focal_tgi, ic:ic, nchan:nchan, chan_name:chan_name, $
        chan_full_name:chan_full_name, polarizer_id:polarizer_id, jd_ref:jd_ref,$
	central_freq_id:central_freq_id }


END
