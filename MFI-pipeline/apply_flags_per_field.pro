;   30/Nov/2017   JARM
;
;   For a given field, it applies the corresponding flagging tables.
;   Identified by inspection of data. 
;
;
;   HISTORY:
;     30/11/2017 - original version. Includes flags for W63, period 2, H3 11Ghz
;     07/12/2017 - specific flags for TAURUS region.
;     23/03/2018 - flags for NOMINAL30 and NOMINAL50, period2
;     29/08/2018 - flags for NOMINAL65, period6. There is a source in (l,b) to be removed. The code
;                  now also requires as input the (l,b) coordinates.
;     14/01/2019 - final flags for W49 and W51N. Also displays information on screen when flagging is applied. Flags by DT.
;     21/02/2019 - new flags from Federica for NOMINAL30 and NOMINAL35.
;      6/03/2019 - New subroutine to flag an elliptical region. New
;                  flags from Federica for all nominal.
;     15/03/2019 - FG: additinal flagging. New ellispes and only a
;                  few AZ flag. AZ flag for H3 11 only from el 60
;                  allperiods: 1,2,6 
;     21/03/2019 - FG: correction of 'az' keyword for ellipse
;                  flagging that was creating problems for AZ flag. 
;                  AZ flag of H3 11 only from el 65 per 1,2,6 and
;                  el 70 per 6
;     26/03/2019 - Correction of the flag_ellipse routine. Removing a
;                  not necessaru az flag in H3 11 el 65 per 1. Adding
;                  az flags for H3 11 from el 35 and 40, and for H3 13
;                  from el 65. 
;-
;     27/03/2019 - Add an ellpse in H4 19 el 65 per 6. Split in two
;                  ellipse an RFI of H3 13 el35 per 6. Change az range
;                  of H3 11 el35 per 6. Change some range of H3 11 el
;65.
;     05/04/2019 - Take out two ellipses of nom 40 per 2 H3 11 and 13
;                  at coord gl, gb = (3.4, 5.4)
;     25/04/2019 - New exclusion region (5 deg) for the QRT1. Also flagging now EL=50, period 6 [JARM]
;     26/04/2019 - Additional azimuth ranges for H3 11 and
;                  13 and H2/H4 17 and 19.
;     29/04/2019 - removing some flags (unnecesary), and reordering of the code.
;     03/05/2019 - Adding more ellipses for period 5.
;     07/05/2019 - Adding more ellipses and changing few az ranges.
;     07/05/2019 - Adding one ellipse in el 35 H4 19 and one az range
;                  in el 50 allperiods.
;     08/07/2019 - Adding AZ flags for H1 11 and 13, period 1 and 2.
;     24/09/2020 - adding a flag KEEPPLANETS, to avoid flagging planets
;     19/01/2021 - updated to check active flags per bit. Only activating when was previously inactive
;     20/01/2021 - new bit is used for planets: pbit=7. bit=8 is reserved for flagging regions.
;     05/03/2021 - updated flag for horn1, 13ghz, intensity
;
;--
;


FUNCTION FLAG_ELLIPSE, X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl, GB=gb, IAZ=iaz, AZ_INPUT = az_input
;Determine if a sample is in the ellipse given in input and return the
;indeces to flag
;Inputs: (X0, Y0) = (gl ,gb) in degrees of the center of the ellipse
;           alpha = rotation angle of the ellipse in degrees
;  rmajor, rminor = major and minor axes of the ellipse
;          gl, gb = longitude and latitude of the TOD of a determined
;          horn
;              az = 1 or 2 if the azimuth to flag is <180 or >180
;Output: indeces of the samples to be flagged.

  ;If we are close to the galactic meridian, force the continuity in gl.
  if ((x0 lt 40.) or (x0 gt 320.)) then begin
     gl_neg = where(gl gt 180.)
     gl[gl_neg] = gl[gl_neg]-360.
     if(x0 gt 320) then x0=x0-360.
  endif

  az   = (az_input + 360.d0*10.d0) mod 360.d0 ; guarantees that AZ is in [0,360] 

;Rotate coordinate system of alpha and
  ;center in the center of the ellipse
  xp = (gl-x0)*cos(gb*!dtor)
  yp = gb - y0
  
  a = alpha * !dtor
  x =  xp*COS(a) - yp*SIN(a)
  y =  xp*SIN(a) + yp*COS(a)
  
  
  ; Ellipse equation
  dist = sqrt(  (x/rmajor)^2. + (y/rminor)^2. )
  if (iaz eq 1) then cut = where(dist le 1 and az le 180.0, ncut)
  if (iaz eq 2) then cut = where(dist le 1 and az ge 180.0, ncut)

  return, cut
end

;=============

pro do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az_input, el, flag, one, mfi=mfi, kmin=kmin, kmax=kmax
if not keyword_set(mfi) then info_quijote,mfi=mfi
if not keyword_set(kmin) then kmin=0
if not keyword_set(kmax) then kmax=3
az   = (az_input + 360.d0*10.d0) mod 360.d0 ; guarantees that AZ is in [0,360]
nrfi = n_elements(el_rfi)
bit  = round(alog(one)/alog(2))+1
if one ne 2^(bit-1L) then stop," ERROR in bit definition. "
sancho_check_activebit_inflags, flag, bit, bitmasc
for j=0, nrfi-1 do begin
   cut = where( az ge azmin_rfi[j] and az le azmax_rfi[j] and abs(el - el_rfi[j]) le 0.5, ncut)
   if ncut ge 1 then begin
         for k=kmin,kmax do begin	
            ichan = mfi.ic[horn-1,freq,k]
            flag[ichan,cut] = flag[ichan,cut] + one*(1-bitmasc[ichan,cut])
         endfor
   endif 
endfor
end

;-----------

PRO apply_flags_per_field, filename, az, el, gl,gb,jd, flag, bit=bit, period=period, keepplanets=keepplanets, pbit=pbit

; Keywords
if not keyword_set(bit) then bit = 8
if not keyword_set(pbit) then pbit = 7 ; for planets and transients
if not keyword_set(period) then period=2
if not keyword_set(keepplanets) then keepplanets=0

; Info
if n_params() lt 6 then begin
   print,""
   print,"  Syntax -- apply_flags_per_field, filename, az, el, gl,gb,jd, flag, bit=bit, pbit=pbit, period=period, /KEEPPLANETS "
   print,""
   print,""
   return
endif

version = "March-12-2020"
version = "January-19-2021"
version = "March-5-2021"

; Params
one  = 2^(bit-1L)
onep = 2^(pbit-1L)
info_quijote,mfi=mfi

; Get info of fields
sancho_extract_rootname, filename, root, modpos, campo


; Flagging applied?
flageado = 0

;-----------------
; W49. Flags by D.Tramonte, 9/Jan/2019. Updated JARM 11/Jan/2019.
;      Final update 14/Jan/2018 by Denis & Alberto
;
if (campo eq "W49") then begin

; Horn3, 11ghz
   horn = 3
   freq = 0
   ;azmin_rfi = [  95.0, 140.0, 250.0]
   ;azmax_rfi = [ 125.0, 170.0, 254.0 ]
   ;el_rfi    = [ 49.94, 65.73, 33.06]
   azmin_rfi = [94.0,   130.0,  250.0, 270.0 ]
   azmax_rfi = [126.0,  165.0,  254.0, 275.0 ]
   el_rfi    = [49.9,   65.7,   33.1,  33.1  ]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi
; Horn3, 13ghz
   horn = 3
   freq = 1
   azmin_rfi = [113.0,  147.0,  250.0  ]
   azmax_rfi = [126.0,  165.0,  253.0  ]
   el_rfi    = [49.9,   65.7,   33.1   ]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi

; Horn4, 17ghz
   horn = 4
   freq = 0
   azmin_rfi = [246.0 ]
   azmax_rfi = [249.0 ]
   el_rfi    = [33.1 ]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi
; Horn4, 19ghz
   horn = 4
   freq = 1
   azmin_rfi = [246.0 ]
   azmax_rfi = [249.0 ]
   el_rfi    = [33.1 ]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi

   flageado = 1

endif


if (campo eq "W51N") then begin

; Horn2, 19ghz
   horn = 2
   freq = 1
   azmin_rfi = [ 255.0,  275.0 ]
   azmax_rfi = [ 267.0,  285.0 ]
   el_rfi    = [ 35.5,   35.5  ]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi

; Horn3, 11ghz
   horn = 3
   freq = 0
   azmin_rfi = [ 92.0,  104.0, 255.0, 180.0 ]
   azmax_rfi = [ 105.0, 120.0, 285.0, 204.0 ]
   el_rfi    = [ 35.2,  51.8,  35.5,  72.3  ]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi

; Horn4, 17ghz
   horn = 4
   freq = 0
   azmin_rfi = [ 256.0 ]
   azmax_rfi = [ 260.0 ]
   el_rfi    = [ 55.2 ]    ; --->> NOTICE: this could actually be a galactic source!
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi

; Horn4, 19ghz
   horn = 4
   freq = 1
   azmin_rfi = [ 256.0 ]
   azmax_rfi = [ 260.0 ]
   el_rfi    = [ 55.2 ]    ; --->> NOTICE: this could actually be a galactic source!
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi

   flageado = 1
endif


;-----------------
; W63. Flags by JARM, 30/Nov/2017
;
if (campo eq "W63") then begin
; Horn3, 11ghz
   horn = 3
   freq = 0
   azmin_rfi = [42.7, 51.1, 43.2, 51.1]
   azmax_rfi = [48.4, 56.9, 48.4, 55.5]
   el_rfi    = [52.0, 52.0, 53.4, 53.4]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi
; Horn2, 19ghz
   horn = 2
   freq = 1
   azmin_rfi = [42.9]
   azmax_rfi = [47.6]
   el_rfi    = [53.4]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi

   flageado = 1
endif

;-----------------
; TAURUS. Flags by JARM, 7/Dec/2017
;
if (campo eq "TAURUS") then begin
; Horn2, 17ghz
   horn = 2
   freq = 0
   azmin_rfi = [84.0,   0.0, 260.0, 283.0, 180.0,   0.0,   0.0]
   azmax_rfi = [93.0, 180.0, 270.0, 295.0, 360.0, 180.0, 180.0]
   el_rfi    = [43.2,  58.2,  55.7,  55.7,  80.1,  34.5,  49.1]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi
; Horn2, 19ghz
   horn = 2
   freq = 1
   azmin_rfi = [273.0, 180.0, 180.0]
   azmax_rfi = [282.0, 360.0, 360.0]
   el_rfi    = [ 35.4,  80.1,  65.2]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi
; Horn3, 11ghz
   horn = 3
   freq = 0
   azmin_rfi = [  0.0, 180.0]
   azmax_rfi = [180.0, 360.0]
   el_rfi    = [ 58.2,  71.2]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi
; Horn3, 13ghz
   horn = 3
   freq = 1
   azmin_rfi = [180.0]
   azmax_rfi = [360.0]
   el_rfi    = [ 71.2]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi
   flageado = 1
endif



;***************************************
; FLAGS FOR MFI WIDE SURVEY - 2/MAY/2019
;***************************************
;
; [1] FLAGGING BY AZ RANGES.


;-----------------
; NOMINAL30. Flags by JARM, 23/Mar/2018.
;
if (campo eq "NOMINAL30" and period eq 2) then begin

   print," Flagging ",campo
; Horn1, 11ghz
   horn = 1
   freq = 0
   azmin_rfi = [ 125.0 ]
   azmax_rfi = [ 140.0 ]
   el_rfi    = [ 30.0 ]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi
   flageado = 1

   print," Flagging ",campo
; Horn1, 13ghz
   horn = 1
   freq = 1
   azmin_rfi = [ 105.0 ]
   azmax_rfi = [ 135.0 ]
   el_rfi    = [ 30.0 ]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi
   flageado = 1

   print," Flagging ",campo
; Horn3, 11ghz
   horn = 3
   freq = 0
   azmin_rfi = [  0.0, 335.0]
   azmax_rfi = [ 25.0, 360.0]
   el_rfi    = [ 30.0,  30.0]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi
   flageado = 1
endif

if (campo eq "NOMINAL30" and period eq 6) then begin
   print," Flagging ",campo
   ; Horn3, 13ghz, cor/uncor
   horn = 3
   freq = 1 
   azmin_rfi = [ 117.0 ] ; Upper side of low dec band (Bright)
   azmax_rfi = [ 122.0 ]
   el_rfi    = [ 30.0 ]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi 
   flageado = 1
endif


;-----------------
; NOMINAL35. AZ flags by FG, 21/02/2019.
;
if (campo eq "NOMINAL35" and period eq 6) then begin
   print," Flagging ",campo
 
; ADDITIONAL AZ FLAGS for H2 19 07/05/19.
   horn = 2
   freq = 1
   azmin_rfi = [ 237.90 ]
   azmax_rfi = [ 240.95 ]
   el_rfi    = [  35.0 ]

   
   ; Horn3
   horn = 3
   freq = 0
   azmin_rfi = [  95.0 ] ; Around upper side of flagged band 
   azmax_rfi = [ 109.0 ]
   el_rfi    = [  35.0 ]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi

   
   ; Horn 4
   horn = 4
   freq = 0
   azmin_rfi = [ 120.0 ]
   azmax_rfi = [ 122.2 ]
   el_rfi    = [  35.0 ]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi

 
; UPDATE AZ FLAGS for H4 19 07/05/19.
   horn = 4
   freq = 1
   azmin_rfi = [ 119.0, 132.0 ]
   azmax_rfi = [ 126.0, 137.0 ]
   el_rfi    = [  35.0, 35.0 ]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi

   ; ADDITIONAL AZ FLAGS for H3 11 and 13 26/04/19. REVISED ON 29/04/2019
   ; Horn3, 11ghz, cor/uncor
   horn = 3
   freq = 0 
   azmin_rfi = [ 122.0 ] ; Upper side of low dec band (Bright)
   azmax_rfi = [ 137.0 ]
   el_rfi    = [  35.0 ]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi

   ; Horn3, 13ghz, cor/uncor  Back to the version of 03/05/19
   horn = 3
   freq = 1
   azmin_rfi = [ 121.8, 130.0 ] ; Upper side of low declination band and thin stripe below
   azmax_rfi = [ 126.8, 132.0 ]
   el_rfi    = [  35.0, 35.0]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi

   
   flageado = 1

endif

;-----------------------------
; NOMINAL40. AZ flags for periods 2 and 5
;
if (campo eq "NOMINAL40" and period ge 2) then begin
   print," Flagging ",campo

; Horn1, 11ghz
   horn = 1
   freq = 0
   azmin_rfi = [ 140.0 ]
   azmax_rfi = [ 150.0 ]
   el_rfi    = [ 40.0 ]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi

; Horn1, 13ghz
   horn = 1
   freq = 1
   azmin_rfi = [ 140.0, 210.0 ]
   azmax_rfi = [ 150.0, 250.0 ]
   el_rfi    = [ 40.0, 40.0 ]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi

   ; Horn3, 11ghz, cor/uncor
   horn = 3
   freq = 0
   azmin_rfi = [  99.0, 128.0 ] ;Around upper side of flagged region; ;Lower side of low dec band (Bright)
   azmax_rfi = [ 112.0, 135.0 ]
   el_rfi    = [  40.0,  40.0 ]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi

   
   ; Horn3, 13ghz, cor/uncor
   horn = 3
   freq = 1 
   azmin_rfi = [ 128.8, 104.0 ] ;Borders of low declinatin band (Bright)
   azmax_rfi = [ 137.3, 112.6 ]
   el_rfi    = [  40.0,  40.0 ]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi
   
   flageado = 1
   
endif

;-----------------
; NOMINAL50. AZ flags by JARM, 8/Mar/2018. New flags by FG merged
;
if (campo eq "NOMINAL50" and period eq 2) then begin 
   print," Flagging ",campo
   
; Horn1, 11ghz
   horn = 1
   freq = 0
   azmin_rfi = [ 50.0, 170.0, 228.0 ]
   azmax_rfi = [ 60.0, 200.0, 235.0 ]
   el_rfi    = [ 50.0, 50.0, 50.0 ]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi
   flageado = 1
   
   print," Flagging ",campo
; Horn1, 13ghz
   horn = 1
   freq = 1
   azmin_rfi = [ 45.0, 170.0 ]
   azmax_rfi = [ 60.0, 200.0 ]
   el_rfi    = [ 50.0, 50.0 ]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi
   flageado = 1
endif

if (campo eq "NOMINAL50") then begin 
   print," Flagging ",campo
   ; Horn 2, 19ghz
   horn = 2
   freq = 1
   azmin_rfi = [ 138.0 ]
   azmax_rfi = [ 140.6 ]
   el_rfi    = [  50.0 ]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi

   ; Horn3, 11ghz, cor/uncor
   horn = 3
   freq = 0 
   azmin_rfi = [ 42.0, 105.0, 148.0 ] ; 1) Satellite crossing Cygnus; 2) Near sidelobes of sats around dec=0 ; 3) Same, but for negative declinations. 
                                ; NOte by JARM: second cut could be
                                ; extended. We did it and it was not
                                ; improving the maps
   azmax_rfi = [ 60.0, 126.0, 154.0 ]
   el_rfi    = [ 50.0,  50.0,  50.0 ]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi

   ; Horn3, 13ghz, cor/uncor
   horn = 3
   freq = 1
   azmin_rfi = [ 190.8 ] ;This is the bright sat in the center of low dec band
   azmax_rfi = [ 199.5 ]
   el_rfi    = [ 50.0 ]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi

; Horn 4, 19ghz ADDED on 09/05/19
   horn = 4
   freq = 1
   azmin_rfi = [ 180.0 ]
   azmax_rfi = [ 200.0 ]
   el_rfi    = [  50.0 ]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi

   flageado = 1 
endif

if (campo eq "NOMINAL50" and period eq 5) then begin
   print," Flagging ",campo

   ; Horn3, 11ghz, cor/uncor
   horn = 3
   freq = 0 
   azmin_rfi = [ 196.0 ] ; adding a small satellite at negative declinations. 
   azmax_rfi = [ 210.0 ]
   el_rfi    = [  50.0 ]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi
   
   flageado = 1   
endif


;-----------------
; NOMINAL60. AZ ranges
;
if (campo eq "NOMINAL60" and period eq 1) then begin
   print," Flagging ",campo

; Horn1, 13ghz
   horn = 1
   freq = 1
   azmin_rfi = [ 190.0 ]
   azmax_rfi = [ 210.0 ]
   el_rfi    = [ 60.0 ]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi

; Horn1, 13ghz. Added by JARubino on March 5th, 2021.
   horn = 1
   freq = 1
   azmin_rfi = [ 145.0 ]
   azmax_rfi = [ 155.0 ]
   el_rfi    = [ 60.0 ]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi

   ; Horn3, 11ghz, cor/uncor
   horn = 3
   freq = 0 
   azmin_rfi = [ 46.5, 102.0, 213.0 ] ;One large ring around the NCP, band with far sidelobes and one thin residual satellite on the right of the flagged part
   azmax_rfi = [ 51.1, 130.0, 218.0 ]
   el_rfi    = [ 60.0,  60.0,  60.0 ]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi

   flageado = 1
endif

if (campo eq "NOMINAL60" and period ge 2) then begin
   print," Flagging ",campo

   print," Flagging ",campo
; Horn1, 11ghz
   horn = 1
   freq = 0
   azmin_rfi = [ 215.0 ]
   azmax_rfi = [ 219.0 ]
   el_rfi    = [ 60.0 ]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi


; Horn1, 13ghz
   horn = 1
   freq = 1
   azmin_rfi = [ 195.0 ]
   azmax_rfi = [ 230.0 ]
   el_rfi    = [ 60.0 ]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi


; Horn3, 11ghz, cor/uncor
   horn = 3
   freq = 0 
   azmin_rfi = [ 102.0, 213.0 ] ;Band with far sidelobes and one thin residual satellite on the right of the flagged part that I see in az stack
   azmax_rfi = [ 132.0, 218.0 ]
   el_rfi    = [  60.0,  60.0 ]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi

   flageado = 1
endif



;-----------------
; NOMINAL65. AZ flags
;
if (campo eq "NOMINAL65" and period eq 1) then begin
   print," Flagging ",campo
   
   ; Horn3, 11ghz, cor/uncor
   horn = 3
   freq = 0 
   azmin_rfi = [ 167.0 ]
   azmax_rfi = [ 220.0 ]
   el_rfi    = [ 65.0 ]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi

   ; Horn3, 13ghz, cor/uncor
   horn = 3
   freq = 1 
   azmin_rfi = [ 190.0 ] 
   azmax_rfi = [ 197.0 ]
   el_rfi    = [ 65.0 ]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi

   flageado = 1
endif


if (campo eq "NOMINAL65" and period eq 2) then begin
   print," Flagging ",campo
   
  ; Horn3, 11ghz, cor/uncor
   horn = 3
   freq = 0 
   azmin_rfi = [ 160.0 ] ;Thin and bright sat in low declination band
   azmax_rfi = [ 200.0 ]
   el_rfi    = [ 65.0 ]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi


   ; Horn3, 13ghz, cor/uncor
   horn = 3
   freq = 1 
   azmin_rfi = [ 190.0 ] ;Bright features in stack az
   azmax_rfi = [ 197.0 ]
   el_rfi    = [  65.0 ]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi

   flageado = 1
endif


if (campo eq "NOMINAL65" and period eq 6) then begin
   print," Flagging ",campo

   ; Horn3, 11ghz, cor/uncor
   horn = 3
   freq = 0 
   azmin_rfi = [ 165.0 ];Thin and bright sat in low declination band
   azmax_rfi = [ 200.0 ]
   el_rfi    = [ 65.0 ]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi

   ; Horn3, 13ghz, cor/uncor
   horn = 3
   freq = 1 
   azmin_rfi = [ 190.0 ] ;Bright reaatures in stack az
   azmax_rfi = [ 196.0 ]
   el_rfi    = [ 65.0 ]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi

   
   ; Horn3, 11ghz, cor/uncor
   horn = 3
   freq = 0 
   azmin_rfi = [ 110.0, 211.0 ] ; RFI that crosses the zone flagged for the QRT1. JARM note: leave the last one.
   azmax_rfi = [ 143.0, 223.0 ]
   el_rfi    = [  65.0,  65.0 ]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi


   ; Additional flags by FG, 2/May/2019
   ; Horn2, 17 and 19ghz, cor/uncor
   horn = 2
   freq = 0 
   azmin_rfi = [ 44.5 ] ;Large Ring around NCP
   azmax_rfi = [ 49.5 ]
   el_rfi    = [ 65.0 ]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi

   horn = 2
   freq = 1
   azmin_rfi = [ 44.5 ] ;Large Ring around NCP
   azmax_rfi = [ 49.5 ]
   el_rfi    = [ 65.0 ]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi


   ; Horn4, 17  cor/uncor
   horn = 4
   freq = 0 
   azmin_rfi = [ 43.0 ] ;Large Ring around NCP
   azmax_rfi = [ 49.0 ]
   el_rfi    = [ 65.0 ]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi


   ; Horn4, 19ghz, cor/uncor   
   horn = 4
   freq = 1
   azmin_rfi = [ 43.0 ] ;Large Ring around NCP
   azmax_rfi = [ 49.0 ]
   el_rfi    = [ 65.0 ]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi

   flageado = 1
endif


;-----------------
; NOMINAL70. 
if (campo eq "NOMINAL70" and period eq 6) then begin
   print," Flagging ",campo
   
   ; Horn3, 11ghz, cor/uncor
   horn = 3
   freq = 0 
   azmin_rfi = [ 135.0 ]
   azmax_rfi = [ 213.0 ]
   el_rfi    = [  70.0 ]
   do_flags_regions, horn, freq, azmin_rfi, azmax_rfi, el_rfi, az, el, flag, one, mfi=mfi

   flageado = 1
endif




;***************************************
; FLAGS FOR MFI WIDE SURVEY - 2/MAY/2019
;***************************************
;
; [2] FLAGGING BY ELLIPSES.

;-----------------
; NOMINAL30, period 2. 
; Adding flags 27/03/2019, FG
if (campo eq "NOMINAL30" and period eq 2) then begin
   print," Flagging ",campo
; Adding H2 19 flags on 07/05/19 
; Horn2, 19ghz, cor/uncor
;(1)
   horn = 2
   freq = 1
   ihorn = horn-1
      
   ;Ellipse parameters
   x0 = 17.53                   ;gl center degrees
   y0 = -1.33                   ;gb center degrees
   alpha = 65.0                 ;orientation of the ellipse, degrees
   rmajor = 8.50                ;semi-major axis
   rminor = 1.5                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
      
   ;Flag
   if n_elements(cut) ge 1 then begin
      sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif
   
;(2)
   horn = 2
   freq = 1
   ihorn = horn-1
      
   ;Ellipse parameters
   x0 = 104.0                   ;gl center degrees
   y0 = -74.1                   ;gb center degrees
   alpha = 10.0                 ;orientation of the ellipse, degrees
   rmajor = 9.0                 ;semi-major axis
   rminor = 2.0                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
      
   ;Flag
   if n_elements(cut) ge 1 then begin
      sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif
   
; Horn4, 19ghz, cor/uncor
   horn = 4
   freq = 1
   ihorn = horn-1
      
   ;Ellipse parameters
   x0 = 34.77                   ;gl center degrees
   y0 = -44.0                   ;gb center degrees
   alpha = 68.0                 ;orientation of the ellipse, degrees
   rmajor = 22.0                ;semi-major axis
   rminor = 1.5                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
      
   ;Flag
   if n_elements(cut) ge 1 then begin
      sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif


   ;(3) ADDED on 13/05/19
      horn = 4
      freq = 1
      ihorn = horn-1
      
                                ;Ellipse parameters
      x0 = 25.38              ;gl center degrees
      y0 = -24.55             ;gb center degrees
      alpha = 70.0             ;orientation of the ellipse, degrees
      rmajor = 3.0              ;semi-major axis
      rminor = 2.0             ;semi-minor axis
      iaz = 1                    ;azimuth minor (1) or major (2) than 180
      
      cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
      
      ;Flag
      if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
         ichan = reform(mfi.ic[ihorn,freq,0:3],4)
         for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
      endif

   flageado = 1
endif



;-----------------
; NOMINAL35, period 6. Flags by FG, 21/02/2019.
;
if (campo eq "NOMINAL35" and period eq 6) then begin
   print," Flagging ",campo

   
; Horn2, 19ghz, cor/uncor

;(1)
      horn = 2
      freq = 1
      ihorn = horn-1
      
                                ;Ellipse parameters
      x0 = 11.0                 ;gl center degrees
      y0 = 18.50                ;gb center degrees
      alpha = 60.0              ;orientation of the ellipse, degrees
      rmajor = 8.0              ;semi-major axis
      rminor = 1.0              ;semi-minor axis
      iaz = 1                    ;azimuth minor (1) or major (2) than 180
      
      cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
      
                                ;Flag
      if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
         ichan = reform(mfi.ic[ihorn,freq,0:3],4)
         for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
      endif



;(2)
      horn = 2
      freq = 1
      ihorn = horn-1
      
                                ;Ellipse parameters
      x0 = 27.14                ;gl center degrees
      y0 = -16.43               ;gb center degrees
      alpha = 65.0              ;orientation of the ellipse, degrees
      rmajor = 14.0             ;semi-major axis
      rminor = 0.8              ;semi-minor axis
      iaz = 1                    ;azimuth minor (1) or major (2) than 180
      
      cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
      
                                ;Flag
      if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
         ichan = reform(mfi.ic[ihorn,freq,0:3],4)
         for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
      endif


;(3)
      horn = 2
      freq = 1
      ihorn = horn-1
      
                                ;Ellipse parameters
      x0 = 74.24                ;gl center degrees
      y0 = -67.80               ;gb center degrees
      alpha = 40.0              ;orientation of the ellipse, degrees
      rmajor = 10.0             ;semi-major axis
      rminor = 1.0              ;semi-minor axis
      iaz = 1                    ;azimuth minor (1) or major (2) than 180
      
      cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
      
                                ;Flag
      if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
         ichan = reform(mfi.ic[ihorn,freq,0:3],4)
         for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
      endif


;(4)
      horn = 2
      freq = 1
      ihorn = horn-1
      
                                ;Ellipse parameters
      x0 = 103.1                ;gl center degrees
      y0 = -70.50               ;gb center degrees
      alpha = 22.0              ;orientation of the ellipse, degrees
      rmajor = 5.00             ;semi-major axis
      rminor = 1.5              ;semi-minor axis
      iaz = 1                    ;azimuth minor (1) or major (2) than 180
      
      cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
      
                                ;Flag
      if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
         ichan = reform(mfi.ic[ihorn,freq,0:3],4)
         for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
      endif



; Horn3, 11ghz, uncor ONLY

;(1)
      horn = 3
      freq = 0
      ihorn = horn-1
      
      ;Ellipse parameters
      x0 = 15.4                 ;gl center degrees
      y0 = -3.3                 ;gb center degrees
      alpha = 65.0              ;orientation of the ellipse, degrees
      rmajor = 4.0              ;semi-major axis
      rminor = 1.0              ;semi-minor axis
      iaz = 1                    ;azimuth minor (1) or major (2) than 180
      cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
      
      ;Flag
      if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
         ichan = reform(mfi.ic[ihorn,freq,2:3],2)
         for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
      endif


;(2)
      horn = 3
      freq = 0
      ihorn = horn-1
      
                                ;Ellipse parameters
      x0 = 189.0                ;gl center degrees
      y0 = -63.0                ;gb center degrees
      alpha = -63.0             ;orientation of the ellipse, degrees
      rmajor = 9.0              ;semi-major axis
      rminor = 1.5              ;semi-minor axis
      iaz = 1                    ;azimuth minor (1) or major (2) than 180
      
      cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
      
                                ;Flag
      if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
         ichan = reform(mfi.ic[ihorn,freq,2:3],2)
         for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
      endif

   

; Horn3, 13ghz, cor/uncor

;(1)
   horn = 3
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 16.0                    ;gl center degrees
   y0 = -5.5                    ;gb center degrees
   alpha = 65.0                 ;orientation of the ellipse, degrees
   rmajor = 4.0                 ;semi-major axis
   rminor = 1.0                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif


;(2)
   horn = 3
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 14.2                    ;gl center degrees
   y0 = 0.3                     ;gb center degrees
   alpha = 65.0                 ;orientation of the ellipse, degrees
   rmajor = 3.0                 ;semi-major axis
   rminor = 1.0                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif



;(3)
   horn = 3
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 188.0                   ;gl center degrees
   y0 = -64.0                   ;gb center degrees
   alpha = -58.0                ;orientation of the ellipse, degrees
   rmajor = 16.0                ;semi-major axis
   rminor = 1.5                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif

   ;(4) Same RFI captured in two ellipses
   horn = 3
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 6.8                    ;gl center degrees
   y0 = 7.8                   ;gb center degrees
   alpha = 62.0                 ;orientation of the ellipse, degrees
   rmajor = 12.0                ;semi-major axis
   rminor = 1.0                 ;semi-minor axis
   iaz = 2                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif
   

   horn = 3
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 17.8                    ;gl center degrees
   y0 = -14.8                   ;gb center degrees
   alpha = 68.0                 ;orientation of the ellipse, degrees
   rmajor = 18.0                ;semi-major axis
   rminor = 1.0                 ;semi-minor axis
   iaz = 2                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif
   


; Horn4, 17ghz, cor/uncor

;(1)
   horn = 4
   freq = 0
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 34.3 		        ;gl center degrees
   y0 = -34.0		        ;gb center degrees
   alpha = 66.0                 ;orientation of the ellipse, degrees
   rmajor = 32.0                ;semi-major axis
   rminor = 1.5                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif
   
   

   
; Horn4, 19ghz, cor/uncor

;UPDATE 07/05/19
;This taken out and replaced by the az range
   ;(1)
;      horn = 4
;      freq = 1
;      ihorn = horn-1
      
                	                ;Ellipse parameters
;      x0 = 32.7 		        ;gl center degrees
;      y0 = -31.1		       ;gb center degrees
;      alpha = 66.0                      ;orientation of the ellipse, degrees
;      rmajor = 34.0                     ;semi-major axis
;      rminor = 1.7              ;semi-minor axis
;      iaz = 1                    ;azimuth minor (1) or major (2) than 180
      
;      cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
      
                                ;Flag
;      if n_elements(cut) ge 1 then begin
;         ichan = reform(mfi.ic[ihorn,freq,0:3],4)
;         for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
;      endif


;(the central long stripe satellites that gets brighter in these coordinates. The effect is mostly in polarization.)

;(2)

      horn = 4
      freq = 1
      ihorn = horn-1
      
                                ;Ellipse parameters
      x0 = 92.5                ;gl center degrees
      y0 = -64.9               ;gb center degrees
      alpha = 22.0             ;orientation of the ellipse, degrees
      rmajor = 15.0             ;semi-major axis
      rminor = 3.0             ;semi-minor axis
      iaz = 1                    ;azimuth minor (1) or major (2) than 180
      
      cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
      
                                ;Flag
      if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
         ichan = reform(mfi.ic[ihorn,freq,0:3],4)
         for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
      endif

;(3)
      horn = 4
      freq = 1
      ihorn = horn-1
      
                                ;Ellipse parameters
      x0 = 7.1                ;gl center degrees
      y0 = 30.0               ;gb center degrees
      alpha = 55.0             ;orientation of the ellipse, degrees
      rmajor = 15.0             ;semi-major axis
      rminor = 2.0             ;semi-minor axis
      iaz = 1                    ;azimuth minor (1) or major (2) than 180
      
      cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
      
                                ;Flag
      if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
         ichan = reform(mfi.ic[ihorn,freq,0:3],4)
         for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
      endif
      
   flageado = 1

endif



;-----------------------------
; NOMINAL40, period 2.

if (campo eq "NOMINAL40" and period eq 2) then begin
   print," Flagging ",campo

; Horn2, 19ghz, cor/uncor

   horn = 2
   freq = 1
   ihorn = horn-1
   
;(1)                                ;Ellipse parameters
   x0 = 312.0                   ;gl center degrees
   y0 = 64.7                    ;gb center degrees
   alpha = 12.0                 ;orientation of the ellipse, degrees
   rmajor = 8.0                 ;semi-major axis
   rminor = 1.6                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif


;(2)   
   horn = 2
   freq = 1
   ihorn = horn-1
   
                               ;Ellipse parameters
   x0 = 197.5                   ;gl center degrees
   y0 = -26.5                   ;gb center degrees
   alpha = -68.0                ;orientation of the ellipse, degrees
   rmajor = 6.0                 ;semi-major axis
   rminor = 1.5                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif


;(3)   
   horn = 2
   freq = 1
   ihorn = horn-1
   
                               ;Ellipse parameters
   x0 = 30.2                   ;gl center degrees
   y0 = -26.0                   ;gb center degrees
   alpha = 70.00                ;orientation of the ellipse, degrees
   rmajor = 6.2                 ;semi-major axis
   rminor = 1.1                 ;semi-minor axis
   iaz = 2                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif



; Horn3, 11ghz, ONLY uncor

   horn = 3
   freq = 0
   ihorn = horn-1
   
;(1)                                ;Ellipse parameters
   x0 = 345.0                   ;gl center degrees
   y0 = 70.0                    ;gb center degrees
   alpha = 45.0                 ;orientation of the ellipse, degrees
   rmajor = 8.0                 ;semi-major axis
   rminor = 2.0                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,2:3],2)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif


;(2)

   horn = 3
   freq = 0
   ihorn = horn-1
   
   x0 = 213.6                   ;gl center degrees
   y0 = 16.0                    ;gb center degrees
   alpha = 105.0                 ;orientation of the ellipse, degrees
   rmajor = 8.0                 ;semi-major axis
   rminor = 1.5                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,2:3],2)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif





; Horn3, 13ghz, cor

   horn = 3
   freq = 1
   ihorn = horn-1
   
;(1)                                ;Ellipse parameters
   x0 = 343.0                   ;gl center degrees
   y0 = 69.0                    ;gb center degrees
   alpha = 43.0                 ;orientation of the ellipse, degrees
   rmajor =14.0                 ;semi-major axis
   rminor = 3.0                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:1],2)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif

;(2)
   horn = 3
   freq = 1
   ihorn = horn-1
   
   x0 = 214.0                   ;gl center degrees
   y0 = 16.5                    ;gb center degrees
   alpha = 105.0                ;orientation of the ellipse, degrees
   rmajor =11.0                 ;semi-major axis
   rminor = 1.5                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:1],2)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif

;(3)
   horn = 3
   freq = 1
   ihorn = horn-1
   
   x0 = 1.0                     ;gl center degrees
   y0 = 26.0                    ;gb center degrees
   alpha = 52.0                 ;orientation of the ellipse, degrees
   rmajor = 6.0                 ;semi-major axis
   rminor = 1.0                 ;semi-minor axis
   iaz = 2                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
   ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:1],2)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif

; Horn3, 13ghz, uncor
   
;(1)
   horn = 3
   freq = 1
   ihorn = horn-1
                                ;Ellipse parameters
   x0 = 343.0                   ;gl center degrees
   y0 = 69.0                    ;gb center degrees
   alpha = 43.0                 ;orientation of the ellipse, degrees
   rmajor =14.0                 ;semi-major axis
   rminor = 3.0                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
   ;Flag
   if n_elements(cut) ge 1 then begin
      sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,2:3],2)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif


;(2)
   horn = 3
   freq = 1
   ihorn = horn-1
   
   x0 = 214.3                   ;gl center degrees
   y0 = 17.5                    ;gb center degrees
   alpha = 105.0                ;orientation of the ellipse, degrees
   rmajor =11.0                 ;semi-major axis
   rminor = 1.5                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
   ;Flag
   if n_elements(cut) ge 1 then begin
      sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,2:3],2)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif


;(3)
   horn = 3
   freq = 1
   ihorn = horn-1
   
   x0 = 1.0                     ;gl center degrees
   y0 = 26.0                    ;gb center degrees
   alpha = 52.0                 ;orientation of the ellipse, degrees
   rmajor = 6.0                 ;semi-major axis
   rminor = 1.0                 ;semi-minor axis
   iaz = 2                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
   ;Flag
   if n_elements(cut) ge 1 then begin
      sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,2:3],2)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif


   flageado = 1
   
endif

;-----------------
; NOMINAL40, period 5.
   if (campo eq "NOMINAL40" and period eq 5) then begin
      print," Flagging ",campo
      
;03/05/19 Additional ellipse 
; Horn3, 13ghz, cor/uncor

   ;(1)
   horn = 3
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 354.19                ;gl center degrees
   y0 = 32.75                  ;gb center degrees
   alpha =  50.0                 ;orientation of the ellipse, degrees
   rmajor = 6.0                 ;semi-major axis
   rminor = 1.0                 ;semi-minor axis
   iaz = 2                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
   ;Flag
   if n_elements(cut) ge 1 then begin
      sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif
;03/05/19 End additional ellipse 

   
; Horn4, 19ghz, cor (the same as uncorr but little bigger)

;(1)
      horn = 4
      freq = 1
      ihorn = horn-1
      
                                ;Ellipse parameters
      x0 = 36.5                ;gl center degrees
      y0 = -50.0               ;gb center degrees
      alpha = 70.0             ;orientation of the ellipse, degrees
      rmajor = 28.0             ;semi-major axis
      rminor = 2.0             ;semi-minor axis
      iaz = 2                    ;azimuth minor (1) or major (2) than 180
      
      cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
      
                                ;Flag
      if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
         ichan = reform(mfi.ic[ihorn,freq,0:1],2)
         for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
      endif


      
; Horn4, 19ghz, uncor (the same as corr but little smaller)

;(1)
      horn = 4
      freq = 1
      ihorn = horn-1
      
                                ;Ellipse parameters
      x0 = 38.0                ;gl center degrees
      y0 = -51.1               ;gb center degrees
      alpha = 70.0             ;orientation of the ellipse, degrees
      rmajor = 25.0             ;semi-major axis
      rminor = 1.5             ;semi-minor axis
      iaz = 2                    ;azimuth minor (1) or major (2) than 180
      
      cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
      
                                ;Flag
      if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
         ichan = reform(mfi.ic[ihorn,freq,2:3],2)
         for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
      endif

   flageado = 1

endif


;---------------------------------------
; NOMINAL50, period 2.
if (campo eq "NOMINAL50" and period eq 2) then begin
   print," Flagging ",campo

   
; Horn2, 19ghz, cor/uncor
   
;(1)
   horn = 2
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 157.6                   ;gl center degrees
   y0 = 23.5                    ;gb center degrees
   alpha = 110.0                ;orientation of the ellipse, degrees
   rmajor = 2.5                 ;semi-major axis
   rminor = 1.0                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif


;(2)
   horn = 2
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 161.6                   ;gl center degrees
   y0 = 38.6                    ;gb center degrees
   alpha = 45.0                 ;orientation of the ellipse, degrees
   rmajor = 4.0                 ;semi-major axis
   rminor = 2.0                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif

;(3)
   horn = 2
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 91.9                    ;gl center degrees
   y0 = 50.5                    ;gb center degrees
   alpha = 150.0                ;orientation of the ellipse, degrees
   rmajor = 3.5                 ;semi-major axis
   rminor = 1.0                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif


;(4)
   horn = 2
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 105.5                   ;gl center degrees
   y0 = 8.5                     ;gb center degrees
   alpha =  50.0                ;orientation of the ellipse, degrees
   rmajor = 3.0                 ;semi-major axis
   rminor = 0.8                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif

 
;(5)
   horn = 2
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 122.4                   ;gl center degrees
   y0 = -5.1                    ;gb center degrees
   alpha =   0.0                ;orientation of the ellipse, degrees
   rmajor = 1.5                 ;semi-major axis
   rminor = 1.5                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif


;(6) 
   horn = 2
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 9.9                     ;gl center degrees
   y0 = 33.8                    ;gb center degrees
   alpha =   60.0               ;orientation of the ellipse, degrees
   rmajor = 5.0                 ;semi-major axis
   rminor = 1.0                 ;semi-minor axis
   iaz = 2                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif


;(7)
   horn = 2
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 207.5                   ;gl center degrees
   y0 = -14.3                   ;gb center degrees
   alpha =   110.0              ;orientation of the ellipse, degrees
   rmajor = 4.5                 ;semi-major axis
   rminor = 1.0                 ;semi-minor axis
   iaz = 2                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif

 ;(8) ADDED on 07/05/19
   horn = 2
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 147.5                   ;gl center degrees
   y0 = 51.4                   ;gb center degrees
   alpha =  0.0                ;orientation of the ellipse, degrees
   rmajor = 2.5                 ;semi-major axis
   rminor = 1.5                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif

   
; Horn4, 19ghz, cor/uncor

 ;(1)
   horn = 4
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 152.76                  ;gl center degrees
   y0 = 20.71                   ;gb center degrees
   alpha =   125.0              ;orientation of the ellipse, degrees
   rmajor = 3.0                 ;semi-major axis
   rminor = 1.0                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif

   
 ;(2)
   horn = 4
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 155.0                   ;gl center degrees
   y0 = 32.2                    ;gb center degrees
   alpha =   100.0              ;orientation of the ellipse, degrees
   rmajor = 2.5                 ;semi-major axis
   rminor = 1.0                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif


 ;(3)
   horn = 4
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 145.5                   ;gl center degrees
   y0 = 34.5                    ;gb center degrees
   alpha = 77.0                 ;orientation of the ellipse, degrees
   rmajor = 6.0                 ;semi-major axis
   rminor = 1.0                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif

 ;(4)
   horn = 4
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 147.4                   ;gl center degrees
   y0 = 44.9                    ;gb center degrees
   alpha = 35.0                 ;orientation of the ellipse, degrees
   rmajor = 3.0                 ;semi-major axis
   rminor = 1.0                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif

   
  ;(5)
   horn = 4
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 106.8                   ;gl center degrees
   y0 = 51.6                    ;gb center degrees
   alpha = -10.0                ;orientation of the ellipse, degrees
   rmajor = 5.0                 ;semi-major axis
   rminor = 2.0                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif


  ;(6)
   horn = 4
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 16.1                    ;gl center degrees
   y0 = 19.3                    ;gb center degrees
   alpha = 60.0                 ;orientation of the ellipse, degrees
   rmajor = 4.5                 ;semi-major axis
   rminor = 1.0                 ;semi-minor axis
   iaz = 2                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif



  ;(7)
   horn = 4
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 202.9                   ;gl center degrees
   y0 = -30.1                   ;gb center degrees
   alpha = 110.0                ;orientation of the ellipse, degrees
   rmajor = 4.0                 ;semi-major axis
   rminor = 0.8                 ;semi-minor axis
   iaz = 2                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif

;03/05/19 Additional ellipses
; Horn4, 19ghz, cor/uncor

    ;(8)
   horn = 4
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 126.8                 ;gl center degrees
   y0 = -0.6                  ;gb center degrees
   alpha = -10                  ;orientation of the ellipse, degrees
   rmajor = 2.5                 ;semi-major axis
   rminor =  2.0                ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif  



    ;(9)
   horn = 4
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 127.14                 ;gl center degrees
   y0 = 6.2                  ;gb center degrees
   alpha = -20                  ;orientation of the ellipse, degrees
   rmajor = 8.0                 ;semi-major axis
   rminor =  1.0                ;semi-minor axis
   iaz = 2                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif
;03/05/19 End additional ellipses
   
   flageado = 1

endif




;-----------------
; NOMINAL50, period 5.

if (campo eq "NOMINAL50" and period eq 5) then begin
   print," Flagging ",campo

; Horn2, 19ghz, cor/uncor

 ;(1)
   horn = 2
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 131.3                   ;gl center degrees
   y0 = -3.5                    ;gb center degrees
   alpha =   0.0                ;orientation of the ellipse, degrees
   rmajor = 3.0                 ;semi-major axis
   rminor = 1.0                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif



 ;(2)
   horn = 2
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 94.1                    ;gl center degrees
   y0 = 16.1                    ;gb center degrees
   alpha =  80.0                ;orientation of the ellipse, degrees
   rmajor = 2.0                 ;semi-major axis
   rminor = 1.0                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif


 ;(3)
   horn = 2
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 144.5                   ;gl center degrees
   y0 = 13.8                    ;gb center degrees
   alpha =  130.0               ;orientation of the ellipse, degrees
   rmajor = 3.0                 ;semi-major axis
   rminor = 0.5                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif


 ;(4)
   horn = 2
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 139.15                  ;gl center degrees
   y0 = -7.75                   ;gb center degrees
   alpha =  130.0               ;orientation of the ellipse, degrees
   rmajor = 4.0                 ;semi-major axis
   rminor = 2.0                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif

   
 ;(5)
   horn = 2
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 128.6                   ;gl center degrees
   y0 = 4.0                     ;gb center degrees
   alpha =  170.0               ;orientation of the ellipse, degrees
   rmajor = 2.4                 ;semi-major axis
   rminor = 0.8                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif

 ;(6)
   horn = 2
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 4.0                     ;gl center degrees
   y0 = 39.7                    ;gb center degrees
   alpha =  55.0                ;orientation of the ellipse, degrees
   rmajor = 3.5                 ;semi-major axis
   rminor = 0.5                 ;semi-minor axis
   iaz = 2                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif

 ;(7)
   horn = 2
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 202.2                   ;gl center degrees
   y0 = -24.2                   ;gb center degrees
   alpha =  115.0               ;orientation of the ellipse, degrees
   rmajor = 3.5                 ;semi-major axis
   rminor = 1.0                 ;semi-minor axis
   iaz = 2                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif

 ;(8)



; Horn3, 11ghz, cor/uncor

 ;(1)
   horn = 3
   freq = 0
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 346.0                   ;gl center degrees
   y0 = 67.8                    ;gb center degrees
   alpha =  50.0                ;orientation of the ellipse, degrees
   rmajor = 8.0                 ;semi-major axis
   rminor = 1.5                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif

;03/05/19 Additional ellipses
 ; Horn3, 11ghz, cor/uncor

;(2)
   horn = 3
   freq = 0
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 263.52                  ;gl center degrees
   y0 = 58.29                  ;gb center degrees
   alpha =  149.0                 ;orientation of the ellipse, degrees
   rmajor = 13.0                 ;semi-major axis
   rminor = 1.5                 ;semi-minor axis
   iaz = 2                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif

;(3)
   horn = 3
   freq = 0
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 211.69                  ;gl center degrees
   y0 = 1.86                  ;gb center degrees
   alpha =  115.0                 ;orientation of the ellipse, degrees
   rmajor = 15.0                 ;semi-major axis
   rminor = 1.5                 ;semi-minor axis
   iaz = 2                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif

;03/05/19 End additional ellipses
   


; Horn3, 13ghz, cor/uncor

;(1)
   horn = 3
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 346.0                   ;gl center degrees
   y0 = 67.8                    ;gb center degrees
   alpha =  48.0                ;orientation of the ellipse, degrees
   rmajor = 13.0                ;semi-major axis
   rminor = 2.0                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif


;(2)
   horn = 3
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 205.8                   ;gl center degrees
   y0 = -2.3                    ;gb center degrees
   alpha =  105.0               ;orientation of the ellipse, degrees
   rmajor = 13.0                ;semi-major axis
   rminor = 2.0                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif
   


   
; Horn4, 19ghz, cor/uncor
;(1)
   horn = 4
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 100.4                   ;gl center degrees
   y0 = 13.8                    ;gb center degrees
   alpha =  50.0                ;orientation of the ellipse, degrees
   rmajor = 2.0                 ;semi-major axis
   rminor = 1.0                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif


;(2)
   horn = 4
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 154.0                   ;gl center degrees
   y0 = 34.8                    ;gb center degrees
   alpha =  90.0                ;orientation of the ellipse, degrees
   rmajor = 3.0                 ;semi-major axis
   rminor = 1.0                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif
   
;03/05/19 Additional ellipses
; Horn4, 19ghz, cor/uncor
  
;(3)
   horn = 4
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 140.50                 ;gl center degrees
   y0 = 1.62                  ;gb center degrees
   alpha = 90.0                  ;orientation of the ellipse, degrees
   rmajor = 3.5                 ;semi-major axis
   rminor =  1.5                ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif


;(4)
   horn = 4
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 125.8                 ;gl center degrees
   y0 = -0.37                  ;gb center degrees
   alpha = 0.0                  ;orientation of the ellipse, degrees
   rmajor = 3.5                 ;semi-major axis
   rminor =  1.7                ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif 


;(5)
   horn = 4
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 10.19                 ;gl center degrees
   y0 = 27.52                  ;gb center degrees
   alpha = 60.0                  ;orientation of the ellipse, degrees
   rmajor = 6.0                 ;semi-major axis
   rminor =  1.5                ;semi-minor axis
   iaz = 2                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif


;(6)
   horn = 4
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 197.04                 ;gl center degrees
   y0 = -39.05                  ;gb center degrees
   alpha = 125.0                  ;orientation of the ellipse, degrees
   rmajor = 6.0                 ;semi-major axis
   rminor =  1.5                ;semi-minor axis
   iaz = 2                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif

   
;(7)
   horn = 4
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 19.63                 ;gl center degrees
   y0 = -2.62                  ;gb center degrees
   alpha = 63.0                  ;orientation of the ellipse, degrees
   rmajor = 10.0                 ;semi-major axis
   rminor =  0.6                ;semi-minor axis
   iaz = 2                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif

   
;(8)
   horn = 4
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 41.27                 ;gl center degrees
   y0 = -44.37                  ;gb center degrees
   alpha = 64.0                  ;orientation of the ellipse, degrees
   rmajor = 10.0                 ;semi-major axis
   rminor =  0.6                ;semi-minor axis
   iaz = 2                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif

;(9)
   horn = 4
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 98.64                 ;gl center degrees
   y0 = 32.09                  ;gb center degrees
   alpha = 90.0                  ;orientation of the ellipse, degrees
   rmajor = 5.0                 ;semi-major axis
   rminor =  1.0                ;semi-minor axis
   iaz = 2                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif

;(10)
   horn = 4
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 114.47                 ;gl center degrees
   y0 = 7.36                  ;gb center degrees
   alpha = 30.0                  ;orientation of the ellipse, degrees
   rmajor = 4.0                 ;semi-major axis
   rminor =  1.0                ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif
   
;03/05/19 End additional ellipses
   
   flageado = 1
   
endif



;-----------------
; NOMINAL50, period 6.

if (campo eq "NOMINAL50" and period eq 6) then begin
   print," Flagging ",campo

; Horn2, 19ghz, cor/uncor
;(1)
   horn = 2
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 158.3                   ;gl center degrees
   y0 = 24.5                    ;gb center degrees
   alpha =  120.0               ;orientation of the ellipse, degrees
   rmajor = 2.0                 ;semi-major axis
   rminor = 1.0                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif


;(2)
   horn = 2
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 122.0                   ;gl center degrees
   y0 = 50.5                    ;gb center degrees
   alpha =  0.0                 ;orientation of the ellipse, degrees
   rmajor = 2.0                 ;semi-major axis
   rminor = 1.0                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif



;(3)
   horn = 2
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 98.20                   ;gl center degrees
   y0 = 40.8                    ;gb center degrees
   alpha =  125.0               ;orientation of the ellipse, degrees
   rmajor = 2.5                 ;semi-major axis
   rminor = 0.5                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif



;(4)
   horn = 2
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 92.30                   ;gl center degrees
   y0 = 25.9                    ;gb center degrees
   alpha =  85.0                ;orientation of the ellipse, degrees
   rmajor = 1.5                 ;semi-major axis
   rminor = 1.0                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif


;(5)
   horn = 2
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 119.5                   ;gl center degrees
   y0 = 3.9                     ;gb center degrees
   alpha =  5.0                 ;orientation of the ellipse, degrees
   rmajor = 2.0                 ;semi-major axis
   rminor = 0.5                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif


;(6)
   horn = 2
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 118.8                   ;gl center degrees
   y0 = -4.0                    ;gb center degrees
   alpha =  40.0                ;orientation of the ellipse, degrees
   rmajor = 2.5                 ;semi-major axis
   rminor = 0.8                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif



;(7)
   horn = 2
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 120.3                   ;gl center degrees
   y0 = -10.6                   ;gb center degrees
   alpha =  120.0               ;orientation of the ellipse, degrees
   rmajor = 2.5                 ;semi-major axis
   rminor = 1.0                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif



;(8)
   horn = 2
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 1.8                     ;gl center degrees
   y0 = 41.6                    ;gb center degrees
   alpha =  50.0                ;orientation of the ellipse, degrees
   rmajor = 3.0                 ;semi-major axis
   rminor = 1.0                 ;semi-minor axis
   iaz = 2                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif




;(9)
   horn = 2
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 200.6                   ;gl center degrees
   y0 = -26.7                   ;gb center degrees
   alpha =  115.0               ;orientation of the ellipse, degrees
   rmajor = 3.8                 ;semi-major axis
   rminor = 1.0                 ;semi-minor axis
   iaz = 2                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif




;(10)
   horn = 2
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 73.4                   ;gl center degrees
   y0 = -66.0                   ;gb center degrees
   alpha =  40.0               ;orientation of the ellipse, degrees
   rmajor = 8.0                 ;semi-major axis
   rminor = 0.5                 ;semi-minor axis
   iaz = 2                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif



; Horn3, 11ghz, cor/uncor
;(1)
   horn = 3
   freq = 0
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 340.0                    ;gl center degrees
   y0 = 69.60                   ;gb center degrees
   alpha =  45.00               ;orientation of the ellipse, degrees
   rmajor = 6.0                 ;semi-major axis
   rminor = 1.5                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif
   
;03/05/19 End additional ellipses
; Horn3, 11ghz, cor/uncor

;(2)
   horn = 3
   freq = 0
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 263.52                  ;gl center degrees
   y0 = 58.29                  ;gb center degrees
   alpha =  149.0                 ;orientation of the ellipse, degrees
   rmajor = 13.0                 ;semi-major axis
   rminor = 1.5                 ;semi-minor axis
   iaz = 2                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif

;(3)
   horn = 3
   freq = 0
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 208.15                  ;gl center degrees
   y0 = -5.04                  ;gb center degrees
   alpha =  115.0                 ;orientation of the ellipse, degrees
   rmajor = 15.0                 ;semi-major axis
   rminor = 1.5                 ;semi-minor axis
   iaz = 2                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif

;03/05/19 End additional ellipses

; Horn3, 13ghz, cor/uncor
;(1)
   horn = 3
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 340.0                    ;gl center degrees
   y0 = 69.60                   ;gb center degrees
   alpha =  42.00               ;orientation of the ellipse, degrees
   rmajor = 13.0                ;semi-major axis
   rminor = 2.5                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif

;(2)
   horn = 3
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 205.0                    ;gl center degrees
   y0 = -3.8                    ;gb center degrees
   alpha =  105.0               ;orientation of the ellipse, degrees
   rmajor = 10.0                ;semi-major axis
   rminor = 1.5                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif





; Horn4, 17ghz, cor/uncor
;(1)
   horn = 4
   freq = 0
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 195.2                    ;gl center degrees
   y0 = -41.5                   ;gb center degrees
   alpha =  120.0               ;orientation of the ellipse, degrees
   rmajor = 4.0                 ;semi-major axis
   rminor = 1.0                 ;semi-minor axis
   iaz = 2                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif



; Horn4, 19ghz, cor/uncor
;(1)
   horn = 4
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 153.2                    ;gl center degrees
   y0 = 21.11                   ;gb center degrees
   alpha =  120.0               ;orientation of the ellipse, degrees
   rmajor = 2.0                 ;semi-major axis
   rminor = 1.0                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif


;(2)
   horn = 4
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 97.5                     ;gl center degrees
   y0 = 21.6                    ;gb center degrees
   alpha =  80.0                ;orientation of the ellipse, degrees
   rmajor = 2.5                 ;semi-major axis
   rminor = 1.5                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif


;(3)
   horn = 4
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 114.3                    ;gl center degrees
   y0 = 1.7                     ;gb center degrees
   alpha =  40.0                ;orientation of the ellipse, degrees
   rmajor = 3.0                 ;semi-major axis
   rminor = 1.0                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif



;(4)
   horn = 4
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 126.7                    ;gl center degrees
   y0 = -3.5                    ;gb center degrees
   alpha =  135.0               ;orientation of the ellipse, degrees
   rmajor = 3.0                 ;semi-major axis
   rminor = 1.5                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif


;(5)
   horn = 4
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 195.1                    ;gl center degrees
   y0 = -41.5                   ;gb center degrees
   alpha =  120.0               ;orientation of the ellipse, degrees
   rmajor = 5.0                 ;semi-major axis
   rminor = 1.0                 ;semi-minor axis
   iaz = 2                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif

   flageado  = 1
   
endif

;-----------------
; NOMINAL60, period 5.

if (campo eq "NOMINAL60" and period eq 5) then begin
   print," Flagging ",campo
   
;03/05/19 Additional ellipses   
; Horn4, 19ghz, cor/uncor

;(1)
   horn = 4
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 87.21                  ;gl center degrees
   y0 = 17.7                  ;gb center degrees
   alpha =  90.0                 ;orientation of the ellipse, degrees
   rmajor = 3.5                 ;semi-major axis
   rminor = 1.7                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif

;(2)
   horn = 4
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 75.6                  ;gl center degrees
   y0 = 21.17                  ;gb center degrees
   alpha =  90.0                 ;orientation of the ellipse, degrees
   rmajor = 3.5                 ;semi-major axis
   rminor = 3.5                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif


;(3)
   horn = 4
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 123.34                  ;gl center degrees
   y0 = 58.43                  ;gb center degrees
   alpha =  180.0                 ;orientation of the ellipse, degrees
   rmajor = 2.5                 ;semi-major axis
   rminor = 1.0                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif

;(4)
   horn = 4
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 156.0                  ;gl center degrees
   y0 = -14.5                  ;gb center degrees
   alpha =  22.0                 ;orientation of the ellipse, degrees
   rmajor = 9.0                 ;semi-major axis
   rminor = 1.2                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif

;(5)
   horn = 4
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 151.16                  ;gl center degrees
   y0 = -16.04                  ;gb center degrees
   alpha =  22.0                 ;orientation of the ellipse, degrees
   rmajor = 9.0                 ;semi-major axis
   rminor = 1.2                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif
;03/05/19 End additional ellipses
   
   flageado  = 1
     
endif

 ;-----------------
; NOMINAL60, period 6.

if (campo eq "NOMINAL60" and period eq 6) then begin
   print," Flagging ",campo
   
; Horn2, 19ghz, cor/uncor
   
;(1)
   horn = 2
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 64.4                    ;gl center degrees
   y0 = 36.0                    ;gb center degrees
   alpha =  0.0                 ;orientation of the ellipse, degrees
   rmajor = 3.5                 ;semi-major axis
   rminor = 3.0                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif

   ;(2)
   horn = 2
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 81.2                    ;gl center degrees
   y0 = 22.1                    ;gb center degrees
   alpha =  90.0                ;orientation of the ellipse, degrees
   rmajor = 3.0                 ;semi-major axis
   rminor = 1.5                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif


      ;(3)
   horn = 2
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 190.8                   ;gl center degrees
   y0 = -36.2                   ;gb center degrees
   alpha =  120.0               ;orientation of the ellipse, degrees
   rmajor = 3.0                 ;semi-major axis
   rminor = 0.5                 ;semi-minor axis
   iaz = 2                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif


   
;(4) Added on 07/05/19
   horn = 2
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 146.85                    ;gl center degrees
   y0 = 57.48                    ;gb center degrees
   alpha =  45.0                 ;orientation of the ellipse, degrees
   rmajor = 3.0                 ;semi-major axis
   rminor = 1.0                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif

; Horn4, 19ghz, cor/uncor
   
;(1)
   horn = 4
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 75.7                    ;gl center degrees
   y0 = 31.9                    ;gb center degrees
   alpha =  100.0                 ;orientation of the ellipse, degrees
   rmajor = 4.0                 ;semi-major axis
   rminor = 3.0                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif

   
;(2)
   horn = 4
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 278.9                   ;gl center degrees
   y0 = 68.0                    ;gb center degrees
   alpha =  160.0               ;orientation of the ellipse, degrees
   rmajor = 13.0                ;semi-major axis
   rminor = 2.5                 ;semi-minor axis
   iaz = 2                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif

;(3)
   horn = 4
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 219.6                   ;gl center degrees
   y0 = 30.4                    ;gb center degrees
   alpha =  112.0               ;orientation of the ellipse, degrees
   rmajor = 13.0                ;semi-major axis
   rminor = 1.5                 ;semi-minor axis
   iaz = 2                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif
   
;03/05/19 Additional ellipses
; Horn4, 19ghz, cor/uncor

   ;(4)
   horn = 4
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 180.2                  ;gl center degrees
   y0 = 41.1                  ;gb center degrees
   alpha =  140.0                 ;orientation of the ellipse, degrees
   rmajor = 8.0                 ;semi-major axis
   rminor = 2.0                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif

   ;(5)
   horn = 4
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 139.24                  ;gl center degrees
   y0 = 56.33                  ;gb center degrees
   alpha =  30.0                 ;orientation of the ellipse, degrees
   rmajor = 3.0                 ;semi-major axis
   rminor = 1.0                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif


   ;(6)
   horn = 4
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 84.32                 ;gl center degrees
   y0 = 25.9                  ;gb center degrees
   alpha =  90.0                 ;orientation of the ellipse, degrees
   rmajor = 4.0                 ;semi-major axis
   rminor = 1.5                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif   
;03/05/19 End additional ellipses

;(7) ADDED on 07/05/19
   horn = 4
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 114.4                    ;gl center degrees
   y0 = 60.0                    ;gb center degrees
   alpha = -40.0                 ;orientation of the ellipse, degrees
   rmajor = 1.5                 ;semi-major axis
   rminor = 1.0                 ;semi-minor axis
   iaz = 2                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif  
   flageado = 1
   
endif


;-----------------
; NOMINAL65, period 6.

if (campo eq "NOMINAL65" and period eq 6) then begin
   print," Flagging ",campo

  ; Horn4, 19ghz, cor/uncor
;(1)
   horn = 4
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 199.6                   ;gl center degrees
   y0 = -19.5                    ;gb center degrees
   alpha =  120.0               ;orientation of the ellipse, degrees
   rmajor = 6.0                ;semi-major axis
   rminor = 0.8                 ;semi-minor axis
   iaz = 2                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif

   
;03/05/19 Additional ellipses
   ;(2)
   horn = 4
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 158.86                  ;gl center degrees
   y0 = 54.12                   ;gb center degrees
   alpha =  70.0                 ;orientation of the ellipse, degrees
   rmajor = 3.5                 ;semi-major axis
   rminor = 1.5                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif

   ;(3)
   horn = 4
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 77.34                  ;gl center degrees
   y0 = 25.70                  ;gb center degrees
   alpha =  100.0                 ;orientation of the ellipse, degrees
   rmajor = 4.0                 ;semi-major axis
   rminor = 1.0                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif
;03/05/19 End additional ellipses

;(4) ADDED on 07/05/19
   horn = 4
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 165.28                  ;gl center degrees
   y0 = 42.32                    ;gb center degrees
   alpha =  65.0                ;orientation of the ellipse, degrees
   rmajor = 2.5                ;semi-major axis
   rminor = 1.0                 ;semi-minor axis
   iaz = 2                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
         sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif
   
      flageado = 1

endif

;-----------------
; NOMINAL70, period 6.
if (campo eq "NOMINAL70" and period eq 6) then begin
   print," Flagging ",campo

;03/05/19 Additional ellipses
; Horn4, 19ghz, cor/uncor

   ;(1)
   horn = 4
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 172.69                   ;gl center degrees
   y0 = 53.87                   ;gb center degrees
   alpha =  65.0                 ;orientation of the ellipse, degrees
   rmajor = 3.5                 ;semi-major axis
   rminor = 1.5                 ;semi-minor axis
   iaz = 1                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
   ;Flag
   if n_elements(cut) ge 1 then begin
      sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif

   ;(2)
   horn = 4
   freq = 1
   ihorn = horn-1
   
                                ;Ellipse parameters
   x0 = 142.34                  ;gl center degrees
   y0 = 69.57                  ;gb center degrees
   alpha =  45.0                 ;orientation of the ellipse, degrees
   rmajor = 3.5                 ;semi-major axis
   rminor = 1.5                 ;semi-minor axis
   iaz = 2                       ;azimuth minor (1) or major (2) than 180
   
   cut = FLAG_ELLIPSE(X0=x0, Y0=y0, ALPHA=alpha, RMAJOR=rmajor, RMINOR=rminor, GL=gl[ihorn,*], GB=gb[ihorn,*], IAZ=iaz, AZ_INPUT = az)
   
                                ;Flag
   if n_elements(cut) ge 1 then begin
      sancho_check_activebit_inflags, flag, bit, bitmasc
      ichan = reform(mfi.ic[ihorn,freq,0:3],4)
      for k=0,n_elements(ichan)-1 do flag[ichan[k],cut] = flag[ichan[k],cut] + one*(1-bitmasc[ichan[k],cut])
   endif
;03/05/19 End additional ellipses

   flageado = 1
   
endif



;******************************************************
; TRANSIENTS and PLANETS. Flags by JARM,BRG 29/Aug/2018. Removing QRT1=VENUS, r1=12deg around Venus. Also flagging horn 1 (July/2019)
;
jd0          = jd[0]
jd1          = 2457793.d0 ; 2017-2-8
jd_inter     = jd1 + 55.d0
jd2          = 2457859.d0 ; 2017-4-15
MINDISTVENUS1 =  5.0 ; deg
MINDISTVENUS2 = 12.0 ; deg
if ( (campo eq "NOMINAL65" or campo eq "NOMINAL70" or campo eq "NOMINAL50") and period eq 6 and jd0 ge jd1 and jd0 le jd2 and keepplanets eq 0) then begin
   MINDISTVENUS = MINDISTVENUS1
   if jd0 ge jd_inter then MINDISTVENUS = MINDISTVENUS2

   print," Flagging ",campo," in period ",period
   print,"   --> for QRT1, using radius = ",MINDISTVENUS 

   sancho_check_activebit_inflags, flag, pbit, bitmasc

   ; Horns 1,2,3,4
   for horn=1,4 do begin
	ihorn = horn-1
	planet_coords, jd0, ravenus, decvenus, planet="venus",/jd
	glactc,ravenus, decvenus, 2000.0, glvenus, gbvenus, 1, /degree
	
	d = distancia(gl[ihorn,*],gb[ihorn,*],glvenus[0],gbvenus[0],/degree)
	cut = where(d le MINDISTVENUS, ncut)

	if ncut ge 1 then begin
		ichan = reform(mfi.ic[ihorn,*,*],8)
		for k=0,7 do flag[ichan[k],cut] = flag[ichan[k],cut] + onep * (1-bitmasc[ichan[k],cut])
	endif
   endfor
   flageado = 1
endif


; TRANSIENT2: Jupiter. It is clearly seen inside period6, nominal35 data. From March 12th, 2020, I flag it always
jd0        = jd[0]
MINDISTJUP = 2.0 ; deg

if keepplanets eq 0 then begin
print," Flagging ",campo," in period ",period
print,"   --> for Jupiter, using radius = ",MINDISTJUP

sancho_check_activebit_inflags, flag, pbit, bitmasc
for horn=1,4 do begin
	ihorn = horn-1
	planet_coords, jd0, ra_jup, dec_jup, planet="jupiter",/jd
	glactc,ra_jup, dec_jup, 2000.0, gl_jup, gb_jup, 1, /degree
	d = distancia(gl[ihorn,*],gb[ihorn,*],gl_jup[0],gb_jup[0],/degree)
	cut = where(d le MINDISTJUP, ncut)
	if ncut ge 1 then begin
		ichan = reform(mfi.ic[ihorn,*,*],8)
		for k=0,7 do flag[ichan[k],cut] = flag[ichan[k],cut] + onep*(1-bitmasc[ichan[k],cut])
	endif
endfor
flageado = 1
endif

; TRANSIENT3: Venus. Added march 12th, 2020
jd0        = jd[0]
MINDISTPLA = 2.0 ; deg

if keepplanets eq 0 then begin
print," Flagging ",campo," in period ",period
print,"   --> for Venus, using radius = ",MINDISTPLA

sancho_check_activebit_inflags, flag, pbit, bitmasc
for horn=1,4 do begin
	ihorn = horn-1
	planet_coords, jd0, ra_pla, dec_pla, planet="venus",/jd
	glactc,ra_pla, dec_pla, 2000.0, gl_pla, gb_pla, 1, /degree
	d = distancia(gl[ihorn,*],gb[ihorn,*],gl_pla[0],gb_pla[0],/degree)
	cut = where(d le MINDISTPLA, ncut)
	if ncut ge 1 then begin
		ichan = reform(mfi.ic[ihorn,*,*],8)
		for k=0,7 do flag[ichan[k],cut] = flag[ichan[k],cut] + onep*(1-bitmasc[ichan[k],cut])
	endif
endfor
flageado = 1
endif


; TRANSIENT4: Mars. Added March 12th, 2020
jd0        = jd[0]
MINDISTPLA = 2.0 ; deg

if keepplanets eq 0 then begin
print," Flagging ",campo," in period ",period
print,"   --> for Mars, using radius = ",MINDISTPLA

sancho_check_activebit_inflags, flag, pbit, bitmasc
for horn=1,4 do begin
	ihorn = horn-1
	planet_coords, jd0, ra_pla, dec_pla, planet="mars",/jd
	glactc,ra_pla, dec_pla, 2000.0, gl_pla, gb_pla, 1, /degree
	d = distancia(gl[ihorn,*],gb[ihorn,*],gl_pla[0],gb_pla[0],/degree)
	cut = where(d le MINDISTPLA, ncut)
	if ncut ge 1 then begin
		ichan = reform(mfi.ic[ihorn,*,*],8)
		for k=0,7 do flag[ichan[k],cut] = flag[ichan[k],cut] + onep*(1-bitmasc[ichan[k],cut])
	endif
endfor
flageado = 1
endif




; End of subroutine
if flageado then begin
	print," APPLY_FLAGS_PER_FIELD. Flagging applied in field ",campo,". Version = ", version
endif else begin
	print," APPLY_FLAGS_PER_FIELD. No flagging applied. "
endelse

END



