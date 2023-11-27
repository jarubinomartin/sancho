;   1/Feb/2013   J.A.Rubino
;
;   For a given jd, AZ and EL vector, and a given horn of the MFI
;   instrument, it returns RA,DEC (or GL,GB)
; 
;   HISTORY:
;      1/Februa/2013 - original version
;      23/March/2013 - implemented two versions for the precession-nutation-aberration corrections: 
;                      one based on IDL routines, and another one
;                      based on TPOINT equations. 
;      30/March/2013 - Includes pointing model corrections, and
;                      apparent-to-J2000 transformation.  The keyword
;                      OLD runs the old version (no corrections). The
;                      pointing  correction is NOT implemented yet in
;                      the parallactic angle computation.
;      04/July/2013  - Focal plane position of each horn is now computed
;                      exactly. A typo has been corrected in the parallactic
;                      angle computation in galactic coordinates. 
;      19/Dec/2013   - include the option of having a different pmodel per horn.
;      27/Jan/2014   - typo corrected in the /standard version of the code. Now
;                      uses fp2sky (correct transformation).
;      24/Jul/2014   - corrects the calculation of parallactic angle. It now
;                      includes the pointing correction also. Checked both for
;                      RADEC and Galactic.
;      10/Jul/2015   - on output, also the (corrected) true local coordinates (AZ,EL) for the
;                      requested horn are given.
;
;
;-
;


PRO sancho_mfi_get_skycoords, jd, AZ, EL, ihorn, RA, DEC, POSANG, AZ_true, EL_true, galactic=galactic, $
                              old=old, standard=standard, verbose=verbose, idl=idl

if not keyword_set(galactic) then galactic = 0  
if not keyword_set(old) then old = 0  
if not keyword_set(standard) then standard = 0  
if not keyword_set(verbose) then verbose = 0  
if not keyword_set(idl) then idl = 0  


; Info
if n_params() lt 4 then begin
   print,""
   print,"  Syntax -- sancho_mfi_get_skycoords, jd, AZ, EL, ihorn, RA, DEC, PARANGLE, AZ_c, EL_c [, galactic=] "
   print,""
   print,"    ihorn = 1, 2, 3 or 4  for the MFI instrument (0=center of the focal plane). "
   print,"    Input AZ and EL values should be given in degrees. Output is RA, DEC (J2000) in degrees. "
   print,"    If /GALACTIC keyword is activated, then the output coordinates are GL, GB (J2000)in degrees. "
   print,"    PARANGLE is the parallactic angle in degrees"
   print,""
   return
endif

; MFI info and constants
info_quijote,mfi=mfi

d2r   = !DPI/180.d0
r2d   = 1.d0 / d2r
h2d   = 15.d0

; Identify horn
if (ihorn le mfi.nhorn and ihorn ge 1) then begin
   ihorn = fix(ihorn)
   print," --> Computing skycoords for MFI horn number ",ihorn
endif else begin
   if (ihorn eq 0) then begin
      ihorn = 0
      print," --> Computing skycoords for the center of MFI focal plane. "
   endif else begin
      print," --> Number of horns of the MFI instrument: ",mfi.nhorn
      return
   endelse
endelse


; Compute offset for that horn.
xpos = 0.0
ypos = 0.0
if (ihorn ge 1 and ihorn le mfi.nhorn) then begin
   xpos     = -1.d0 * mfi.xpos[ihorn-1]/mfi.focal*r2d
   ypos     = mfi.ypos[ihorn-1]/mfi.focal*r2d
   ;delta_EL = ypos
   ;delta_AZ = xpos / COS( reform(EL+delta_EL)*d2r )
endif


; Pointing model
sancho_mfi_pmodel, pmodel, ihorn=ihorn


; Compute LMST (local mean sidereal time)
izana,longitude=LON,latitud=LAT,altitude=ALTITUDE,/silent
; Use sancho_jd_to_lst, jd, lmst_hours, caso="" for the old 
; definition, CT2LST. That was equivalent to use
;  CT2LST, lmst_hours, LON, 0, jd
sancho_jd_to_lst, jd, lmst_hours, caso="AOKI"



; Standard correction, using HOR2EQ and POSANG
if (standard) then begin

   ; Correct (AZ,EL) to the position of the horn
   AZ0 = AZ
   EL0 = EL
   if (xpos ne 0.0 and ypos ne 0.0) then fp2sky, AZ, EL, xpos, ypos, AZ0, EL0


   ; RA-DEC. 
   ; In old versions I used: HOR2EQ, el+delta_EL, az+delta_AZ, jd, ...
   HOR2EQ, EL0, AZ0, jd, RA, DEC, HA, LAT=LAT , LON=LON, PRECESS_=1, NUTATE_= 1, ALTITUDE=ALTITUDE
   
   ; PAR-ANGLE
   posang = parangle(lmst_hours*h2d - RA, DEC, LAT, /degree)

   ; For this case, there is no pointing model correction for (AZ,EL):
   AZ_true = AZ0
   EL_true = EL0


endif else begin


; My conversion
   n       = n_elements(AZ)
   ra      = fltarr(n)
   dec     = fltarr(n)
   posang  = fltarr(n)
   AZ_true = fltarr(n)
   EL_true = fltarr(n)
   FOR i=0L, n-1L DO BEGIN

      ; Change to (AZ,EL) for the selected horn
      ;AZ0 = AZ[i] + delta_AZ[i]
      ;EL0 = EL[i] + delta_EL
      AZ0 = AZ[i]
      EL0 = EL[i]
      if (xpos ne 0.0 and ypos ne 0.0) then fp2sky, AZ[i], EL[i], xpos, ypos, AZ0, EL0


      ; Pointing model corrections.
      if (old) then begin
         AZ_c = AZ0
         EL_c = EL0
      endif else begin
         if (i eq 0) then print," --> Performing also the pointing model correction. Pmodel = ",pmodel, mfi.focal
         ;azel2vec, AZ0, EL0, x_in, /degree 
         ;pointing_model, x_in, x_out, pmodel=pmodel,direction=1 ; enc2sky
         ;vec2azel, x_out, AZ_c, EL_c, /degree

         pointing_model_denis, pmodel, AZ0, EL0, AZ_c, EL_c, direction=1
         
      endelse
      AZ_true[i] = AZ_c
      EL_true[i] = EL_c
      

      ; Display info
      if (verbose) then begin
         print," (*) Information: "
         print,"     AZ, EL input values    = ",AZ[i],EL[i] 
         print,"     AZ, EL center of horn  = ",AZ0, EL0
         print,"     AZ, EL after pmodel    = ",AZ_c,EL_c
      endif


      ; Coordinates and parallactic angle
      lst  = lmst_hours[i]*h2d  ; deg
      jd_i = jd[i]
      if (old) then jd_i = 0    ; the apparent to J2000 correction is not done
      denis = 1 ; this was default for azel2skycoord.pro
      if (idl) then denis = 0
      azel2skycoord, AZ_c, EL_c, lst, LAT, RA2, DEC2, parang, /degree, galactic=galactic, jd=jd_i, $
                     AZ0=AZ0,EL0=EL0,verbose=verbose, denis=denis, idl=idl
      ra[i]     = RA2
      dec[i]    = DEC2
      posang[i] = parang
      ;print,parang


   ENDFOR

endelse




END
