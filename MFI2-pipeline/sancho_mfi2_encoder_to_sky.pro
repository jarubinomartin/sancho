;   December/2016  D.Tramonte
;
;   New routine for obtaining the pointing model by Denis. It changes the time definition
;   to CT2LST from IDL, and uses HOR2EQ for the horizontal to equatorial transformation.
;
;   PURPOSE: 
;   For given jd, AZ_enc and EL_enc input vectors, and a given horn of the MFI
;   instrument, it returns RA,DEC (or GL,GB).
;
;   HISTORY:
;       Dec/2016 - original version [DT]
;    12/Nov/2017 - name changed to sancho_mfi_encoder_to_sky.pro. Now the code computes the parangle. Output adapted [JARM]
;    28/May/2024 - adapted to MFI2
;
;-
;

function scalarprod,x,y
return,total(x*y)
end

function vectorprod,A,B
  A1 = A[0]
  A2 = A[1]
  A3 = A[2]
  B1 = B[0]
  B2 = B[1]
  B3 = B[2]
return,[ (A2*B3 - A3*B2), (A3*B1 - A1*B3), (A1*B2 - A2*B1) ]
end

;============

pro sancho_mfi2_encoder_to_sky, jd, az_enc, el_enc, horn, coord1, coord2, parang, APPARENT=apparent, GALACTIC=galactic, $
                               AZ_TEL=az_tel, EL_TEL=el_tel, AZ_HORN=az_horn, EL_HORN=el_horn, $
                               RA_APP=ra_app, DEC_APP=dec_app, RA_J2000=ra_j2000, DEC_J2000=dec_j2000, $ 
                               GL=gl, GB=gb, VERBOSE=verbose, POSANG=posang

; Description message
if n_params() lt 1 then begin 
   print, ""
   print, " ============================================================================================="
   print, ""
   print, " -- Syntax: "
   print, ""
   print, ""
   print, "   sancho_encoder_to_sky, jd, az_enc, el_enc, horn, APPARENT=apparent, GALACTIC=galactic, "
   print, "                          AZ_TEL=az_tel, EL_TEL=el_tel, AZ_HORN=az_horn, EL_HORN=el_horn, "
   print, "                          RA_APP=ra_app, DEC_APP=dec_app, RA_J2000=ra_j2000, DEC_J2000=dec_j2000, "
   print, "                          GL=gl, GB=gb, VERBOSE=verbose, POSANG=posang"
   print, ""
   print, ""
   print, "   *** ALL ANGULAR COORDINATES (INPUT AND OUTPUT) IN DEGREES ***"
   print, ""
   print, ""
   print, "                       ----  INPUT  ----    "
   print, ""
   print, "   jd : julian date (scalar or vector)"
   print, "   (az_enc, el_enc) : telescope encoder coordinates (scalars or vectors)"
   print, "   horn : {0,1,2,3,4}. Set 0 for focal plane center"
   print, "   /APPARENT : if set, computation of apparent equatorial coordinates is performed"
   print, "   /GALACTIC : if set, computation of galactic coordinates is performed"
   print, ""
   print, ""
   print, "                       ----  OUTPUT ----    "
   print, ""
   print, "   (AZ_TRUE, EL_TRUE): telescope horizon coordinates after pointing correction"
   print, "   (AZ_HORN, EL_HORN): horn horizon coordinates after pointing correction"
   print, "   (RA_APP, DEC_APP): apparent equatorial coordinates"
   print, "   (RA_J2000, DEC_J2000): equatorial coordinates for epoch J2000"
   print, "   (GL, GB): galactic coordinates"
   print, "   (posang): paralactic angle, in Equatorial or Galactic coordinates, as requested"
   print, ""
   print, " ============================================================================================="
   print, ""
   return
endif


;---------------------------------------------------------------------------------------------------------
; Input settings
appt = 0
glct = 0
if keyword_set(APPARENT) then appt = 1
if keyword_set(GALACTIC) then glct = 1
d2r   = !dpi/180.d
r2d   = 1.d / d2r
h2d   = 15.d0

;check input
njd = n_elements(jd)
naz = n_elements(az_enc)
nel = n_elements(el_enc)
if (nel ne naz) then begin 
   print, ""
   print, " -- Input arrays do not have the same length. Returning."
   print, ""
   return
endif
if (njd eq 1) and (naz gt 1) then jd = dblarr(naz) + jd
if (njd gt 1) and (njd ne naz) then begin 
   print, ""
   print, " -- Input arrays do not have the same length. Returning."
   print, ""
   return
endif

;output vectors
nps = n_elements(jd)
az_tel = dblarr(nps) 
el_tel = dblarr(nps)
az_horn = dblarr(nps) 
el_horn = dblarr(nps)
ra_app = dblarr(nps)
dec_app = dblarr(nps)
ra_j2000 = dblarr(nps)
dec_j2000 = dblarr(nps)
gl = dblarr(nps)
gb = dblarr(nps)
posang = dblarr(nps)

;Get horn position. Specific for MFI2.
info_quijote, newmfi=mfi2
if ((horn gt mfi2.nhorn) or (horn lt 0)) then begin 
   print, ""
   print, " -- You must specify a horn number between 0 and "+strtrim(mfi2.nhorn,2)+". Returning. "
   print, ""
   return
endif
x_horn = 0.d
y_horn = 0.d
if (horn ge 1) then begin 
   x_horn = -1.d0 * mfi2.xpos[horn-1]/mfi2.focal*r2d
   y_horn = mfi2.ypos[horn-1]/mfi2.focal*r2d
endif

;Load pointing model, referred to center of focal plane 
print, ""
sancho_mfi2_pmodel, pmodel
print, ""
   
; Observatory coordinates
izana, longitude=lon, latitud=lat, altitude=altitude, /silent

;---------------------------------------------------------------------------------------------------------
;loop over (eventual) input jd vector
print," Doing Pmodel"
pointing_model_denis, pmodel, az_enc, el_enc, az_tel, el_tel, direction=1


;from focal plane center to horn position
print," Doing fp2sky"
fp2sky, az_tel, el_tel, x_horn, y_horn, az_horn, el_horn
;; Use this only for parangle computation
fp2sky, az_enc, el_enc, x_horn, y_horn, az_horn_nocorr, el_horn_nocorr


;from horizon coordinates to J2000 equatorial
print," Doing hor2eq"
hor2eq, el_horn, az_horn, jd, ra_j2000, dec_j2000, lat=lat, lon=lon, $ 
        precess_=1, nutate_=1, refract_=0, aberration_=1, ALTITUDE=altitude

; computing LST
sancho_jd_to_lst, jd, lmst_hours, CASO=""
lst = lmst_hours*h2d 

;if specified, computes apparent equatorial coordinates
if appt then begin 
   hor2eq, el_horn, az_horn, jd, ra_app, dec_app, lat=lat, lon=lon, $ 
           precess_=0, nutate_=0, refract_=0, aberration_=0, altitude=altitude
endif


; Computing parangle
print," Doing parangle"
if (glct) then print," Doing galct"
;posang2 = parangle(lst - ra_j2000, dec_j2000, lat, /degree)
for ip=0L, nps-1L do begin

     ra_j  = ra_j2000[ip]
     dec_j = dec_j2000[ip]

     radec2vec, ra_j, dec_j, x_j2000, /degree
     x_out = x_j2000
     R     = rot_hor2eq( lst[ip], lat, /degree )

     if (glct) then begin 
	R_gal = rot_eq2gal()      ; Computed for J2000
        x_gal = Matrix_multiply(R_gal, x_j2000)
        vec2radec, x_gal, glon, glat, /degree
        gl[ip] = glon
        gb[ip] = glat
        R_1 = R
        R   = Matrix_multiply(R_gal,R_1)
	x_out = x_gal
     endif

     ; Parangle
     ;azel2orientvec, az_tel[ip], el_tel[ip], x_EL, x_AZ, /degree ; option 1.
     ;azel2orientvec, az_horn_nocorr[ip], el_horn_nocorr[ip], x_EL, x_AZ, /degree ; option 2
     azel2orientvec, az_horn[ip], el_horn[ip], x_EL, x_AZ, /degree ; option 3
     x_EL_out = Matrix_multiply(R, x_EL)
     z        = [0.0,0.0,1.0]      ; unitary vector in Z direction. 
     parang   = ATAN( scalarprod(x_EL_out,vectorprod(z,x_out)), $
                      scalarprod(x_EL_out,vectorprod(x_out,vectorprod(z,x_out))) ) / d2r
     posang[ip] = parang 

endfor


;-------------------
; Output
coord1 = ra_j2000
coord2 = dec_j2000
parang = posang
if (glct) then begin 
	coord1 = gl
	coord2 = gb
endif   

;---------------------------------------------------------------------------------------------------------
;print info
if keyword_set(VERBOSE) then begin 
   index=0
   print, "   ========================================================================"
   print, ""
   print, "      ENCODER             (AZ,EL) = ", az_enc[index], el_enc[index]
   print, ""
   print, "      HORIZON (TEL.)      (AZ,EL) =  ", az_tel[index], el_tel[index]
   print, ""
   print, "      HORIZON (HORN)      (AZ,EL) =  ", az_horn[index], el_horn[index]
   if appt then begin 
      print, ""
      print, "      EQUATORIAL APP.     (RA,DEC) =  ", ra_app[index], dec_app[index]
   endif
   print, ""
   print, "      EQUATORIAL J2000    (RA,DEC) =  ", ra_j2000[index], dec_j2000[index]
   if glct then begin 
      print, ""
      print, "      GALACTIC            (GL, GB) =  ", gl[index], gb[index]
   endif
   print, ""
   print, "   ========================================================================"
   print, ""
endif


return
end
