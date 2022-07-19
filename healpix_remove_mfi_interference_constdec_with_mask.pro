;   23/05/2017   JARM
;
;   Removes MFI interference by substracting a constant pattern in
;   declination (FDEC). 
;   Input is a healpix map (RING ordering is assumed).
;
;-
;

PRO get_function_constdec_for_onemap, mapanomask, mapa, dec, map_out, fdec=fdec

npix   = n_elements(mapa)

; 1) Stack in declination
map_out = mapanomask
cutno0  = where(map_out ne 0.0, nonzeros)
cutno0m = where(mapa ne 0.0, nonzerosm)
decmin  = min(dec[cutno0m])
decmax  = max(dec[cutno0m])
deltad_arcmin  = 20.0 ; in arcmin. Tested for 1deg beam
deltad  = deltad_arcmin/60.0
ndecs   = round( (decmax-decmin)/deltad ) 
if ndecs eq 0 then ndecs = 300
;ndecs   = 300 ; this was equivalent to approx 20 arcmin
;deltad  = (decmax-decmin)/float(ndecs)
vdec    = findgen(ndecs)*deltad + decmin
  
stack  = fltarr(ndecs)
for j=0,ndecs-1 do begin
   cut = where(dec ge (vdec[j]-deltad/2.0) and dec lt (vdec[j]+deltad/2.0), ncut)
   if ncut ge 2 then begin
      pixset = mapa[cut]
      ind = where(pixset ne 0, count)
      if count ge 1 then pixset = pixset[ind]
      if n_elements(pixset) ge 2 then stack[j] = median(pixset)
   endif
endfor
         
; 1.1 Project template back onto sky
template = fltarr(npix)
template2 = fltarr(npix)
  
if (nonzeros ge 1) then begin
   ii = (dec - decmin)/deltad
   cut = where(ii ge 0 and ii lt ndecs, ncut)
   if ncut ge 1 then template[cut] = stack[ii[cut]]
     
   X  = vdec
   Y  = stack
   Y2 = SPL_INIT(X, Y)
   cut = where(dec ge decmin and dec le decmax, ncut)
   X2  = dec[cut] 
   res = SPL_INTERP(X, Y, Y2, X2)
     
   if ncut ge 1 then template2[cut] = res[*]
endif
         
; 1.2 Substract template
if nonzeros ge 1 then map_out[cutno0] = mapanomask[cutno0] - template2[cutno0]

; Output fdec
fdec = {dec:vdec, fdec:stack}

END

;========================================

PRO healpix_remove_mfi_interference_constdec_with_mask, map_in, map2, order_in=order_in, nside=nside, $
   dipole=dipole, gal_cut=gal_cut,pol=pol, maskI=maskI, maskP=maskP, $
   Ifdec=fdec_I, Qfdec=fdec_Qr, Ufdec=fdec_Ur

; Info
if n_params() eq 0 then begin
   print,""
   print,"   Syntax -- healpix_remove_mfi_interference_constdec_with_mask, map_in, map_out, order_in=order_in, nside=nside, [, /DIPOLE, /POL, maskI=maskI, maskP=maskP, Ifdec=, Qfdec=, Ufdec= ] "

   print,""
   print,"   Default ordering is RING. Output ordering is RING"
   print,"   If set /dipole, it also removes a dipole from each map. "
   print,"   Default maskI is a map of ones everywhere. Default maskP is maskI."
   print,"   If /POL is active, then the FDEC is applied to polarization maps also. "
   print," "
   print,"   INPUT: "
   print,"        map_in, either (npix) or (npix,3). "
   print,"   OUTPUT: "
   print,"        map_out after filtering with FDEC. If requested, the fitted function of declination is "
   print,"        also returned as Ifdec=fdec_I, Qfdec=fdec_Qr, Ufdec=fdec_Ur"
   print,""
   print,"   Example:  "
   print,"       IDL> maskI = annular_mask(10.0,-10.0,nside) "
   print,"       IDL> healpix_remove_mfi_interference_constdec_with_mask, map_in, map_out, maskI=maskI"
   print,""
   return
endif

; Keywords
if not keyword_set(order_in) then order_in="RING"
if not keyword_set(nside) then nside = 0
if not keyword_set(dipole) then dipole = 0
if not keyword_set(gal_cut) then gal_cut = 0.0
if not keyword_set(pol) then pol=0

; Compute npix and nside
npix = n_elements(map_in[*,0])
if nside eq 0 then begin
   nside = npix2nside(npix)
endif 

; Compute nmaps (1 for intensity, 3 for int and polarization).
nmaps = n_elements(map_in[0,*])

; additional masks
if maskI eq !null then begin
   dmaskI = replicate(1,npix)
endif else begin
   nside_maskI = npix2nside(n_elements(maskI))
   if nside_maskI ne nside $
      then ud_grade, maskI, dmaskI, nside_out=nside, order_in=order_in $
      else dmaskI = maskI
endelse
if maskP eq !null then begin
   dmaskP = dmaskI
endif else begin
   nside_maskP = npix2nside(n_elements(maskP))
   if nside_maskP ne nside $
      then ud_grade, maskP, dmaskP, nside_out=nside, order_in=order_in $
      else dmaskP = maskP
endelse


; Ordering
map1nomask = map_in
if (strupcase(order_in) ne "RING") then begin
   map1nomask = reorder(map_in, /N2R)
   dmaskI = reorder(dmaskI, /N2R)
   dmaskP = reorder(dmaskP, /N2R)
endif

map1 = map1nomask
map1[*,0] *= dmaskI
if nmaps ge 3 then begin
   for p=1,2 do map1[*,p] *= dmaskP
endif

; Output map
map2 = map1nomask
print," (*) Removing interference pattern (assuming constant declination) with mask"

npix = nside2npix(nside*1L)
ipix = lindgen(npix)
PIX2ANG_RING, nside, ipix, theta, phi
d2r  = !DPI/180.0
gl   = phi / d2r
gb   = 90.0 - theta/d2r 
glactc,ra, dec, 2000.0, gl, gb, 2, /DEGREE
  
; Correcting intensity maps. itype=0  
print,"     --> Intensity maps"
itype = 0
         
mapanomask = reform(map1nomask[*,0]) 
mapa       = reform(map1[*,0]) 
         
; 1) Substract the stack in declination, fitting a function f(dec)
get_function_constdec_for_onemap, mapanomask, mapa, dec, map_out, fdec=fdec_I

; 2) Removing dipole if requested
if (dipole ne 0) then begin
     cut0 = where(map_out eq 0.0, n0)
     remove_dipole, map_out, nside=nside, ordering="RING", /silent, gal_cut=gal_cut
     if n0 ge 1 then map_out[cut0] = 0.0
endif
	
; 3) Save map
map2[*,0] = map_out[*]
            

; Correcting polarization maps. itype=1,2
if (pol ne 0) then begin

   print,"     --> Polarization maps"
   mapaQnomask = reform(map1nomask[*,1]) 
   mapaUnomask = reform(map1nomask[*,2]) 
   mapaQ = reform(map1[*,1]) 
   mapaU = reform(map1[*,2]) 

   ; 1) Rotate maps to RA-DEC parangle for interference correction
   rotate_stokes_healpixmap, mapaQnomask, mapaUnomask, mapaQnomaskr, mapaUnomaskr, /gal2eq, /silent
   rotate_stokes_healpixmap, mapaQ, mapaU, mapaQr, mapaUr, /gal2eq, /silent
            
   ; 2) Substract the stack in declination, fitting a function f(dec)
   get_function_constdec_for_onemap, mapaQnomaskr, mapaQr, dec, mapQ_out, fdec=fdec_Qr
   get_function_constdec_for_onemap, mapaUnomaskr, mapaUr, dec, mapU_out, fdec=fdec_Ur

   ; 3) Rotate back the maps to GAL parangle after the correction
   rotate_stokes_healpixmap, mapQ_out, mapU_out, Qfinal,Ufinal, /eq2gal, /silent

   ; 2) Save maps
   map2[*,1] = Qfinal[*]
   map2[*,2] = Ufinal[*]
            
endif else begin
   print,"     --> Polarization maps: not corrected."
endelse


end
