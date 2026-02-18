;    28/05/2024   JARM
;
;    PURPOSE:
;    Obtains the MFI2 pointing model. Based on the MFI code.
;
;    HISTORY:
;    28/05/2024 - original version
;
;-
;

PRO sancho_mfi2_get_pointingmodel, jd, az, el, gl, gb, par, azhorn, elhorn, version=version, nhorns=nhorns

; Keywords
if not keyword_set(version) then version="Version1"
if not keyword_set(nhorns) then nhorns=4 ; using old MFI DAS


; Info
if n_params() lt 3 then begin
	print,""
	print,"  Syntax -- sancho_mfi2_get_pointingmodel, jd, az, el, gl, gb, par, azhorn, elhorn, NHORNS= "
	print,""
	print,"  INPUTS:"
	print,"     jd, az, el (vectors), NHORNS=, VERSION= 'Version1' "
	print,""
	print,"  OUTPUTS:"
	print,"     gl, gb, par, azhorn, elhorn. All are arrays with dimensions (4,ndata). "
	print,""
	print,"  VERSIONS:"
	print,"     Version1. Model fitted in May 2024 by MFT."
	print,""
	return
endif

; Main code
n      = n_elements(jd)
gl     = fltarr(nhorns,n)
gb     = fltarr(nhorns,n)
par    = fltarr(nhorns,n)
azhorn = fltarr(nhorns,n)
elhorn = fltarr(nhorns,n)
for ihorn=0,nhorns-1 do begin
   horn = ihorn+1
   sancho_mfi2_encoder_to_sky, jd, az, el, horn, gl=gl1, gb=gb1, posang=par1, $
                               AZ_HORN=az_horn1, EL_HORN=el_horn1, /galactic  
   gl[ihorn,*]  = gl1
   gb[ihorn,*]  = gb1
   par[ihorn,*] = par1
   azhorn[ihorn,*] = az_horn1
   elhorn[ihorn,*] = el_horn1
endfor


END
