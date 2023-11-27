;    22/03/2018   JARM
;
;    PURPOSE:
;    Obtains the MFI pointing model. It has two versions.1=until 2017, 2=after march 2018
;
;    HISTORY:
;    22/03/2018 - original version
;
;-
;

PRO sancho_mfi_get_pointingmodel, jd, az, el, gl, gb, par, azhorn, elhorn, version=version

; Keywords
if not keyword_set(version) then version="Version2"

; Info
if n_params() lt 3 then begin
	print,""
	print,"  Syntax -- sancho_mfi_get_pointingmodel, jd, az, el, gl, gb, par, azhorn, elhorn, version=version "
	print,""
	print,"  INPUTS:"
	print,"     jd, az, el (vectors), VERSION= 'Version1' or 'Version2'"
	print,""
	print,"  OUTPUTS:"
	print,"     gl, gb, par, azhorn, elhorn. All are arrays with dimensions (4,ndata). "
	print,""
	print,"  VERSIONS:"
	print,"     Version1. Model generated in 2013, and used until Feb 2018."
	print,"     Version2. Updated model used from Mar 2018. Corrected precession, and now referred to center of focal plane."
	print,""
	return
endif

; Main code
case version OF
	"Version1": begin
        	sancho_mfi_get_skycoords, jd, az, el, 1, gl1, gb1, par1, az_horn1, el_horn1, /galactic
        	sancho_mfi_get_skycoords, jd, az, el, 2, gl2, gb2, par2, az_horn2, el_horn2, /galactic
        	sancho_mfi_get_skycoords, jd, az, el, 3, gl3, gb3, par3, az_horn3, el_horn3, /galactic
        	sancho_mfi_get_skycoords, jd, az, el, 4, gl4, gb4, par4, az_horn4, el_horn4, /galactic
	end
	"Version2": begin
		sancho_mfi_encoder_to_sky, jd, az, el, 1, gl=gl1, gb=gb1, posang=par1, AZ_HORN=az_horn1, EL_HORN=el_horn1, /galactic
		sancho_mfi_encoder_to_sky, jd, az, el, 2, gl=gl2, gb=gb2, posang=par2, AZ_HORN=az_horn2, EL_HORN=el_horn2, /galactic
		sancho_mfi_encoder_to_sky, jd, az, el, 3, gl=gl3, gb=gb3, posang=par3, AZ_HORN=az_horn3, EL_HORN=el_horn3, /galactic
		sancho_mfi_encoder_to_sky, jd, az, el, 4, gl=gl4, gb=gb4, posang=par4, AZ_HORN=az_horn4, EL_HORN=el_horn4, /galactic
	end
	else: stop," SANCHO_MFI_GET_POINTINGMODEL. Undefined pointing model: ",version
endcase

; Outputs
n = n_elements(jd)

gl       = fltarr(4,n)
gl[0,*]  = gl1
gl[1,*]  = gl2
gl[2,*]  = gl3
gl[3,*]  = gl4

gb       = fltarr(4,n)
gb[0,*]  = gb1
gb[1,*]  = gb2
gb[2,*]  = gb3
gb[3,*]  = gb4

par      = fltarr(4,n)
par[0,*] = par1
par[1,*] = par2
par[2,*] = par3
par[3,*] = par4

azhorn      = fltarr(4,n)
azhorn[0,*] = az_horn1
azhorn[1,*] = az_horn2
azhorn[2,*] = az_horn3
azhorn[3,*] = az_horn4

elhorn      = fltarr(4,n)
elhorn[0,*] = el_horn1
elhorn[1,*] = el_horn2
elhorn[2,*] = el_horn3
elhorn[3,*] = el_horn4

END
