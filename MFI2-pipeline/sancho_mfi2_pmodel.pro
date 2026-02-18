;  28/May/2024   J.A.Rubino
;
;  Pointing model for QT1 and MFI2
;
;  HISTORY:
;    28/May/2024 : original version
;  


FUNCTION pars2pmodel_structure, a
  pmodel = { P_f:a[0], P_x:a[1], P_y:a[2], P_c:a[3], P_n:a[4], P_a:a[5], P_b:a[6] } 
  return,pmodel
END

;-----------------

PRO sancho_mfi2_pmodel, pmodel, ihorn=ihorn, center=center, zero=zero, pointingmodel=pointingmodel

; Keywords
if not keyword_set(ihorn) then ihorn=0
if not keyword_set(center) then center=0
if not keyword_set(zero) then zero=0
if not keyword_set(pointingmodel) then pointingmodel=''

arsec2rad  = !DPI/3600.0/180.0 

; Initializing the structure to zero.
get_params_pointing_model, pmodel
txt = "Pointing model initialized to zero. "
pointingmodel = 'Zero'

; Computed pointing models
; [1] Model for central horn. Computed by M.Fernandez with CASS, CRAB observations. May-2024
;  Estos son; por orden: Pa, Pb, Pc, Pn, Pf, Px, Py
mateo_orig = [0.84, -0.72, 0.33, 0.04, 0.07, 0.03, -0.03]*3600.0 * arsec2rad
mateo_h1   = mateo_orig[[5,6,7,3,4,1,2]-1]
pmodel     = pars2pmodel_structure(mateo_h1)
txt    = "Model for QT1 telescope (central pixel, Mateo Fernandez, May/2024)"
pointingmodel = 'MFT-May2024'


;==========================
; Set pmodel to zero if requested
if zero then begin
	get_params_pointing_model, pmodel
	txt =  "     > Using null pointing model (all parameters to zero)."
endif
print,txt


END
