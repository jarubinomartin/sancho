;  4/May/2017   JARM
;  For a given vector of AZ values, it checks if it corresponds to NOMINAL data
;

pro check_if_nominal, az, is_nominal

if n_params() eq 0 then begin
   print,""
   print," Syntax -- check_if_nominal, az, is_nominal "
   print,""
   print,""
   return
endif

  offsets   = [0.0, 90.0, 180.0, 270.0]
  noffs     = n_elements(offsets)
  deltaaz_v = fltarr(noffs)
  
  for i=0,noffs-1 do begin
     azmin   = min( (az+offsets[i]) mod 360.0)
     azmax   = max( (az+offsets[i]) mod 360.0)
     deltaaz = azmax-azmin
     deltaaz_v[i] = deltaaz
  endfor
  
  is_nominal = 0
  deltaaz = min(deltaaz_v)
  if (deltaaz ge 355.0) then is_nominal = 1

; Display info
resp = ["NO","YES"]
if n_params() eq 1 then print," > Is nominal? ",resp[is_nominal]

end

