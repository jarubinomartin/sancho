;   24/05/2024   J.A.Rubino
;
;   Read MFI2 data from TODs. Based on sancho_read_raster_scan.pro for MFI.
;


PRO sancho_read_mfi2_tod, ff_fits, jd, az, el, data, time, cal, doplot=doplot, str=str

; Keywords
if not keyword_set(moon) then moon=0
if not keyword_set(doplot) then doplot=0
if not keyword_set(str) then str=0


; Params
if n_params() lt 1 then begin
   print,""
   print,"  Syntax -- sancho_read_mfi2_tod, filename, jd, az, el, data, time, cal"
   print,""
   print,""
   return
endif


; Read TOD. Assuming fits format
if does_exists(ff_fits) then begin

   print," > Reading ",ff_fits

   ; Read file
   a     = mrdfits(ff_fits,1,hdr)
   names = tag_names(a)

   ; Data
   az      = a.az
   el      = a.el
   data    = a.data
   n       = n_elements(a.az)
   nchan   = n_elements(data[*,0])

   
   ; Out of range values for elevation.
   cut = where(el ge 90.0, ncut)
   if ncut ge 1 then el[cut] = 90.0

   ; Calibration signal. CAL=1 means CAL signal is on.
   ; As for MFI. Added on Nov-2024 for MFI2
   cut3 = where(names eq "CAL", n3)
   if (n3 eq 1) then begin
      cal = a.cal
   endif else begin
      cal = intarr(n)
   endelse

   ; JD. Referenced to the first day of the QUIJOTE Commissioning (13/Nov/2012, 0.00h)
   ; Use only in raw data.
   jd_ref = 2456244.5d
   jd     = a.jd
   if (min(jd) lt jd_ref) then jd = a.jd + jd_ref

   ; jd  = jd1970 + time/factor
   factor = 24.d * 60.d * 60.d  ; secs per day.
   jd1970 = 2440587.5d
   time   = (jd - jd1970)*factor 
   
   ; Output information
   delvar = a
   str    = {jd:jd, az:az, el:el, data:data, time:time, cal:cal }

endif else begin

   print," > File not found! "
   stop

endelse


; Plot full raster
if (doplot) then begin
   !p.charsize=1.5
   window,0,title="AZ-EL trajectory"
   !p.multi=0
   cut = where(az ne 0 and el ne 0)
   plot,az[cut],el[cut],/xs,/ys,xtitle="AZ (deg)",ytitle="EL(deg)"
endif

; Info
print," (*) SANCHO_READ_MFI2_TOD. File read. "

return
END
