;   23/12/2012   J.A.Rubino
;
;   Read data from TODs.
;


PRO sancho_read_raster_scan, ff_fits, jd, az, el, data, time, mod_ang, cal, $
                             moon=moon, doplot=doplot, str=str

; Keywords
if not keyword_set(moon) then moon=0
if not keyword_set(doplot) then doplot=0
if not keyword_set(str) then str=0


; Params
if n_params() lt 1 then begin
   print,""
   print,"  Syntax -- sancho_read_raster_scan, filename, jd, az, el, data, time, mod_ang, cal [,/MOON, /DOPLOT]"
   print,""
   print,""
   return
endif


; Read TOD. Assuming fits format
if does_exists(ff_fits) then begin

   print," > Reading ",ff_fits

   ; read file
   a  = mrdfits(ff_fits,1,hdr)

   names = tag_names(a)

   ; data
   az      = a.az
   el      = a.el
   data    = a.data
   n       = n_elements(a.az)
   nchan   = n_elements(data[*,0])
   
   ; Mod angle
   cut1 = where(names eq "MOD_ANGLE", n1)
   cut2 = where(names eq "MOD_ANG", n2)
   if (n1 eq 1) then mod_ang = a.mod_angle
   if (n2 eq 1) then mod_ang = a.mod_ang


   ; Out of range values for elevation.
   cut = where(el ge 90.0, ncut)
   if ncut ge 1 then el[cut] = 90.0

   ; Calibration signal. CAL=1 means CAL signal is on.
   cut3 = where(names eq "CAL", n3)
   if (n3 eq 1) then begin
      cal = a.cal
   endif else begin
      cal = intarr(n)
   endelse

   ; JD. Referenced to the first day of the Commissioning (13/Nov/2012, 0.00h)
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
   str    = {az:az, el:el, jd:jd, time:time, data:data, mod_angle:mod_ang, cal:cal }

endif else begin

   print," > File not found! "
   stop

endelse



; Compute Moon location
if (moon) then begin
   h2d = 15.0
   moonpos,jd, ra_moon, dec_moon, d_moon
   
   izana,longitude=LON,latitud=LAT,ALTITUDE=ALTITUDE
   CT2LST, lst, LON, 0, jd
;   HA = lst*h2d - ra_moon
;   HADEC2ALTAZ, HA, dec_moon, LAT, el_moon, az_moon
   EQ2HOR, ra_moon, dec_moon, jd, el_moon, az_moon, HA, $
           LAT=LAT , LON=LON, PRECESS_=1, $
           NUTATE_= 1, ALTITUDE=ALTITUDE
   print," Moon location computed. "
endif



; Plot full raster
if (doplot) then begin
	cfp
        !p.charsize=1.5

        window,0,title="AZ-EL trajectory"
        !p.multi=0
	cut = where(az ne 0 and el ne 0)
	plot,az[cut],el[cut],/xs,/ys,xtitle="AZ (deg)",ytitle="EL(deg)"
	if (moon) then oplot,az_moon,el_moon,psym=3,col=2

        window,1,title="Modulator"
        !p.multi=[0,1,4]
	cut = where(az ne 0 and el ne 0)
        for i=0,3 do begin
           plot,time-time[0],mod_ang[i,*],/xs,/ys,xtitle="Time (s)",$
                ytitle="Modulator (deg)",xr=[0,40]
        endfor

endif

; Info
print," (*) SANCHO_READ_RASTER_SCAN. File read. "

return
END
