;   11/Oct/2013   J.A.Rubino
;
;   HISTORY:
;     11/10/2013 - original version, based on sancho_average_data.pro. This code
;                  is used for binning data, accounting properly for gaps of the
;                  CAL diode.
;     02/12/2013 - The sig_avg vector now is divided by sqrt(number of points).
;     06/02/2015 - Computation of navg uses round() instead of long().
;     14/03/2016 - option for nulltests of using only even or odd index numbers. Nulltest=1,2 selects the even and odd parts.
;
;
;-
;


PRO sancho_bin_data, data, jd, el, az, mod_ang, wei, $
                     data_avg, jd_avg, el_avg, az_avg, mod_ang_avg, sig_avg, $
                     msbin=msbin, nulltest=nulltest

; Keywords
if not keyword_set(msbin) then msbin = 100.0 ; bin size in ms.
if not keyword_set(nulltest) then nulltest = 0 ; carry out nulltest splitting of data. 


; Info
if n_params() lt 4 then begin
   print,""
   print,"  Syntax -- sancho_bin_data, data, jd, el, az, mod_ang, wei, data_avg, jd_avg, el_avg, az_avg, $"
   print,"                                 mod_ang_avg, rms_avg [,MSBIN=, NULLTEST=]"
   print,""
   print,""
   return
endif


; Params
d2r   = !DPI/180.0
nchan = n_elements(data[*,0])
n     = n_elements(az)
nq    = n_elements(mod_ang[*,0])

; Check nulltest values
nulltest = fix(nulltest)
if nulltest gt 2 then nulltest=0
if nulltest lt 0 then nulltest=0


; Averaging
IF (msbin gt 1.0) THEN BEGIN

   printf,-1," (*) SANCHO_BIN_DATA: averaging in bins of ",msbin," ms.",format='(a47,i5,a4)'


   t_avg = msbin/1e3            ; averaging over this number of secs. 
   fsamp = 1.d3                 ; sampling at this rate (Hz). Default value.
   navg  = round(t_avg * fsamp) ; number of samples of the average.
   MINSAMP = 5                  ; minimum no of samples to compute averages
   
   ns          = round(n/navg)
   data_avg    = fltarr(nchan, ns)
   sig_avg     = fltarr(nchan,ns)
   mod_ang_avg = fltarr(nq, ns)
   jd_avg      = dblarr(ns)
   el_avg      = fltarr(ns)
   az_avg      = fltarr(ns)


   ; data
   print,"      Binning data..."
   print,"      Info: t_avg, fsamp, navg, ns, nulltest = ", t_avg, fsamp, navg, ns, nulltest
   for i=0,nchan-1 do begin
      for j=0L,ns-1L do begin
         i1 = j * navg
         i2 = (j+1L) * navg - 1L
         im = (i1+i2)/2

         wei_d = wei[i1:i2]
         ddd   = data[i,i1:i2]*wei_d
	 ; Re-select data for the nulltest case. Added on March 14th, 2018
	 if nulltest ne 0 then begin
		indice  = indgen(n_elements(ddd))
		evenodd = nulltest-1 ; values are 0 or 1
		if (evenodd ne 0 and evenodd ne 1) then stop," Wrong value of nulltest"
         	cut     = where(ddd ne 0.0 and wei_d ne 0.0 and ((indice mod 2) eq evenodd), ncut)
	 endif else begin
         	cut = where(ddd ne 0.0 and wei_d ne 0.0, ncut)
	 endelse
         if (ncut ge MINSAMP) then begin
            data_avg[i,j] = mean( ddd[cut] )
            sig_avg[i,j]  = stdev( ddd[cut] ) / sqrt(ncut*1.d0)
         endif
         jd_avg[j]     = jd[im]
         az_avg[j]     = az[im]
         el_avg[j]     = el[im]
      endfor
   endfor

   ; modulator angles
   print,"      Binning modulator angles..."
   for i=0,nq-1 do begin
      for j=0L,ns-1L do begin
         i1 = j * navg
         i2 = (j+1L) * navg - 1L
         im = (i1+i2)/2
;         mod_ang_avg[i,j] = mean( mod_ang[i,i1:i2] ) ; mod_ang[i,im] 
         mod_ang_avg[i,j] = mod_ang[i,im] 
      endfor
   endfor


   

ENDIF ELSE BEGIN
   ; Data stay at their original binning
   data_avg    = data
   jd_avg      = jd
   az_avg      = az
   el_avg      = el
   sig_avg     = fltarr(nchan,n) + 1.0
   mod_ang_avg = mod_ang

   cut = where(wei eq 0, ncut)
   for ichan=0,nchan-1 do begin
      sig_avg[ichan,cut] = 0.0
   endfor

   printf,-1," (*) SANCHO_BIN_DATA: no averaging. Keeping bins of 1 ms."

ENDELSE


END
