;   24/May/2024   J.A.Rubino
;
;   HISTORY:
;     24/05/2024 - original version, based on sancho_bin_data_mfi.pro
;     12/11/2024 - adding wei from CAL signal
;-
;

PRO sancho_bin_data_mfi2, data, jd, el, az, wei, data_avg, jd_avg, el_avg, az_avg, $
                          sig_avg, sigsum_avg, sigdiff_avg, $
                          cov_avg, msbin=msbin, nulltest=nulltest, newmfi=mfi2

; Keywords
if not keyword_set(msbin) then msbin = 40.0  ; default bin size in ms.
if not keyword_set(nulltest) then nulltest = 0 ; carry out nulltest splitting of data. 
if not keyword_set(newmfi) then info_quijote, newmfi=mfi2

; Info
if n_params() lt 4 then begin
   print,""
   print,"  Syntax -- sancho_bin_data_mfi2, data, jd, el, az, wei, data_avg, jd_avg, el_avg, az_avg, $"
   print,"                              rms_avg, sigsum_avg, sigdiff_avg, cov_avg [,MSBIN=, NULLTEST=, NEWMFI=]"
   print,""
   print,""
   return
endif


; Params
d2r   = !DPI/180.0
nchan = n_elements(data[*,0])
n     = n_elements(az)
nhorn = 4                       ; Not using mfi2.NHORN, but OLD MFI values, OLD DAS
nfreq = mfi2.NFREQ_PER_HORN
ncorr = 2                       ; old notation of MFI.

; Check nulltest values
nulltest = fix(nulltest)
if nulltest gt 2 then nulltest=0
if nulltest lt 0 then nulltest=0


; Averaging
IF (msbin gt 1.0) THEN BEGIN

   printf,-1," (*) SANCHO_BIN_DATA_MFI2: averaging in bins of ",msbin," ms.",format='(a47,i5,a4)'

   t_avg = msbin/1e3            ; averaging over this number of secs. 
   fsamp = 1.d3                 ; sampling at this rate (Hz). Default value.
   navg  = round(t_avg * fsamp) ; number of samples of the average.
   MINSAMP = 5                  ; minimum no of samples to compute averages
   
   ns          = round(n/navg)
   data_avg    = fltarr(nchan, ns)
   sig_avg     = fltarr(nchan,ns)
   jd_avg      = dblarr(ns)
   el_avg      = fltarr(ns)
   az_avg      = fltarr(ns)
   ; New variables, JARM, 7/7/2019
   sigsum_avg  = fltarr(nhorn,nfreq,ncorr,ns)
   sigdiff_avg = fltarr(nhorn,nfreq,ncorr,ns)


   ; data
   print,"      Binning data..."
   print,"      Info: t_avg, fsamp, navg, ns, nulltest = ", t_avg, fsamp, navg, ns, nulltest

   ; Alternative way of running over all database, keeping info of sum and 
   ; difference of channels. 7/7/2019, JARM.
   for ihorn=0,nhorn-1 do begin
      for ifreq=0,nfreq-1 do begin
         for icorr=0,ncorr-1 do begin

            for j=0L,ns-1L do begin
               i1 = j * navg
               i2 = (j+1L) * navg - 1L
               im = (i1+i2)/2
               
	       ichan_s = mfi2.ic[ihorn,ifreq,icorr*2]
	       ichan_d = mfi2.ic[ihorn,ifreq,icorr*2+1]

               wei_d = wei[i1:i2]
               ddd_s = data[ichan_s,i1:i2]*wei_d
               ddd_d = data[ichan_d,i1:i2]*wei_d

	       ; Re-select data for the nulltest case. Added on March 14th, 2018
               ; Note on 7/7/2019. It forces to have same selection in sum and diff
	       if nulltest ne 0 then begin
	  	  indice  = indgen(n_elements(ddd_s))
		  evenodd = nulltest-1 ; values are 0 or 1
		  if (evenodd ne 0 and evenodd ne 1) then stop," Wrong value of nulltest"
         	  cut = where(ddd_s ne 0.0 and ddd_d ne 0.0 and wei_d ne 0.0 and ((indice mod 2) eq evenodd), ncut)
	       endif else begin
         	  cut = where(ddd_s ne 0.0 and ddd_d ne 0.0 and wei_d ne 0.0, ncut)
	       endelse
               if (ncut ge MINSAMP) then begin
                  ; ichan_s
                  data_avg[ichan_s,j] = mean( ddd_s[cut] )
                  sig_avg[ichan_s,j]  = stddev( ddd_s[cut] ) / sqrt(ncut*1.d0)
		  ; ichan_d
                  data_avg[ichan_d,j] = mean( ddd_d[cut] )
                  sig_avg[ichan_d,j]  = stddev( ddd_d[cut] ) / sqrt(ncut*1.d0)
                  ; stddev for sum and diff. New JARM, 7/7/2019. 
                  sigsum_avg[ihorn,ifreq,icorr,j]  = stddev( ddd_s[cut] + ddd_d[cut]) / sqrt(ncut*1.d0)
                  sigdiff_avg[ihorn,ifreq,icorr,j] = stddev( ddd_s[cut] - ddd_d[cut]) / sqrt(ncut*1.d0)
               endif
               jd_avg[j]     = jd[im]
               az_avg[j]     = az[im]
               el_avg[j]     = el[im]
               
            endfor
            
         endfor
      endfor
   endfor


ENDIF ELSE BEGIN
   ; Data stay at their original binning
   data_avg    = data
   jd_avg      = jd
   az_avg      = az
   el_avg      = el
   sig_avg     = fltarr(nchan,n) + 1.0
   sigsum_avg  = fltarr(nhorn,nfreq,ncorr,n)+1.0
   sigdiff_avg = fltarr(nhorn,nfreq,ncorr,n)+1.0

   cut = where(wei eq 0, ncut)
   for ichan=0,nchan-1 do begin
      sig_avg[ichan,cut] = 0.0
   endfor

   printf,-1," (*) SANCHO_BIN_DATA_MFI2: no averaging. Keeping bins of 1 ms."

ENDELSE

; Generate cov 
cov_avg = (sigsum_avg^2. - sigdiff_avg^2.)/4.0


END
