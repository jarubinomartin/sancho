;   2/09/2014  J.A.Rubino
;
;   PURPOSE:
;   Generates .ctod files for QUIJOTE MFI
;
;
;   HISTORY: 
;   02/09/2014 - original version
;   04/05/2015 - calibration is done separately for A,B,C,D positions.
;   23/02/2016 - using new factors from Cass derived by RGS. Calibrating on Y
;                channel.
;   10/05/2016 - calibration now based on x+y for uncorrelated, and (x+y)+(x-y)
;                for correlated channels.
;    9/02/2017 - Now using mfi_params.pro to obtain gains and diode
;                amplitudes. Old code is still usable (with select_values_for_gain_correction_old.pro) 
;   07/03/2017 - Now using .gain files as inputs.
;   20/10/2017 - Updating wei after the recalibration. Weight units are now K^(-2).
;   12/01/2018 - change paths from /net/nas4/quijote/ to /net/nas/proyectos/quijote/
;   12/02/2018 - flagging from Diego's code implemented by JD instead of index of pixels.
;   05/03/2018 - typo corrected. Now writing wei after the recalibration.
;   19/03/2018 - new options added to remove long-term drifts (baseline model) and f(az) for RFI.
;   22/03/2018 - Added on output the pointingmodel info, and the AZhorn, ELhorn values
;   10/07/2019 - now including sigmacov at the pixel level. Already incorporated in the BTOD files
;   19/01/2021 - apply flag table revisited. Now checks if the bit is activated,
;                before adding it. Also apply_flags_transients and apply_flags_transients_byjd.
;
;
;-
; 

PRO global_ctod_params,par

  outpath    = "/net/nas/proyectos/quijote/ctod/"  ;raw/"  ;"crab_nofred_as2/"  ;"dipole2/" ; "gain-all-sum/"
  ;outpath    = "/net/nas/proyectos/quijote/ctod/test_baselines/"  ;"crab_nofred_as2/"  ;"dipole2/" ; "gain-all-sum/"
  ;outpath    = "/net/nas/proyectos/quijote/ctod/gain-sumc/"
  ;outpath    = "/net/nas/proyectos/quijote/ctod/gain-y/"

  ; Gain model: Y, X, SC (sum of corr), SNC (sum of uncorr),
  ;             AS (all referenced to the sum)
  gain_model = "AS2" ; "AS" ; "SC" ; "SNC" ; "Y" ; "AS"   

  
  ; Fit method. Default calibration method used by RGS.
  fit_method = 6
  
  ; Diode calibration. D(theta).
  ff = paths_for_quijote() + "etc/TAB_OF_AVERAGED_DIODE_FIT_POL_PAR_DRIFT_CORRECTED_RESULTS.txt" 
  tab_dtheta = ddread(ff)


  ; Diode temperature, T_0
  ff = paths_for_quijote() + "etc/TAB_OF_DIODE_TEMP_ESTIMATES_all.txt" 
  tab_t0 = ddread(ff)


  ; Modulators calibration, M(theta).
  ff = paths_for_quijote() + "etc/polarmod_transmission_allhorns.dat"
  tab_mtheta = ddread(ff)
  
  ; Diode calibration. CASS factors. Added on May 2015. 
  ff = paths_for_quijote() + "etc/raw_cal_factors_cass_from_1407.dat"
  tab_diode_cass = ddread(ff)


  ; Cass-A and diode factors. Added on May 2016
  ; Using: MED_GAINS_CASS  FLOAT     = Array[3, 5, 4, 2, 4] = [i,j,k,l,m]
  ; i=0,1,2 para los tres periodos que comentamos esta mañana
  ; j=0,1,2,3 para A,B,C y D, j=4 tomando todas las observaciones
  ; [4,2,4] es mi formato habitual.
  ;ff = "/home/jalberto/quijote/idl/sancho/gains_cass_amp_diode.sav"
  path_to_sav = "/net/nas/proyectos/quijote/etc/gains/" ;"/home/jalberto/quijote/idl/sancho/"
  ff = path_to_sav + "gains_cass_amp_diode_y.sav"
  restore,ff
  tab_med_amp_diode_y  = MED_AMP_DIODE
  tab_med_gains_cass_y = MED_GAINS_CASS
  ff = path_to_sav + "gains_cass_amp_diode_x.sav"
  restore,ff
  tab_med_amp_diode_x  = MED_AMP_DIODE
  tab_med_gains_cass_x = MED_GAINS_CASS
  ff = path_to_sav + "gains_cass_amp_diode_sumnc.sav"
  restore,ff
  tab_med_amp_diode_sumnc  = MED_AMP_DIODE
  tab_med_gains_cass_sumnc = MED_GAINS_CASS
  ff = path_to_sav + "gains_cass_amp_diode_sumc.sav"
  restore,ff
  tab_med_amp_diode_sumc  = MED_AMP_DIODE
  tab_med_gains_cass_sumc = MED_GAINS_CASS

  
  ; Maximum-minimum correction allowed
  MINMAX_CORRECTION = [0.8,1.2] ;[0.5,1.5]

  
  ; CAL files
  path_cal = "/net/nas/proyectos/quijote/cal/"


  ; GAIN files
  path_gain = "/net/nas/proyectos/quijote/gain/"

  ; BASE files
  path_base = "/net/nas/proyectos/quijote/base/"

  ; Transients directory (by index). Old version
  path_transients = "/net/nas/proyectos/quijote/btod_transient/stackXcheck_55/" ; 3rd version   - Sept 2017. It does not work well
  path_transients = "/net/nas/proyectos/quijote/btod_transient/with_signif/"    ; First version - June 2017. 
  path_transients = "/net/nas/proyectos/quijote/btod_transient/nogal/"          ; 2nd version - July 2017. Keep this as default

  ; Transients directory (by JD). New version
  path_byjd_transients = "/net/nas/proyectos/quijote/btod_transient/byjd_40ms/"

  
  ; structure
  par     = { outpath:outpath, tab_dtheta:tab_dtheta, tab_t0:tab_t0, tab_mtheta:tab_mtheta,$
              tab_diode_cass:tab_diode_cass, path_cal:path_cal, path_gain:path_gain, MINMAX_CORRECTION:MINMAX_CORRECTION, $
              tab_med_amp_diode_y:tab_med_amp_diode_y,         tab_med_gains_cass_y:tab_med_gains_cass_y, $
              tab_med_amp_diode_x:tab_med_amp_diode_x,         tab_med_gains_cass_x:tab_med_gains_cass_x, $
              tab_med_amp_diode_sumnc:tab_med_amp_diode_sumnc, tab_med_gains_cass_sumnc:tab_med_gains_cass_sumnc, $
              tab_med_amp_diode_sumc:tab_med_amp_diode_sumc,   tab_med_gains_cass_sumc:tab_med_gains_cass_sumc, $
              gain_model:gain_model, fit_method:fit_method, path_transients:path_transients, path_byjd_transients:path_byjd_transients,$
	      path_base:path_base }
  

END



;=====================
; Input angles are in degrees
function mfi_diode_cal, theta, ichan, tab_dtheta=tab_dtheta
  if not keyword_set(tab_dtheta) then begin
     global_ctod_params,par
     tab_dtheta = par.tab_dtheta
  endif

  if (ichan lt 0 or ichan gt 31) then begin
     salida = theta * 0.0
  endif else begin
     d2r    = !DPI/180.0
     AA     = tab_dtheta[1,ichan]
     BB     = tab_dtheta[3,ichan]
     salida = 1.0 + AA * cos(4.*theta*d2r) + BB * sin(4.*theta*d2r) 
  endelse

return,salida
end


;==========================

function mfi_diode_t0, ichan, tab_t0=tab_t0
  if not keyword_set(tab_t0) then begin
     global_ctod_params,par
     tab_t0 = par.tab_t0
  endif

  if (ichan lt 0 or ichan gt 31) then begin
     salida = 0.0
  endif else begin
     salida = tab_t0[1,ichan]   ; K
  endelse

return,salida
end

;=====================
; Input angles are in degrees
function mfi_polar_transmission, theta, ichan, tab_mtheta=tab_mtheta
  if not keyword_set(tab_mtheta) then begin
     global_ctod_params,par
     tab_mtheta = par.tab_mtheta
  endif

  if (ichan lt 0 or ichan gt 31) then begin
     salida = theta * 0.0
  endif else begin
     d2r    = !DPI/180.0
     AA     = tab_mtheta[1,ichan]
     BB     = tab_mtheta[3,ichan]
     salida = 1.0 + AA * cos(2.*theta*d2r) + BB * sin(2.*theta*d2r) 
  endelse

return,salida
end

;=====================

pro smooth_gain_old,gain,flag,gain_smooth

KERNEL_CAL_SMTH_MIN = 24.0 ;min. Kernel for smoothing the CAL data
TIME_ONE_CAL        = 31.0 ;sec. 

nx = round(KERNEL_CAL_SMTH_MIN/TIME_ONE_CAL*60.0) 
n  = n_elements(gain)
gain_smooth = gain

if (n gt 2*nx) then begin
   for i=0,n-1 do begin
      i1 = i-nx/2
      if (i1 lt 0) then i1=0
      i2 = i+nx/2
      if (i2 gt (n-1)) then i2=n-1
      datos = gain[i1:i2]
      peso  = flag[i1:i2]
      xxx   = datos[where(peso eq 0)]
      gain_smooth[i] = median(xxx)
   endfor
endif else begin
   cut = where(flag eq 0)
   gain_smooth[*] = median(gain[cut])
endelse
end

;---

pro smooth_gain, jd, gain, flag, gain_smooth, gain_stdev, tkernel=tkernel
  
if not keyword_set(tkernel) then tkernel = 120.0 ; minutes
  
d2r = !DPI/180.0

dt_min = tkernel                ; minutes
dt     = dt_min/60.0/24.0       ; jd  

; initialize to zero
n  = n_elements(gain)
gain_smooth = fltarr(n)
gain_stdev  = fltarr(n)

; loop
for i=0L,n-1L do begin
   cut = where(flag eq 0 and abs(jd-jd[i]) le dt/2.0, ncut)
   if ncut ge 10 then begin
      gain_smooth[i] = median(gain[cut])
      gain_stdev[i]  = stddev(gain[cut])
   endif
endfor


end


;============

pro read_mfi_gainmodel, root, jd, gains, corr, flag, par=par, JD_BTOD=JD_BTOD, MIN_TOBS = MIN_TOBS  

if not keyword_set(par) then global_ctod_params,par

  ; obtain GAIN file
  nlen   = strlen(root)
  root2  = ""
  if nlen ge 5 then root2  = strmid(root,0,nlen-4) ; removes -001 if present

  ff1 = par.path_gain + root  + ".gain" 
  ff2 = par.path_gain + root2 + ".gain" 
  if does_exists(ff2) then ff1 = ff2
  print," --> GAIN file: ",ff1


  ; Obtain smoothed gains only if they exist 
  if does_exists(ff1) then begin
     a     = mrdfits(ff1,1)
     jd    = a.jd
     n     = n_elements(jd)
     nchan = n_elements(a.gainmodel[*,0])
     gains = a.gainmodel
     corr  = a.oneplusdeltaG
     MIN_TOBS   = a.MIN_TOBS_HOURS 
     model      = a.model
     fit_method = a.fit_method
     flag       = a.flag
  endif else begin
     print, "No GAIN file . "
     if not keyword_set(JD_BTOD) then stop," No GAIN file. Please supply JD_BTOD. "

     ; Returns gains=0 in this case.
     jd    = minmax(JD_BTOD)
     n     = n_elements(jd)
     nchan = 32
     gains = fltarr(nchan,n)
     flag  = intarr(nchan,n) + 1
     corr  = gains
     MIN_TOBS = 2.0 ; hours 
  endelse
end

;=======================

pro get_mfi_cal_gains, root_cal, jd, gains, gains_orig, par=par, JD_BTOD=JD_BTOD

if not keyword_set(par) then global_ctod_params,par

  ; obtain cal file
  nlen   = strlen(root_cal)
  root2  = ""
  if nlen ge 5 then root2  = strmid(root_cal,0,nlen-4) ; removes -001 if present

  ffcal  = par.path_cal + root_cal + ".cal" 
  ffcal2 = par.path_cal + root2    + ".cal" 
  if does_exists(ffcal2) then ffcal = ffcal2
  print," --> CAL file: ",ffcal


  ; Obtain smoothed gains only if they exist 
  if does_exists(ffcal) then begin
     a     = mrdfits(ffcal,1)
     jd    = a.jd
     n     = n_elements(jd)
     nchan = n_elements(a.gain[*,0])
     gains = fltarr(nchan,n)
     for ichan=0,nchan-1 do begin

        ;smooth_gain_old, a.gain[ichan,*], a.flag[ichan,*], gains_v2
        smooth_gain, jd, a.gain[ichan,*], a.flag[ichan,*], gains_v, gains_stdev_v

        gains[ichan,*] = gains_v[*]

     endfor
     gains_orig = a.gain

     ; Check ranges of GAINs. If they deviate too much from typical values,
     ; then exclude data.
     ; not done yet


  endif else begin
     print, "No CAL file. "
     if not keyword_set(JD_BTOD) then stop," No CAL file. Please supply JD_BTOD. "

     ; Returns gains=0 in this case.
     jd    = minmax(JD_BTOD)
     n     = n_elements(jd)
     nchan = 32
     gains = fltarr(nchan,n)
     gains_orig = gains
     
  endelse
  
end

;=================

pro mfi_interpolate_gains, jd, g, jd2, g2, flag, MIN_TOBS = MIN_TOBS

if not keyword_set(MIN_TOBS) then MIN_TOBS = 2.0 ; hr

jdref = jd2[0]

nx = n_elements(jd)
x  = jd  - jdref

nx2 = n_elements(jd2) 
x2  = jd2 - jdref

nchan = n_elements(g[*,0])

g2 = fltarr(nchan, nx2)

; Short observations: use fit to a straight line.
t_obs    = (max(jd) - min(jd) )*24.0

IF (t_obs le MIN_TOBS) THEN BEGIN

   ; Short observation. Fit to straight line.
   for i=0,nchan-1 do begin
      if n_params() eq 5 then begin
         cut = where(g[i,*] gt 0.0 and flag[i,*] eq 0, ncut)
      endif else begin
         cut  = indgen(nx)
         ncut = nx
      endelse
      if (ncut ge 2) then begin
         ;help,jd,g,flag
         y    = reform(g[i,*])
         coef = POLY_FIT( x[cut], y[cut], 1)
         g2[i,*] = coef[0] +  coef[1]*x2
      endif
   endfor
   print," INFO MFI_INTERPOLATE_GAINS: short observation, with t(hrs) < ",MIN_TOBS
   
ENDIF ELSE BEGIN

   ; Long observation. Spline fitting. Not using flag
   for i=0,nchan-1 do begin
      Y  = reform(g[i,*])
      if (total(Y) ne 0.0) then begin
         Y2 = SPL_INIT(x, Y)
         sp = SPL_INTERP(x, Y, Y2, x2)
         g2[i,*] = sp          
      endif
   endfor
   print," INFO MFI_INTERPOLATE_GAINS: long observation. T_obs(hr) = ",t_obs

ENDELSE

end

;==========

pro identify_mfi_mod_position, ang_v, k

d2r   = !DPI/180.0
casos = [0.0, 22.5, 45.0, 67.5]

; NOTE. Before, I was taking ang = median(ang_v), but for the
; case of A position, it could produce errors if we have values
; like 0.1deg and 359.9 deg.
ang = ang_v[0]
ang = (ang + 360.0) mod 360.0

k   = -1
for j=0,3 do begin
   dist = ACOS( COS( (ang - casos[j])*d2r ) )/d2r
   if ( abs(dist) le 1.0) then k=j
endfor

end


;=======================

function nlines,file,help=help
on_error,2
if keyword_set(help) then begin & doc_library,'nlines' & return,0 & endif
tmp = ' '
nl = 0L
on_ioerror,NOASCII
if n_elements(file) eq 0 then file = pickfile()
openr,lun,file,/get_lun
while not eof(lun) do begin
  readf,lun,tmp
  nl = nl + 1L
  endwhile
close,lun
free_lun,lun
NOASCII:
return,nl
end

;=====================

pro apply_flag_table, filename, jd, flag, bit=bit
  if not keyword_set(bit) then bit = 2

  ; name of the flag table file
  kk1 = strsplit(filename,"/",/extract)
  kk2 = kk1[n_elements(kk1)-1]
  kk3 = strsplit(kk2,"-",/extract)

  if n_elements(kk3) ge 3 then begin
     root = kk3[0] + "-" + kk3[1] + "-" + kk3[2]
     ff_flag = "/net/nas/proyectos/quijote/flag/" + root + ".flag"
  endif else begin
     print," --> Error in the name format! "
     return
  endelse

  ; Read flagging table if exists, and if containing lines
  existe = does_exists(ff_flag)
  nlines = 0
  if (existe) then nlines = nlines(ff_flag)
  if nlines eq 0 then existe = 0
  
  if (existe) then begin
     print," > Reading flagging table ",ff_flag
     readcol,ff_flag,kk,ichan,jdmin,jdmax,format="A,I,D,D"
     if n_elements(kk) ge 1 then existe = 1
  endif

  ; Apply
  if (existe) then begin
     one = 2^(bit-1)
     sancho_check_activebit_inflags, flag, bit, bitmasc
     for i=0,n_elements(ichan)-1 do begin
        cut = where(jd ge jdmin[i] and jd le jdmax[i], ncut)
        ; single channel
        if (ichan[i] ge  0 and ncut ge 1) then flag[ichan[i],cut] = flag[ichan[i],cut] + one * (1-bitmasc[ichan[i],cut])
        ; all channels
        if (ichan[i] eq -1 and ncut ge 1) then flag[*,cut] = flag[*,cut] + one * (1-bitmasc[*,cut])
     endfor
  endif

end


;=====================
; Applies the flags from Diego's code, generated by Flavien
;

pro apply_flags_transients, filename, jd, flag, bit=bit, path=path, mfi=mfi
  if not keyword_set(bit) then bit = 2
  if not keyword_set(path) then begin
     global_ctod_params, par
     path = par.path_transients
  endif
  if not keyword_set(mfi) then info_quijote, mfi=mfi

  ; name of the transients file
  kk1 = strsplit(filename,"/",/extract)
  kk2 = kk1[n_elements(kk1)-1]
  kk3 = strsplit(kk2,".",/extract)

  root = kk3[0]
  ff   = path + root + "_transient.fits"

  
  ; Read flagging table if exists, and if containing lines
  existe = does_exists(ff)
  vacio  = 0 
  
  if (existe) then begin
     print," > Reading transients file ",ff
     hdr    = headfits(ff,ext=1)
     naxis1 = SXPAR( hdr, "NAXIS1") 
     if naxis1 eq 0 then begin
        print,"     This file is empty: ",ff
        return
     endif else begin
        a     = mrdfits(ff,1)
        names = TAG_NAMES(a)
     endelse
  endif else begin
     print,"     This file is absent: ",ff
     return
  endelse

  ; apply
  if (existe) then begin
     one = 2^(bit-1)
     sancho_check_activebit_inflags, flag, bit, bitmasc
     for ichan=0,mfi.nchan-1 do begin
        if ichan le 9 then fullname_ichan = "CHANNEL0" + string(ichan,format="(i1)")
        if ichan gt 9 then fullname_ichan = "CHANNEL"  + string(ichan,format="(i2)")

        j = where(names eq fullname_ichan, ncut)
        
        if ncut eq 1 then begin
           ii  = j[0]
           idx = a.(ii)
           if n_elements(idx) ge 1 then flag[ichan,idx] = flag[ichan,idx] + one * (1-bitmasc[ichan,idx]) ; new. Only if bit is not active
        endif
        
     endfor
  endif

end


;=====================
; Applies the flags from Diego's code BY JD (new version 12/Feb/2018)
;
; Format: JD_start   JD_end   JD_peak    amplitude   GL(deg)   GB(deg)
; 

pro apply_flags_transients_byjd, filename, jd, flag, bit=bit, path=path, mfi=mfi, noflag=noflag
  if not keyword_set(bit) then bit = 2
  if not keyword_set(path) then begin
     global_ctod_params, par
     path = par.path_byjd_transients
  endif
  if not keyword_set(mfi) then info_quijote, mfi=mfi

  ; name of the transients file
  kk1 = strsplit(filename,"/",/extract)
  kk2 = kk1[n_elements(kk1)-1]
  kk3 = strsplit(kk2,".",/extract)

  root = kk3[0]
  ff   = path + root + "_transient.fits"

  
  ; Read flagging table if exists, and if containing lines
  existe = does_exists(ff)
  vacio  = 0 
  
  if (existe) then begin
     print," > Reading transients file by JD ",ff
     hdr    = headfits(ff,ext=1)
     naxis1 = SXPAR( hdr, "NAXIS1") 
     if naxis1 eq 0 then begin
        print,"     This file is empty: ",ff
        return
     endif else begin
        a     = mrdfits(ff,1)
        names = TAG_NAMES(a)
     endelse
  endif else begin
     print,"     This file is absent: ",ff
     return
  endelse

  ; apply flags
  if (existe) then begin
     one = 2^(bit-1)
     sancho_check_activebit_inflags, flag, bit, bitmasc
     for ichan=0,mfi.nchan-1 do begin
        if ichan le 9 then fullname_ichan = "CHANNEL0" + string(ichan,format="(i1)")
        if ichan gt 9 then fullname_ichan = "CHANNEL"  + string(ichan,format="(i2)")

        j = where(names eq fullname_ichan, ncut)
        
        if ncut eq 1 then begin
           ii   = j[0]
	   data = a.(ii) ; Format: JD_start   JD_end   JD_peak    amplitude   GL(deg)   GB(deg)
           jdv  = reform(data[0:1,*]) ; vector of JD ranges
	   njd  = n_elements(jdv[0,*])
	   glv  = reform(data[4,*])
	   gbv  = reform(data[5,*])
           check_if_in_mask, glv,gbv, noflag, status
	   for k=0,njd-1 do begin
		cutjd = where(jd ge jdv[0,k] AND jd le jdv[1,k], ncutjd)
		if (ncutjd ge 1  and status[k] eq 0) then flag[ichan,cutjd] = flag[ichan,cutjd] + one * (1-bitmasc[ichan,cutjd])
	   endfor
        endif
        
     endfor
  endif
end
;---
pro check_if_in_mask, glv,gbv, noflag, status
nside = 128
d2r   = !DPI/180.0
theta = (90.0-gbv)*d2r
phi   = glv*d2r
ANG2PIX_RING, nside, theta, phi, ipix
status = noflag[ipix]
end

;================

pro select_values_for_gain_correction_old, caso, mfi, params, jperiod, imodpos, ihorn, ifreq, itype, gains_ic, amp_diode, $
                                       amp_diode_cal, iref_gain, ref_gain, G_cass, iref_gain_d, rfactor

  CASE caso OF
     "Y": begin 
        iref_gain = mfi.ic[ihorn,ifreq,3] ; Y 
        ref_gain  = params.tab_med_amp_diode_y[ jperiod, imodpos, ihorn, ifreq, 3]
        G_cass    = params.tab_med_gains_cass_y[jperiod,       4, ihorn, ifreq, itype] ; avg in A,B,C,D
     end
     "Y-mod": begin
        iref_gain = mfi.ic[ihorn,ifreq,3] ; Y 
        ref_gain  = params.tab_med_amp_diode_y[ jperiod, imodpos, ihorn, ifreq, 3]
        G_cass    = params.tab_med_gains_cass_y[jperiod, imodpos, ihorn, ifreq, itype] ; separates in A,B,C,D
     end
     "X": begin 
        iref_gain = mfi.ic[ihorn,ifreq,2] ; X 
        ref_gain  = params.tab_med_amp_diode_x[ jperiod, imodpos, ihorn, ifreq, 2]
        G_cass    = params.tab_med_gains_cass_x[jperiod,       4, ihorn, ifreq, itype] ; avg in A,B,C,D
     end
     "SC": begin                                   ; All referenced to the SUM of corr 
        iref_gain   = mfi.ic[ihorn,ifreq,0]        ; X+Y
        iref_gain_d = mfi.ic[ihorn,ifreq,1]        ; X-Y
        rfactor     = params.tab_med_gains_cass_sumc[jperiod, 4, ihorn, ifreq, 1] / $ 
                      params.tab_med_gains_cass_sumc[jperiod, 4, ihorn, ifreq, 0] ; d/s 
        ref_gain    = params.tab_med_amp_diode_sumc[ jperiod, imodpos, ihorn, ifreq, 0] + rfactor * $
                      params.tab_med_amp_diode_sumc[ jperiod, imodpos, ihorn, ifreq, 1]  
        G_cass      = params.tab_med_gains_cass_sumc[jperiod,       4, ihorn, ifreq, itype] ; avg in A,B,C,D
     end
     "SNC": begin                                  ; All referenced to the SUM of non corr 
        iref_gain   = mfi.ic[ihorn,ifreq,2]        ; X
        iref_gain_d = mfi.ic[ihorn,ifreq,3]        ; Y
        rfactor     = params.tab_med_gains_cass_sumnc[jperiod, 4, ihorn, ifreq, 3] / $ 
                      params.tab_med_gains_cass_sumnc[jperiod, 4, ihorn, ifreq, 2] ; y/x 
        ref_gain    = params.tab_med_amp_diode_sumnc[ jperiod, imodpos, ihorn, ifreq, 2] + rfactor * $
                      params.tab_med_amp_diode_sumnc[ jperiod, imodpos, ihorn, ifreq, 3]  
        G_cass      = params.tab_med_gains_cass_sumnc[jperiod,       4, ihorn, ifreq, itype] ; avg in A,B,C,D
        
     end

     "AS": begin                ; All referenced to the SUM, corr and uncorr with their respective cases.
        if (itype eq 0 or itype eq 1) then begin
           iref_gain   = mfi.ic[ihorn,ifreq,0]        ; X+Y
           iref_gain_d = mfi.ic[ihorn,ifreq,1]        ; X-Y
           rfactor     = params.tab_med_gains_cass_sumc[jperiod, 4, ihorn, ifreq, 1] / $ 
                         params.tab_med_gains_cass_sumc[jperiod, 4, ihorn, ifreq, 0] ; d/s 
           ref_gain    = params.tab_med_amp_diode_sumc[ jperiod, imodpos, ihorn, ifreq, 0] + rfactor * $
                         params.tab_med_amp_diode_sumc[ jperiod, imodpos, ihorn, ifreq, 1]  
           G_cass      = params.tab_med_gains_cass_sumc[jperiod,       4, ihorn, ifreq, itype] ; avg in A,B,C,D
        endif 
        if (itype eq 2 or itype eq 3) then begin
           iref_gain   = mfi.ic[ihorn,ifreq,2]        ; X
           iref_gain_d = mfi.ic[ihorn,ifreq,3]        ; Y
           rfactor     = params.tab_med_gains_cass_sumnc[jperiod, 4, ihorn, ifreq, 3] / $ 
                         params.tab_med_gains_cass_sumnc[jperiod, 4, ihorn, ifreq, 2] ; y/x 
           ref_gain    = params.tab_med_amp_diode_sumnc[ jperiod, imodpos, ihorn, ifreq, 2] + rfactor * $
                         params.tab_med_amp_diode_sumnc[ jperiod, imodpos, ihorn, ifreq, 3]  
           G_cass      = params.tab_med_gains_cass_sumnc[jperiod,       4, ihorn, ifreq, itype] ; avg in A,B,C,D
        endif 
     end            
     "AS2": begin               ; All referenced to the SUM, corr and uncorr with their respective cases. ref_gain to avg of ABCD
        if (itype eq 0 or itype eq 1) then begin
           iref_gain   = mfi.ic[ihorn,ifreq,0]        ; X+Y
           iref_gain_d = mfi.ic[ihorn,ifreq,1]        ; X-Y
           rfactor     = params.tab_med_gains_cass_sumc[jperiod, 4, ihorn, ifreq, 1] / $ 
                         params.tab_med_gains_cass_sumc[jperiod, 4, ihorn, ifreq, 0] ; d/s 
           ref_gain    = params.tab_med_amp_diode_sumc[ jperiod, 4, ihorn, ifreq, 0] + rfactor * $
                         params.tab_med_amp_diode_sumc[ jperiod, 4, ihorn, ifreq, 1]  
           G_cass      = params.tab_med_gains_cass_sumc[jperiod, 4, ihorn, ifreq, itype] ; avg in A,B,C,D
        endif 
        if (itype eq 2 or itype eq 3) then begin
           iref_gain   = mfi.ic[ihorn,ifreq,2]        ; X
           iref_gain_d = mfi.ic[ihorn,ifreq,3]        ; Y
           rfactor     = params.tab_med_gains_cass_sumnc[jperiod, 4, ihorn, ifreq, 3] / $ 
                         params.tab_med_gains_cass_sumnc[jperiod, 4, ihorn, ifreq, 2] ; y/x 
           ref_gain    = params.tab_med_amp_diode_sumnc[ jperiod, 4, ihorn, ifreq, 2] + rfactor * $
                         params.tab_med_amp_diode_sumnc[ jperiod, 4, ihorn, ifreq, 3]  
           G_cass      = params.tab_med_gains_cass_sumnc[jperiod, 4, ihorn, ifreq, itype] ; avg in A,B,C,D
        endif 
     end  
     
  ENDCASE

end


; New Code, based on mfi_params.pro
pro select_values_for_gain_correction, caso, mfi, params, jperiod, imodpos, ihorn, ifreq, itype, gains_ic, amp_diode, $
                                       amp_diode_cal, iref_gain, ref_gain, G_cass, iref_gain_d, rfactor

  CASE caso OF
     "AS": begin                ; All referenced to the SUM, corr and uncorr with their respective cases.
        if (itype eq 0 or itype eq 1) then begin
           iref_gain   = mfi.ic[ihorn,ifreq,0]                                 ; X+Y
           iref_gain_d = mfi.ic[ihorn,ifreq,1]                                 ; X-Y
           rfactor     = gains_ic[ihorn, ifreq, 1] / gains_ic[ihorn, ifreq, 0] ; d/s 
           ref_gain    = amp_diode_cal[ imodpos, ihorn, ifreq, 0] + rfactor * amp_diode_cal[ imodpos, ihorn, ifreq, 1]  
           G_cass      = gains_ic[ihorn, ifreq, itype] ; avg in A,B,C,D
        endif 
        if (itype eq 2 or itype eq 3) then begin
           iref_gain   = mfi.ic[ihorn,ifreq,2]                                 ; X
           iref_gain_d = mfi.ic[ihorn,ifreq,3]                                 ; Y
           rfactor     = gains_ic[ihorn, ifreq, 3] / gains_ic[ihorn, ifreq, 2] ; y/x 
           ref_gain    = amp_diode_cal[ imodpos, ihorn, ifreq, 2] + rfactor * amp_diode_cal[ imodpos, ihorn, ifreq, 3]  
           G_cass      = gains_ic[ihorn, ifreq, itype] ; avg in A,B,C,D
        endif 
     end            
     "AS2": begin               ; All referenced to the SUM, corr and uncorr with their respective cases. ref_gain to avg of ABCD
        if (itype eq 0 or itype eq 1) then begin
           iref_gain   = mfi.ic[ihorn,ifreq,0]                                 ; X+Y
           iref_gain_d = mfi.ic[ihorn,ifreq,1]                                 ; X-Y
           rfactor     = gains_ic[ihorn, ifreq, 1] / gains_ic[ihorn, ifreq, 0] ; d/s 
           ref_gain    = amp_diode[ ihorn, ifreq, 0] + rfactor * amp_diode[ ihorn, ifreq, 1]  
           G_cass      = gains_ic[ihorn, ifreq, itype] ; avg in A,B,C,D
        endif 
        if (itype eq 2 or itype eq 3) then begin
           iref_gain   = mfi.ic[ihorn,ifreq,2]                                 ; X
           iref_gain_d = mfi.ic[ihorn,ifreq,3]                                 ; Y
           rfactor     = gains_ic[ihorn, ifreq, 3] / gains_ic[ihorn, ifreq, 2] ; y/x 
           ref_gain    = amp_diode[ ihorn, ifreq, 2] + rfactor * amp_diode[ ihorn, ifreq, 3]  
           G_cass      = gains_ic[ihorn, ifreq, itype] ; avg in A,B,C,D
        endif 
     end     
  ENDCASE

end


;=================

pro read_mfi_baselinemodel, root, jd, bases, flag, par=par, JD_BTOD=JD_BTOD, MIN_TOBS = MIN_TOBS  

if not keyword_set(par) then global_ctod_params,par

  ; obtain base file
  nlen   = strlen(root)
  root2  = ""
  if nlen ge 5 then root2  = strmid(root,0,nlen-4) ; removes -001 if present

  ff1 = par.path_base + root  + ".base" 
  ff2 = par.path_base + root2 + ".base" 
  if does_exists(ff2) then ff1 = ff2
  print," --> BASE file: ",ff1


  ; Obtain smoothed baselines only if they exist 
  if does_exists(ff1) then begin
     a     = mrdfits(ff1,1)
     jd    = a.jd + jd_btod[0]
     n     = n_elements(jd)
     nchan = n_elements(a.baseline[*,0])
     bases = a.baseline
     flag  = fltarr(nchan,n)
     cut0  = where(bases eq 0,n0)
     if n0 ge 1 then flag[cut0] = 1
     MIN_TOBS   = a.MIN_TOBS_HOURS 
  endif else begin
     print, "No BASE file . "
     if not keyword_set(JD_BTOD) then stop," No BASE file. Please supply JD_BTOD. "

     ; Returns gains=0 in this case.
     jd    = minmax(JD_BTOD)
     n     = n_elements(jd)
     nchan = 32
     bases = fltarr(nchan,n)
     flag  = fltarr(nchan,n)
     MIN_TOBS = 2.0 ; hours 
  endelse
end



;=================
;=================
;=================

PRO prepare_ctod_for_mfi, root, overwrite=overwrite, path=path, mfi=mfi, params=params, rms_flag=rms_flag, verbose=verbose, $
                          nomedian=nomedian, display=display, int_rms_flag=int_rms_flag, noflag=noflag, baseline=baseline, rfifaz=rfifaz, outpath=outpath,$
				keeptransients=keeptransients

; Keyword
if not keyword_set(overwrite) then overwrite = 0
if not keyword_set(path) then path = ""
if not keyword_set(rms_flag) then rms_flag = 0
if not keyword_set(int_rms_flag) then int_rms_flag = 0
if not keyword_set(verbose) then verbose = 0
if not keyword_set(nomedian) then nomedian = 0
if not keyword_set(display) then display = 0
if not keyword_set(noflag) then noflag = fltarr(nside2npix(128))
if not keyword_set(baseline) then baseline = 0
if not keyword_set(rfifaz) then rfifaz = 0
if not keyword_set(keeptransients) then keeptransients = 0


; Get global parameters
if not keyword_set(params) then global_ctod_params, params
if not keyword_set(mfi) then info_quijote, mfi=mfi
jd_ref    = mfi.jd_ref
jd_140411 = 2456758.5d0
jd_151201 = 2457357.5d0
jd_160301 = 2457448.5d0
jd_160501 = 2457509.5d0
jd_160528 = 2457536.5d0
jd_161015 = 2457676.5d0


; Obtain all the CTOD files associated
print," (*) MAIN PATH to btod: ",path
print,"     ROOT name: ",root
spawn,"ls -1 " + path + root + "*.btod", ff_list, error

if error ne "" then return


; Count TOD maps
ntod = n_elements(ff_list)
print," (*) Total number of BTOD files: ",ntod



; MAIN LOOP over blocks
cfp
if not keyword_set(outpath) then outpath = params.outpath
FOR itod=0,ntod-1L DO BEGIN

   filename = ff_list[itod]

   ; Output ctod file. I do a test to the file name
   kk    = strsplit(filename,"/",/extract)
   kk2   = kk[n_elements(kk)-1] 
   root2 = strtrim(strmid(kk2,0,strlen(kk2)-5),2)
   
   txt       = strtrim(string(itod),2)
   case strlen(txt) of
      1: tail = "-00"+txt 
      2: tail = "-0"+txt
      3: tail = "-"+txt
   endcase
   if (ntod eq 1) then tail = ""
   root3 = strtrim(root,2) + tail

   if (root2 eq root3) then begin
      ff_out   = outpath + root3 + ".ctod"
   endif else begin
      print," Problem with the naming of the file. "
      stop
   endelse

   
   ; Check nominal mod position
   kk3      = strsplit(root2,"-",/extract)
   campo    = strtrim(kk3[0],2)
   txtmod   = strmid(campo,0,1,/reverse) 
   casos    = ["A","B","C","D","x"]
   imodnomi = where(casos eq txtmod) 


   ; Check if the output file already exists
   existe  = does_exists(ff_out)
   goahead = 1
   if (existe) then goahead = 0
   if (existe and overwrite) then goahead = 1


   ; Main part of code
   if (goahead) then begin

      print," --> Reading btod: ",filename

      ;-----------------------
      ; [1] Read data
      a     = mrdfits(filename,1,/silent)
      names = tag_names(a)
      data  = a.data
      wei   = a.wei
      flag  = a.flag
      jd    = a.jd
      njd   = n_elements(jd)
      msbin = a.msbin
      t     = (jd - jd[0])*24.0 ; hrs

      ; Get AZ_horn, EL_horn if defined
      if total(names eq "AZHORN") gt 0 then begin
      	az_horn = a.azhorn
      	el_horn = a.elhorn
      endif else begin
	az_horn = fltarr(4,njd)
	el_horn = fltarr(4,njd)
      endelse       

      ; Get pointing model 
      pointingmodel = "Version1"
      if total(names eq "POINTINGMODEL") gt 0 then pointingmodel = a.pointingmodel

      ; Get sigmacov. Added on 10/07/2019.
      if total(names eq "SIGMACOV") gt 0 then begin
	sigmacov = a.sigmacov
      endif else begin
	print," No SIGMACOV was found. Setting SIGMACOV=0."
	sigmacov = fltarr(4,2,2,njd)
      endelse

      
      ; Update flags from Fred code
      if (rms_flag ne 0) then begin
         print,"     Updating flag with Fred's code and rms_flag = ",rms_flag
;         mfi_flag_stats,filename=root,f_data=1,flag=flag,jd=jd,bin=rms_flag
         mfi_flag_stats2,filename=root2,f_diff=1,flag=flag,jd=jd,bin=rms_flag
      endif else begin
         print,"     Not using flags from Fred's code. Keyword rms_flag = ",rms_flag
      endelse

      if (int_rms_flag ne 0) then begin
         print,"     Updating flag with Fred's code for intensity and int_rms_flag = ",int_rms_flag
         mfi_flag_stats,filename=root2,f_data=1,flag=flag,jd=jd,bin=int_rms_flag
      endif else begin
         print,"     Not using flags from Fred's code for intensity. Keyword rms_flag = ",int_rms_flag
      endelse


      
      ; Global gains and diode cal factors 
      fit_method = params.fit_method
      mfi_params, MOD_REF_ANGLE=theta0, JD_OBS=jd[0], gains_ic=gains_ic_all, /vector_ic, fit_method=fit_method
      mfi_params, JD_OBS=jd[0], CHAN_REF=1, FIT_METHOD=fit_method, GAINS_IC=gains_ic, GAINS_MOD_POS_IC=gic_permod,$
                  AMP_DIODE=amp_diode, AMP_MOD_POS_DIODE=amp_diode_cal 

      
      ; [1.1] Derive the angle for each measurement: 4*(mod_ang-theta0) + 2*par
      ; NOTE on Jan/2019: with theta0 depending on frequency, there are two phi values in reality per horn.
      ; Thus, phi is obsolete. It should not be used. Indeed, it is not written
      ;;phi   = 4.0 * (a.mod_ang - replicate(1,njd)##theta0[*,0]) + 2.0*a.par 
      
      ; [1.2] Load flagging tables and additional flags
      print," (*) Loading and applying flagging tables (eye inspection). "
      apply_flag_table, filename, jd, flag, bit=8 ; bit=8 is defined by the user.
     
      ; [1.3] Flags from Diego's code. Only for old binning at the moment. 
      ;if (msbin gt 30.0) then begin
      	;apply_flags_transients, filename, jd, flag, bit=8, path=params.path_transients ; bit=8 is defined by the user.
      ;endif
      if (keeptransients eq 0) then begin
      	print," (*) Loading and applying flags from Diegos code (compact objects). "
      	apply_flags_transients_byjd, filename, jd, flag, bit=8, noflag=noflag, path=params.path_byjd_transients ; bit=8 is defined by the user.
      endif else begin
      	print," (*) KEEPTRANSIENTS activated --> Flags from Diegos code (compact objects) are NOT applied. "
      endelse

      
      ;-----------------------
      ; [2] Gain correction from .gain files
      ;get_mfi_cal_gains, root, jd_cal, gains, gains_orig, par=params, JD_BTOD=jd
      read_mfi_gainmodel, root, jd_cal, gains, oneplusdeltaG_cal, flag_cal, par=params, JD_BTOD=jd, MIN_TOBS = MIN_TOBS
      mfi_interpolate_gains, jd_cal, gains, jd, gains2, flag_cal, MIN_TOBS = MIN_TOBS
      mfi_interpolate_gains, jd_cal, oneplusdeltaG_cal,  jd, oneplusdeltaG_matrix, flag_cal, MIN_TOBS = MIN_TOBS

      
      ; Display gain model correction
      if (display) then begin
         window,0,xs=1400,ys=900
         !p.multi=[0,1,4]
         !p.charsize=1.7
         for ihorn=0,3 do begin
            txth = ", horn "+string(ihorn+1,format="(i1)") 
            i1   = ihorn*8
            i2   = i1+3
            yran = [0.8,1.2] 
            xran = minmax(t)
            for i=i1,i2 do begin
               cut = where(flag[i,*] eq 0, ncut)
               if ncut eq 0 then cut = indgen(n_elements(t))
               if (i eq i1) then plot,t[cut],oneplusdeltaG_matrix[i,cut],yr=yran,xr=xran, /xsty, $
                                      title=root+txth,xtitle="time (h)",ytitle=" 1 + delta_G",psym=3,charsize=1.7
               if (i ne i1) then oplot,t[cut],oneplusdeltaG_matrix[i,cut],col=i-i1+2,psym=3
               if (i eq i1) then oplot,(jd_cal- jd[0])*24.0,oneplusdeltaG_cal[i,*],linestyle=2,thick=3
               if (i ne i1) then oplot,(jd_cal- jd[0])*24.0,oneplusdeltaG_cal[i,*],linestyle=2,col=i-i1+2,thick=3
            endfor
         endfor
      endif
   

      ; [2.1] Calculate gain period. Obsolete. Not used now. 
      jperiod = 0
      if (jd[0] ge jd_140411) then jperiod = 1
      if (jd[0] ge jd_151201) then jperiod = 2
      ; New computation of period. 
      find_mfi_period,jd[0],periodo
      
          
      ; [2.2] Amply and calibrate in temperature.
      ctod      = data
      print,"*  ICHAN IHORN IFREQ  MOD      minG         maxG         medG        G_cass     ichan"
      
      contador_sigcov = intarr(4,2,2) 
      for ichan=0,mfi.nchan-1 do begin

         ; Get ihorn, ifreq and itype
         ihorn = mfi.polarizer_id[ichan] - 1
         freq  = mfi.central_freq_id[ichan]
         ifreq = 0                                        ; low-freq
         if (freq eq 13.0 or freq eq 19.0) then ifreq = 1 ; high-freq
         canales = ["s","d","x","y"]
         cut     = where(canales eq mfi.chan_id[ichan])
         itype   = cut[0]
	 icorr   = 0                ; 0 and 1 are correlated.
	 if itype ge 2 then icorr=1 ; 2 and 3 are uncorrelated
         contador_sigcov[ihorn,ifreq,icorr] = contador_sigcov[ihorn,ifreq,icorr] + 1
         
         ; Check Mod position
         identify_mfi_mod_position, a.mod_ang[ihorn,*], imodpos
         
         if (ihorn ge 1 and imodnomi ne imodpos) then print," FILE name and imodpos do not agree. Using imodpos from mod_angle. "
         if (ihorn eq 0 and jd[0] ge jd_140411 and imodpos ne 0) then begin
            print," -- WARNING -- Error in horn1"
            imodpos = 0
         endif
         if (imodpos eq -1) then print," -- WARNING -- imodpos is different from A,B,C,D. "


         ; G_cass factor. Note that it is insensitive to the value of imodpos
         G_cass = gains_ic[ihorn, ifreq, itype] ; avg in A,B,C,D
         
         
         ; Compute gain correction
         oneplusdeltaG = reform( oneplusdeltaG_matrix[ichan,*] )
         minG = min(oneplusdeltaG)
         maxG = max(oneplusdeltaG)
         medG = median(oneplusdeltaG)


         printf,-1, ichan, ihorn, ifreq, imodpos, minG, maxG, medG, G_cass, mfi.ic[ihorn,ifreq,itype], format='(4i6,4f13.5,i6)'

         
         ; Check range of correction.
         do_correction = 1
;         if (medG lt params.MINMAX_CORRECTION[0] or medG gt params.MINMAX_CORRECTION[1] ) then do_correction = 0
         if (minG lt params.MINMAX_CORRECTION[0] or maxG gt params.MINMAX_CORRECTION[1] ) then do_correction = 0
         if do_correction eq 0 then begin
            print," * Correction out of range. Data are set to zero. "
            oneplusdeltaG = fltarr(njd)
         endif

         
         ; Correcting data and weights
         cut  = where(oneplusdeltaG ne 0.0,ncut)
         if ncut ge 1 then begin
		ctod[ichan,cut] = ctod[ichan,cut] * (G_cass / oneplusdeltaG[cut])    ; K
		wei[ichan,cut]  = wei[ichan,cut]  * (oneplusdeltaG[cut] / G_cass)^2. ; K^(-2)  
		sigmacov[ihorn,ifreq,icorr,cut] = sigmacov[ihorn,ifreq,icorr,cut] * (G_cass / oneplusdeltaG[cut])  ; K^2, after multiplying the two terms
	 endif
         cut0 = where(oneplusdeltaG eq 0.0,ncut0)
         if ncut0 ge 1 then begin
		ctod[ichan,cut0] = 0.0
		wei[ichan,cut0]  = 0.0
		sigmacov[ihorn,ifreq,icorr,cut0] = 0.0 
         endif
                  
         ; tests
         ;factor = 0.0
         ;if (imodpos ge 0 and imodpos le 3) then factor = params.tab_diode_cass[imodpos,ichan]
         ;ctod[ichan,*] = ctod[ichan,*] * factor

      endfor
      ; Health check. Important for the correct normalization of sigmacov. Added on 10/07/2019
      if total(contador_sigcov-2) ne 0 then stop," Error in running over all possible channels. " 

      
      ;------------
      ; [3] Interference correction
      ;
      sancho_get_mfi_subscans, root, a.az, a.el, k_ini, k_end, obsmode=obsmode
      check_if_nominal, a.az, is_nominal
      if (is_nominal and obsmode ne 5) then stop,"* Error in obsmode"
      if (obsmode eq 5 and is_nominal eq 0) then stop,"* Error in obsmode"

      if rfifaz then begin
	print," > Doing RFI correction by substracting a function of AZ, f(AZ). "
	sancho_get_rfifazmodel, root, periodo, obsmode, xaz, faz
      endif else begin
      	print," > No interference correction at this stage. "
      endelse

      

      ;------------
      ; [4] Baseline correction. 

      if (baseline and is_nominal) then begin
	; Removing baseline as a function of time
        read_mfi_baselinemodel, root, jd_base, bases, flag_base, par=params, JD_BTOD=jd, MIN_TOBS = MIN_TOBS
        mfi_interpolate_gains, jd_base, bases, jd, bases2, flag_base, MIN_TOBS = MIN_TOBS

        cut  = where(bases2 ne 0.0 and ctod ne 0.0, ncut)
	ctod[cut] = ctod[cut] - bases2[cut] 

	cut0 = where(ctod eq 0 and bases2 eq 0.0, ncut0)
	ctod[cut0] = 0.0 
      endif


      ; 4.1 Based on the median in each sub-scan. 
      nscans       = n_elements(k_ini)
      GBCUT        = 20.0

      jd_per_scan  = fltarr(nscans) 
      med_per_scan = fltarr(mfi.nchan, nscans) 
      gbmat = fltarr(mfi.nchan, njd)
      for i= 0, 7 do gbmat[i,*] = a.gb[0,*]
      for i= 8,15 do gbmat[i,*] = a.gb[1,*]
      for i=16,23 do gbmat[i,*] = a.gb[2,*]
      for i=24,31 do gbmat[i,*] = a.gb[3,*]

      ctod_nomed = ctod * 0.0
      ctod_orig  = ctod
      
      ; 4.1.1 Removes the median. For all database unless user specifies the opposite
      dothis = 1
      if (nomedian ne 0) then dothis = 0
      if (dothis) then begin
         for iscan=0L,nscans-1L do begin
            i1     = k_ini[iscan] 
            i2     = k_end[iscan]
            
            jd_per_scan[iscan] = mean( jd[i1:i2] - jd[0], /double ) 
            
            for ichan=0,mfi.nchan-1 do begin
               azs     = a.az[i1:i2] 
               datos   = ctod[ichan,i1:i2]
               flags   = flag[ichan,i1:i2]
               g_bs    = gbmat[ichan,i1:i2]
               if (obsmode eq 5) then begin
                  cut  = where(datos ne 0.0 and flags eq 0 and abs(g_bs) ge GBCUT, ncut) ; for nominal, also excludes galactic plane
               endif else begin
                  cut  = where(datos ne 0.0 and flags eq 0, ncut)
               endelse
               cut2    = where(datos ne 0.0, ncut2)
               
               if ncut ge 6 then begin
                  mediana                   = median( datos[cut] )
                  med_per_scan[ichan,iscan] = mediana
                  datos[cut2]               = datos[cut2] - mediana ; remove median in all non-zero values
               endif else begin
                  datos[*]   = 0.0
               endelse
               
               ctod_nomed[ichan,i1:i2] = datos[*]         
               
            endfor
         endfor

      endif else begin
         ctod_nomed = ctod
      endelse

      
      ; 4.1.2 Removes a dipole/quadrupole for nominal data only
      dothis = 0
      if (dothis and obsmode eq 5) then begin
         NUMSCAN    = 50        ; doing the fit in blocks of NUMSCAN scans

         ctod_model_dipole = ctod*0.0
         ctod_nomed        = ctod*0.0
         factor_cmb        = 0.0
         
         ; First computes the CMB dipole
         cmbdipole = fltarr(4,njd)
         if (factor_cmb ne 0.0) then begin
            for ihorn=0,3 do begin
               get_dipole_prediction, a.gl[ihorn,*], a.gb[ihorn,*], ihdip
               cmbdipole[ihorn,*] = ihdip[*]
            endfor
         endif
         
         for iscan=0L,nscans-1L, NUMSCAN do begin
            iscan1 = iscan
            iscan2 = iscan + NUMSCAN -1L
            if (iscan2 gt (nscans-1L)) then iscan2 = nscans-1L
            i1     = k_ini[iscan1] 
            i2     = k_end[iscan2]
            
            if (verbose) then print,iscan,iscan1,iscan2, i1, i2
            
            for ichan=0,mfi.nchan-1 do begin
               azs   = a.az[i1:i2] 
               datos = ctod[ichan,i1:i2]
               flags = flag[ichan,i1:i2]
               g_bs  = gbmat[ichan,i1:i2]
               cut   = where(datos ne 0.0 and flags eq 0 and abs(g_bs) ge GBCUT, ncut) ; for nominal, also excludes galactic plane
               cut2  = where(datos ne 0.0, ncut2)

               ihorn = mfi.polarizer_id[ichan] - 1
               
               if ncut ge 6 then begin

                  mediana = median( datos[cut] )
                  std     = stddev( datos[cut] )
                  X       = azs[cut]*!DPI/180.0
                  Y       = datos[cut] - factor_cmb*cmbdipole[ihorn,cut]
                  weights = wei[cut] 
                  ;weights = fltarr(ncut) + 1.0
                  coefs   = [mediana, std/4.0, 0.0]
                  yfit    = CURVEFIT(X, Y, weights, coefs, FUNCTION_NAME='fseries', ITMAX=50, status=status)
                  if ichan eq 0 and iscan eq 0 then begin
                     print," --> Removing a dipole (1) or a quadrupole (2) = ",(n_elements(coefs)-1)/2
                     print,"     and using an average over NSCANS          = ",NUMSCAN
                     print,"     The CMB dipole is also removed. "
                  endif


                  if (status eq 2) then begin
                     datos[*]  = 0.0
                     ndip      = i2-i1+1
                     dipolo    = fltarr(ndip)
                  endif else begin
                     fseries, azs*!DPI/180.0, coefs, dipolo
                     datos[cut2] = datos[cut2] - dipolo[cut2] - factor_cmb*cmbdipole[ihorn,cut2]
                  endelse

               endif else begin
                  ndip       = i2-i1+1
                  dipolo     = fltarr(ndip)
                  datos[*]   = 0.0
               endelse
               
               ctod_nomed[ichan,i1:i2]        = datos[*]            
               ctod_model_dipole[ichan,i1:i2] = dipolo
               
            endfor
         endfor
      endif


      
      ; 4.2 Smoothing the median in each sub-scan and interpolating
      dothis = 0
      if (dothis) then begin
         base  = fltarr(mfi.nchan, njd)
         ctod2 = fltarr(mfi.nchan, njd)
         if (nscans gt 1) then begin
            base_per_scan = fltarr(mfi.nchan,nscans)
            for ichan=0,mfi.nchan-1 do begin
               smooth_gain, jd_per_scan, med_per_scan[ichan,*], intarr(nscans), med_smth, med_stdev, tkernel=30.0 ; 30min
               base_per_scan[ichan,*] = med_smth[*]
            endfor
            jdnew = jd - jd[0]
            mfi_interpolate_gains, jd_per_scan, base_per_scan, jdnew, 1, base ; baseline levels
            
            cut2  = where(ctod ne 0.0, ncut2)
            if ncut2 ge 1 then ctod2[cut2] = ctod[cut2] - base[cut2]
         endif
      endif


      ; 4.3 Ring stack
      dothis = 0
      if (dothis) then begin
         masc  = intarr(mfi.nchan, njd)
         cut   = where(ctod ne 0  and flag eq 0 and abs(gbmat) le 20.0, ncut)
         if ncut ge 1 then masc[cut] = 1
         sancho_ring_stack, a.az, ctod_nomed, AZ_v, stack, masc=masc
      endif


      ; 4.4 Sustraer mediana global
      dothis = 0
      if (nomedian ne 0) then dothis = 1
      if (dothis) then begin
	 print," > Removing global median in the CTOD file. "
         ctod3 = ctod * 0.0
         for ichan=0,mfi.nchan-1 do begin
            cut     = where(ctod[ichan,*] ne 0.0 and flag[ichan,*] eq 0 and abs(gbmat[ichan,*]) ge GBCUT, ncut) ; also excludes gal plane
            cut2    = where(ctod[ichan,*] ne 0.0, ncut2)
            if ncut ge 5 then begin
               mediana           = median( ctod[ichan,cut] )
               ctod3[ichan,cut2] = ctod[ichan,cut2] - mediana ; remove median in all non-zero values
            endif 
         endfor
	ctod_nomed = ctod3
      endif


      ; 4.5 Final choice
      ;ctod = ctod2
      ctod = ctod_nomed
      ;ctod = ctod3

      
;      for ichan=0,mfi.nchan-1 do begin
;      
;         mediana       = 0.0
;         cut           = where( flag[ichan,*] eq 0 , ncut)
;         if ncut ge 10 then mediana = median( ctod[ichan,cut] )
;
;         di            = where(ctod[ichan,*] ne 0.0, ndi)
;         if ndi ge 1 then ctod[ichan,di] = ctod[ichan,di] - mediana
;
;      endfor
      



      ;----------------------
      ; [5] Write data
      ; Originally, we proposed this format: b = {jd:a.jd, gl:a.gl, gb:a.gb, phi:phi, data:ctod, wei:wei, flag:flag }
      ; But in practice, we use the same format as btod files. 
      ; jflag, beaflag flags are zero for the moment. They will be updated in sancho_reprocess_ctod.pro
      print," --> and writing ",ff_out
      b = { jd:a.jd, az:a.az, el:a.el, gl:a.gl, gb:a.gb, par:a.par, $
            data:ctod, mod_ang:a.mod_ang, wei:wei, flag:flag, $
            gains:gains_ic_all, mod_ref_angle:theta0, msbin:msbin, gainmodel:params.GAIN_MODEL, $
            fit_method:fit_method, rms_flag:rms_flag, int_rms_flag:int_rms_flag, transients:params.path_byjd_transients,$
	    jflag:0, beaflag:0, rfifaz:0, baseline:baseline, pointingmodel:pointingmodel, azhorn:az_horn, elhorn:el_horn, sigmacov:sigmacov }      
      mwrfits, b, ff_out, /create
         
   endif else begin
      print," --> File already exists: ",ff_out
   endelse
   
   print,""

ENDFOR
         

end

;=====================

PRO sancho_genera_ctod, filename, path=path, rms_flag=rms_flag, verbose=verbose, overwrite=overwrite, nomedian=nomedian, $
                        display=display, int_rms_flag=int_rms_flag, baseline=baseline, rfifaz=rfifaz, outpath=outpath, $
			keeptransients=keeptransients

; Keywords
if not keyword_set(rms_flag) then rms_flag = 0
if not keyword_set(int_rms_flag) then int_rms_flag = 0
if not keyword_set(verbose) then verbose = 0
if not keyword_set(nomedian) then nomedian = 0
if not keyword_set(display) then display = 0
if not keyword_set(baseline) then baseline = 0
if not keyword_set(rfifaz) then rfifaz = 0
if not keyword_set(keeptransients) then keeptransients = 0

; Info
if n_params() eq 0 then begin
   print,""
   print,"   Syntax --  sancho_genera_ctod, txt [,/OVERWRITE, /NOMEDIAN, /VERBOSE, RMS_FLAG=, PATH=, OUTPATH=]"
   print,""
   print,"  txt could be filename or directly the list of files"
   print,""
   print,"  Default options:  IDL> sancho_genera_ctod, txt, rms_flag=6  "
   print,""
   return
endif

; Keywords
if not keyword_set(path) then path = "/net/nas/proyectos/quijote/btod/"


; Identify what is filename
if size(filename,/type) ne 7 then stop," FILENAME should be of type STRING. "
caso = size(filename,/n_dimension)

CASE caso OF
   0: readcol, filename, lista, format="A"
   1: lista = filename
   else: stop," Enter filename or list of files."
ENDCASE

nlist = n_elements(lista)
print," (*) Total number of input ROOT files: ",nlist

; Params
global_ctod_params, params
if not keyword_set(outpath) then outpath=params.outpath
info_quijote, mfi=mfi

; Mask to avoid flags, Nside=128
read_fits_map,"/net/nas/proyectos/quijote/btod_transient/nogal/dontflag.fits",masc,nside=nside,ordering=ordering,coordsys=coordsys


; Display info
print,"  > Path to BTOD files : ",path
print,"  > Output path (CTODS): ",outpath
print,"  > GAIN MODEL         : ",params.gain_model
wait,2

; Run main code
for i=0L,nlist-1L do begin
   prepare_ctod_for_mfi, lista[i], path=path, mfi=mfi, params=params, overwrite=overwrite, rms_flag=rms_flag, verbose=verbose, $
                         nomedian=nomedian, display=display, int_rms_flag=int_rms_flag, noflag=masc, baseline=baseline, rfifaz=rfifaz, outpath=outpath,$
			 keeptransients=keeptransients
endfor


; END message
print," CODE ended succesfully. "

END
