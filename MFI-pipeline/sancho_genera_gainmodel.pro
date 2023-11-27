;   21/02/2017  J.A.Rubino
;
;   PURPOSE:
;   Generates .gain files for QUIJOTE MFI
;
;
;   HISTORY: 
;   21/02/2017 - original version. based on the previous version of sancho_genera_ctod.pro
;   19/05/2017 - health check. If we have data points with
;                             (1+delta)=0, then the full correction is set to zero.
;   12/01/2018 - change paths from /net/nas4/quijote/ to /net/nas/proyectos/quijote/
;
;
;-
; 

PRO global_gain_params,par

  outpath = "/net/nas/proyectos/quijote/gain/" 

  ; Kernel for smoothing gains. Four values, one per Horn
  tk_lowfreq  =  30.0 ; min
  tk_highfreq = 120.0 ; min
  tkernel_min = [tk_lowfreq, tk_highfreq, tk_lowfreq, tk_highfreq ] ;min
 
  ; Mininum gain allowed. Below this value, gain is assumed to be zero
  MIN_GAIN_PER_HORN  = [ 8d-4, 1d-4, 1d-4, 1d-4 ] ; volts 
  
  ; Maximum Chi2 allowed
  MAXCHI2 = 3.0
  
  ; Maximum-minimum correction allowed
  MINMAX_CORRECTION = [0.8,1.2] ;[0.5,1.5]

  ; Minimum duration of observation for doing linear fit to gain
  MIN_TOBS = 2.0                ; hr
  

   ; Gain model: Y, X, SC (sum of corr), SNC (sum of uncorr),
  ;             AS (all referenced to the sum)
  gain_model = "AS2"  
    
  ; Fit method. Default calibration method used by RGS.
  fit_method = 6
  
  ; CAL files
  path_cal = "/net/nas/proyectos/quijote/cal/"

  ; structure
  par     = { outpath:outpath, tkernel_min:tkernel_min, maxchi2:MAXCHI2, $
              path_cal:path_cal, gain_model:gain_model, fit_method:fit_method, $
              MINMAX_CORRECTION:MINMAX_CORRECTION, MIN_TOBS:MIN_TOBS, MIN_GAIN_PER_HORN:MIN_GAIN_PER_HORN }
  
END




;=====================

pro smooth_gain_and_reflag, jd, gain, flag, gain_smooth2, gain_stdev2, tkernel=tkernel, reflag=reflag, ihorn=ihorn

if not keyword_set(ihorn) then ihorn=0

global_gain_params, par
if not keyword_set(tkernel) then begin
   tkernel = par.tkernel_min[ihorn] ; minutes
endif
   
dt_min = tkernel                ; minutes
dt     = dt_min/60.0/24.0       ; jd  


; Min gain
MIN_GAIN = par.MIN_GAIN_PER_HORN[ihorn]

; initialize to zero
n  = n_elements(gain)
gain_smooth = fltarr(n)
gain_stdev  = fltarr(n)


; first iteration
for i=0L,n-1L do begin
   cut = where(gain gt 0 and flag eq 0 and abs(jd-jd[i]) le dt/2.0, ncut)
   if ncut ge 10 then begin
      gain_smooth[i] = median(gain[cut])
      gain_stdev[i]  = stdev(gain[cut])
   endif
endfor


; reflag using a N-sigma clipping around the smoothed gain.
sigma  = median(gain_stdev)
clipp  = 6.0
reflag = flag
cut    = where( abs(gain-gain_smooth) ge clipp*sigma or gain le MIN_GAIN, n_reflagged)
if n_reflagged ge 1 then reflag[cut] = 1
   
; second iteration
gain_smooth2 = gain_smooth
gain_stdev2  = gain_stdev
if n_reflagged ge 1 then begin
   for i=0L,n-1L do begin
      cut = where( reflag eq 0 and abs(jd-jd[i]) le dt/2.0, ncut)
      if ncut ge 10 then begin
         gain_smooth2[i] = median(gain[cut])
         gain_stdev2[i]  = stdev(gain[cut])
      endif
   endfor
endif



; Smooth and reflag, main part
use_new_implementation = 1

IF (use_new_implementation) THEN BEGIN
   ; first iteration
   sancho_smooth_gain, jd, gain, flag, dt, gain_smooth, gain_stdev_ignore
   ; reflag using a N-sigma clipping around the smoothed gain.
   sigma  = median(gain_stdev)
   clipp  = 6.0
   reflag = flag
   cut    = where( abs(gain-gain_smooth) ge clipp*sigma or gain le MIN_GAIN, n_reflagged)
   if n_reflagged ge 1 then reflag[cut] = 1
   ; second iteration
   sancho_smooth_gain, jd, gain, reflag, dt, gain_smooth2, gain_stdev_ignore2
ENDIF


; Filling gaps
gain_smooth3 = sancho_lininterpol_gain(jd-jd[0], gain_smooth2, reflag)
gain_smooth2 = gain_smooth3


end


;=================

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

;==========================
; New Code, based on mfi_params.pro
pro select_values_for_gain_correction, caso, mfi, params, imodpos, ihorn, ifreq, itype, gains_ic, amp_diode, $
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
;=================
;=================

PRO prepare_gain_for_mfi, root, overwrite=overwrite, path=path, mfi=mfi, params=params, verbose=verbose, $
                          display=display, interactive=interactive


; Keyword
if not keyword_set(overwrite) then overwrite = 0
if not keyword_set(path) then path = ""
if not keyword_set(verbose) then verbose = 0
if not keyword_set(display) then display = 0
if not keyword_set(interactive) then interactive = 0

; Get global parameters
if not keyword_set(params) then global_gain_params, params
if not keyword_set(mfi) then info_quijote, mfi=mfi

; Params
tkernel_v = params.tkernel_min
jd_140411 = 2456758.5d0
jd_151201 = 2457357.5d0
jd_160301 = 2457448.5d0
jd_160501 = 2457509.5d0
jd_160528 = 2457536.5d0
jd_161015 = 2457676.5d0



; CAL file
print," (*) MAIN PATH to cal: ",path
print,"     ROOT name: ",root
print,"     Smooth kernel (min): ",tkernel_v
nlen   = strlen(root)
root2  = ""
if nlen ge 5 then root2  = strmid(root,0,nlen-4) ; removes -001 if present

ffcal  = path + root  + ".cal" 
ffcal2 = path + root2 + ".cal" 
if does_exists(ffcal2) then ffcal = ffcal2

; Output file
outpath = params.outpath
ff_out   = outpath + root + ".gain"


; Check nominal mod position
kk3      = strsplit(root,"-",/extract)
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
   
   print," --> Reading CAL: ",ffcal

   existe_cal = does_exists(ffcal)
   if existe_cal eq 0 then begin
	print,"**** CAL file does no exists ***"
	wait,3
	return
   endif

   a     = mrdfits(ffcal,1)
   jd    = a.jd
   flag  = a.flag
   t     = (jd - jd[0])*24.0 ; hrs
   
   n     = n_elements(jd)
   nchan = mfi.nchan
   gains = fltarr(nchan,n)
   gains_std = fltarr(nchan,n)
   chi2  = fltarr(nchan)

   TRANGE_HRS = (max(jd) - min(jd))*24.0 ; Trange in hrs
   
   ; Smooth and reflag
   for ichan=0,nchan-1 do begin

      ihorn   = mfi.polarizer_id[ichan]-1
      tkernel = tkernel_v[ihorn]

      smooth_gain_and_reflag, jd, a.gain[ichan,*], flag[ichan,*], gains_v, gains_stdev_v, tkernel=tkernel, reflag=reflag,ihorn=ihorn
      flag[ichan,*] = reflag[*]
      
      ; chi2 estimate
      cut = where(flag[ichan,*] eq 0, ncut)
      chi2[ichan] = total( (gains_v[cut] - a.gain[ichan,cut])^2. / gains_stdev_v[cut]^2. ) / float(ncut)
      
      ; health check with the chi2
      factor = 1.0
      if (chi2[ichan] gt params.maxchi2) then factor = 0.0
     
      ; final gain
      gains[ichan,*] = gains_v[*]*factor
      gains_std[ichan,*] = gains_stdev_v[*]
      
   endfor
   gains_orig = a.gain

   
   ; Additional smoothing in real space (2) or low-pass filter (1), or no change (0). --> gains2
   method = 1                   ; default for long observations
   if (TRANGE_HRS le params.MIN_TOBS) then method = 0
   CASE method OF
      0: gains2 = gains
      1: low_pass_gain_flag, jd, gains, gains2, flag, min_period=tkernel_v/60.0
      2: gaussian_smoothing_gain, jd, gains, gains2
   ENDCASE

   ; Display smoothed gains
   if (display) then begin
      cfp
      wset,0
      !p.multi=[0,1,4]
      !p.charsize=1.7
      for ihorn=0,3 do begin
         txth = ", horn "+string(ihorn+1,format="(i1)") 
         i1   = ihorn*8
         i2   = i1+7
         yran = minmax(gains[i1:i2,*])*1d3*1.05
         xran = minmax(t)
         for i=i1,i2 do begin
            cut = where(flag[i,*] eq 0, ncut)
            if ncut eq 0 then cut = indgen(n_elements(t))
            if (i eq i1) then plot,t[cut],gains_orig[i,cut]*1d3,yr=yran,xr=xran, /xsty, title=root+txth,xtitle="time (h)",ytitle="Gain (mV)"
            if (i ne i1) then oplot,t[cut],gains_orig[i,cut]*1d3
            oplot,t,gains[i,*]*1d3,col=2 ; red
            oplot,t,gains2[i,*]*1d3,col=3,linestyle=2 ; green
            ;oplot,t,(gains2[i,*]+gains_std[i,*])*1d3,col=4,linestyle=2
         endfor
      endfor
   endif


   
   ; Building gain model
   fit_method = params.fit_method
   mfi_params,JD_OBS=jd[0],CHAN_REF=1,FIT_METHOD=fit_method,GAINS_IC=gains_ic,$
              GAINS_MOD_POS_IC=gic_permod,AMP_DIODE=amp_diode,AMP_MOD_POS_DIODE=amp_diode_cal
   njd        = n_elements(jd)
   oneplusdeltaG = fltarr(mfi.nchan, njd)
   corr_orig     = fltarr(mfi.nchan, njd)
   Phase_diode   = fltarr(mfi.nchan, njd)
   for ichan=0,mfi.nchan-1 do begin

      ; Get ihorn, ifreq and itype
      ihorn = mfi.polarizer_id[ichan] - 1
      freq  = mfi.central_freq_id[ichan]
      ifreq = 0                                        ; low-freq
      if (freq eq 13.0 or freq eq 19.0) then ifreq = 1 ; high-freq
      canales = ["s","d","x","y"]
      cut     = where(canales eq mfi.chan_id[ichan])
      itype   = cut[0]
      
         
      ; Check Mod position
      identify_mfi_mod_position, a.mod_ang[ihorn,*], imodpos
      
      if (ihorn ge 1 and imodnomi ne imodpos) then print," FILE name and imodpos do not agree. Using imodpos from mod_angle. "
      if (ihorn eq 0 and jd[0] ge jd_140411 and imodpos ne 0) then begin
         print," -- WARNING -- Error in horn1"
         imodpos = 0
      endif

         
      ; Select reference channels for gain correction 
      caso        = params.gain_model
      rfactor     = 0.0
      ref_gain    = 0.0
      iref_gain_d = 0
      select_values_for_gain_correction, caso, mfi, params, imodpos, ihorn, ifreq, itype, gains_ic, amp_diode, $
                                         amp_diode_cal, iref_gain, ref_gain, G_cass, iref_gain_d, rfactor

      ; Computes (1 + delta_G) for the smoothed gains:
      func = 0.0
      if (ref_gain ne 0.0) then func = ( gains2[iref_gain,*] + rfactor * gains2[iref_gain_d,*] ) / ref_gain
      minG = min(func)
      maxG = max(func)
      medG = median(func)
      if (minG lt params.MINMAX_CORRECTION[0] or maxG gt params.MINMAX_CORRECTION[1] ) then func = 0.0
      oneplusdeltaG[ichan,*] = func

      ; and also computes (1 + delta_G) at the original resolution
      func = 0.0
      if (ref_gain ne 0.0) then func = ( gains_orig[iref_gain,*] + rfactor * gains_orig[iref_gain_d,*] ) / ref_gain
      if (minG lt params.MINMAX_CORRECTION[0] or maxG gt params.MINMAX_CORRECTION[1] ) then func = 0.0
      corr_orig[ichan,*] = func

      
      ; print values on screen
      printf,-1, ichan, ihorn, ifreq, imodpos, iref_gain, ref_gain, median(gains2[iref_gain,*]), minG, maxG, medG, $
             G_cass, mfi.ic[ihorn,ifreq,itype], caso, format='(5i6,2f13.6,4f13.5,i6,a8)'


      ; Computes also the phase of the diode
      denom = ( gains2[iref_gain,*] + rfactor * gains2[iref_gain_d,*] )
      cutd  = where(denom ne 0.0, ncutd)         
      if ncutd ge 1 then Phase_diode[ichan,cutd] = ( gains2[iref_gain,cutd] - rfactor * gains2[iref_gain_d,cutd] ) / denom[cutd]
                            
      
   endfor


   ; Display gain model correction
   if (display) then begin
      window,1
      !p.multi=[0,1,4]
      !p.charsize=1.7
      for ihorn=0,3 do begin
         txth = ", horn "+string(ihorn+1,format="(i1)") 
         i1   = ihorn*8
         i2   = i1+3
         yran = [0.8,1.2] ;params.MINMAX_CORRECTION
         xran = minmax(t)
         for i=i1,i2 do begin
            cut = where(flag[i,*] eq 0, ncut)
            if ncut eq 0 then cut = indgen(n_elements(t))
            if (i eq i1) then plot,t[cut],oneplusdeltaG[i,cut],yr=yran,xr=xran, /xsty, title=root+txth,xtitle="time (h)",ytitle=" 1 + delta_G"
            if (i ne i1) then oplot,t[cut],oneplusdeltaG[i,cut],col=i-i1+1
            oplot,t[cut],corr_orig[i,cut],psym=3
         endfor
      endfor
   endif


   ; Display phase of the diode
   if (display) then begin
      window,2
      !p.multi=[0,1,4]
      !p.charsize=1.7
      for ihorn=0,3 do begin
         txth = ", horn "+string(ihorn+1,format="(i1)") 
         i1   = ihorn*8
         i2   = i1+3
         yran = [-1.0,1.0]      ;params.MINMAX_CORRECTION
         xran = minmax(t)
         for i=i1,i2 do begin
            cut = where(flag[i,*] eq 0, ncut)
            if ncut eq 0 then cut = indgen(n_elements(t))
            if (i eq i1) then plot,t[cut],phase_diode[i,cut],yr=yran,xr=xran, /xsty, title=root+txth, $
                                   xtitle="time (h)",ytitle="Diode pol (Q or U)"
            if (i ne i1) then oplot,t[cut],phase_diode[i,cut],col=i-i1+1
         endfor
      endfor
   endif

      
   ; Write file   
   print," --> and writing ",ff_out
   b = {jd:jd, gainmodel:gains2, flag:flag, mod_ang:a.mod_ang, chi2:chi2, tkernel_min:tkernel_v, oneplusdeltaG:oneplusdeltaG, $
        corr_orig:corr_orig, model:params.gain_model, fit_method:params.fit_method, MIN_TOBS_HOURS:params.MIN_TOBS,$
        phase_diode:Phase_diode }
   spawn,"rm -f "+ff_out
   mwrfits, b, ff_out

endif else begin
   print," --> File already exists: ",ff_out
endelse

print,""

; Interactive
if (interactive) then begin
   gofasser = ""
   read,gofasser,prompt=" Should I go fasser (y/n) [y]? "
   if gofasser eq "" then gofasser = "y"
   if gofasser ne "y" then stop, " STOP -- inside prepare_gain_for_mfi() subroutine"
endif

end

;=====================

PRO sancho_genera_gainmodel, filename, path=path, verbose=verbose, overwrite=overwrite, display=display, $
                             dops=dops, interactive=interactive

; Keywords
if not keyword_set(verbose) then verbose = 0
if not keyword_set(display) then display = 0
if not keyword_set(dops) then dops = 0
if not keyword_set(interactive) then interactive = 0

; Info
if n_params() eq 0 then begin
   print,""
   print,"   Syntax --  sancho_genera_gainmodel, txt [,/OVERWRITE, /DISPLAY]"
   print,""
   print,"  txt could be filename or directly the list of files"
   print,""
   return
endif


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
global_gain_params, params
info_quijote, mfi=mfi

; Keywords
if not keyword_set(path) then path = params.path_cal


; Display info
print,"  > Output path: ",params.outpath

; Output PS or screen
cfp
if (dops) then begin
   abre_ps,"gainmodel.ps",/landscape,/color
   display = 1
endif else begin
   if (display) then window,0,xs=1400,ys=900
endelse

; Run main code
for i=0L,nlist-1L do begin
   prepare_gain_for_mfi, lista[i], path=path, mfi=mfi, params=params, overwrite=overwrite, verbose=verbose,$
                         display=display, interactive=interactive
endfor

if dops then cierra_ps

END
