;   4/05/2017  J.A.Rubino
;
;   PURPOSE:
;   Reads existing .ctod files for QUIJOTE MFI and carries out reprocessing. 
;
;
;   HISTORY: 
;   04/05/2017 - original version
;   23/10/2017 - now includes flagging tables based on RMS (/JFLAG) and Beatriz inspection(/beaflag)
;   12/01/2018 - update paths from nas4/quijote to nas/proyectos/quijote
;   12/07/2019 - now including sigmacov at the pixel level. Already incorporated in the BTOD and CTOD files
;   19/01/2021 - updated flagging routines, to check for activated bit before adding. 
;   05/03/2021 - keep baseline keyword.
;
;-
; 

PRO global_ctod_params,par

  outpath    = "/net/nas/proyectos/quijote/ctod/final/"
  ;outpath    = "/net/nas/proyectos/quijote/ctod/w49w51ic443/"
  ;outpath    = "/net/nas/proyectos/quijote/ctod/test_rfictod/"

  ; Gain model: Y, X, SC (sum of corr), SNC (sum of uncorr),
  ;             AS (all referenced to the sum)
  gain_model = "AS2"            ; "AS" ; "SC" ; "SNC" ; "Y" ; "AS"   

  
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


  ; Transients directory by index. Old version
  ;path_transients = "/net/nas/proyectos/quijote/btod_transient/with_signif/"
  path_transients = "/net/nas/proyectos/quijote/btod_transient/nogal/"

  ; Transients directory (by JD). New version
  path_byjd_transients = "/net/nas/proyectos/quijote/btod_transient/byjd_40ms/"

 
  ; structure
  par     = { outpath:outpath, tab_dtheta:tab_dtheta, tab_t0:tab_t0, tab_mtheta:tab_mtheta,$
              tab_diode_cass:tab_diode_cass, path_cal:path_cal, path_gain:path_gain, MINMAX_CORRECTION:MINMAX_CORRECTION, $
              tab_med_amp_diode_y:tab_med_amp_diode_y,         tab_med_gains_cass_y:tab_med_gains_cass_y, $
              tab_med_amp_diode_x:tab_med_amp_diode_x,         tab_med_gains_cass_x:tab_med_gains_cass_x, $
              tab_med_amp_diode_sumnc:tab_med_amp_diode_sumnc, tab_med_gains_cass_sumnc:tab_med_gains_cass_sumnc, $
              tab_med_amp_diode_sumc:tab_med_amp_diode_sumc,   tab_med_gains_cass_sumc:tab_med_gains_cass_sumc, $
              gain_model:gain_model, fit_method:fit_method, path_transients:path_transients, path_byjd_transients:path_byjd_transients  }

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
      gain_stdev[i]  = stdev(gain[cut])
   endif
endfor


end

;==================

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

pro apply_flag_table, filename, jd, flag, bit=bit, path=path
  if not keyword_set(bit) then bit = 2
  if not keyword_set(path) then path = "/net/nas/proyectos/quijote/flag/"

  ; name of the flag table file
  kk1 = strsplit(filename,"/",/extract)
  kk2 = kk1[n_elements(kk1)-1]
  kk3 = strsplit(kk2,"-",/extract)

  if n_elements(kk3) ge 3 then begin
     root = kk3[0] + "-" + kk3[1] + "-" + kk3[2]
     ff_flag = path + root + ".flag"
  endif else begin
     print," --> Error in the name format! "
     return
  endelse

  ; Read flagging table if exists, and if containing lines
  existe = does_exists(ff_flag)
  if (existe) then begin
     nlines = nlines(ff_flag)
  endif else begin
     nlines = 0
     print," No file found: ",ff_flag
  endelse
  if nlines eq 0 then existe = 0
  
  if (existe) then begin
     print," > Reading flagging table ",ff_flag
     readcol,ff_flag,kk,ichan,jdmin,jdmax,format="A,I,D,D"
     if n_elements(kk) ge 1 then existe = 1
  endif 

  ; apply
  if (existe) then begin
     one = 2^(bit-1)
     sancho_check_activebit_inflags, flag, bit, bitmasc
     for i=0,n_elements(ichan)-1 do begin
        cut = where(jd ge jdmin[i] and jd le jdmax[i])
        if (ichan[i] ge  0) then flag[ichan[i],cut] = flag[ichan[i],cut] + one * (1-bitmasc[ichan[i],cut])
        if (ichan[i] eq -1) then flag[*,cut]        = flag[*,cut]        + one * (1-bitmasc[*,cut])
     endfor
  endif else begin
     print,"     No correction applied to flags. "
  endelse
  print,""
  
end


;=====================
; Applies the flags from Diego's code

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
           if n_elements(idx) ge 1 then flag[ichan,idx] = flag[ichan,idx] + one * (1-bitmasc[ichan,idx])
        endif
        
     endfor
  endif
  print,""
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


;=================
;=================
;=================

PRO reprocess_ctod_for_mfi, root, overwrite=overwrite, path=path, mfi=mfi, params=params, rms_flag=rms_flag, $
                            int_rms_flag=int_rms_flag, verbose=verbose, rmdipole=rmdipole, beaflag=beaflag, jflag=jflag, eflag=eflag,$
				noflag=noflag, updatepointing=updatepointing, rfifaz=rfifaz, rfictod=rfictod, outpath=outpath, $
				keepplanets=keepplanets, keeptransients=keeptransients, keepbaseline=keepbaseline

; Keyword
if not keyword_set(overwrite) then overwrite = 0
if not keyword_set(path) then path = ""
if not keyword_set(rms_flag) then rms_flag = 0
if not keyword_set(int_rms_flag) then int_rms_flag = 0
if not keyword_set(verbose) then verbose = 0
if not keyword_set(rmdipole) then rmdipole = 0
if not keyword_set(beaflag) then beaflag = 0
if not keyword_set(jflag) then jflag = 0
if not keyword_set(eflag) then eflag = 0
if not keyword_set(noflag) then noflag = fltarr(nside2npix(128))
if not keyword_set(updatepointing) then updatepointing = 0
if not keyword_set(rfifaz) then rfifaz = 0
if not keyword_set(rfictod) then rfictod = 0
if not keyword_set(keepplanets) then keepplanets = 0
if not keyword_set(keeptransients) then keeptransients = 0
if not keyword_set(keepbaseline) then keepbaseline = 0


; Get global parameters
if not keyword_set(params) then global_ctod_params, params
if not keyword_set(mfi) then info_quijote, mfi=mfi
if not keyword_set(outpath) then outpath = params.outpath
jd_ref    = mfi.jd_ref
jd_140411 = 2456758.5d0
jd_151201 = 2457357.5d0
jd_160301 = 2457448.5d0
jd_160501 = 2457509.5d0
jd_160528 = 2457536.5d0
jd_161015 = 2457676.5d0


; Obtain all the CTOD files associated
print,"==================================================================="
print," (*) MAIN PATH to ctod: ",path
print,"     ROOT name: ",root
spawn,"ls -1 " + path + root + "*.ctod", ff_list, error

if error ne "" then return


; Count TOD maps
ntod = n_elements(ff_list)
print," (*) Total number of CTOD files: ",ntod



; MAIN LOOP over blocks
cfp
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

   
   ; Check if the output file already exists
   existe  = does_exists(ff_out)
   goahead = 1
   if (existe) then goahead = 0
   if (existe and overwrite) then goahead = 1

   ; Main part of code
   if (goahead) then begin

      print," --> Reading ctod: ",filename

      ;-----------------------
      ; [1] Read data
      a     = mrdfits(filename,1,/silent)
      names = tag_names(a)
      ctod  = a.data
      wei   = a.wei
      flag  = a.flag
      jd    = a.jd
      az    = a.az
      el    = a.el
      gl    = a.gl
      gb    = a.gb
      par   = a.par
      njd   = n_elements(jd)
      msbin = a.msbin

      find_mfi_period,jd[0],periodo

      ; Get AZ_horn, EL_horn if defined
      if total(names eq "AZHORN") gt 0 then begin
      	azhorn = a.azhorn
      	elhorn = a.elhorn
      endif else begin
	azhorn = fltarr(4,njd)
	elhorn = fltarr(4,njd)
      endelse       

      ; Get pointing model 
      pointingmodel = "Version1"
      if total(names eq "POINTINGMODEL") gt 0 then pointingmodel = a.pointingmodel

      ; Get sigmacov. Added on 12/07/2019.
      if total(names eq "SIGMACOV") gt 0 then begin
	sigmacov = a.sigmacov
      endif else begin
	print," No SIGMACOV was found. Setting SIGMACOV=0."
	sigmacov = fltarr(4,2,2,njd)
      endelse


      ; Update pointing model if requested
      if (updatepointing and (pointingmodel eq "Version1")) then begin
        print," UPDATING POINTING MODEL from Version1 to Version2."
	pointingmodel = "Version2"
        sancho_mfi_get_pointingmodel, jd, az, el, gl, gb, par, azhorn, elhorn, version= pointingmodel
      endif
      
      ; Update flags from Fred code. bit=9
      if (rms_flag ne 0) then begin
         print,"     Updating flag with Fred's code for polarization and rms_flag = ",rms_flag
         mfi_flag_stats2,filename=root2,f_diff=1,flag=flag,jd=jd,bin=rms_flag
      endif else begin
         print,"     Not using flags from Fred's code for polarization. Keyword rms_flag = ",rms_flag
      endelse

      if (int_rms_flag ne 0) then begin
         print,"     Updating flag with Fred's code for intensity and int_rms_flag = ",int_rms_flag
         mfi_flag_stats,filename=root2,f_data=1,flag=flag,jd=jd,bin=int_rms_flag
      endif else begin
         print,"     Not using flags from Fred's code for intensity. Keyword rms_flag = ",int_rms_flag
      endelse


      ; [1.2] Load flagging tables and additional flags
      ;       bit=8 is defined by the user.
      print," (*) Loading and applying flagging tables (eye inspection of the gain). "
      apply_flag_table, filename, jd, flag, bit=8, path="/net/nas/proyectos/quijote/flag/" 
     
      ;if (msbin gt 30.0) then begin
      ;	apply_flags_transients, filename, jd, flag, bit=8, path=params.path_transients ; bit=8 is defined by the user.
      ;endif
      if (keeptransients eq 0) then begin
      	print," (*) Loading and applying flags from Diegos code (compact objects). "
      	apply_flags_transients_byjd, filename, jd, flag, bit=8, noflag=noflag, path=params.path_byjd_transients ; bit=8 is defined by the user.
      endif else begin
      	print," (*) KEEPTRANSIENTS activated --> Flags from Diegos code (compact objects) are NOT applied. "
      endelse



      ; [1.3] Loading flags from eye inspection (BRG)
      if (beaflag) then begin
         print," (*) Loading and applying flagging tables (eye inspection of the btod, BRG). "
         apply_flag_table, filename, jd, flag, bit=8, path="/net/nas/proyectos/quijote/flag_data/btod/" 
         
         print," (*) Loading and applying flagging tables (eye inspection of the ctod, BRG). "
         apply_flag_table, filename, jd, flag, bit=8, path="/net/nas/proyectos/quijote/flag_data/ctod/test2/" 
      endif
      
      ; [1.4] Loading flags from RMS statistics (JARM)
      if (jflag) then begin
         print," (*) Loading and applying flagging tables (rms of the CTOD, JARM). "
         ;;;sancho_mfi_flag_rfi, root, jd, az, el, gl, gb, ctod, flag, reflag, bit=8, noflag=noflag, sky=2, perelev=0
         ;sancho_mfi_flag_rfi, root, jd, az, el, gl, gb, ctod, flag, reflag, bit=8, noflag=0, sky=2, perelev=0
         sancho_mfi_flag_rfi, root, jd, az, el, gl, gb, ctod, flag, reflag, bit=8, noflag=0, sky=2, perelev=0, $
	 		outpath="/net/nas/proyectos/quijote/flag_stats_ctod/jalberto/withcov2/", withcov=1
         flag = reflag
      endif

      ; [1.5] Loading specific flags for regions (JARM, FG). bit=8 for flags, bit=7 for planets (pbit)
      if (eflag) then begin
	 print," (*) Loading and applying flagging tables (specific for each field, JARM/FG/DT, /eflag, Planets). "
	 apply_flags_per_field, filename, az, el, gl, gb, jd, flag, bit=8, pbit=7, period=periodo, keepplanets=keepplanets
      endif 

      
      ;------------
      ; [3] Interference correction
      ;
      sancho_get_mfi_subscans, root, a.az, a.el, k_ini, k_end, obsmode=obsmode
      check_if_nominal, a.az, is_nominal
      if (is_nominal and obsmode ne 5) then stop,"* Error in obsmode"
      if (obsmode eq 5 and is_nominal eq 0) then stop,"* Error in obsmode"

      file_rfifaz = ""
      if (rfifaz ne 0) then begin
	print," > Doing RFI correction by substracting a function of AZ, f(AZ). "
        ;txt_rfifaz = "azstack_sims_sky_dipole_oofnoise_may2019" ; use this for sims
        ;txt_rfifaz = "azstack_apr2019" ; use this for withcov
        txt_rfifaz = "azstack_feb2021" ; use this for withcov2
        ;txt_rfifaz = "azstack_mar2021" ; use this for withcov3
	sancho_get_rfifazmodel, root, periodo, obsmode, a.mod_ang, az, ctod, newctod, file_rfifaz, xaz=xaz, faz=faz, txtname=txt_rfifaz 
	ctod = newctod
      endif else begin
	if rfictod then begin
		print," > Doing RFI correction by substracting a function of AZ, f(AZ) PER INDIVIDUAL CTOD file. "
		sancho_calc_rfifazmodel_4onectod, root, periodo, obsmode, a.mod_ang, az, ctod, flag, newctod, xaz=xaz, faz=faz
		ctod = newctod
	endif else begin
      		print," > No interference correction at this stage. "
      	endelse
      endelse


      ;------------
      ; [4] Baseline correction. 

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
      ctod_orig  = ctod          ; original ctod, no baseline substraction

      ; 4.1.1 Removes the median. Done for all database.
      nomedian = 0
      dothis   = 1
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
;            if (obsmode eq 5) then begin
;               cut  = where(datos ne 0.0 and flags eq 0 and abs(g_bs) ge GBCUT, ncut) ; for nominal, also excludes galactic plane
;            endif else begin
;            cut  = where(datos ne 0.0 and flags eq 0, ncut)
;            endelse
            cut  = where(datos ne 0.0 and flags eq 0, ncut)
            cut2 = where(datos ne 0.0, ncut2)
            
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
      dothis = rmdipole
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
;               cut   = where(datos ne 0.0 and flags eq 0 and abs(g_bs) ge GBCUT, ncut) ; for nominal, also excludes galactic plane
               cut   = where(datos ne 0.0 and flags eq 0 , ncut) ; DO not exclude galactic plane
               cut2  = where(datos ne 0.0, ncut2)

               ihorn = mfi.polarizer_id[ichan] - 1
               
               if ncut ge 6 then begin

                  mediana = median( datos[cut] )
                  std     = stdev( datos[cut] )
                  X       = azs[cut]*!DPI/180.0
                  Y       = datos[cut] - factor_cmb*cmbdipole[ihorn,cut]
                  weights = fltarr(ncut) + 1.0  ; do not use wei[cut]               
                  ;coefs   = [mediana, 0.0, 0.0, 0.0, 0.0]
                  coefs   = [mediana, 0.0, 0.0]
                  status  = 0
                  yfit    = CURVEFIT(X, Y, weights, coefs, FUNCTION_NAME='fseries', ITMAX=50, status=status)
                  if ichan eq 0 and iscan eq 0 then begin
                     print," --> Removing a dipole (1) or a quadrupole (2) = ",(n_elements(coefs)-1)/2
                     print,"     and using an average over NSCANS          = ",NUMSCAN
                     if (factor_cmb ne 0.0) then print,"     The CMB dipole is also removed. "
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


      ; spin-synchronous correction based on an atmospheric model
      dothis = 0
      if dothis && (obsmode eq 5) then begin
         ss_method = 'fit'      ;'ilc'
         case ss_method of
            'fit': begin
               ss_time_range = 1 ; minutes
               nmode = 2         ; 1:dipole, 2:dipole+quadrupole, ...
            end
            'ilc': begin
               ss_time_range = 5 ; minutes
               nmode = 0
            end
         endcase
         print," (*) Spin-synchronous correction base on atmospheric model: ",ss_method
         ctod_nomed = suppress_spinsynch(a.jd, ctod, flag, a.el, a.az, method=ss_method, nmode=nmode, time_range=ss_time_range)
      endif


      ; 4.5 Final choice
      ;ctod = ctod2
      ctod = ctod_nomed
      ;ctod = ctod3
      if keepbaseline then begin
	ctod=ctod_orig
	print," NO BASELINE substraction applied. "
      endif
      

      ;----------------------
      ; [5] Write data
      print," --> and writing ",ff_out
      baseline = 0
      names = tag_names(a)
      cut1  = where(names eq "BASELINE", n1)
      if n1 eq 1 then baseline = a.baseline
      b = { jd:a.jd, az:a.az, el:a.el, gl:gl, gb:gb, par:par, $
            data:ctod, mod_ang:a.mod_ang, wei:wei, flag:flag, $
            gains:a.gains, mod_ref_angle:a.mod_ref_angle, msbin:a.msbin, gainmodel:a.gainmodel, $
            fit_method:a.fit_method, rms_flag:a.rms_flag, int_rms_flag:int_rms_flag, transients:params.path_byjd_transients,$
            jflag:jflag, beaflag:beaflag, eflag:eflag, baseline:baseline, rfifaz:rfifaz, pointingmodel:pointingmodel, azhorn:azhorn, elhorn:elhorn, sigmacov:sigmacov, rfictod:rfictod, file_rfifaz:file_rfifaz   }
      mwrfits, b, ff_out, /create
      
   endif else begin
      print," --> File already exists: ",ff_out
   endelse
   
   print,""

ENDFOR
         

end

;=====================

PRO sancho_reprocess_ctod, filename, path=path, rms_flag=rms_flag, verbose=verbose, overwrite=overwrite, int_rms_flag=int_rms_flag,$
                           rmdipole=rmdipole, beaflag=beaflag, jflag=jflag, eflag=eflag, updatepointing=updatepointing, rfifaz=rfifaz, rfictod=rfictod, outpath=outpath, $
				keepplanets=keepplanets, keeptransients=keeptransients,keepbaseline=keepbaseline

; Keywords
if not keyword_set(rms_flag) then rms_flag = 0
if not keyword_set(int_rms_flag) then int_rms_flag = 0
if not keyword_set(verbose) then verbose = 0
if not keyword_set(rmdipole) then rmdipole = 0
if not keyword_set(beaflag) then beaflag = 0
if not keyword_set(jflag) then jflag = 0
if not keyword_set(eflag) then eflag = 0
if not keyword_set(updatepointing) then updatepointing = 0
if not keyword_set(rfifaz) then rfifaz = 0
if not keyword_set(rfictod) then rfictod = 0
if not keyword_set(keepplanets) then keepplanets = 0
if not keyword_set(keeptransients) then keeptransients = 0
if not keyword_set(keepbaseline) then keepbaseline = 0



; Info
if n_params() eq 0 then begin
   print,""
   print,"   Syntax --  sancho_reprocess_ctod, txt [/overwrite, /JFLAG, /BEAFLAG, /EFLAG, rms_flag=, int_rms_flag=]"
   print,""
   print,"  txt could be filename or directly the list of files"
   print,""
   return
endif

; Keywords with paths
if not keyword_set(path) then path = "/net/nas/proyectos/quijote/ctod/"


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
info_quijote, mfi=mfi

; Output path
if not keyword_set(outpath) then outpath = params.outpath

; Mask to avoid flags, Nside=128
read_fits_map,"/net/nas/proyectos/quijote/btod_transient/nogal/dontflag.fits",masc,nside=nside,ordering=ordering,coordsys=coordsys


; Display info
print,"  > Input path : ",path
print,"  > Output path: ",outpath
print,"  > GAIN MODEL : ",params.gain_model
wait,2

; Run main code
for i=0L,nlist-1L do begin
   reprocess_ctod_for_mfi, lista[i], path=path, mfi=mfi, params=params, overwrite=overwrite, $
                           rms_flag=rms_flag, int_rms_flag=int_rms_flag, verbose=verbose, rmdipole=rmdipole, $
                           beaflag=beaflag, jflag=jflag, eflag=eflag,  noflag=masc, updatepointing=updatepointing, $
			   rfifaz=rfifaz, rfictod=rfictod, outpath=outpath, keepplanets=keepplanets, keeptransients=keeptransients, $
			   keepbaseline=keepbaseline
endfor


; END message
print," CODE ended succesfully. "

END
