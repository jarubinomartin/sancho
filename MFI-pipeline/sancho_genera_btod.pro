;   22/07/2014  J.A.Rubino
;
;   PURPOSE:
;   Generates .btod files for QUIJOTE MFI
;
;   HISTORY: 
;   22/07/2014 - original version
;   31/07/2014 - BTOD are now grouped in blocks of 20 TOD files.
;   19/10/2014 - flag_bad_data also removes now the negative voltage values (Volt_min)
;   03/03/2015 - flag_hk now also uses the TEMP_BEM values. Also a REFLAG mode
;                implemented. Data are not rebined in this mode, just the flags are changed.
;   07/05/2015 - if the data corresponds to NOMINAL mode with vel(AZ)>10 deg/s,
;                then uses 40ms binning. See routine update_msbin_from_hk().
;   20/10/2016 - updated upper limit for pressure. Originally it was 2d-6. On March 21st 2016,
;                we set 4d-6. Now we put 1d-5.
;   12/01/2018 - paths changed from /net/nas4/quijote/ to /net/nas/proyectos/quijote/
;   22/03/2018 - New pointing model included in the computation. Also AZ_horn and EL_horn are added
;                at the end of the new .btod file.
;   07/09/2018 - Some temperature sensors are broken in 2018. We change HK limits using "update_hk_temp_limits".
;   23/01/2019 - all nominal observations now are binned to 40ms, including those with 6deg/s.
;   07/07/2019 - using new sancho_bin_data_mfi.pro, which also returns the sigmacov. Output now includes it.
;   16/11/2020 - flags 217x and 217y channels for any observation in the period 1/5/2016 - 28/5/2016, during 
;                which channel 217x underwent a strong gain variation. [RGS]
;   19/01/2021 - flag routines revisited. Now, they check if the corresponding bit is activated before adding it. In this
;                way, we can always trace back the status of flags. Also there was a typo for flags in data.
;-
; 

PRO global_btod_params,par

  msbin      = 40.0;60.0             ; ms. Bin size
  msbin_nomi = 40.0             ; ms. Bin size used for nominal mode

  outpath  = "/net/nas/proyectos/quijote/btod/"

  volt_max = 3.0 ; Volts. Maximum voltage allowed 
  volt_min = 0.0 ; Volts. Minimum voltage allowed 

  ntod_per_btod = 20 ; Number of TOD files in a BTOD.


  ; radii around sources for flagging
  radius_sat  =  5.0            ; deg
  radius_nsl  =  5.0            ; deg. nsl=near side lobe
  radius_sun  = 10.0            ; deg
  radius_moon = 10.0            ; deg


  ; files for HK
  path_hk           = "/net/nas/proyectos/quijote/tod/"
  MAX_PRESS         = 1d-4      ; bar. Updated. Before we used 2d-6 and 4d-6 from March 21st 2016.
  MAX_20K_SHIELD    = 30.0      ; K
  MAX_TL_20K_BASE   = 18.0      ; K
  MAX_CCC_1ST       = 40.0      ; K
  MAX_DIFF_TEMP_BEM =  2.5      ; Celsius. Chosen by eye by JARM.



  ; files for gain flagging
  filegain = "/net/nas/proyectos/quijote/etc/median_stdev_G_stats_listALL_NOM_NOSHIELD_avcal.txt"
  path_cal = "/net/nas/proyectos/quijote/cal/"


  ; Pointing model
  POINTINGMODEL = "Version2" ; "Version1" from 2013. "Version2" from March2018

  ; structure
  par     = { msbin:msbin, outpath:outpath, volt_max:volt_max, volt_min:volt_min, $
              ntod_per_btod:ntod_per_btod, nominal_msbin:msbin_nomi, $
              radius_sat:radius_sat,radius_nsl:radius_nsl, radius_sun:radius_sun, $
              radius_moon:radius_moon, path_hk:path_hk, filegain:filegain, path_cal:path_cal, $
              MAX_PRESS:MAX_PRESS, MAX_20K_SHIELD:MAX_20K_SHIELD, $
              MAX_TL_20K_BASE:MAX_TL_20K_BASE, MAX_CCC_1ST:MAX_CCC_1ST, $
              MAX_DIFF_TEMP_BEM:MAX_DIFF_TEMP_BEM, POINTINGMODEL:POINTINGMODEL }

END


;===================

pro update_hk_temp_limits, jd0,  MAX_20K_SHIELD, MAX_TL_20K_BASE, MAX_CCC_1ST 

  global_btod_params, par
  MAX_20K_SHIELD    = par.MAX_20K_SHIELD
  MAX_TL_20K_BASE   = par.MAX_TL_20K_BASE
  MAX_CCC_1ST       = par.MAX_CCC_1ST

  ; after March/2018, we change temperature settings
  jd_180301         = 2458178.5d
  if (jd0 ge jd_180301) then begin
     MAX_20K_SHIELD    = 1.d9      ; K. Sensor broken. Effectively not using this bound
     MAX_TL_20K_BASE   = 18.0      ; K. Same value as before
     MAX_CCC_1ST       = 45.0      ; K. Increased to 45K.
     print,"   ---> Updated values for temperature settings (after march/2018)."
  endif
end

;===================

pro update_msbin_from_hk, root_hk, msbin

  global_btod_params, par
  msbin = par.msbin             ; default value

  ffhk  = par.path_hk + root_hk + ".hk" 
  if does_exists(ffhk) then begin
     a     = mrdfits(ffhk,1)
     ;if (a.obs_params[0] eq 5 AND a.obs_params[1] ge 10) then begin
     if (a.obs_params[0] eq 5 AND a.obs_params[1] ge 5) then begin
        msbin = par.nominal_msbin ; New bin size.
        print," "
        printf,-1," WARNING. ",root_hk," corresponds to a NOMINAL mode observation", format='(a,a,a)'
        printf,-1,"          with vel(AZ)=",a.obs_params[1]," deg/s. Using msbin (ms) = ",$
               msbin,format='(a,f6.1,a,f5.1)'
        print," "
     endif
  endif
  
end


;===================

pro flag_mfi_hk, root_hk, jd, flag, bit=bit
  if not keyword_set(bit) then bit = 1
  
  global_btod_params,par
  info_quijote,mfi=mfi
  
  ; variables
  nchan = 32
  one   = 2^(bit-1)

  ; obtain HK file
  ffhk  = par.path_hk + root_hk + ".hk" 
  print,"  --> HK file: ",ffhk


  ; compute flags only if HK file exists
  if does_exists(ffhk) then begin

     a     = mrdfits(ffhk,1)
     jd_hk = a.jd + mfi.jd_ref 
     cut   = where(jd_hk ge min(jd) and jd_hk le max(jd), ncut)

     ; TEMPERATURE LIMITS
     update_hk_temp_limits, jd_hk[0],  MAX_20K_SHIELD, MAX_TL_20K_BASE, MAX_CCC_1ST 

     ; Apply here the different criteria
     doflag = 0
     if (ncut ge 1) then begin
        ; Internal Temperatures. Any temperature in the tod.
        ; See example of NOMINAL65D-131230-1034.hk with spikes due to pb in one compresor
        if ( max(a.temp[11,cut]) gt MAX_20K_SHIELD  ) then doflag = 1
        if ( max(a.temp[13,cut]) gt MAX_TL_20K_BASE ) then doflag = 1
        if ( max(a.temp[14,cut]) gt MAX_CCC_1ST     ) then doflag = 1
        ; BEM rack temperatures. Difference between maximum and mininum
        diff_temp_bem = max(a.TEMP_BEM[0,cut]) - min(a.TEMP_BEM[0,cut]) 
        if ( diff_temp_bem gt par.MAX_DIFF_TEMP_BEM ) then doflag = 1
        ; PRESSURE. Note that .hk file contains press in mbar.
        if ( max(a.press*1d-3) gt par.MAX_PRESS ) then doflag = 1
     endif

     ; Flag only where the HK flag is not active.
     if (doflag) then begin
        sancho_check_activebit_inflags, flag, bit, bitmasc
        flag = flag + one * (1-bitmasc) 
     endif

  endif else begin
     print,"     HK file not found. " 
  endelse

end


PRO flag_bad_data, data, flag, bit=bit
  if not keyword_set(bit) then bit = 2
  global_btod_params,par
  one = 2^(bit-1)

  ; Remove high and low voltage data
  VOLT_MAX = par.volt_max
  VOLT_MIN = par.volt_min
  sancho_check_activebit_inflags, flag, bit, bitmasc
  cut = where( (data ge VOLT_MAX OR data le VOLT_MIN) AND bitmasc eq 0, ncut)
  if ncut ge 1 then flag[cut] = flag[cut] + one

END

;---

PRO flag_sun_moon, jd, gl, gb, flag

  d2r         = !DPI/180.d0
  global_btod_params,par
  radius_sun  = par.radius_sun
  radius_moon = par.radius_moon
  
  flag        = intarr(n_elements(jd))

  glactc, ra,dec, 2000.0, gl,gb, 2, /degree

  ; SUN
  sunpos, jd, ras, decs
  d_sun = distancia(ra,dec,ras,decs,/degree)
  cut   = where(d_sun le radius_sun, ncut)
  if (ncut ge 1) then flag[cut] = 1

  ; Moon 
  moonpos, jd, ram, decm
  d_moon = distancia(ra,dec,ram,decm,/degree)
  cut    = where(d_moon le radius_moon, ncut)
  if (ncut ge 1) then flag[cut] = 1

END

;---

PRO flag_mfi_satellites, az, el, flag, bit=bit
  if not keyword_set(bit) then bit = 4

  global_btod_params,par
  radius_sat = par.radius_sat
  one        = 2^(bit-1)

  sancho_check_activebit_inflags, flag, bit, bitmasc
  
  for ihorn=1,4 do begin
     ichan1 = 8*(ihorn-1) 
     ichan2 = ichan1+7

     txt      = string(ihorn,format="(i1)")
     filename = paths_for_quijote()+"satellites/final_master_catalog_Horn"+txt+".dat" 
     readcol,filename,x,y,format="F,F" ; AZ, EL
     ns = n_elements(x)
     for j=0,ns-1 do begin
        dist = distancia(az,el,x[j],y[j],/degree)
        cut  = where(dist le radius_sat, ncut) 
        if (ncut ge 1) then flag[ichan1:ichan2,cut] = flag[ichan1:ichan2,cut] + one * (1-bitmasc[ichan1:ichan2,cut])
     endfor
  endfor

END


PRO flag_mfi_sat_nsl, az, el, flag, bit=bit
  if not keyword_set(bit) then bit = 4

  global_btod_params,par
  radius_nsl = par.radius_nsl
  one        = 2^(bit-1)

  sancho_check_activebit_inflags, flag, bit, bitmasc
  
  for ihorn=1,4 do begin
     ichan1 = 8*(ihorn-1) 
     ichan2 = ichan1+7

     txt      = string(ihorn,format="(i1)")
     filename = paths_for_quijote()+"satellites/visual_catalog_nearlobes_Horn"+txt+".dat" 
     if does_exists(filename) then begin
        readcol,filename,x,y,ampl,format="F,F,F" ; AZ, EL, Amplitude
        ns = n_elements(x)
        for j=0,ns-1 do begin
           dist = distancia(az,el,x[j],y[j],/degree)
           cut  = where(dist le radius_nsl, ncut) 
           if (ncut ge 1) then flag[ichan1:ichan2,cut] = flag[ichan1:ichan2,cut] + one * (1-bitmasc[ichan1:ichan2,cut])
        endfor
     endif
  endfor

END


;---------

pro flag_mfi_gain, root_cal, jd, flag, bit=bit
  if not keyword_set(bit) then bit = 6

  global_btod_params,par

  ; variables
  nchan = 32
  one   = 2^(bit-1)

  sancho_check_activebit_inflags, flag, bit, bitmasc
  
  ; obtain cal file
  ffcal    = par.path_cal + root_cal + ".cal" 
  print," --> CAL file: ",ffcal

  ; compute flags only if exists
  if does_exists(ffcal) then begin
     
     a   = mrdfits(ffcal,1)
     
     readcol,par.filegain, id, THRES, format="I,F"
     
     ; Only HORN3 at the moment
     ;for ichan=0,nchan-1 do begin
     for ichan=16,23 do begin
        cut  = where(a.jd ge min(jd) and a.jd le max(jd) and a.flag[ichan,*] eq 0, ncut)
        stdv = 0.0
        if ncut ge 2 then begin
           gain = a.gain[ichan,cut]
           stdv = STDEV(gain)
        endif
        if (stdv ge 2.0*THRES[ichan] ) then flag[ichan,*] = flag[ichan,*] + one * (1-bitmasc[ichan,*])
        
     endfor
     
  endif else begin
     print,"     CAL file not found. " 
  endelse
  
  ; Flag 217x and 217y channels in period 1/5/2016 - 28/5/2016, where the 217x channel
  ; underwent a strong gain variation
  jd_obs = median(jd)
  jd_min = 2457509.5d0          ; JD corresponding to 1/5/2016, 00:00
  jd_max = 2457537.5d0          ; JD corresponding to 29/5/2016, 00:00
  
  if jd_obs ge jd_min and jd_obs le jd_max then begin
     flag[11,*] = flag[11,*] + one * (1-bitmasc[11,*]) ; Channel 217x
     flag[13,*] = flag[13,*] + one * (1-bitmasc[13,*]) ; Channel 217y
  endif

  
end


;=================

PRO prepare_btod_for_mfi, root, overwrite=overwrite, path=path, reflag=reflag, outpath=outpath

; Keyword
if not keyword_set(overwrite) then overwrite = 0
if not keyword_set(path) then path = ""
if not keyword_set(reflag) then reflag = 0


; Get global parameters
global_btod_params,params
info_quijote,mfi=mfi

; Binning
msbin = params.msbin

; Obtain all the TOD files associated
print," (*) MAIN PATH: ",path
print,"     ROOT name: ",root
spawn,"ls -1 " + path + root + "*.tod", ff_list, error

if error ne "" then begin
   print," > No file found with root name ",root
   return
endif

; Count TOD maps
ntod = n_elements(ff_list)
print," (*) Total number of TOD files: ",ntod


; Number of blocks
ntod_per_btod  = params.ntod_per_btod
one            = 1
if ((ntod mod ntod_per_btod) eq 0 ) then one = 0
nblocks        = max( [1, ntod/ntod_per_btod+one] )
print," (*) Number of blocks: ",nblocks



; Prepare output names and root names
if not keyword_set(outpath) then outpath = params.outpath
ff_out  = strarr(nblocks)
for i=0,nblocks-1 do begin
   txt       = strtrim(string(i),2)
   case strlen(txt) of
      1: tail = "-00"+txt 
      2: tail = "-0"+txt
      3: tail = "-"+txt
      else: stop," Max number of nblocks allowed is 1000"
   endcase
   if (nblocks eq 1) then tail = ""
   ff_out[i] = outpath + strtrim(root,2) + tail + ".btod"
endfor


; MAIN LOOP over blocks
;cfp
FOR iblock=0,nblocks-1 DO BEGIN

   ; range of TOD files
   itod_min = ntod_per_btod*iblock
   itod_max = min( [ntod-1, itod_min + ntod_per_btod -1 ] )

   print," (*) Block: ",iblock 
   print,"     Range of TOD files: ",itod_min,itod_max

   ; check if the output file already exists
   existe  = does_exists(ff_out[iblock])
   goahead = 1
   rebinea = 1
   if (existe) then goahead = 0
   if (existe and overwrite) then goahead = 1

   ; Reflag mode: reads btod from scratch, and only changes weights
   if (reflag eq 1) then begin
      goahead = 1
      rebinea = 0
      if (existe eq 0) then return 
      ; read from scratch
      btod_scratch = mrdfits(ff_out[iblock],1)
      print,"    REFLAG mode. Reading btod from ",ff_out[iblock]
   endif



   ; Main part of code
   if (goahead) then begin

      ; Update msbin value. Only need for nominal data is vel(AZ) = 12 deg/s
      update_msbin_from_hk, root, msbin


      primero  = 1
      FOR itod=itod_min,itod_max DO BEGIN
         

         filename = ff_list[itod]
         print," --> Reading ",filename
         print,"     File i within ntod = ",itod,itod_max


         if (rebinea eq 1) then begin


            ;----------------------
            ; [1] Read map
            sancho_read_raster_scan, filename, jd_all, az_all, el_all, data_all, time_all, mod_ang_all, cal
            n_all = n_elements(jd_all)
            nchan = n_elements(data_all[*,0])
         
            cut = where(az_all ne 0 and el_all ne 0)
            ;plot,az_all[cut],el_all[cut],/xs,/ys,xtitle="AZ (deg)",ytitle="EL(deg)"




            ;----------------------
            ; [2] Bin data
            ; Exclude blocks of CAL data
            wei_all = fltarr(n_all) + 1.0
            cut     = where(cal eq 1,ncut)
            wei_all[cut]    = 0.0
            wei_all[cut-1L] = 0.0
            wei_all[cut-2L] = 0.0
            wei_all[cut-3L] = 0.0
            


            ; Bin data. Reference data vectors are:
            ;    data, jd, el, az, mod_ang, sigma
            sancho_bin_data_mfi, data_all, jd_all, el_all, az_all, mod_ang_all, wei_all, $
                             data, jd, el, az, mod_ang, sigma, sigmasum, sigmadiff, sigmacov, msbin=msbin
            
            n = n_elements(jd)


            ;----------------------
            ; [3] Obtaining the pointing for each horn
            sancho_mfi_get_pointingmodel, jd, az, el, gl, gb, par, azhorn, elhorn, version= params.POINTINGMODEL 
      
      
            ;----------------------
            ; [4] Define final weights
            wei = fltarr(nchan,n)
            cut = where(sigma ne 0.0, ncut)
            if (ncut ge 1) then wei[cut] = 1.0/sigma[cut]^2.

         endif else begin

            ; [1-4] Reads from pre-existing btod file.
            ; It only reads the JD array.
            ;fxbopen,unit,filename, 1
            ;fxbread,unit,jd_all,"JD"
            ;fxbclose,unit
            sancho_read_raster_scan, filename, jd_all, az_all, el_all, data_all, time_all, mod_ang_all, cal
            n_all = n_elements(jd_all)
            nchan = 32
         
            cut = where(btod_scratch.jd ge min(jd_all) and btod_scratch.jd le max(jd_all) )
            n   = n_elements(cut)

            jd = btod_scratch.jd[cut]
            az = btod_scratch.az[cut]
            el = btod_scratch.el[cut]

            gl   = fltarr(4,n)
            gb   = fltarr(4,n)
            par  = fltarr(4,n)
            gl   = btod_scratch.gl[*,cut]
            gb   = btod_scratch.gb[*,cut]
            par  = btod_scratch.par[*,cut]

            data = fltarr(32,n)
            data = btod_scratch.data[*,cut]

            mod_ang = fltarr(4,n)
            mod_ang = btod_scratch.mod_ang[*,cut]

            wei = fltarr(4,n)
            wei = btod_scratch.wei[*,cut]

	    azhorn = btod_scratch.azhorn[cut]
	    elhorn = btod_scratch.elhorn[cut]

	    sigmacov = fltarr(4,2,2,n)
	    sigmacov = btod_scratch.sigmacov[*,*,*,cut]

         endelse

         

         ;----------------------
         ; [5] Flagging. 
	 ; Initialize to zero values.
         flag = intarr(nchan,n)
         
      
         ; [5.1] HK data
         print," > Flagging using HK data."
         flag_mfi_hk, root, jd, flag, bit=1

      
         ; [5.2] Bad data
         print," > Flagging bad data."
         flag_bad_data, data, flag, bit=2


         ; [5.3] Sun, Moon. Do this only once. Otherwise,
         ; the bit will be activated twice.
         print," > Flagging data due to Sun and Moon."
         bit = 3
         one =  2^(bit-1)
         flag_sun_moon, jd, gl[0,*], gb[0,*], flag1
         flag_sun_moon, jd, gl[1,*], gb[1,*], flag2
         flag_sun_moon, jd, gl[2,*], gb[2,*], flag3
         flag_sun_moon, jd, gl[3,*], gb[3,*], flag4
         
         for ichan=0,nchan-1 do begin
            ihorn = mfi.polarizer_id[ichan]
            if (ihorn eq 1) then flag[ichan,*] = flag[ichan,*] + flag1[*]*one
            if (ihorn eq 2) then flag[ichan,*] = flag[ichan,*] + flag2[*]*one
            if (ihorn eq 3) then flag[ichan,*] = flag[ichan,*] + flag3[*]*one
            if (ihorn eq 4) then flag[ichan,*] = flag[ichan,*] + flag4[*]*one
         endfor


         ; [5.4] Satellites.
         print," > Flagging satellites."
         flag_mfi_satellites, az,el, flag, bit=4
         flag_mfi_sat_nsl, az,el, flag, bit=4


         ; [5.5] Bad weather
         print," > No flag based on bad weather yet."
         bit = 5
      

         ; [5.6] Gain drifts
         print," > Flagging according to gain."
         flag_mfi_gain, root, jd, flag, bit=6 

         
         ; [5.7] Flag for planets: Mars, Venus, Jupiter
         bit = 7 
         
         ; [5.8] Free. To be defined by the user.
         bit = 8

         ; [5.9] Flag based on RMS (from Fred)
         bit = 9 
         
         
         ; [6] Store data 
         if (primero eq 1) then begin
            jd_btod      = jd
            az_btod      = az
            el_btod      = el
            gl_btod      = gl
            gb_btod      = gb
            par_btod     = par
            data_btod    = data
            mod_ang_btod = mod_ang
            wei_btod     = wei
            flag_btod    = flag
            azhorn_btod  = azhorn
            elhorn_btod  = elhorn
	    sigmacov_btod= sigmacov
            primero      = 0
         endif else begin
            jd_btod      = [jd_btod, jd]
            az_btod      = [az_btod, az]
            el_btod      = [el_btod, el]
            gl_btod      = [[gl_btod], [gl]]
            gb_btod      = [[gb_btod], [gb]]
            par_btod     = [[par_btod], [par]]
            data_btod    = [[data_btod], [data]]
            mod_ang_btod = [[mod_ang_btod], [mod_ang]]
            wei_btod     = [[wei_btod], [wei]]
            flag_btod    = [[flag_btod], [flag]]
            azhorn_btod  = [[azhorn_btod], [azhorn]]
            elhorn_btod  = [[elhorn_btod], [elhorn]]

	    n1 = n_elements(sigmacov_btod[0,0,0,*])
	    n2 = n_elements(sigmacov[0,0,0,*])
	    sigmacov_btod= reform([[reform(sigmacov_btod,16,n1)],[reform(sigmacov,16,n2)] ],4,2,2,n1+n2)
         endelse


      ENDFOR

      
      ; [6.1] Get the gains for these files. Uses also the info from this code
      mfi_params, JD_OBS=jd_btod[0], GAINS_IC=gains_ic, mod_ref_angle=mod_ref_angle, /vector
      
         
      ;----------------------
      ; [7] Write data
      a = {jd:jd_btod, az:az_btod, el:el_btod, gl:gl_btod, gb:gb_btod, par:par_btod, $
           data:data_btod, mod_ang:mod_ang_btod, wei:wei_btod, flag:flag_btod, $
           gains:gains_ic, mod_ref_angle:mod_ref_angle, msbin:msbin, pointingmodel:params.POINTINGMODEL, $
	   azhorn:azhorn_btod, elhorn:elhorn_btod, sigmacov:sigmacov_btod }
      mwrfits, a, ff_out[iblock], /create


   endif else begin
      print," --> File already exists: ",ff_out[iblock]
   endelse

ENDFOR

end

;=====================

PRO sancho_genera_btod, filename, path=path, reflag=reflag, overwrite=overwrite, outpath=outpath

; Keywords
if not keyword_set(path) then path = "/net/nas/proyectos/quijote/tod/"
if not keyword_set(reflag) then reflag = 0
if not keyword_set(overwrite) then overwrite = 0

; Info
global_btod_params,params
if not keyword_set(outpath) then outpath=params.outpath
if n_params() eq 0 then begin
   print,""
   print,"   Syntax --  sancho_genera_btod, filename  [OVERWRITE=, PATH=, REFLAG=]"
   print,""
   print,"   Default path to TOD files is PATH = ",path
   print,"   Default output path of BTOD files = ",outpath
   print,"   Default pointing model            = ",params.POINTINGMODEL
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

print," (*) Total number of TOD files: ",nlist


; Run main code
for i=0L,nlist-1L do begin
   prepare_btod_for_mfi, lista[i], path=path, overwrite=overwrite, reflag=reflag, outpath=outpath
endfor

END
