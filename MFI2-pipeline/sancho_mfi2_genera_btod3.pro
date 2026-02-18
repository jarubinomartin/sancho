;   24/05/2024  J.A.Rubino
;
;   PURPOSE:
;   Generates .btod3 files for QUIJOTE MFI2 data, using the old MFI DAS
;
;   HISTORY: 
;   24/05/2024 - original version. Based on sancho_genera_btod.pro
;   12/11/2024 - including CAL signal.
;
;  

PRO global_btod3_params, par

  msbin      = 20.0             ;60.0             ; ms. Bin size
  msbin_nomi = 20.0             ; ms. Bin size used for nominal mode
  
  outpath  = "/net/calp-nas/proyectos/quijote3/btod/"

  volt_max = 3.0 ; Volts. Maximum voltage allowed 
  volt_min = 0.0 ; Volts. Minimum voltage allowed 

  ntod_per_btod = 20 ; Number of TOD3 files in a BTOD3.

  ; Radii around sources for flagging. Using MFI values
  radius_sat  =  5.0            ; deg
  radius_nsl  =  5.0            ; deg. nsl=near side lobe
  radius_sun  = 10.0            ; deg
  radius_moon = 10.0            ; deg


  ; Files for HK. Not updated yet. Using MFI values
  path_hk           = "/net/calp-nas/proyectos/quijote3/tod/"
  MAX_PRESS         = 1d-4 ; bar. Updated. Before we used 2d-6 and 4d-6 from March 21st 2016.
  MAX_20K_SHIELD    = 30.0      ; K
  MAX_TL_20K_BASE   = 18.0      ; K
  MAX_CCC_1ST       = 40.0      ; K
  MAX_DIFF_TEMP_BEM =  2.5      ; Celsius. Chosen by eye by JARM.


  ; Pointing model
  POINTINGMODEL = "Version1" ; "Version1" from May-2024. 

  ; structure
  par     = { msbin:msbin, outpath:outpath, volt_max:volt_max, volt_min:volt_min, $
              ntod_per_btod:ntod_per_btod, nominal_msbin:msbin_nomi, $
              radius_sat:radius_sat, radius_nsl:radius_nsl, radius_sun:radius_sun, $
              radius_moon:radius_moon, path_hk:path_hk, $
              MAX_PRESS:MAX_PRESS, MAX_20K_SHIELD:MAX_20K_SHIELD, $
              MAX_TL_20K_BASE:MAX_TL_20K_BASE, MAX_CCC_1ST:MAX_CCC_1ST, $
              MAX_DIFF_TEMP_BEM:MAX_DIFF_TEMP_BEM, POINTINGMODEL:POINTINGMODEL }

END


;===================

pro update_hk_temp_limits, jd0,  MAX_20K_SHIELD, MAX_TL_20K_BASE, MAX_CCC_1ST 

  global_btod3_params, par
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

pro flag_mfi2_hk, root_hk, jd, flag, bit=bit
  if not keyword_set(bit) then bit = 1
  
  global_btod3_params,par
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


PRO flag_mfi2_bad_data, data, flag, bit=bit
  if not keyword_set(bit) then bit = 2
  global_btod3_params,par
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
  global_btod3_params,par
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

PRO flag_mfi2_satellites, az, el, flag, bit=bit
  if not keyword_set(bit) then bit = 4

  global_btod3_params,par
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


PRO flag_mfi2_sat_nsl, az, el, flag, bit=bit
  if not keyword_set(bit) then bit = 4

  global_btod3_params,par
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


;=================

PRO prepare_btod3_for_mfi2, root, overwrite=overwrite, path=path, reflag=reflag, outpath=outpath, nhorns=nhorns

; Keyword
if not keyword_set(overwrite) then overwrite = 0
if not keyword_set(path) then path = ""
if not keyword_set(reflag) then reflag = 0
if not keyword_set(nhorns) then nhorns=4 ; using old MFI DAS

; Get global parameters
global_btod3_params,params
info_quijote,mfi=mfi

; Binning
msbin = params.msbin

; Obtain all the TOD files associated
print," (*) MAIN PATH: ",path
print,"     ROOT name: ",root
spawn,"ls -1 " + path + root + "*.tod3", ff_list, error

if error ne "" then begin
   print," > No file found with root name ",root
   return
endif

; Count TOD maps
ntod = n_elements(ff_list)
print," (*) Total number of TOD3 files: ",ntod


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
   ff_out[i] = outpath + strtrim(root,2) + tail + ".btod3"
endfor


; MAIN LOOP over blocks
;cfp
FOR iblock=0,nblocks-1 DO BEGIN

   ; range of TOD files
   itod_min = ntod_per_btod*iblock
   itod_max = min( [ntod-1, itod_min + ntod_per_btod -1 ] )

   print," (*) Block: ",iblock 
   print,"     Range of TOD3 files: ",itod_min,itod_max

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
      print,"    REFLAG mode. Reading btod3 from ",ff_out[iblock]
   endif


   ; Main part of code
   if (goahead) then begin

      primero  = 1
      FOR itod=itod_min,itod_max DO BEGIN
         
         filename = ff_list[itod]
         print," --> Reading ",filename
         print,"     File i within ntod = ",itod,itod_max


         if (rebinea eq 1) then begin

            ;----------------------
            ; [1] Read map
            sancho_read_mfi2_tod, filename, jd_all, az_all, el_all, data_all, time_all, cal
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

            
            sancho_bin_data_mfi2, data_all, jd_all, el_all, az_all, wei_all, $
                                  data, jd, el, az, sigma, sigmasum, sigmadiff, sigmacov, msbin=msbin
            
            n = n_elements(jd)


            ;----------------------
            ; [3] Obtaining the pointing for each horn
            sancho_mfi2_get_pointingmodel, jd, az, el, gl, gb, par, azhorn, elhorn, $
                                           version= params.POINTINGMODEL, NHORNS=nhorns ; using old MFI DAS 
      
      
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
            sancho_read_mfi2_tod, filename, jd_all, az_all, el_all, data_all, time_all
            n_all = n_elements(jd_all)
            nchan = 32
         
            cut = where(btod_scratch.jd ge min(jd_all) and btod_scratch.jd le max(jd_all) )
            n   = n_elements(cut)

            jd = btod_scratch.jd[cut]
            az = btod_scratch.az[cut]
            el = btod_scratch.el[cut]

            gl   = fltarr(nhorns,n)
            gb   = fltarr(nhorns,n)
            par  = fltarr(nhorns,n)
            gl   = btod_scratch.gl[*,cut]
            gb   = btod_scratch.gb[*,cut]
            par  = btod_scratch.par[*,cut]

            data = fltarr(nchan,n)
            data = btod_scratch.data[*,cut]

            wei = fltarr(nhorns,n)
            wei = btod_scratch.wei[*,cut]

	    azhorn = btod_scratch.azhorn[cut]
	    elhorn = btod_scratch.elhorn[cut]

	    sigmacov = fltarr(nhorns,2,2,n)
	    sigmacov = btod_scratch.sigmacov[*,*,*,cut]

         endelse

         

         ;----------------------
         ; [5] Flagging. 
	 ; Initialize to zero values.
         flag = intarr(nchan,n)
         
         ; [5.2] Bad data
         print," > Flagging bad data."
         flag_mfi2_bad_data, data, flag, bit=2


         ; [5.3] Sun, Moon. Do this only once. Otherwise,
         ; the bit will be activated twice.
         print," > Flagging data due to Sun and Moon."
         bit = 3
         one =  2^(bit-1)

         for ihorn=0,nhorns-1 do begin
            flag_sun_moon, jd, gl[ihorn,*], gb[ihorn,*], flag1
            ichan1 = ihorn*8
            ichan2 = ichan1 + 7
            for ichan=ichan1,ichan2 do begin
               flag[ichan,*] = flag[ichan,*] + flag1[*]*one
            endfor
         endfor

         ; [5.4] Satellites.
         print," > No flagging satellites yet."
         bit = 4

         ; [5.5] Bad weather
         print," > No flag based on bad weather yet."
         bit = 5
      

         ; [5.6] Gain drifts
         print," > No flagging according to gain drift yet."
         bit = 6
         
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
            wei_btod     = [[wei_btod], [wei]]
            flag_btod    = [[flag_btod], [flag]]
            azhorn_btod  = [[azhorn_btod], [azhorn]]
            elhorn_btod  = [[elhorn_btod], [elhorn]]

	    n1 = n_elements(sigmacov_btod[0,0,0,*])
	    n2 = n_elements(sigmacov[0,0,0,*])
	    sigmacov_btod= reform([[reform(sigmacov_btod,16,n1)],[reform(sigmacov,16,n2)] ],4,2,2,n1+n2)
         endelse


      ENDFOR

      
      ;----------------------
      ; [7] Write data
      a = {nhorns:nhorns, jd:jd_btod, az:az_btod, el:el_btod, gl:gl_btod, gb:gb_btod, par:par_btod, $
           data:data_btod, wei:wei_btod, flag:flag_btod, msbin:msbin, pointingmodel:params.POINTINGMODEL, $
	   azhorn:azhorn_btod, elhorn:elhorn_btod, sigmacov:sigmacov_btod }
      mwrfits, a, ff_out[iblock], /create


   endif else begin
      print," --> File already exists: ",ff_out[iblock]
   endelse

ENDFOR

end

;=====================

PRO sancho_mfi2_genera_btod3, filename, path=path, reflag=reflag, overwrite=overwrite, outpath=outpath, nhorns=nhorns

; Keywords
if not keyword_set(path) then path = "/net/calp-nas/proyectos/quijote3/tod/"
if not keyword_set(reflag) then reflag = 0
if not keyword_set(overwrite) then overwrite = 0
if not keyword_set(nhorns) then nhorns=4 ; using old MFI DAS

; Info
global_btod3_params,params
if not keyword_set(outpath) then outpath=params.outpath
if n_params() eq 0 then begin
   print,""
   print,"   Syntax --  sancho_mfi2_genera_btod3, filename  [OVERWRITE=, PATH=, REFLAG=]"
   print,""
   print,"   Default path to TOD3 files is PATH = ",path
   print,"   Default output path of BTOD3 files = ",outpath
   print,"   Default pointing model             = ",params.POINTINGMODEL
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

print," (*) Total number of TOD3 files: ",nlist


; Run main code
for i=0L,nlist-1L do begin
   prepare_btod3_for_mfi2, lista[i], path=path, overwrite=overwrite, reflag=reflag, outpath=outpath, nhorns=nhorns
endfor

END
