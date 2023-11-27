;   21/April/2016   JARM
;
;   Reads HK data for MFI
;   
;   HISTORY:
;      21/04/2016 - original version
;      31/10/2017 - updated version. Now checks the correct name of HK.
;                   It returns the JD corrected, and obsmode and AZ_vel as keywords
;      19/02/2018 - updated path to HK directory (from quijote3 to proyectos/quijote/)
;      01/07/2020 - aditional parameters extracted on output: jd_hk, press, tbem1/2, 
;-
;


PRO sancho_read_mfi_hk, root, a, path=path, obsmode=obsmode, az_vel=az_vel,verbose=verbose, el=el, $
	press=press, tbem1=tbem1, tbem2=tbem2, jd_hk=jd_hk, t_out1=t_out1, t_out2=t_out2


; Keywords
if not keyword_set(path) then path = "/net/nas/proyectos/quijote/tod/"  
if not keyword_set(verbose) then verbose = 0

; Info
if n_params() eq 0 then begin
   print,""
   print,"   Syntax -- sancho_read_mfi_hk, root, output_structure [, path=, obsmode=, az_vel=, /verbose] "
   print,""
   print,"   Default path for HK data is ",path
   print,""
   return
endif

; Global variables
info_quijote,mfi=mfi

; Check root name
kk1     = strsplit(root,"/",/extract)
kk2     = kk1[n_elements(kk1)-1]
kk3     = strsplit(kk2,".",/extract)
root_hk = strtrim(kk3[0],2)

; HK file
nlen  = strlen(root_hk)
root2 = strmid(root_hk,0,nlen-4) ; removes -001 if present, for nominal data

ffhk  = path + root_hk + ".hk" 
ffhk2 = path + root2 + ".hk" 
if does_exists(ffhk2) then ffhk = ffhk2

print,"  --> HK file: ",ffhk

; Reading
a = {jd:0, temp:fltarr(20,1), temp_bem:fltarr(2,1), press:0 }
obsmode = 0
az_vel  = 0
el      = 0
sancho_extract_parameters_mfi_hk, a, jd_hk, T_MOTOR1, T_MOTOR2, T_MOTOR3, T_MOTOR4, T_OUT1, T_OUT2, T_20K_BasePla, T_OMT1, $
			T_80K_cryobr, T_OMT3, T_OMT2, T_20K_SHIELD, T_TL_20K_BASE, T_CCC_1ST, TBEM1,TBEM2, PRESS
if does_exists(ffhk) then begin
   
   a = mrdfits(ffhk,1)
   if (min(a.jd) lt 1d4) then a.jd = a.jd + mfi.jd_ref
   
; Parameters
   sancho_extract_parameters_mfi_hk, a, jd_hk, T_MOTOR1, T_MOTOR2, T_MOTOR3, T_MOTOR4, T_OUT1, T_OUT2, T_20K_BasePla, T_OMT1, $
			T_80K_cryobr, T_OMT3, T_OMT2, T_20K_SHIELD, T_TL_20K_BASE, T_CCC_1ST, TBEM1,TBEM2, PRESS

; Obsmode: 10=raster, 11=sky raster, 5=nominal, 4=Homed
   obsmode       = fix(a.obs_params[0])
; AZ velocity (deg/s)
   if (obsmode eq 5) then az_vel = a.obs_params[1]
   if (obsmode ne 5) then az_vel = a.obs_params[12]

; EL (deg)
   if (obsmode eq 5) then el = a.obs_params[2]  ; EL for nominal
   if (obsmode ne 5) then el = a.obs_params[10] ; Commanded initial EL value.

; Verbose
   if (verbose ne 0) then begin
      txtmodes = ["OFF","Standby","Homing","Homed","Nominal","Pointing",$
                  "tracking","blind spot", "parking","raster", "sky raster", $
                  "absolute moving", "relative moving" ]
      
      print," (*) Basic information: "
      print,"     Root         = ",root
      print,"     Obsmode      = ",obsmode
      print,"     Txt_obsmode  = ",txtmodes[obsmode-1]
      print,"     v_AZ (deg/s) = ",az_vel 
      print,"     EL (deg)     = ",el
      print,""
      print," (*) Minmax values: "
      print,"     Press (mbar)      = ",minmax(PRESS)
      print,"     TBEM1 (C)         = ",minmax(TBEM1)
      print,"     TBEM2 (C)         = ",minmax(TBEM2)
      print,"     T_20K_SHIELD (K)  = ",minmax(T_20K_SHIELD)
      print,"     T_TL_20K_BASE (K) = ",minmax(T_TL_20K_BASE)
      print,"     T_CCC_1ST (K)     = ",minmax(T_CCC_1ST)
      print,""
   endif

endif

END  
