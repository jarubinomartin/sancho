;   21/Nov/2017   JARM
;
;   Reads .cal file for MFI
;   
;   HISTORY:
;      21/10/2017 - original version. Based on sancho_read_mfi_hk.pro
;      12/01/2018 - update paths from nas4 to nas      
;
;-
;


PRO sancho_read_mfi_cal, root, a, jd, gain, path=path, verbose=verbose,status=status


; Keywords
if not keyword_set(path) then path = "/net/nas/proyectos/quijote/cal/"  
if not keyword_set(verbose) then verbose = 0
if not keyword_set(status) then status = 0

; Info
if n_params() eq 0 then begin
   print,""
   print,"   Syntax -- sancho_read_mfi_cal, root, output_structure, jd, gain [, path=, status=, /verbose] "
   print,""
   print,"   Default path for CAL data is ",path
   print,"   On output, status=0 if everything is ok. status=1 if no file read."
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

ffhk  = path + root_hk + ".cal" 
ffhk2 = path + root2 + ".cal" 
if does_exists(ffhk2) then ffhk = ffhk2

print,"  --> CAL file: ",ffhk

; Reading
a = 0
status = 1
if does_exists(ffhk) then begin
   
   a = mrdfits(ffhk,1)
   if (min(a.jd) lt 1d4) then a.jd = a.jd + mfi.jd_ref
   status = 0
   
; Parameters
   jd      = a.jd 
   gain    = a.gain
   rain    = a.rain
   base    = a.base
   acti    = a.acti
   flag    = a.flag
   mod_ang = a.mod_ang
   rms     = a.rms

; Verbose
   if (verbose ne 0) then begin
      
      print," (*) Basic information: "
      print,"     Root         = ",root
      print,""
      print," (*) Minmax values: "
      print,"     gain  = ",minmax(gain)
      print,"     base  = ",minmax(base)
      print,"     acti  = ",minmax(acti)
      print,""
   endif

endif

END  
