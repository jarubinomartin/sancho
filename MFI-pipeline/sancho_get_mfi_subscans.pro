;   17/March/2016   JARM
;
;   For a given MFI btod file, it returns the scans.
;
;  HISTORY:
;     17/03/2016 - original version
;     04/05/2017 - updated version. If HK file is missing, it checks if
;                  nominal. Otherwise, it assumes that it is raster mode.
;     12/01/2018 - updated paths from nas4 to nas
;     19/02/2018 - updated path to hk directory
;
;-
;

PRO sancho_get_mfi_subscans, root, az_in, el_in, k_ini, k_end, path_hk=path_hk, path_btod = path_btod, doplot=doplot,$
                             obsmode=obsmode

; Keywords
if not keyword_set(path_hk) then path_hk = "/net/nas/proyectos/quijote/tod/"
if not keyword_set(path_btod) then path_btod = "/net/nas/proyectos/quijote/btod/"
if not keyword_set(doplot) then doplot = 0


; Global params
info_quijote,mfi=mfi


; Info
if n_params() eq 0 then begin
   print,""
   print,"  Syntax --  sancho_get_mfi_subscans, root, az_in, el_in "
   print,""
   print,"  > Examples: IDL> sancho_get_mfi_subscans, 'W63B-150418-0353'"
   print,""
   root = "W63B-150418-0353"    ; "CASSA-151218-1801"
   return
endif


; Read BTOD file if no (AZ,EL) is provided on entry.
if n_params() eq 1 then begin
   ff = path_btod + root + ".btod"
   a     = mrdfits(ff,1)
   az_in = a.az
   el_in = a.el
endif


; Read HK file and obtain information
nlen  = strlen(root)
root2 = strmid(root,0,nlen-4) ; removes -001 if present

ffhk   = path_hk + root  + ".hk"
ffhk2  = path_hk + root2 + ".hk" 
if does_exists(ffhk2) then ffhk = ffhk2

check_if_nominal, az_in, is_nominal

txtmodes = ["OFF","Standby","Homing","Homed","Nominal","Pointing",$
            "tracking","blind spot", "parking","raster", "sky raster", $
            "absolute moving", "relative moving" ]

if does_exists(ffhk) then begin
   print,"  --> HK file: ",ffhk
   b     = mrdfits(ffhk,1)
   jd_hk = b.jd + mfi.jd_ref
   obsmode = b.obs_params[0]    ; 10=raster, 11=sky raster, 5=nominal, 4=Homed
   if (obsmode eq 5) then az_vel  = b.obs_params[1]
   if (obsmode ne 5) then az_vel  = b.obs_params[12]
endif else begin
   ; default values
   obsmode = 10                 ; raster
   if (is_nominal) then obsmode = 5 
   az_vel  = 1.0
   print,"  --> No HK file read. Assuming ",txtmodes[obsmode-1]
endelse



; [1] NOMINAL MODE
status = 0
IF (obsmode eq 5) THEN BEGIN

   status = 1
   
   ; 1.1 Identify changes in the AZ angle
   az = az_in
   n  = n_elements(az)
   di = az-shift(az,1)
   id = intarr(n)
   id[where(di ge 0)] =  1.0
   id[where(di lt 0)] = -1.0

   uuu   = where(id lt 0)

   
ENDIF

; [2] RASTER MODE
IF (obsmode eq 10 or obsmode eq 11) THEN BEGIN

   status = 1
   
   ; 2.1 First, it checks the range of AZ. Shifts in case of close to 0 or 360.0
   az   = az_in
   sancho_get_azel_radec_ranges, az, el_in, az0, el0, azmin, azmax, elmin,elmax
   dist = azmax - azmin
   if (dist ge 355.0) then az = (az_in + 180.0) mod 360.0


   ; 2.2 Smooths the AZ distribution in scales of 1-2s.
   n  = n_elements(az)
   az = smooth(az,40)

   ; 2.3 Identify changes in the scan direction
   di = az-shift(az,1)
   id = intarr(n)
   id[where(di ge 0)] =  1.0
   id[where(di lt 0)] = -1.0

   uuu   = uniq(id)

ENDIF

; Now computes ini and end of each subscan.
if (status eq 1) then begin
   if uuu[0] eq 0 then begin
      ns    = n_elements(uuu)
      k_ini = [0, uuu[1:ns-1L]+1L]
      k_end = [uuu[1:ns-1L],n-1L ]
   endif else begin
      k_ini = [0, uuu+1L]
      k_end = [uuu,n-1L ]
   endelse
   
   if max(k_ini) gt (n-1) then begin
      ns    = n_elements(k_ini)
      k_ini = k_ini[0:ns-2]
      k_end = k_end[0:ns-2]
   endif
   
   nscans = n_elements(k_ini)-1


   ; Plots
   if (doplot) then begin
      cfp
      window,0,xs=1800
      nmax = n
      if (nmax gt 50000L) then nmax = 50000L
      plot,az_in,xrange=[0,nmax],yr=[0,360],ystyle=1,xtitle="No of sample",ytitle="AZ (deg)"
      for i=0,nscans do vline,k_ini[i],col=2,linestyle=1
   endif
   
endif


; If no subscans where selected, then uses full range
if status eq 0 then begin
   nscans = 1
   k_ini  = [0L]
   k_end  = [n_elements(az_in)-1L]
endif

; Display info on screen
print,"(*) Derived number of scans = ", nscans
print,"    Average lenght of scans = ", round(median(k_end-k_ini))
print,"    Mean lenght of scans    = ", round(mean(k_end-k_ini))
print,"    Observing mode          = ", round(obsmode),", ",txtmodes[obsmode-1]

;print,k_ini
;print,k_end
;print,n
;print,n_elements(k_ini)


END
