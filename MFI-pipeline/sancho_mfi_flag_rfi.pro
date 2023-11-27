;  17/10/2017   JARM
;
;  Flags MFI calibrated data based on the RMS
;
;  HISTORY:
;     17/10/2017 - original version
;     12/01/2018 - updated paths (nas4-->nas/proyectos)
;     19/01/2021 - updated to check active bits
;
;-
;
 

PRO sancho_mfi_flag_rfi, root, jd, az, el, gl, gb, data, flag, reflag, bit=bit, display=display, overwrite=overwrite, path=path, outpath=outpath, noflag=noflag, sky=sky, $
	perelev=perelev, withcov=withcov

; Keywords
if not keyword_set(path) then path = "/net/nas/proyectos/quijote/ctod/"
if not keyword_set(outpath) then outpath = "/net/nas/proyectos/quijote/flag_stats_ctod/jalberto/"
if not keyword_set(display) then display = 0
if not keyword_set(overwrite) then overwrite=0
if not keyword_set(bit) then bit = 8
if not keyword_set(noflag) then noflag = fltarr(nside2npix(128))
if not keyword_set(perelev) then perelev = 0
if not keyword_set(withcov) then withcov = 0

one = 2^(bit-1)


; Info
if n_params() eq 0 and not keyword_set(ffctod) then begin
	print,""
	print,"  Syntax -- sancho_mfi_flag_rfi, root, jd, az, el, gl, gb, data, flag, reflag, bit=bit [/OVERWRITE, noflag=, sky=] "
	print,""
	print,"  Note: if only the root name is provide, it reads the file from disk. "
	print,""
	print,'  Example:  sancho_mfi_flag_rfi, "NOMINAL30A-160907-1635-001" '
	print,""
	return
endif

; Read the ctod file if not proving all fields
;  Example.  root   = "NOMINAL30A-160907-1635-001"
ffctod = path + root + ".ctod" 
if n_params() lt 8 then begin
	print," ==> Reading CTOD from file ",ffctod
	a=mrdfits(ffctod,1)
	jd = a.jd
	az = a.az
	el = a.el
	gl = a.gl
	gb = a.gb
	data = a.data
	flag = a.flag
endif
sancho_check_activebit_inflags, flag, bit, bitmasc

; Parameters
info_quijote,mfi=mfi
nchan  = n_elements(data[*,0])
n      = n_elements(data[0,*])
djd    = median(shift(jd,-1)-jd)
factor = 24.d * 60.d * 60.d  ; secs per JDday
dt     = djd * factor        ; secs
tkernel_sec = 30.0           ; secs. 
if not keyword_set(sky) then sky = 2 ; Sky model


; Identify if nominal data. obsmode: 5=nominal, 10=raster
sancho_read_mfi_hk, root, str_hk, obsmode=obsmode, el=el_nominal, /verbose
check_if_nominal, az, is_nominal
if (is_nominal) then obsmode = 5
txtfield = "_nominal_sky"+string(sky,format="(i1)")
if obsmode ne 5 then txtfield = ""

; Identify period. Based on RGS definition of periods
find_mfi_period,jd[0],periodo


; Identify elevation
elev = round(median(el)) 
if elev ne el_nominal then print," SANCHO_MFI_FLAG_RFI: inconsistent value of elevation. "
txtelev = string(elev,format="(i2)")


; Check if a mask is used, and define "extra" text
domask = max(noflag) gt 0.0 
extra    = "_nomask"
if domask then extra = "_mask"  


; FILENAME for thresholds. For nominal files, this is thres_nominal_periodN.fits
;get_filename_for_thresholds, periodo, sky, ffthres, obsmode=obsmode, txtextra=extra
ffthres  = "/home/jalberto/quijote/idl/sancho/thres.fits"
;if periodo ge 1 then ffthres = "/home/jalberto/quijote/idl/sancho/thres"+txtfield+"_period"+string(periodo,format="(i1)")+"_20ms.fits"
;if periodo ge 1 then ffthres = "/home/jalberto/quijote/idl/sancho/thres"+txtfield+"_period"+string(periodo,format="(i1)")+".fits"
if periodo ge 1 then ffthres = "/home/jalberto/quijote/idl/sancho/thres"+txtfield+"_period"+string(periodo,format="(i1)")+extra+".fits"
if periodo ge 1 and withcov eq 1 then ffthres = "/home/jalberto/quijote/idl/sancho/thres"+txtfield+"_period"+string(periodo,format="(i1)")+extra+"_withcov.fits"
if periodo ge 1 and perelev ne 0 then ffthres = "/home/jalberto/quijote/idl/sancho/thres"+txtfield+"_period"+string(periodo,format="(i1)")+"_elev"+txtelev+extra+".fits" 
if periodo ge 1 and withcov eq 1 then ffthres = "/home/jalberto/quijote/idl/sancho/thres"+txtfield+"_period"+string(periodo,format="(i1)")+extra+"_withcov2.fits"
print," > THRESHOLD file: ", ffthres



; Output file with the statistics (separate mask and nomask)
outpath_sky = outpath
if sky ne 2 then outpath_sky = outpath + "sky"+string(sky,format="(i1)")+"/"
ffout = outpath_sky + root +"_stats"+extra+".fits"


; CHECK IF ffout EXISTS
existe = does_exists(ffout)

; READ INFO IF STAT FILE EXISTS
saltacalc = 0
if (existe and (overwrite eq 0)) then begin
	print,"  ==> Reading STAT from file ",ffout
	b           = mrdfits(ffout,1)
	tkernel_sec = b.tkernel_sec
	NSIGUP      = b.NSIGUP
	NSIGLOW     = b.NSIGLOW
	sigsum      = b.sigsum
	sigdiff     = b.sigdiff
	sigma       = b.sigma
	sky         = b.sky
	saltacalc   = 1
endif
if (existe eq 0) then overwrite=1


; Bins in time
ndk         = round( tkernel_sec / dt)
ntot        = n/ndk
if (ntot*ndk lt n) then ntot=ntot+1

imin = lindgen(ntot)*ndk
imax = imin+ndk-1
dime = where(imax gt n-1, ndime)
if ndime ge 1 then imax[dime] = n-1


; SKIP computations if info already read
IF (saltacalc eq 0) THEN BEGIN

; Get sky model
sancho_mfi_skymodel, gl, gb, data, model, sky=sky

; Define mask for statistics. 
masc  = fltarr(4,n)
if domask then begin
	d2r   = !DPI/180.0
	for ihorn=0,3 do begin
		theta = (90.0-gb[ihorn,*])*d2r
		phi   = gl[ihorn,*]*d2r
		ANG2PIX_RING, 128, theta, phi, ipix
		masc[ihorn,*] = noflag[ipix[*]]
	endfor	
endif

; 1. Using the rms of each channel.
sigma = fltarr(nchan,ntot)
for ichan=0,nchan-1 do begin
	ihorn = mfi.POLARIZER_ID[ichan]-1
	for j=0,ntot-1 do begin
		i1    = imin[j] 
		i2    = imax[j]
		datos = data[ichan,i1:i2]-model[ichan,i1:i2]
		flags = flag[ichan,i1:i2]
		masca = masc[ihorn,i1:i2]

		cut   = where(flags eq 0 and datos ne 0.0 and masca eq 0, ncut)
		if ncut ge 3 then begin
			sigma[ichan,j] = stddev(datos[cut])
		endif
	endfor
endfor

; 2. Using the rms of the sum and difference
sigsum  = fltarr(4,2,2,ntot)
sigdiff = fltarr(4,2,2,ntot)
for ihorn=0,3 do begin
	for ifreq=0,1 do begin
		for icorr=0,1 do begin
			
			ichan_s = mfi.ic[ihorn,ifreq,0] ; X+Y
			ichan_d = mfi.ic[ihorn,ifreq,1] ; X-Y
			if icorr eq 1 then begin
				ichan_s = mfi.ic[ihorn,ifreq,2] ; X
				ichan_d = mfi.ic[ihorn,ifreq,3] ; Y
			endif
			
			;print,ihorn,ifreq,icorr,ichan_s,ichan_d," ", mfi.chan_name[ichan_s]," ",mfi.chan_name[ichan_d]
			
			for j=0,ntot-1 do begin
				i1    = imin[j] 
				i2    = imax[j]
				
				datos_SUM   = data[ichan_s,i1:i2]  + data[ichan_d,i1:i2]
				modelo_SUM  = model[ichan_s,i1:i2] + model[ichan_d,i1:i2]
				datos_DIFF  = data[ichan_s,i1:i2]  - data[ichan_d,i1:i2]
				modelo_DIFF = model[ichan_s,i1:i2] - model[ichan_d,i1:i2]
				flags       = flag[ichan_s,i1:i2]  + flag[ichan_d,i1:i2]
				masca       = masc[ihorn,i1:i2]
			
				cut   = where(flags eq 0 and datos_SUM ne 0.0 and masca eq 0, ncut)
				if ncut ge 3 then begin
					sigsum[ihorn,ifreq,icorr,j]  = stddev(datos_SUM[cut]-modelo_SUM[cut])
					sigdiff[ihorn,ifreq,icorr,j] = stddev(datos_DIFF[cut]-modelo_DIFF[cut])
				endif
			endfor
		endfor
	endfor
endfor

ENDIF

; Define provisional threshold cuts
NSIGUP       = 1.7
thres_chan_up = median(sigma,dim=2)*NSIGUP
thres_sum_up  = median(sigsum,dim=4)*NSIGUP
thres_diff_up = median(sigdiff,dim=4)*NSIGUP
NSIGLOW      = 0.5
thres_chan_lo = median(sigma,dim=2)*NSIGLOW
thres_sum_lo  = median(sigsum,dim=4)*NSIGLOW
thres_diff_lo = median(sigdiff,dim=4)*NSIGLOW


; Read actual threshold cuts from file:
existe     = does_exists(ffthres)
if (existe) then begin
	print," ==> Reading threshold cuts from ",ffthres
	print,"     SKY = ",sky 
	str        = mrdfits(ffthres,1)
	NSIGUP     = str.NSIGUP
	NSIGLOW    = str.NSIGLOW
	thres_chan_up = str.thres_chan * NSIGUP
	thres_sum_up  = str.thres_sum  * NSIGUP
	thres_diff_up = str.thres_diff * NSIGUP
	thres_chan_lo = str.thres_chan * NSIGLOW
	thres_sum_lo  = str.thres_sum  * NSIGLOW
	thres_diff_lo = str.thres_diff * NSIGLOW

	; update values for sum_lo, which are changed
	thres_sum_lo  = str.thres_sum_lo
endif

; Print thresholds
print," THRESHOLDS: "
print,"Chan: ",thres_chan_up*1000.0
print,"Suma: ",thres_sum_up*1000.0
print,"Diff: ",thres_diff_up*1000.0



; Reflagging
reflag = flag

caso   = 2

; Case 1. Based on rms of each channel
IF caso eq 1 THEN BEGIN
print," --> Flagging RFI based on the rms of each channel "
for ichan=0,nchan-1 do begin
	for j=0,ntot-1 do begin
		i1    = imin[j] 
		i2    = imax[j]
		if (sigma[ichan,j] eq 0 OR sigma[ichan,j] gt thres_chan_up[ichan] OR sigma[ichan,j] lt thres_chan_lo[ichan] ) then reflag[ichan,i1:i2] = flag[ichan,i1:i2] + one*(1-bitmasc[ichan,i1:i2])
	endfor
endfor
ENDIF


; Case 2. Based on rms of the sum and difference
IF CASO eq 2 THEN BEGIN
print," --> Flagging RFI based on the rms of the difference"
for ihorn=0,3 do begin
	for ifreq=0,1 do begin
		for icorr=0,1 do begin
			
			ichan_s = mfi.ic[ihorn,ifreq,0] ; X+Y
			ichan_d = mfi.ic[ihorn,ifreq,1] ; X-Y
			if icorr eq 1 then begin
				ichan_s = mfi.ic[ihorn,ifreq,2] ; X
				ichan_d = mfi.ic[ihorn,ifreq,3] ; Y
			endif
			
			for j=0,ntot-1 do begin
				i1    = imin[j] 
				i2    = imax[j]
				if (sigsum[ihorn,ifreq,icorr,j] eq 0 OR sigsum[ihorn,ifreq,icorr,j] gt thres_sum_up[ihorn,ifreq,icorr] OR sigsum[ihorn,ifreq,icorr,j] lt thres_sum_lo[ihorn,ifreq,icorr] ) then begin
					reflag[ichan_s,i1:i2] = flag[ichan_s,i1:i2] + one*(1-bitmasc[ichan_s,i1:i2])
					reflag[ichan_d,i1:i2] = flag[ichan_d,i1:i2] + one*(1-bitmasc[ichan_d,i1:i2])
				endif	
				if (sigdiff[ihorn,ifreq,icorr,j] eq 0 OR sigdiff[ihorn,ifreq,icorr,j] gt thres_diff_up[ihorn,ifreq,icorr] OR sigdiff[ihorn,ifreq,icorr,j] lt thres_diff_lo[ihorn,ifreq,icorr] ) then begin
					reflag[ichan_s,i1:i2] = flag[ichan_s,i1:i2] + one*(1-bitmasc[ichan_s,i1:i2])
					reflag[ichan_d,i1:i2] = flag[ichan_d,i1:i2] + one*(1-bitmasc[ichan_d,i1:i2])
				endif
			endfor
		endfor
	endfor
endfor
ENDIF

; Save values
str = {root:root, tkernel_sec:tkernel_sec, NSIGUP:NSIGUP, NSIGLOW:NSIGLOW, sigsum:sigsum, sigdiff:sigdiff, sigma:sigma, sky:sky, elev:elev, periodo:periodo, domask:domask }
if (overwrite ne 0) then begin
	print," ==> Writing STATS file to ",ffout
	spawn,"rm -f "+ffout
	mwrfits, str, ffout
endif else begin
	print," ==> STATS file already exists. FILE NOT OVERWRITEN. "
endelse

; Display test
if (display) then begin
	ichan=0 											    
	data2 = data*0  										    
	data2[ichan,where(reflag[ichan,*] eq 0) ] =data[ichan,where(reflag[ichan,*] eq 0)]		    

	cfp
	plot,data[ichan,where(flag[ichan,*] eq 0)]							    
	oplot,data2[ichan,where(flag[ichan,*] eq 0)],col=2  
	oplot,model[ichan,where(flag[ichan,*] eq 0)],col=3

	cut=where(flag ne reflag)
	print,"Fraccion de diferencia = ",n_elements(cut)*1.0/(n_elements(flag)*1.0)

	stop
endif
print,"   --> SANCHO_MFI_FLAG_RFI done. "
print,"---------------------------------"

END
