;  18/Oct/2017   JARM
;
;  This code generates the basic statistics for a set of ctod files
;  This information is used to generate the criterium for flagging RFI signals
;
;-
;

PRO sancho_genera_mfi_rfi_stat, filename, overwrite=overwrite, ctodpath=ctodpath, noflag=noflag, sky=sky, outpath=outpath

; Info
if n_params() eq 0 then begin
	print,""
	print,"  Syntax -- sancho_genera_mfi_rfi_stat, txt"
	print,""
	print,"  txt could be filename or directly the list of files"
	print,""
	return
endif

; Keywords
if not keyword_set(overwrite) then overwrite=0
if not keyword_set(ctodpath) then ctodpath="/net/nas/proyectos/quijote/ctod/"
if not keyword_set(noflag) then noflag = fltarr(nside2npix(128))
if not keyword_set(sky) then sky = 2
if not keyword_set(outpath) then outpath = "/net/nas/proyectos/quijote/flag_stats_ctod/jalberto/"

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

; Run main code
for i=0L,nlist-1L do begin
   ;ffctod = "/net/nas/proyectos/quijote/ctod/20ms/"+lista[i]+".ctod"
   ffctod = ctodpath+lista[i]+".ctod"
   print," Index and filename: ",i,nlist," ", ffctod
   sancho_mfi_flag_rfi, lista[i], path=ctodpath, overwrite=overwrite, noflag=noflag, sky=sky, outpath=outpath
endfor


end

