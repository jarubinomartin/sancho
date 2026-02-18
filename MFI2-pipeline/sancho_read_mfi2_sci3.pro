;   18/Mar/2024   JARM
;
;   Reads SCI3 data for MFI2
;
;-
;

PRO sancho_read_mfi2_sci3,ff, data, verbose=verbose, path=path

if not keyword_set(verbose) then verbose=0
if not keyword_set(path) then path = '/net/calp-nas/proyectos/quijote3/data/2024-03/'

if n_params() eq 0 then begin
   ff   = 'MOON-24-03-15-15-29-03-0000.sci3'
   ff   = 'MOON-24-03-15-17-59-04-0005.sci3'
;ff   = 'MOON-24-03-15-15-29-03-0002.sci3'
;ff   = 'TEST0-24-03-12-16-30-43-0000.sci3'
endif

print,' Reading MFI2 data in SCI3 mode: ',ff

nchan = 32
nsamp = 1000L
nlen  = 60L*5L ; duration of each file: 5min

; data types
uint16 = UINT(0)
uint8  = BYTE(0)
uint32 = ULONG(0)
doble  = DOUBLE(0)
int16  = FIX(0)
int32  = LONG(0)

data    = fltarr(nchan,nsamp*nlen)
dataAux = LONARR(nchan,nsamp)


close,1
openr,1,path+ff,/swap_endian
idx=0L
for k=1L,nlen do begin
   for i=0,nchan-1 do begin
      campo1 = uint16
      campo2 = uint16
      campo3 = uint16
      readu,1,campo1
      readu,1,campo2
      readu,1,campo3
      
      campo4 = uint8
      campo5 = uint8
      campo6 = uint8
      campo7 = uint8
      readu,1,campo4
      readu,1,campo5
      readu,1,campo6
      readu,1,campo7
      
      Tstamp = uint32
      Time1  = uint32
      readu,1,Tstamp
      readu,1,Time1
      
      SID        = uint16
      AqErr      = uint8
      channel_ID = uint8
      readu,1,SID
      readu,1,AqErr
      readu,1,channel_ID
      
      cal_flag  = uint8
      cal_sw    = uint8
      mod_mode  = uint8
      repuesto1 = uint8
      readu,1,cal_flag
      readu,1,cal_sw
      readu,1,mod_mode
      readu,1,repuesto1
      
      pcounter1 = uint16
      samprate  = uint16
      readu,1,pcounter1
      readu,1,samprate
      
      trigger   = doble
      gap       = int16
      Sc_factor = doble
      readu,1,trigger
      readu,1,gap
      readu,1,Sc_factor
      
      nsamples = uint32
      readu,1,nsamples
      
      dataAux1 = LONARR(nsamp)
      readu,1,dataAux1
      dataAux[i,*] = dataAux1[*]
      
; Print info	
      if verbose then begin
         print,'** Info ** '
         print,'campo1,2,3   = ',campo1,campo2,campo3
         print,'campo4,5,6,7 = ',campo4,campo5,campo6,campo6
         print,'Tstamp, Time = ',Tstamp,Time1
         print,'SID, AqErr   = ',SID,AqErr
         print,'channel_ID   = ',channel_ID
         print,'cal_flag, cal_sw, mod_mode, repuesto1 = ',cal_flag,cal_sw,mod_mode,repuesto1
         print,'pcounter1,samprate    = ',pcounter1,samprate
         print,'trigger,gap,Sc_factor = ',trigger,gap,Sc_factor
         print,'nsamples     = ',nsamples
      endif	
   endfor
   data[*,idx:idx+999L] = double(dataAux)*Sc_factor
   idx += 1000L
endfor
close,1

END

