;   4/July/2013   J.A. Rubino
;
;   Given a central coordinate (alpha,delta) on sky, and a
;   certain location in the focal plane in real degrees (xpos, ypos)
;   it derives the coordinate of the point on sky.
;   It uses the transformations for gnomonic projection, from Kovalevsky
;   (1985). See appendix B of my PhD thesis.
;
;   HISTORY:
; 
;      4/07/2013 - original version
;      1/08/2013 - new version includes the case near the poles.
;                  To be tested in more detail. 
;
;-
; 

PRO fp2sky, a, d, xpos, ypos, anew, dnew

; Info
if n_params() lt 4 then begin
   print,""
   print,"  Syntax -- fp2sky, alpha, delta, xpos, ypos, alpha1, delta1"
   print,""
   print,"     Alpha and delta are vectors. (xpos,ypos) are numbers. "
   print,"     All inputs in degrees. "
   print,""
   return
endif

; Variables
d2r = !DPI/180.0
x   = xpos*d2r
y   = ypos*d2r
n   = n_elements(a)


; Null case:
if (x eq 0.0 and y eq 0.0) then begin
   anew = a
   dnew = d
   return
endif

; Denominator
deno = ( 1.0 - y*tan(d*d2r) )
neg  = where( deno lt 0.0, nneg )


; Exact computation for delta-alpha
tanda = fltarr(n)
da    = fltarr(n)

if (x eq 0.0) then begin
   da[*]       = 0.0                   ; deg
   if (nneg ge 1) then da[neg] = 180.0 ; deg
endif else begin

   cut   = where(d ne 0.0, ncut)
   if (ncut ge 1) then begin
      tanda[cut] = (1.0/tan(d[cut]*d2r) + tan(d[cut]*d2r)) * x * sin(d[cut]*d2r) * cos(d[cut]*d2r) / $
                   ( cos(d[cut]*d2r) - y*sin(d[cut]*d2r) )
   endif

   ; Case of dec=0
   cut   = where(d eq 0.0, ncut)
   if (ncut ge 1) then tanda[cut] = x
   
   ; Compute da
   da    = ATAN(tanda)/d2r      ; deg

   ; Check if change of sign
   if (nneg ge 1) then da[neg] = da[neg] + 180.0

endelse


; Exact computation for delta-delta
tandd = fltarr(n)

if (x eq 0.0) then begin
   tandd = cos(da*d2r) * ( y + tan(d*d2r) )  / ( 1.0 - y*tan(d*d2r) )
endif else begin
   cut   = where(d ne 0.0, ncut)
   if (ncut ge 1) then tandd[cut] = ( sin(da[cut]*d2r)/x - cos(d[cut]*d2r)*cos(da[cut]*d2r) ) / sin(d[cut]*d2r)
   
   cut   = where(d eq 0.0, ncut)
   if (ncut ge 1) then tandd[cut] = y * cos(da[cut]*d2r)
endelse
dd   = ATAN(tandd)/d2r          ; deg


; Problems in the transformation. This occurs near the poles, when (d+ypos) ge 90.0
npole = where(dd gt 90.0, nnpole)
if nnpole ge 1 then dd = 180.0 - dd

; OLD VERSION: I simply was stopping the program
;cut = where( (d+ypos) ge 90.0, ncut)
;if (ncut ge 1) then stop, "FP2SKY - problem in the transformation."




; New coordinates
anew = (a + da) mod 360.0
dnew = dd


END

