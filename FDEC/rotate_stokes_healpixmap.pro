;   8/11/2017   JARM
;
;   Routine to transform a healpix map in GALACTIC coordinates, from galactic to equatorial parangles, rotating 
;   the Stokes parameters (Q,U). The /eq2gal option makes the inverse operation. 
;


PRO rotate_stokes_healpixmap, Q_in, U_in, Q_out, U_out, gal2eq=gal2eq, eq2gal=eq2gal, nested=nested, silent=silent

; Keywords
if not keyword_set(gal2eq) then gal2eq = 0
if not keyword_set(eq2gal) then eq2gal = 0
if not keyword_set(nested) then nested = 0
if not keyword_set(silent) then silent = 0

; Info
if n_params() lt 2 then begin
   print,""
   print,"  Syntax -- rotate_stokes_healpixmap, Q_in, U_in, Q_out, U_out, /gal2eq [, /eq2gal, /nested]"
   print,""
   print,""
   return
endif


d2r = !DPI/180.0

; Ordering
ordering = "RING"
if (nested) then ordering = "NESTED"

; Npix and nside. Expecting map_in = [npix, 3]
npix  = n_elements(Q_in[*,0])
nside = npix2nside(npix)

; Direction of the transformation
inco   = ""
outco  = ""
status = 0
txt    = "No transformation applied. Use /gal2eq or /eq2gal."
euler  = transpose(rot_gal2eq())
if (gal2eq) then begin
   inco   = "G" 
   outco  = "Q"
   status = 1
   txt    = "Galactic --> Equatorial."
endif
if (eq2gal) then begin
   inco   = "Q" 
   outco  = "G"
   status = 1
   txt    = "Equatorial --> Galactic. "
endif

; Output map
Q_out = Q_in
U_out = U_in

; Transformation
if (status) then begin

   ipix = lindgen(npix)
   if ordering eq "RING" then pix2vec_ring, nside, ipix, vector
   if ordering ne "RING" then pix2vec_nest, nside, ipix, vector
   stokes     = [ [Q_in], [U_in] ]                                               ; Q,U Stokes parameters
   vector_out = rotate_coord(vector,euler_matrix=euler,stokes_parameters=stokes) ; GAL2EQ
   if (eq2gal) then begin
      stokes     = [ [Q_in], [U_in] ]                                                               ; Q,U Stokes parameters
      vector_out2 = rotate_coord(vector_out,euler_matrix=invert(euler),stokes_parameters=stokes)    ; EQ2GAL
   endif

   Q_out[*]   = stokes[*,0]
   U_out[*]   = stokes[*,1]
   
endif

if silent eq 0 then print," CODE ended succesfully. "+txt

END
