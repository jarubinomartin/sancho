;--------------------------------------------------
;-- annular mask between theta1 et theta2, both
;-- latitudes in degrees (0 between theta1 and theta2).
function annular_mask, theta1, theta2, nside, deltat

  if theta2 gt theta1 then begin
     print, 'ERROR in annular mask: theta2 must be lower or equal to theta1'
     return, -1
  endif else begin

     npix = nside2npix(nside)
     pix = lindgen(npix)
     pix2ang_ring, nside, pix, theta, phi
     lat = (!pi/2-theta)*180./!pi

     mask = dblarr(npix)

     if theta1 eq theta2 then begin

        return, replicate(1,npix)

     endif else begin

        if theta1 eq 90. then begin
     
           if theta2 eq -90 then begin

              return, mask

           endif else begin

              region = where(lat ge theta2)
              mask(region) = 1
              if ~(deltat eq !null) then begin
              if   deltat ne   0    then begin
                 aporeg = where(lat lt theta2 and lat gt theta2-deltat)
                 mask(aporeg) = (1+cos((theta2-lat(aporeg))/deltat*!pi))/2.
              endif
              endif
              return, 1-mask

           endelse

        endif else begin

           if theta2 eq -90 then begin

              region = where(lat le theta1)
              mask(region) = 1
              if ~(deltat eq !null) then begin
              if   deltat ne   0    then begin
                 aporeg = where(lat gt theta1 and lat lt theta1+deltat)
                 mask(aporeg) = (1+cos((theta1-lat(aporeg))/deltat*!pi))/2.
              endif
              endif
              return, 1-mask

           endif else begin

              region = where(lat le theta1 and lat gt theta2)
              mask(region) = 1
              if ~(deltat eq !null) then begin
              if   deltat ne   0    then begin
                 aporeg_high = where(lat gt theta1 and lat lt theta1+deltat)
                 aporeg_low  = where(lat lt theta2 and lat gt theta2-deltat)
                 mask(aporeg_high) = (1+cos((theta1-lat(aporeg_high))/deltat*!pi))/2.
                 mask(aporeg_low)  = (1+cos((theta2-lat(aporeg_low ))/deltat*!pi))/2.
              endif
              endif
              return, 1-mask

           endelse

        endelse

     endelse

  endelse

end
