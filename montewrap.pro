pro montewrap, filename, monteiter, props = mad_props $
               , gal = gal, show = show $
               , fl = forcelin, sigdiscont = sigdiscont $
               , nredun = nredun, clfind = clfind, eclump = eclump $
               , minpeak = minpeak, thresh = thresh, friends = friends $
               , edge = grow_thresh, nodecomp = nodecomp $
               , minarea = minarea, mask = mask, nonuniform = nonuniform $
               , fscale = fscale, deltapeak = deltapeak, minpix = minpixin $
               , assignment = newmask, noextrap = noextrap $
               , makemask = makemask, kindist = kindist $
               , bmfriends = bmfriends, specfriends = specfriends $
               , dist = dist

;+
;
; NAME:
;
; PURPOSE:
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; SOME PROCEDURES THAT AREN'T IDL/GODDARD THAT WE USE:
;
; MODIFICATION HISTORY:
;
;-
  
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; I. THE ORIGINAL FITS TO PROPS CALL
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

  fits2props, filename, props = mad_props, gal = gal, show = show $
    , fl = forcelin, sigdiscont = sigdiscont $
    , nredun = nredun, clfind = clfind, eclump = eclump $
    , minpeak = minpeak, thresh = thresh, friends = friends $
    , edge = grow_thresh, nodecomp = nodecomp $
    , minarea = minarea, mask = mask, nonuniform = nonuniform $
    , fscale = fscale, deltapeak = deltapeak, minpix = minpixin $
    , assignment = newmask, noextrap = noextrap $
    , makemask = makemask, kindist = kindist $
    , bmfriends = bmfriends, specfriends = specfriends $
    , dist = dist  

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; II. ADD NOISE AND REPEAT TO (OVER)ESTIMATE THE ERRORS THOROUGHLY
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

  for i = 0, monteiter-1 do begin
;     **** AL: I'D LIKE TO SIMPLIFY THIS CALL A BIT. THOUGHTS? ****
;     CAREFUL! ONLY PASS A MASK IF THE USER IS ALREADY DOING THAT
    if (makemask le 0) then montemask = mask else $
      print, 'Monte Carlo Iteration ', i+1, ' remaking mask in monte'
;     CALL FITS2PROPS EXACTLY AS IT WAS CALLED, *EXCEPT* DON'T FLAG
;     "MONTEITER" AND *DO* FLAG THE "ADDCORRNOISE" FLAG
    fits2props, filename, props = monteprops, gal = gal, show = show $
      , fl = forcelin, sigdiscont = sigdiscont $
      , nredun = nredun, clfind = clfind, eclump = eclump $
      , minpeak = minpeak, thresh = thresh, friends = friends $
      , edge = grow_thresh, nodecomp = nodecomp $
      , minarea = minarea, nonuniform = nonuniform $
      , fscale = fscale, deltapeak = deltapeak, minpix = minpixin $
      , noextrap = noextrap, makemask = makemask, kindist = kindist $
      , bmfriends = bmfriends, specfriends = specfriends $
      , /addcorrnoise, mask = montemask
    if n_elements(montepropra) eq 0 then begin
      montepropra = monteprops
    endif else begin
      montepropra = [montepropra, monteprops]
    endelse
  endfor

;   NOW ASSIGN EACH MONTE CARLO'D CLOUD TO THE NEAREST CLOUD IN THE
;   ORIGINAL DATASET
  closestcloudra = lonarr(n_elements(montepropra))
  disttoclosestcloud = lonarr(n_elements(montepropra))
  for i = 0, n_elements(montepropra)-1 do begin
;     CRAPPY, UNWEIGHTED, 3-D DISTANCE
    disttoorig = ((montepropra[i].xctr_pix - mad_props.xctr_pix)^2 + $
                  (montepropra[i].yctr_pix - mad_props.yctr_pix)^2 + $
                  (montepropra[i].xctr_pix - mad_props.xctr_pix)^2)^0.5
    mindist = min(disttoorig, minind, /nan)
;     ONLY COUNT CLOUDS THAT ACTUALLY TOUCH THE ORIGINAL CLOUD (MORE
;     OR LESS) TOWARDS THE UNCERTAINTY
    touch = $
      ((abs(montepropra[i].xctr_pix - mad_props[minind].xctr_pix) le $
        2.*mad_props[minind].momxpix) AND $
       (abs(montepropra[i].yctr_pix - mad_props[minind].yctr_pix) le $
        2.*mad_props[minind].momypix) AND $
       (abs(montepropra[i].vctr_pix - mad_props[minind].vctr_pix) le $
        2.*mad_props[minind].momvpix))
    if (touch) then begin
      closestcloudra[i] = minind
      disttoclosestcloud[i] = mindist
    endif else begin
      closestcloudra[i] = -1
      disttoclosestcloud[i] = mindist
    endelse
  endfor

;  THE FACTOR TO SCALE THE ERRORS BY TO ACCOUNT FOR CORRELATION IN
;  DATA DUE TO THE OVERSAMPLED BEAM
;  indfac = (sqrt(mad_props.ppbeam))[0]
;  DON'T BUY IT ANY MORE. NOT FOR THIS TEST.
   indfac = 1.0

;   OKAY, NOW FILL IN THE UNCERTAINTIES FOR EACH ORIGINAL CLOUD BASED
;   ON THE SCATTER OF THE MONTE CARLO CLOUDS
  for i = 0, n_elements(mad_props)-1 do begin
    monteind = where(closestcloudra eq i, ct)
    if (ct gt 0) then begin
      mra = montepropra[monteind]
      mad_props[i].momxpix_uc = indfac*mad(mra.momxpix)/median(mra.momxpix)
      mad_props[i].momypix_uc = indfac*mad(mra.momypix)/median(mra.momypix)
      mad_props[i].momvpix_uc = indfac*mad(mra.momvpix)/median(mra.momvpix)
      mad_props[i].momxpix_noex_uc = $
        indfac*mad(mra.momxpix_noex)/median(mra.momxpix_noex)
      mad_props[i].momypix_noex_uc = $
        indfac*mad(mra.momypix_noex)/median(mra.momypix_noex)
      mad_props[i].momvpix_noex_uc = $
        indfac*mad(mra.momvpix_noex)/median(mra.momvpix_noex)
      mad_props[i].rad_uc = indfac*mad(mra.rad_pc)/median(mra.rad_pc)
      mad_props[i].rad_nodc_uc = indfac*mad(mra.rad_nodc)/median(mra.rad_nodc)
      mad_props[i].rad_noex_uc = indfac*mad(mra.rad_noex)/median(mra.rad_noex)
      mad_props[i].rad_nodc_noex_uc = $
        indfac*mad(mra.rad_nodc_noex)/median(mra.rad_nodc_noex)
      mad_props[i].sigv_uc = indfac*mad(mra.sigv_kms)/median(mra.sigv_kms)
      mad_props[i].sigv_nodc_uc = indfac*mad(mra.sigv_nodc)/median(mra.sigv_nodc)
      mad_props[i].sigv_noex_uc = indfac*mad(mra.sigv_noex)/median(mra.sigv_noex)
      mad_props[i].sigv_nodc_noex_uc = $
        indfac*mad(mra.sigv_nodc_noex)/median(mra.sigv_nodc_noex)
      mad_props[i].flux_uc = indfac*mad(mra.flux_kkms_pc2) $
        /median(mra.flux_kkms_pc2)
      mad_props[i].flux_noex_uc = indfac*mad(mra.flux_noex) $
        /median(mra.flux_noex)
      mad_props[i].mlum_uc = indfac*mad(mra.mlum_msun)/median(mra.mlum_msun)
      mad_props[i].mvir_uc = indfac*mad(mra.mvir_msun)/median(mra.mvir_msun)
    endif else begin
      print, 'Not enough data to monte carlo cloud ', i
    endelse
  endfor

end
