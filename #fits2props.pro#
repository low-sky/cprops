pro fits2props, filename, props = mad_props, gal = gal, show = show $
                , fl = forcelin, sigdiscont = sigdiscont $
                , nredun = nredun, clfind = clfind, eclump = eclump $
                , minpeak = minpeak, thresh = thresh, friends = friends $
                , edge = grow_thresh, nodecomp = nodecomp $
                , minarea = minarea, mask = mask, nonuniform = nonuniform $
                , fscale = fscale, deltapeak = deltapeak, minpix = minpixin $
                , assignment = newmask, noextrap = noextrap $
                , makemask = makemask, kindist = kindist $
                , bmfriends = bmfriends, specfriends = specfriends $
                , adduncorrnoise = adduncorrnoise $
                , addcorrnoise = addcorrnoise $
                , bootstrap = bootstrap, dist = dist, bclip = bclip $
                , physical = physical, xco = xco, minvchan = minvchan $
                , smwidth = smwidth, rms = rmsin, ecube = ecube $
                , near = near, far = far, r0 = r0, v0 = v0, zero2nan = zero2nan
  

;+
;
; NAME:
;   FITS2PROPS
;
; PURPOSE:
;
;   Turn a FITS cube or an IDL save file into cloud properties. This
;   is the outermost wrapper for the cloud properties package.
;
; CALLING SEQUENCE:
;
;   fits2props, filename, props = props
;
; INPUTS:
;
;   filename : the name of the input file. If it ends in '.fits' it
;   will be read in and parsed using READFITS and the header commands,
;   otherwise it will be restored and the program looks for a GAL
;   structure.
;
;   props : the main output of the program, an array of structures
;   containing cloud properties.
;
; KEYWORD PARAMETERS:
;
;
; OUTPUTS:
;
; SOME PROCEDURES THAT AREN'T IDL/GODDARD THAT WE USE:
;
; MAD; RDHD; SIGMA_CUBE; CIRCLE
;
; MODIFICATION HISTORY:
;
;-


; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; I. DEFAULTS, DEFINITIONS, AND INPUTS
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

; MAKE A FILLED CIRCLE THE CUSTOM POINT
  circle, /fill

; DOES FILE EXIST?
  if not file_test(filename, /read) then begin
    message, filename+' is not accessible.', /con
    return
  endif

; PICK BETWEEN THE 'GAL' AND '.fits' FOR CONVENTION. IF THE FILE ENDS
; IN 'fits' THEN WE WILL TREAT IT AS A FITS FILE. OTHERWISE WE RESTORE
; AN IDL SAVE FILE AND TRUST THAT AL's GALAXY INFORMATION STRUCTURE
; COMES ALONG WITH IT.
  if stregex(filename, 'fits', /bool) then begin
    data = readfits(filename, hd)
    rdhd, hd, s = h
  endif else begin
    restore, filename
    data = cube  
  endelse
  if keyword_set(zero2nan) then begin
    badind = where(data eq 0.0, ct)
    if ct gt 0 then data[badind] = !values.f_nan
  endif

; CALCULATE THE MEAN ABSOLUTE DEVIATION OF THE DATA. THIS COULD BEAR
; SOME IMPROVEMENT IN CASE OF PATHOLOGICAL DATASETS WITH EMPTY
; REGIONS.
  sigma = keyword_set(rmsin) ? rmsin : mad(data, /finite)

; ... COMPARE WITH A PRIORI NOISE IF KNOWN/CALCULATED ONLY DO THIS FOR
; THE 'GAL' CONVENTION NOT THE '.fits' CONVENTION.
  if (NOT stregex(filename, 'fits', /bool) AND keyword_set(show)) then $
    print, 'MAD: ', sigma, '; GAL noise: ', gal.noise
  
; CHECK IF A MASK ALREADY EXISTS. IF NOT, WE WILL MAKE ONE BELOW.
  if (n_elements(makemask) eq 0) then $
    makemask = (n_elements(mask) ne n_elements(data))

; FIGURE PIXELS PER BEAM
  if (n_elements(h)) then ppbeam = h.ppbeam else begin
    deltax = abs(gal.daxis[1] - gal.daxis[0])*3600. ; PIX SIZE, AS
    ppbeam = !pi*(gal.beamfwhm/2.)^2./alog(2.)/deltax^2. ; PIX/BEAM

  endelse
  beaminpix = sqrt(ppbeam/!pi*alog(2.))*2./2.354 ; SIGMA

; *******************************************************
; **** MONTE CARLO TOOLS: ADD NOISE TO THE DATA CUBE ****
; *******************************************************

; IF THE ADDNOISE KEYWORD IS SET JUST ADD AN APPROPRIATE AMOUNT OF
; NOISE IN THE DUMBEST WAY POSSIBLE
  if (keyword_set(adduncorrnoise)) then begin
    addnoise, data, sigma
  endif

; IF THE ADDCORRNOISE KEYWORD IS SET, ADD THE APPROPRIATE LEVEL OF
; *CORRELATED* NOISE (MAKE THE NOISE, CONVOLVE WITH A BEAM AND THEN
; ADD).
  if (keyword_set(addcorrnoise)) then begin
    addnoise, data, sigma, /corrnoise, beaminpix = beaminpix
  endif

;   RECALC THE NOISE (IT MAY HAVE GONE UP IF WE ADDED NOISE)
  sigma = keyword_set(rmsin) ? rmsin : mad(data, /finite)

; *******************************************************
; **** FINE TUNING PARAMETERS FOR IDENTIFYING CLOUDS ****
; *******************************************************

; Setting the /PHYSICAL flag sets the parameters of the decomposition
; to be uniform in physical rather than data units.  This _REQUIRES_ a
; well formed FITS header and a characteristic distance DIST.

  if keyword_set(physical) then begin
    message, 'Setting DEFAULTS for decomposition to uniform physical priors.', /con
    if keyword_set(gal) then dist = gal.dist
    if keyword_set(nonuniform) then begin
       if n_elements(rmsin) ne n_elements(data) then begin
          em = errmap_rob(data)
          good_ind = where(em gt 0, ctr)
          if ctr gt 0 then rms = median(em[good_ind]) else $
             rms = keyword_set(rmsin) ? rmsin : mad(data, /finite)
       endif else rms = rmsin
    endif else rms = keyword_set(rmsin) ? rmsin : mad(data, /finite)
    
    if n_elements(hd) gt 0  then begin
      bunit = sxpar(hd, 'BUNIT')
      beamfwhm = sqrt(h.bmaj*h.bmin)
      ppbeam = h.ppbeam
      if stregex(bunit, 'JY/BEAM', /bool) then begin
        rms = rms*h.jypb2k      ; Make this Kelvins.
        if n_elements(bclip) eq 0 then bclip = 2.5*h.k2jypb      
; Scale to Kelvins
      endif
      if stregex(bunit, 'K.KM/S', /bool) or $
        stregex(bunit, 'K km/s', /bool) then begin 
        rms = rms/(h.cdelt[2]/1e3)[0] ; MAKE IT KELVIN!!
        if n_elements(bclip) eq 0 then bclip = 2.5*(h.cdelt[2]/1e3) ;  SCALE to KELVIN!!
      endif
      deltav = h.cdelt[2]
    endif else begin
      if keyword_set(gal) then begin
        beamfwhm = gal.beamfwhm ; ARCSEC
        vaxis = gal.vaxis       ; KM/S
        deltav =abs(vaxis[1]-vaxis[0])*1e3
        rms = rms/(deltav/1e3)
        if n_elements(bclip) eq 0 then bclip = 2.5*deltav/1e3 
;  SCALE to KELVIN!
      endif
    endelse
    beamsz = beamfwhm/206265.*dist

; Set DELTAPEAK = 1 K convolved from 15 pc resolution, but not less than 2 sigma.
    if n_elements(deltapeak) eq 0 then begin
      deltapeak = ((1.0/rms)*(15/beamsz < 1)^2) > 2.0
      if deltapeak eq 2.0 then message, 'Warning: Noise level too high for fully accurate decomposition', /con
    endif

; Set beamfriends to 15 pc scale or 0.5 depending of which is less
    if n_elements(bmfriends) eq 0 then begin
      bmfriends = (15/beamsz)
      if bmfriends lt 0.5 then message, 'Warning: Spatial resolution too coarse for fully accurate decomposition', /con
    endif

; Set specfriends to 2 km/s or 1 
    if n_elements(specfriends) eq 0 then begin
      specfriends = (2e3/deltav) > 1.0
      if specfriends eq 1.0 then message, 'Warning: Velocity resolution too coarse for fully accurate decomposition', /con
    endif
  endif


; MINIMUM AREA OF A CLOUD TO MEASURE IN UNITS OF BEAM AREAS.
  if (n_elements(minarea) eq 0) then $
    minarea = 1.0

; IF THE MINIMUM PIXELS AREN'T SPECIFIED USE THE AREA CUT HERE TOO.
  if (n_elements(minpixin) eq 0) then $
    minpix = minarea*ppbeam else $
    minpix = minpixin*ppbeam

; MINIMUM AMPLITUDE OF A PEAK, IN "SIGMAS" FOR US TO CONSIDER THAT
; CLOUD.
  if (n_elements(minpeak) eq 0) then $
    minpeak = 4.0

; REQUIRED CONTRAST BETWEEN PEAK AND EDGE OF A CLOUD FOR US TO TRY
; MEASURING ITS PROPERTIES.
  if (n_elements(deltapeak) eq 0) then $
    deltapeak = 2.0

; REQUIRE TWO VELOCITY CHANNELS. WITHOUT THIS SOME CALCULATIONS WILL
; RETURN GARBAGE.
  if (n_elements(minvchan) eq 0) then $
    minvchan = 2

; Explore BMFRIENDS beams in every direction if BMFRIEND is set for
; initial guess at local maxima.
  if n_elements(bmfriends) gt 0 then friends = $
    ceil(sqrt(ppbeam/!pi > 1)*bmfriends/2)

  if n_elements(friends) eq 0 then friends = 1

; THRESHOLD FOR MASKING (REQUIRE ADJACENT CHANNELS ABOVE THIS TO START
; THE MASK).
  if (n_elements(thresh) eq 0) then $
    thresh = 4  

; EDGE OF THE MASKING THRESHOLD (MASK GROWS OUT TO THIS SIGNIFICANCE
; CONTOUR).
  if (n_elements(grow_thresh) eq 0) then $
    grow_thresh = 2

; The BCLIP keyword ignores structure above given level in brightness,
; setting it equal to BCLIP.  The values are returned for the property
; analysis

  if n_elements(bclip) gt 0 then begin
    clip_index = where(data gt bclip, clip_count)
    if clip_count gt 0 then orig_data_clip = data[clip_index]
  endif else clip_count = 0


; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; II. MASK AND THEN IDENTIFY INTERESTING LOCAL MAXIMA
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

  if keyword_set(rmsin) then sigma = rmsin

  if keyword_set(nonuniform) then begin
; If the noise in the data is nonuniform, first estimate the noise and
; then move data into units of significance for decomposition:
    message, 'Estimating variance at every position from data', /con
    if n_elements(ecube) eq 0 then $ 
       ecube = sigma_cube(data, width = smwidth) ; THIS IS THE "ERROR" CUBE
    if clip_count gt 0 then data[clip_index] =  $
      (bclip*(1+atan(data[clip_index]/bclip-1)))
    data = data/ecube
    sigma = 1.0
  endif else begin
    if clip_count gt 0 then data[clip_index] = $
      (bclip*(1+atan(data[clip_index]/bclip-1)))
  endelse


  if makemask gt 0 then begin
; PRETTY RUDIMENTARY MASKING, BUT SIMPLE IS GOOD ... PICK 2 CHANNELS
; AT "thresh" SIGMA OR ABOVE AND THEN GROW THE MASK OUT TO THE
; "grow_thresh" SIGMA CONTOUR LEVEL.
    print, 'Masking -- Threshold: ', thresh, '; Edge: ', grow_thresh
    mask = data gt thresh*sigma
    mask = (mask*(shift(mask, 0, 0, 1)+shift(mask, 0, 0, -1)) gt 0)
    constraint = data gt grow_thresh*sigma
    constraint = (constraint*(shift(constraint, 0, 0, 1)+$
                              shift(constraint, 0, 0, -1)) gt 0)
    mask = dilate_mask(mask, constraint = constraint)
  endif



; END OF MASKING ...
;
; ... BEGINNING OF "CLOUD" IDENTIFICATION: HERE WE PICK OUT A SET OF
; CONTIGUOUS REGIONS OF EMISSION THAT WE WILL ANALYZE. THESE REGIONS
; WILL BE SEARCHED FOR INTERESTING LOCAL MAXIMA BELOW AND THERE MAY BE
; MORE THAN ONE MAXIMA PER CLOUD.

; PARE THE MASK, REJECTING CLOUDS THAT ARE TOO SMALL OR LACK
; SIGNIFICANT PEAKS.
  if (keyword_set(show)) then begin
    origmask2d_z = total(mask, 3, /NAN) ge 1
    origmask2d_y = total(mask, 2, /NAN) ge 1
    origmask2d_x = total(mask, 1, /NAN) ge 1
  endif

; IF THE MASK IS TOTALLY EMPTY THEN WE HAVE NO CLOUDS TO
; ANALYZE. RETURN AN EMPTY STRUCTURE.
  if total(mask) eq 0 then begin
    message, 'No clouds worth our effort in this cube.', /con
    mad_props = cloudalyze(data, mask, gal = gal, hd = hd, dist = dist, $
                           near = near, far = far, r0 = r0, v0 = v0)
    return
  endif

; REJECT CLOUDS THAT DO NOT HAVE THE REQUIRED AREA, PEAK VALUE, OR
; CONTRAST
  mask = reject_region(data, mask, minarea = minarea*ppbeam, $
                       minpksig = sigma*minpeak, mindeltav = minvchan)
  
; IF THE "show" KEYWORD IS SET, OPEN THREE WINDOWS TO DISPLAY THE MASK.
  if (keyword_set(show)) then begin

    window, 0, title = 'Z-COMPRESSED'
    newmask2d_z = total(mask, 3, /NAN) ge 1
    mask2d_z = origmask2d_z + newmask2d_z
    loadct, 3
    disp, mask2d_z, min = 0, max = 2

    window, 1, title = 'Y-COMPRESSED'
    newmask2d_y = total(mask, 2, /NAN) ge 1
    mask2d_y = origmask2d_y + newmask2d_y
    loadct, 3
    disp, mask2d_y, min = 0, max = 2

    window, 2, title = 'X-COMPRESSED'
    newmask2d_x = total(mask, 1, /NAN) ge 1
    mask2d_x = origmask2d_x + newmask2d_x
    loadct, 3
    disp, mask2d_x, min = 0, max = 2

    print, 'Showing MASK along all 3 axes. Hit a key to continue with decomposition ...'
    ch = get_kbrd(1)
    wdelete, 0
    wdelete, 1
    wdelete, 2    
  endif

; AGAIN, DOUBLE CHECK FOR AN EMPTY MASK AND RETURN AN EMPTY STRUCTURE
; IF APPROPRIATE
  if total(mask) eq 0 then begin
    message, 'No clouds worth our effort in this cube.', /con
    mad_props = cloudalyze(data, mask, gal = gal, hd = hd, dist = dist, $
                           near = near, far = far, r0 = r0, v0 = v0)
    return
  endif

; IF THE USER REQUEST NO DECOMPOSTION THEN WE TREAT EACH CLOUD AS AN
; INDEPENDENT GMC. THERE'S NO NEED TO SEARCH FOR LOCAL MAXIMA, SO WE
; MEASURE THE PROPERTIES OF THE CLOUDS HERE.
  if keyword_set(nodecomp) then begin
    newmask = label_region(mask)
  endif else begin
    
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; III. NOW DEAL WITH DECOMPOSITION. IDENTIFY CLOUDS SUITABLE FOR
; DECOMPOSITION AND THEN FEED THOSE CLOUDS INTO THE DECOMPOSITION
; WRAPPER.
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

; DECIDE WHICH CLOUDS MIGHT BE DECOMPOSABLE INTO SUBCLOUDS. INPUT IS
; THE DATA AND THE MASK OF SIGNIFICANT CLOUDS, OUTPUT IS THE DECOMP
; ARRAY
    decomp_select, data, mask $
      , x = x, y = y, v = v, t = t, index = index $
      , decompose = decomp, assignment = asgn $
      ,  minpix = minpix

; DO THE DECOMPOSITION. CALLS THE DECOMPOSITION WRAPPER, WHICH TAKES
; THE SELECTED DECOMPOSITION MEASURE AND ASSIGNS EACH PIXEL TO A
; SUBCLOUD.

    decomp_wrap, x, y, v, t, asgn, decomp = decomp $
      , subcloud = subcloud, minpix = minpix $
      , sigdiscont = sigdiscont $
      , nredun = nredun, clfind = clfind $
      , sigma = sigma, kernels = kernels $
      , eclump = eclump, newkernels = newkernels $
      , nodecomp = nodecomp, fscale = fscale $
      , noextrap = noextrap $
      , friends = friends, specfriends = specfriends $
      , delta = deltapeak, ppbeam = ppbeam

; IF THE "show" KEYWORD IS TURNED ON, THEN OPEN THREE WINDOWS AND SHOW
; THE TWO-D VERSION OF THE MASK ALONG EACH OF THE THREE AXES, WITH
; REJECTED REGIONS IN (BLOOD) RED. PAUSE FOR USER TO CONTINUE ...
    if (keyword_set(show)) then begin  
      circle, /fill
      window, 0, title = 'Z-COMPRESSED'
      newmask2d_z = total(mask, 3, /NAN) ge 1
      mask2d_z = origmask2d_z + newmask2d_z
      loadct, 3
      disp, mask2d_z
      oplot, x[kernels]+0.5, y[kernels]+0.5 $
        , psym = 8, symsize = 2., thick = 5 $
        , color = getcolor('blue')
      oplot, x[newkernels]+0.5, y[newkernels]+0.5 $
        , psym = 7, symsize = 1.0, thick = 5 $
        , color = getcolor('green')

      window, 1, title = 'Y-COMPRESSED'
      newmask2d_y = total(mask, 2, /NAN) ge 1
      mask2d_y = origmask2d_y + newmask2d_y
      loadct, 3
      disp, mask2d_y
      oplot, x[kernels]+0.5, v[kernels]+0.5 $
        , psym = 8, symsize = 2., thick = 5 $
        , color = getcolor('blue')
      oplot, x[newkernels]+0.5, v[newkernels]+0.5 $
        , psym = 7, symsize = 1.0, thick = 5 $
        , color = getcolor('green')

      window, 2, title = 'X-COMPRESSED'
      newmask2d_x = total(mask, 1, /NAN) ge 1
      mask2d_x = origmask2d_x + newmask2d_x
      loadct, 3
      disp, mask2d_x
      oplot, y[kernels]+0.5, v[kernels]+0.5 $
        , psym = 8, symsize = 2., thick = 5 $
        , color = getcolor('blue')   
      oplot, y[newkernels]+0.5, v[newkernels]+0.5 $
        , psym = 7, symsize = 1.0, thick = 5 $
        , color = getcolor('green')
      
      print, 'Showing MASK and KERNELS along all 3 axes. Hit a key to continue...'
      ch = get_kbrd(1)
      wdelete, 0
      wdelete, 1
      wdelete, 2
    endif

; Clouds that are NOT decomposed into anything have a SUBCLOUD
; assignment of 0.  Emission in the watershed (shared/contested region
; between clouds) has a SUBCLOUD assignment of 1 and valid subclouds
; have an assignment of 2 or more.  Thus, selecting where the subcloud
; NE 1 gets all non-decomposed clouds and all subclouds that are not
; contested.

    interesting_pixels = where(subcloud ne 1, ct)
    if ct eq 0 then begin
      message, 'There is nothing left.  I am panicking and exiting.', /con
      return
    endif
    x = x[interesting_pixels]
    y = y[interesting_pixels]
    v = v[interesting_pixels]
    t = t[interesting_pixels]
    asgn = asgn[interesting_pixels]
    subcloud = subcloud[interesting_pixels]
    
; This gives each cloud a unique assignment value.  The RECAT function
; maps the N unique cloud assignments to 1,2,..., N.

    new_asgn = recat(asgn*(max(subcloud)+1)+subcloud)

;   EXPAND ASSIGNMENTS BACK UP INTO A CUBE
    sz = size(data)
    newmask = intarr(sz[1], sz[2], sz[3])
    newmask[x, y, v] = new_asgn
  endelse  

; Return data to real units if it's in significance units
  if keyword_set(nonuniform) then $
    data = data*ecube

; Restore clipped data as well
  if clip_count gt 0 then data[clip_index] = orig_data_clip
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; IV. FINALLY, ANALYZE THE CLOUDS USING THE EXTRAPOLATION/MEASUREMENT
; ALGORITHM.
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

; THIS GETS RETURNED IN THE props KEYWORD
  mad_props = cloudalyze(data, newmask, gal = gal, hd = hd $
                         , noextrap = noextrap $
                         , kindist = kindist, filename = filename $
                         , bootstrap = bootstrap, dist = dist $
                         , ppbeam = ppbeam, sigma = sigma, ecube = ecube $
                         , xco = xco, near = near, far = far, r0 = r0, v0 = v0)
  
  return
end
