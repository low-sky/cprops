function are_kernels_disc, xin, yin, vin, tin, lmaxin $
  , shared_contour = shared_contour $
  , levels = levels, sigdiscont = sigdiscont $
  , nredun = nredun, fscale = fscale 
;+
; NAME:
;   ARE_KERNELS_DISC
; PURPOSE:
;   To check if two kernels are discontinuous by comparing the moments
;   of each kernel above and below a merger.
;
; CALLING SEQUENCE:
;   flag = ARE_KERNELS_DISC(x, y, v, t, kernels [, shared_contour =
;   shared_contour, levels = levels, sigdiscont = sigdiscont, nredun =
;   nredun, fscale = fscale, minpix = minpix])
;
; INPUTS:
;   X,Y,V,T -- Vectors containing the data (from VECTORIFY.pro)
;   KERNELS -- Indices of the two kernels (in X, Y, V, T)
;   SHARED_CONTOUR -- T value where kernels merge.
;   
; KEYWORD PARAMETERS:
;   LEVELS -- The levels with which to contour the data.
;   SIGDISCONT -- The fractional change in the moments to register a
;                 discontinuity, a scalar or a 4-elt. vector
;                 corresponding to the jumps in [sigma_x, sigma_y,
;                 sigma_v and int_flux]
;   NREDUN -- The number of moments that need to jump to count the
;             kernels as being discontinuous.
;   FSCALE -- The scale factor to multiply the fractional change in
;             integrated flux by to compare with the same value for
;             the second moments.
;
; OUTPUTS:
;   FLAG -- A binary flag: 1 if the kernels can be regarded as
;           distinct, 0 otherwise.
;
; MODIFICATION HISTORY:
;
;       Documented -- Fri Sep 2 16:24:59 2005, Erik Rosolowsky
;       <erosolow@asgard.cfa.harvard.edu>
;-
  
; WRITE SOME ERROR CHECKING (ONLY 2 LMAXIN, ETC.) -- ALSO, DEBUG!

; EXTRAPOLATE TO ZERO KELVIN
  targett = 0.0

; DEFINE THE LEVEL OF A SIGNIFICANT DISCONTINUITY

; Set scale of Flux relative to 2nd moments.
  if n_elements(fscale) eq 0 then $
    fscale = 2.0

; MINIMUM SIGNIFICANCE FOR A LOG. DERIV JUMP
  if (n_elements(sigdiscont) lt 2) then begin
    if n_elements(sigdiscont) eq 1 then $
    sigdiscont = sigdiscont*[1, 1, 1, fscale] else $
      sigdiscont = [1, 1, 1, fscale]
  endif

  if (n_elements(sigdiscont2) eq 0) then begin
    sigdiscont2 = sigdiscont*2.
  endif

  
; NUMBER OF PROPERTIES THAT MUST DISPLAY JUMP
  if (n_elements(nredun) eq 0) then $
    nredun = 2                    

  if (n_elements(nredun2) eq 0) then $
    nredun2 = 1

; CREATE A SMALL STEP IN INTENSITY AND DEFINE CONTOURS JUST
; ABOVE/BELOW THE SHARED CONTOUR ... NOT CURRENTLY USED.
  if (n_elements(levels) eq 0) then begin
    eps = abs(max(tin[lmaxin], /nan) - shared_contour)/100. 
    above = shared_contour + eps
    below = shared_contour - eps   
  endif else begin
    dummy = min(abs(levels - shared_contour), thislevel, /nan)
    above = levels[(thislevel - 1) > 0]
    below = levels[(thislevel + 1) < (n_elements(levels)-1)]
  endelse
  
; INSTEAD, USE THE CONTOUR LEVELS AT AND JUST ABOVE THE SHARED CONTOUR
; TO LOOK FOR A JUMP.
; Test in case level is 1
  testlevs = levels[((thislevel-1) > 0):(thislevel > 1)]

; THE X,Y,V COORDINATES OF THE TWO LOCAL MAXES
  x1 = xin[lmaxin[0]] & y1 = yin[lmaxin[0]] & v1 = vin[lmaxin[0]]
  x2 = xin[lmaxin[1]] & y2 = yin[lmaxin[1]] & v2 = vin[lmaxin[1]]  

; GET THE MOMENTS OF THE CLOUD FOR EACH LOCAL MAX JUST ABOVE AND JUST
; BELOW THE SHARED CONTOUR.
  momra_1 =  contour_prop(xin, yin, vin, tin, lmaxin[0], $
                          levels = testlevs, /noextrap)
  momra_2 =  contour_prop(xin, yin, vin, tin, lmaxin[1], $
                          levels = testlevs, /noextrap)

; CALCULATE LOGARITHMIC d ln x / d ln y DERIVATIVES OF THE FOUR
; PROPERTIES
  dx1 = (momra_1[1].rmsx-momra_1[0].rmsx)/momra_1[0].rmsx
  dx2 = (momra_2[1].rmsx-momra_2[0].rmsx)/momra_2[0].rmsx
  dy1 = (momra_1[1].rmsy-momra_1[0].rmsy)/momra_1[0].rmsy
  dy2 = (momra_2[1].rmsy-momra_2[0].rmsy)/momra_2[0].rmsy
  dv1 = (momra_1[1].rmsv-momra_1[0].rmsv)/momra_1[0].rmsv
  dv2 = (momra_2[1].rmsv-momra_2[0].rmsv)/momra_2[0].rmsv
  df1 = (momra_1[1].flux-momra_1[0].flux)/momra_1[0].flux
  df2 = (momra_2[1].flux-momra_2[0].flux)/momra_2[0].flux

  d1 = [dx1, dy1, dv1, df1]
  d2 = [dx2, dy2, dv2, df2]
; IDENTIFY DISCONTINUITIES BIG JUMPS IN THE DERIVATIVE
  is_lmax_1_disc = (total(d1 gt sigdiscont) ge nredun) OR $
    (total(d1 gt sigdiscont2) ge nredun2)
  is_lmax_2_disc = (total(d2 gt sigdiscont) ge nredun) OR $
    (total(d2 gt sigdiscont2) ge nredun2)

; RETURN 1 IF *EITHER* SHOW DISCONTINUITY, 0 OTHERWISE
  return, (is_lmax_1_disc OR is_lmax_2_disc)
  
end                             ; OF ARE_KERNELS_DISC


