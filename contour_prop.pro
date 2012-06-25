function contour_prop, xin, yin, vin, tin, kernel, levels = levels_in $
                       , noextrap = noextrap, all_neighbors = all_neighbors, $
                       remove_min = remove_min, forcelin = forcelin, $
                       rotate = rotate, err = err, pavec = pavec
;+
; NAME:
;   CONTOUR_PROP
;
; PURPOSE:
;   To calculate the moments of emission around a fixed local maximum.
;
; CALLING SEQUENCE:
;   moments = CONTOUR_PROP(X,Y,V,T, KERNELS [levels = levels, noextrap
;   = noextrap])
;
; INPUTS:
;   X,Y,V,T -- Vectors containing the data (from VECTORIFY.pro)
;   KERNELS -- Indices of the two kernels (in X, Y, V, T)
;
; KEYWORD PARAMETERS:
;   LEVELS -- The levels with which to contour the data.
;   NOEXTRAP -- Do not extrapolate moments to their values at 0K
;   PAVEC -- Vector of derived rotation angle at all contour levels
;
; OUTPUTS:
;   MOMENTS -- An array of moments with the same number of elements as
;              KERNEL with tags definied by  CLOUDMOM.pro
;
; MODIFICATION HISTORY:
;
;	Fri Mar 24 11:36:22 2006, Erik 
;       Added all_neighbors keyword.
;
;       Documented -- Fri Sep 2 16:29:02 2005, Erik Rosolowsky
;                     <erosolow@asgard.cfa.harvard.edu>
;		
;
;-


  if n_elements(levels_in) eq 0 then levels = contour_values(tin) else $
    levels = levels_in

  nlmax = n_elements(lmaxin)
  x0 = xin[kernel] & y0 = yin[kernel] & v0 = vin[kernel]  

;   TO AVOID EXCESS WORK, IF WE HAVE A MINIMUM CONTOUR LEVEL, DROP
;   DATA BELOW THAT AT THE BEGINNING OF THE LOOP.

  contours = contourcloud(xin, yin, vin, tin, x0 = x0, y0 = y0, $
                          v0 = v0, clev = levels, all_neighbors = all_neighbors)
  use_contours = levels
  pavec = fltarr(n_elements(use_contours))+!values.f_nan 
;   MEASURE THE PROPERTIES AT EACH CONTOUR, THIS CAN BE QUITE TIME
;   INTENSIVE. THE OUTPUT IS AN ARRAY OF MOMENTS CALLED "momra." DO
;   THIS BY CALLING THE "cloudmom" ROUTINE TO CALCULATE THE CUMULATIVE
;   MOMENTS AND THEN EXTRAPOLATE THEM TO THE TARGET CONTOUR, "targett"
  mom = {rmsx:!values.d_nan, rmsy:!values.d_nan, $
         rmsv:!values.d_nan, flux:!values.d_nan, $
         ermsx:!values.d_nan, ermsy:!values.d_nan, $
         ermsv:!values.d_nan, eflux:!values.d_nan, $
         covar:!values.d_nan, number:0L, $
         xcen:!values.d_nan, ycen:!values.d_nan, $
         vcen:!values.d_nan}

  momra = replicate(mom, n_elements(levels))

  for j = 0, n_elements(use_contours) - 1 do begin
    useind = where(contours ge use_contours[j], num)
    if num eq 0 then continue
    xuse = xin[useind]
    yuse = yin[useind]
    vuse = vin[useind]
    tuse = tin[useind]      
    if n_elements(err) gt 0 then euse = err[useind]
    if keyword_set(rotate) then begin
      pa = pa_moment(xuse, yuse, (tuse>0))
      xrot = xuse*cos(pa)+yuse*sin(pa)
      yrot = -xuse*sin(pa)+yuse*cos(pa)    
      xuse = xrot
      yuse = yrot
      pavec[j] = pa
    endif
;   CALL THE "cloudmom" ROUTINE TO CALCULATE THE CUMULATIVE MOMENTS
;   AND THEN EXTRAPOLATE THEM TO THE TARGET CONTOUR, "targett"
    mom = cloudmom(xuse, yuse, vuse, tuse-min(tuse)*$
                   keyword_set(remove_min), mom0t = mom0t $
                   , targett = targett, noextrap = noextrap, $
                   forcelin = forcelin, err = euse)
    momra[j] = mom
  endfor
;  stop
  return, momra
end
