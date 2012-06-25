function contourcloud, xin, yin, vin, tin, x0 = x0, y0 = y0, $
  v0 = v0, clev = clev_in, count = count, all_neighbors = all_neighbors
;+
; NAME:
;   CONTOURCLOUD
;
; PURPOSE:
;   Given input vectors returns a vector of lowest contour levels such
;   that the corresponding element in the input is >= the contour
;   level.
;
; CALLING SEQUENCE:
;   t_cont = CONTOURCLOUD(x, y, v, t, x0 = x0, y0 = y0, v0 = v0,
;   [clev = clev, count = count])
;
; INPUTS:
;   X,Y,V,T -- Vectors containing the data (from VECTORIFY.pro)
;   X0, Y0, V0 -- Position of the kernel
; KEYWORD PARAMETERS:
;   CLEV -- Contour levels to use;
;   COUNT -- Number of pixels >= that contour level.
;  
; OUTPUTS:
;   T_CONT -  a vector of lowest contour levels such
;   that the corresponding element in the input is >= the contour
;   level.
;
; MODIFICATION HISTORY:
;
;	Fri Mar 24 11:35:49 2006, Erik 
;       Added all_neighbors compatibility
;
;       Documented -- Fri Sep 2 16:36:25 2005, Erik Rosolowsky
;                     <erosolow@asgard.cfa.harvard.edu>
;-


; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; DEFINITIONS AND INPUT
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

; COPY ARRAYS SO WE DON'T FUCK ANYTHING UP
  x = xin
  y = yin
  v = vin
  t = tin

; SORT BY DECREASING ANTENNA TEMPERATURE AND REORDER ALL OF THE
; ARRAYS. THE BRIGHTEST PIXEL IS NOW ELEMENT 0.
  sort_t = reverse(sort(t))
  x = x[sort_t]
  y = y[sort_t]
  v = v[sort_t]
  t = t[sort_t]

; MAKE A 3-D CUBE THE SIZE OF THE CLOUD
  cubify,  x, y, v, t, cube = cube, indcube = indcube, /pad, $
           location = location
; GET THE INDEX OF THE SEED POINT
  seedind = location[where((x eq x0) AND (y eq y0) AND (v eq v0))]
; CREATE THE OUTPUT ARRAY AND MAKE IT "NaN" TO START WITH
  contours = t*!values.f_nan

; CREATE THE ARRAY OF POSSIBLE CONTOURS TO LOOP OVER AND MAKE SURE
; THAT THEY ARE IN DESCENDING ORDER
  if (n_elements(clev_in) eq 0) then $
    clev = t else $
      clev = clev_in[reverse(sort(clev_in))]

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; LOOP OVER THE CONTOUR LEVELS IN ORDER TO ASSIGN EACH CLOUD PIXEL TO
; THE RIGHT ONE.
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
;  if keyword_set(halt) then stop

  count = lonarr(n_elements(clev))
  
  for i = 0, n_elements(clev)-1 do begin
    if n_elements(clev) gt 3 then $
      counter, i, n_elements(clev), ' out of '
;   FLAG ALL VALUES ABOVE THE CURRENT CONTOUR
    mask = cube ge clev[i]
; Count number above this contour level
    count[i] = total(mask)

;   IDENTIFY THE REGIONS OF CONTIGUOUS EMISSION
    regions = label_region(mask, /ULONG, all_neighbors = all_neighbors)
;   PICK OUT THE PIXELS IN THE SAME REGION AS THE SEED POINT THAT HAVE
;   NOT YET BEEN ASSIGNED A HIGHER CONTOUR VALUE
    new_in_contour = where((regions[location] eq (regions[seedind])[0]) AND $
                           (finite(contours) eq 0) AND $
                           ((regions[seedind])[0] ne 0), num)

;   IF WE HAVE ANY NEW PIXELS, ASSIGN THEM TO THIS CONTOUR LEVEL
    if (num gt 0) then $
      contours[new_in_contour] = clev[i]

  endfor
; JUGGLE THE ORDER BACK TO MATCH THE INPUT ARRAYS
  dummy = contours
  contours[sort_t] = dummy

  return, contours
  
end                             ; OF CONTOURCLOUD


