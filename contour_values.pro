function contour_values, cube, nlevels = nlevels, $
  all = all, uniform = uniform, $
  minimum = minimum, maximum = maximum
;+
; NAME:
;
;   CONTOUR_VALUES
;
; PURPOSE:
;
;   Accepts cube or vector and returns suggested contours.  Defaults
;   to 250 levels or the number of unique data, whichever is less.
;   Highest contour is first in the vector.
;   
; CALLING SEQUENCE:
;
;   contours = CONTOUR_VALUES ( data [ nlevels = nlevels, /ALL, /UNIFORM ])
;
; INPUTS:
;
;   DATA -- A datacube or vector.
;
; KEYWORD PARAMETERS:
;
;   NLEVELS -- maximum number of levels to include.  Default = 250.
;   ALL -- Return all unique contour levels in the data.
;   UNIFORM -- Controls sampling of contour level vectors.  The
;              default is to sample contour levels weighted by the
;              distribution of the data.  Setting the UNIFORM keyword
;              forces the levels to be quasi-uniformly distributed
;              between the minimum and the maximum (though always
;              taking values present in the data). 
;   MINIMUM, MAXIMUM -- Requires contours to be >= MINIMUM and
;                       <= MAXIMUM
; OUTPUTS:
;
;   CONTOURS -- Reverse sorted vector of contour levels.
;
; MODIFICATION HISTORY:
;
;       Wed Dec 8 09:17:39 2004, Erik Rosolowsky <eros@cosmic>
;		Adapted from about 634 different ways of writing this.
;
;-

; DEFAULTS
  if n_elements(nlevels) eq 0 then nlevels = 250
  if n_elements(minimum) eq 0 then minimum = min(cube, /nan)
  if n_elements(maximum) eq 0 then maximum = max(cube, /nan)

; FIND UNIQUE SORTED LIST OF VALUES IN DATACUBE
  scube = cube[sort(cube)]
  goodind = where(finite(scube), ct)
  if ct eq 0 then return, !values.f_nan
  scube = scube[goodind] > minimum < maximum
  levels = (scube[uniq(scube)])
  if keyword_set(all) then return, reverse(levels)
; IF THERE ARE MORE VALUES THAN NLEVELS THEN DECIMATE ARRAY BY EITHER
  if n_elements(levels) gt nlevels then begin
; SELECTING PSEUDO-UNIFORM LEVELS WITHIN THE DATA SET
    if keyword_set(uniform) then begin
      pseudo_levs = findgen(nlevels)*$
                    (levels[0]-levels[n_elements(levels)-1])/nlevels+$
                    levels[n_elements(levels)-1]
      level_index = interpol(findgen(n_elements(levels)), levels, pseudo_levs)
      levels = levels[level_index]
      levels = levels[sort(levels)]
      levels = levels[uniq(levels)]
      levels = reverse(levels)
    endif else begin
; OR JUST EVENLY SAMPLE THE VALUES THAT WE HAVE.
      level_index = findgen(nlevels)/nlevels*n_elements(levels) 
      levels = levels[level_index]
      levels = reverse(levels)
    endelse
  endif

; GIVE IT TO THEM
  return, levels
end
; STOP SHOUTING AT THE SCREEN.
