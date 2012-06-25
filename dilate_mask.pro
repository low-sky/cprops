function dilate_mask, mask_in, constraint = constraint, $
  all_neighbors = all_neighbors
;+
; NAME:
;   DILATE_MASK
;
; PURPOSE:
;   Starting with a core mask, the mask is expanded to contain all
;   regions that have some overlap with the original mask.
;
; CALLING SEQUENCE:
;   new_mask = DILATE_MASK(mask, constraint = constraint)
;
; INPUTS:
;   MASK -- A binary mask for a data cube.
;   CONSTRAINT -- A binary mask (presumably a superset of the original
;                 mask)
;
; KEYWORD PARAMETERS:
;   ALL_NEIGHBORS -- Keyword passed to LABEL_REGIONS().
;
; OUTPUTS:
;   NEW_MASK -- A dilated mask.
;
; MODIFICATION HISTORY:
;
;       Fri Sep 2 14:07:52 2005, Erik Rosolowsky
;       <erosolow@asgard.cfa.harvard.edu>
;
;		Written
;
;-

; Generate a mask with a 2 everywhere there's an overlap and a 1
; otherwise.  

  mask = (mask_in*constraint)+constraint
  mask_out = mask
  mask_out[*] = 0b
  regions = label_region(mask gt 0, all_neighbors = all_neighbors, /ulong)
  if total(regions) eq 0 then return, !values.f_nan
  h = histogram(regions, binsize = 1, min = 1, rev = ri)
  for k = 0L, n_elements(h)-1 do begin
    inds = ri[ri[k]:(ri[k+1])-1]
    if total(mask[inds] eq 2) gt 1 then mask_out[inds] = 1b
  endfor

  return, mask_out
end
