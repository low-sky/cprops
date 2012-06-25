pro cube_partition, im, level = level, astr = astr, $
                    vwt = vwt, ppbeam = ppbeam, $
                    kernels = kernels, terminate = terminate
;+
; NAME:
;  CUBE_PARTITION
;
; PURPOSE:
;  To decompose a cube of emission into clouds given kernels.  The
;  algorithm works by expanding down from the kernels in set steps of
;  brightness temperature and discrepancies
;  in the assignment are broken by distance to nearest kernel.
;
; CALLING SEQUENCE:
;   CUBE_PARTITION, cube, kernels = kernels, astr = astr [, level =
;   level, vwt = vwt, ppbeam = ppbeam, terminate = terminate]
;
; INPUTS:
;   CUBE -- a data cube
;   KERNELS -- The indices within the datacube that are local maxima
;   
;
; KEYWORD PARAMETERS:
;   LEVEL -- The size of the step to take downward in brightness at
;            every level (in units of the data cube).  Assumed cube is
;            in significance units and defaults to 0.5
;   PPBEAM -- pixels per beam
;   TERMINATE -- The intensity level at which to stop decomposing
;                emission.  
;
; OUTPUTS:
;   ASTR -- A structure containing the assignment information for
;           clouds in a cube.
;
;
; MODIFICATION HISTORY:
;       Documented.
;	Sat Sep  3 16:02:17 2005, Erik
;
;-

  if not keyword_set(level) then level = 2
  if n_elements(vwt) eq 0 then vwt = 1.0
  sz = size(im)

  badvals = where(im ne im, badct)
  if badct gt 0 then im[badvals] = 0.
  ind = where(im gt 0, goodct)

  if n_elements(terminate) eq 0 and goodct gt 0 $
    then terminate = min(im[ind], /nan)
  thresh = max(im, /nan)
  ind = where(im gt 0, ct)
  assignment = intarr(sz[1], sz[2], sz[3])

  thresh = max(im, /nan)
; Drop thorugh cube and dilate!
  repeat begin
    thresh = thresh-level
    mask = im gt thresh
    assignment[kernels] = indgen(n_elements(kernels))+1
    assignment = dilator3d(assignment, [-1, kernels], $
                           constraint = mask, vwt = vwt, ppbeam = ppbeam)
  endrep until (thresh lt terminate)

; Clip to region of interest
  assignment = assignment*(im gt terminate)

; Back out again.
  index = where(assignment gt 0)
  x = index mod sz[1]
  y = (index mod (sz[1]*sz[2]))/sz[1]
  v = index/(sz[1]*sz[2])
  assignment = assignment[index]-1
  assignment = recat(assignment)-1
  counts = histogram(assignment, min = 0) ;, reverse_indices = ri)
  astr = {index:index, x:x, y:y, v:v, $
          assignment:assignment, counts:counts}


  return
end
