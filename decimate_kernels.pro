function decimate_kernels, k_in,  cube, ALL_NEIGHBORS = all_neighbors, $
  sigma = sigma, delta = delta, merger = merger, minpix = minpix;, $
;  levels = levels

;+
; NAME:
;   DECIMATE_KERNELS
; PURPOSE:
;   To eliminate kernels based on area uniquely associated with them
;   and their contrast with the remainder of the emission.
;
; CALLING SEQUENCE:
;   new_kernels = DECIMATE_KERNELS(kernels, cube [,/ALL_NEIGHBORS =
;   all_neighbors, DELTA = delta, SIGMA = sigma]) 
;
; INPUTS:
;   KERNELS -- The potential kernels as indices in the array CUBE
;   CUBE -- The data cube in question.
;
; KEYWORD PARAMETERS:
;   ALL_NEIGHBORS -- Set this keyword to make pixels have 26 neighbors
;                    and not 6.
;
;   DELTA -- The height that a maximum must be (in units of SIGMA)
;            above a saddle point to be considered a unique kernel.
;
;   SIGMA -- The units of error in the data cube.  Defaults to the
;            median, absolute deviation.
;
; OUTPUTS:
;   NEW_KERNELS -- Kernels decimated by requiring a local maximum to
;                  be N-sigma above the valley where it merges with
;                  another cloud.  
; REQUIRES: 
;   MERGEFIND
;
; MODIFICATION HISTORY:
;
;       Documented the dramatic paring down that happened recently.
;       Fri Sep 2 16:05:11 2005, Erik Rosolowsky
;       <erosolow@asgard.cfa.harvard.edu>
;
;       Mon Dec 6 13:11:01 2004, Erik Rosolowsky <eros@cosmic>
;		Removed legacy indexing bug.
;
;       Fri Dec 3 13:32:22 2004, Erik Rosolowsky <eros@cosmic>
;		Brought to life again...
;
;-

; COPY KERNELS TO ALLOW EDITTING/AVOID UNFORTUNATE ACCIDENTS
  kernels = k_in

; CONTRAST CRITERIA
  if n_elements(delta) eq 0 then $
    delta = 2.
  area_rejects = 0
  delta_rejects = 0

; MEAN ABSOLUTE DEVIATION (CONVERTED TO SIGMA) OF DATA IN THE CUBE
  if n_elements(sigma) eq 0 then $
    sigma = mad(cube)

; Establish order of decimation: LOWEST -> HIGHEST
  kernel_value = cube[kernels]
  order = sort(kernel_value)
; Set a vector which tracks good kernels to all kernels good.
  valid_kernel = intarr(n_elements(kernels))+1
  for i = 0, n_elements(order)-1 do begin
    valid_kernel_ind = where(valid_kernel, valid_ct)
    newkernels = kernels[valid_kernel_ind]
    if valid_ct gt 1 then begin
; For this kernel, calculate the levels where it merges will all
; remaining kernels.
      merge_level = minimergefind(cube, kernels[order[i]], $
                                  newkernels, $; levels = levels, $
                                  npixels = npixels, $
                                  all_neighbors = all_neighbors)
; If there aren't enough pixels associated with the kernel it's bad. 
      if npixels lt minpix then valid_kernel[order[i]] = 0
      if npixels lt minpix then area_rejects = area_rejects+1
      
; If it's not high enough above the merge level, then it's bad.
      if (kernel_value[order[i]]-merge_level)/sigma lt delta then begin
        delta_rejects = delta_rejects+1
        valid_kernel[order[i]] = 0
      endif
    endif
  endfor

  message, 'Kernels rejected for area: '+string(area_rejects), /con
  message, 'Kernels rejected for contrast: '+string(delta_rejects), /con
  valid_kernel_ind = where(valid_kernel, valid_ct)
  newkernels = kernels[valid_kernel_ind]
  return, newkernels
end

