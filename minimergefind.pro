function minimergefind, minicube, primary_kernel, kernels $
  , terminate = terminate, everypoint = everypoint $
  , all_neighbors = all_neighbors, uniform = uniform $
  , levels = levels, npixels = npixels

;+
; NAME:
;    MINIMERGEFIND
;
; PURPOSE:
;    Calculates the largest contour level that contains only 1 kernel
;    and the number of pixels above that contour level.
;
; CALLING SEQUENCE:
;    merge_vector = MINIMERGEFIND(minicube, primary_kernel, kernels)
;
; INPUTS:
;    MINICUBE -- A data cube
;    PRIMARY_KERNEL -- Index in MINICUBE of the kernel under consideration.
;    KERNELS -- Indices of all kernels in data cube including the
;               PRIMARY_KERNEL 
;
; KEYWORD PARAMETERS:
;    LEVELS -- Contour levels to consider.
;    ALL_NEIGHBORS -- each pixel has 26 vs. 6 neighbors
;    TERMINATE -- Minimum contour level
;
; OUTPUTS:
;    MERGE_VECTOR -- A vector corresponding to KERNELS with the level
;                    where the PRIMARY_KERNEL merges with all other
;                    kernels. 
;
; MODIFICATION HISTORY:
;
;       Fri Jul 8 14:41:24 2005, Erik Rosolowsky
;       <erosolow@transit.cfa.harvard.edu>
;		Written.
;
;-
  
; Initialize npixels
  npixels = 0

  kernel_value = minicube[primary_kernel]

  if (n_elements(levels) eq 0) then $
    levels = contour_values(minicube, nlevels = nlevels, $
                            minimum = terminate, $
                            all = everypoint, uniform = uniform, $
                            maximum = kernel_value)

; First check if there's only 1 cloud in the kernel.  
  m = minicube ge min(levels)
  asgn = label_region(m, all_neighbors = all_neighbors, /ulong)
  primary_asgn = asgn[primary_kernel]
  inds = where(asgn eq primary_asgn, npixels)
  kernel_int = intersection(inds, kernels, ct)
; If there's only 1, then stop workin' so hard.
  if ct eq 1 then return, min(levels)

; Pick where to start...
  start_k = min(where(levels le kernel_value))

  for k = start_k, n_elements(levels)-1 do begin
    thresh = levels[k]
    m = minicube ge thresh
    asgn = label_region(m, all_neighbors = all_neighbors, /ulong)
    primary_asgn = asgn[primary_kernel]
    inds = where(asgn eq primary_asgn, npixels)
    kernel_int = intersection(inds, kernels, ct)
    if (ct gt 1) then begin
; If we find two+ kernels, step back up a level and return.
      if k gt 0 then begin
        k = k-1
        thresh = levels[k]
        m = minicube ge thresh 
        asgn = label_region(m, all_neighbors = all_neighbors, /ulong)
        primary_asgn = asgn[primary_kernel]
        inds = where(asgn eq primary_asgn, npixels)
        return, thresh
      endif else begin
; However, if there are 2 kernels above the highest level, then search
; for value where kernel is distinct.  
        sublevels = minicube[where(minicube gt levels[0])]
        sublevels = sublevels[reverse(sort(sublevels))]
        start_kk = min(where(sublevels le kernel_value, counter))
        if counter gt 0 then begin
          for kk = start_kk, n_elements(sublevels)-1 do begin

            thresh = sublevels[kk]
            m = minicube ge thresh
            asgn = label_region(m, all_neighbors = all_neighbors, /ulong)
            primary_asgn = asgn[primary_kernel]
            inds = where(asgn eq primary_asgn, npixels)
            kernel_int = intersection(inds, kernels, ct)
            if (ct gt 1) then begin
; If we find two+ kernels, step back up a level and return.
              kk = kk-1
              thresh = sublevels[kk]
              m = minicube ge thresh 
              asgn = label_region(m, all_neighbors = all_neighbors, /ulong)
              primary_asgn = asgn[primary_kernel]
              inds = where(asgn eq primary_asgn, npixels)
              return, thresh
            endif
          endfor
        endif 
      endelse
    endif
  endfor
; If you get here, something is wrong.
  npixels = !values.f_nan
  return, !values.f_nan
end



