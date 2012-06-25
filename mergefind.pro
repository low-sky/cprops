function mergefind, cube, kernels, nlevels = nlevels $
                    , terminate = terminate, everypoint = everypoint $
                    , all_neighbors = all_neighbors, uniform = uniform $
                    , levels = levels, augmented_levels=lvs

; ATTEMPT TO CREATE AN EFFICIENT ALGORITHM TO POPULATE THE MERGER
; MATRIX WITH INFINITE RESOLUTION


; CHECK THAT WE HAVE SOME KERNELS TO TRY MERGING.
  kernel_ct = n_elements(kernels)
  if kernel_ct eq 0 then $
    return, !values.f_nan


; INITIALIZE THE RESULTING DATA PRODUCTS
  merger = fltarr(kernel_ct, kernel_ct)+!values.f_nan
;  merger[0, 1:*] = cube[kernels]
;  merger[1:*, 0] = cube[kernels]
  merger[indgen(kernel_ct), indgen(kernel_ct)] = cube[kernels]
  area = fltarr(kernel_ct, kernel_ct)+!values.f_nan

  kval = cube[kernels]
  maxvalue = max(kval, /nan)
  minvalue = min(kval, /nan)

;  merger_upper =  fltarr(kernel_ct, kernel_ct)+maxvalue
  merger_lower = fltarr(kernel_ct, kernel_ct)+!values.f_nan

  if n_elements(levels) eq 0 then begin
    if n_elements(nlevels) eq 0 then  nlevels = (n_elements(x)/50 > 250) < 500
    levels = contour_values(cube, nlevels = nlevels, minimum = terminate, $
                            all = everypoint, uniform = uniform)

  endif
  lvs = levels[sort(levels)]
  mld = abs(median(lvs-shift(lvs,1)))

  nuniqold = 0
  doneflag = 0b
  z = 0L
  repeat begin 
; Mask and label a data cube for a given contour value
    testvalue = lvs[z]
    l = label_region(cube ge testvalue, all_neighbors = all_neighbors, /ulong)
; Figure out what the kernels are assigned to
    asgns = l[kernels]
; Figure out which assignments show up more than once
    h = histogram(asgns, min = 1, max = max(l), binsize = 1, reverse = ri)
    multiples = where(h gt 1, n_multi)
    nuniq = total(h ge 1)
    count = intarr(kernel_ct)
; Some data cubes really can't be split so this requires a minimum
; resolution of 1e-3 times the median level difference.
    if (nuniq-nuniqold) gt 1 and z gt 0 and ((lvs[z]-lvs[z-1 > 0]) ge $
      1e-3*mld) then begin
; If we get a non-binary split then
      lvs = [lvs[0:(z-1)], (lvs[z-1]+lvs[z])*0.5, lvs[z:*]]
    endif else begin
; If we have a duplicate then ....
      if n_multi gt 0 then begin
        for k = 0, n_multi-1 do begin
; Figure out which kernels are sharing an assignment
          multikern = ri[ri[multiples[k]]:(ri[multiples[k]+1]-1)]
          multikern = multikern#(intarr(h[multiples[k]])+1)        
; Write this into the merger matrix.
          merger_lower[multikern, transpose(multikern)] = testvalue        
        endfor  
      endif
; When a pair of kernels stops merging it won't write into the merger matrix.
      nuniqold = nuniq
      z = z+1
    endelse
    if z eq n_elements(lvs) then doneflag = 1b
  endrep until (doneflag eq 1)

 
  merger_lower[indgen(kernel_ct), indgen(kernel_ct)] = cube[kernels]
; If the top hasn't been resolved, do it.
  ind = where(merger_lower eq max(lvs), ct)
  if ct gt 0 and ct mod 2 eq 0 then begin
    k1 = ind mod n_elements(kernels)
    k2 = ind / n_elements(kernels)
    kern1 = kernels[k1]
    kern2 = kernels[k2]
    for i = 0, n_elements(ind)-1 do begin
      if kern1[i] eq kern2[i] then continue
      merger_lower[k1[i], k2[i]] = $
         minimergefind($
         cube, kern1[i], [kern1[i], kern2[i]], $
         all_neighbors = all_neighbors)
    endfor
 endif
  merger_lower = merger_lower > transpose(merger_lower)
  return, merger_lower
end
