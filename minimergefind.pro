function minimergefind, minicube, primary_kernel, kernels $
               , all_neighbors = all_neighbors $
               , npixels = npixels, tolerance = tolerance, _extra = ex

; For the decimation, try a very quick run through the data to return
; the lowest contour containing only that pixel 

; Require all blank pixels to carry NaNs 

  if n_elements(tolerance) eq 0 then tolerance = 1e-5


; Determine intial clouds.  
  clouds = label_region(minicube eq minicube, all_neighbors = all_neighbors, /ulong)
  primary_asgn = clouds[primary_kernel]
  asgns = clouds[kernels]

  maxvalue = minicube[primary_kernel]
  minvalue =  min(minicube[where(clouds eq primary_asgn, ct)], /nan)
  mld = (maxvalue-minvalue)
  samecloud = kernels[where(asgns eq primary_asgn, ct)]

  cloudmask = (clouds eq primary_asgn)
; If there's only 1 kernel in the cloud return the minimum value in
; the cloud -- Do we want to do this?

  if ct eq 1 then begin
    npixels = total(cloudmask)
    return, minvalue
  endif
; Begin refinement

  blankedcube = minicube*cloudmask
  lastsolo = minvalue
  repeat begin
    testvalue = sqrt(minvalue*maxvalue)
    if minvalue lt 0 or testvalue ne testvalue then testvalue = (minvalue+maxvalue)*0.5
    l = label_region(blankedcube ge testvalue, all_neighbors = all_neighbors, /ulong)
    primary_asgn = l[primary_kernel]
    asgns = l[samecloud]
    sameclump = where(asgns eq primary_asgn, ct)
    if ct eq 1 then begin
      maxvalue = testvalue 
      lastsolo = testvalue
    endif else begin
      minvalue = testvalue
    endelse
  endrep until abs(maxvalue-minvalue) lt (tolerance*mld) > 1e-6
; This finds a value not necessarily in the data cube.  We need to
; find the lowest contour _IN THE CLUMP_ that lassos only the clump.  Alas

  l = label_region(blankedcube ge lastsolo, all_neighbors = all_neighbors, /ulong)
  value = min(blankedcube[where(l eq l[primary_kernel], npixels)],/nan)
  return, value
;  return, !values.f_nan
end



