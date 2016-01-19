function mergefind, cube, kernels, nlevels = nlevels $
                    , terminate = terminate, everypoint = everypoint $
                    , all_neighbors = all_neighbors, uniform = uniform $
                    , levels = levels, area = area
;+
; NAME:
;   MERGEFIND
; PURPOSE:
;   To isolate the contour level at which two maxima share contours.
;
; CALLING SEQUENCE:
;   merge_matrix = MERGEFIND(cube, index_maxima [,/ALL_NEIGHBORS ,
;   /EVERYPOINT, NLEVELS=nlevels])
;
; INPUTS:
;   CUBE -- A data cube
;   INDEX_MAXIMA -- The indices of the maxima in the data cube.
;
; KEYWORD PARAMETERS:
;   NLEVELS -- The number of levels to divide the range of the
;              datacube into.
;   EVERYPOINT -- Use every point in the data cube to achieve maximum
;                 accuracy in the contour levels.
;   ALL_NEIGHBORS -- Set this keyword to make pixels have 26 neighbors
;                    and not 6.
;
; OUTPUTS:
;   MERGE_MATRIX -- For N kernels, an N X N matrix with the j,k
;                   element being the contour level where the jth
;                   kernel merges with the kth kernel. 
;                
; MODIFICATION HISTORY:
;
;       Wed Dec 8 09:50:54 2004, Erik Rosolowsky <eros@cosmic>
;		Added CONTOUR_VALUES() call to get the same contour
;		values we always want.
;
;       Thu Dec 2 17:23:31 2004, Erik Rosolowsky <eros@cosmic>
;		Adapted for television.
;
;-

; CHECK THAT WE HAVE SOME KERNELS TO TRY MERGING.
  kernel_ct = n_elements(kernels)
  if kernel_ct eq 0 then $
    return, !values.f_nan

; TERMINATE IS THE MINIMUM LEVEL TO LOOK FOR SHARED CONTOURS AT. IT
; DEFAULTS TO THE SMALLEST DATA VALUE IN THE CUBE, AT WHICH EVERYTHING
; WILL MERGE.
  
  if (n_elements(levels) eq 0) then $
    levels = contour_values(cube, nlevels = nlevels, minimum = terminate, $
                            all = everypoint, uniform = uniform)


; MAKE THE MERGER MATRIX, WHICH CONTAINS THE DATA VALUES ON THE DIAGONAL
; AND SHARED CONTOUR INDICES BETWEEN KERNEL ELEMENTS.
  merger = fltarr(kernel_ct+1, kernel_ct+1)+!values.f_nan
;  merger[0, 1:*] = cube[kernels]
;  merger[1:*, 0] = cube[kernels]
  merger[indgen(kernel_ct)+1, indgen(kernel_ct)+1] = cube[kernels]

  area = fltarr(kernel_ct+1, kernel_ct+1)+!values.f_nan

  for k = 0, n_elements(levels)-1 do begin

    thresh = levels[k]
    m = cube ge thresh
    asgn = label_region(m, all_neighbors = all_neighbors, /ulong)

    for q = 1, max(asgn) do begin
      inds = where(asgn eq q, pixct)
; Find kernels in each distinct region
      kernel_int = intersection(inds, kernels, ct)
      if ct eq 1 then begin
        kernel_index = where(kernels eq kernel_int[0])+1
        area[kernel_index, kernel_index] = $
          (area[kernel_index, kernel_index] > pixct)        
      endif
; If there are more than 1 kernel then ...
      if ct gt 1 then begin
        kernel_index = lonarr(n_elements(kernel_int))
        for lame = 0, n_elements(kernel_int)-1 do $
          kernel_index[lame] = where(kernels eq kernel_int[lame])
        xmatrix = kernel_index#(intarr(n_elements(kernel_int))+1)+1
        ymatrix = transpose(xmatrix)
        merger[xmatrix, ymatrix] = $
          (merger[xmatrix, ymatrix] > thresh)        
        area[xmatrix, ymatrix] = $
          (float(pixct) < area[xmatrix, ymatrix]) 
      endif
    endfor
    if total(merger ne merger) le 1 then return,  merger[1:*, 1:*]
  endfor
  area = area[1:*, 1:*]
  return, merger[1:*, 1:*]
end
