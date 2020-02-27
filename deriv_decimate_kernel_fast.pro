function deriv_decimate_kernel_fast, cube, merger, lmaxin, levels = levels $
  , sigdiscont = sigdiscont, nredun = nredun $
  , fscale = fscale
;+
; NAME:
;   DERIV_DECIMATE_KERNEL
; PURPOSE:
;   To decimate a kernel set based on the jumps in the moments of the
;   emission before and after merging.
;
; CALLING SEQUENCE:
;   new_kernels = DERIV_DECIMATE_KERNEL( cube, merger_matrix, kernels 
;                 [ levels = levels, sigdistcont = sigdistcont, $
;                   nredun = nredun, fscale = fscale, noextrap = noextrap])
;
; INPUTS:
;   CUBE -- A data cube
;   MERGER_MATRIX -- A matrix generated by the MERGEFIND program that
;                    tracks the levels where various kernels merge.
;   KERNELS -- The indices of local maxima within CUBE
;
; KEYWORD PARAMETERS:
;   LEVELS -- The levels with which to contour the data.
;   SIGDISCONT -- The fractional change in the moments to register a
;                 discontinuity, a scalar or a 4-elt. vector
;                 corresponding to the jumps in [sigma_x, sigma_y,
;                 sigma_v and int_flux]
;   NREDUN -- The number of moments that need to jump to count the
;             kernels as being discontinuous.
;   FSCALE -- The scale factor to multiply the fractional change in
;             integrated flux by to compare with the same value for
;             the second moments.
;
; OUTPUTS:
;
;
; MODIFICATION HISTORY:
;
;       Tue Jun 18 21:05:40 2013, erosolo <erosolo@>
;
;		Added ability to actually skip derivative decimation for
;		performance in large data cubes.  
;
;       Documented -- Fri Sep 2 16:20:12 2005, Erik Rosolowsky
;                     <erosolow@asgard.cfa.harvard.edu>
;-


  if n_elements(sigdiscont) eq 0 then sigdiscont = 0.5

  if total(sigdiscont) eq 0 then begin
     message,'Skipping Derivative Decimation',/con
     return, lmaxin

  endif

  vectorify, cube, mask = (cube eq cube) $
    , x = x, y = y, v = v, t = t, id = id $
    , sz = sz, indvec = indvec

  lmaxes = lonarr(n_elements(lmaxin))
  for k = 0, n_elements(lmaxin)-1 do $
    lmaxes[k] = where(lmaxin[k] eq indvec)

  kernel_value = t[lmaxes]
  order = reverse(sort(kernel_value))
  valid_kernel = intarr(n_elements(order)) + 1 
  checked = intarr(n_elements(order)) 
  merger[indgen(n_elements(order)), $
         indgen(n_elements(order))] = !values.f_nan 
  
  for i = 0, n_elements(order)-1 do begin
     if valid_kernel[order[i]] eq 0 then continue
     if checked[order[i]] eq 1 then continue
     highest_merge = max(merger[order[i],*], other_kernel, /nan)
     kernels = [lmaxes[order[i]], lmaxes[other_kernel]]
     keep = are_kernels_disc(x, y, v, t, kernels $
                             , shared_contour = highest_merge $
                             , levels = levels $
                             , sigdiscont = sigdiscont, nredun = nredun $
                             , fscale = fscale)  
     if keep eq 0 then begin
        valid_kernel[other_kernel] = 0
        merger[other_kernel, *] = !values.f_nan
        merger[*, other_kernel] = !values.f_nan
     endif else begin
        checked[order[i]] = 1
        checked[other_kernel] = 1
     endelse
  endfor
  
  newlmax = lmaxin[where(valid_kernel eq 1)]
  kernel_ct = strcompress(string(n_elements(newlmax)))
  message, /con,  'Number of Kernels after derivative decimation: '+kernel_ct
  return, newlmax

end                             ; OF DERIV_DECIMATE_KERNELS


;;   xgrid = (intarr(n_elements(lmaxes))+1) ## indgen(n_elements(lmaxes))
;;   ygrid = indgen(n_elements(lmaxes)) ## (intarr(n_elements(lmaxes))+1)
;;   mask = xgrid gt ygrid  
;;   uselmax = lonarr(n_elements(lmaxes)) + 1
  
;;   ind_of_interest = where(mask eq 1, nindofint)

;;   while (nindofint gt 0) do begin
;;     highest_merge = max(merger*mask, /NAN, maxind)  
;;     kernel1_ind = xgrid[maxind]
;;     kernel2_ind = ygrid[maxind]

;;     kernels = [lmaxes[kernel1_ind], lmaxes[kernel2_ind]]

;; ;   CHECK IF THE TWO KERNELS ARE DISCONTINUOUS IN PROPERTIES AT THE
;; ;   MERGER. IF SO, WE WILL KEEP BOTH.
;;     keep = are_kernels_disc(x, y, v, t, kernels $
;;                             , shared_contour = highest_merge $
;;                             , levels = levels $
;;                             , sigdiscont = sigdiscont, nredun = nredun $
;;                             , fscale = fscale)  

;; ;   IF NOT, WE KEEP ONLY THE HIGHER OF THE TWO KERNELS
;;     if (keep eq 0) then begin
;;       if (merger[kernel1_ind, kernel1_ind] gt merger[kernel2_ind, kernel2_ind]) then begin
;;         uselmax[kernel2_ind] = 0
;;         mask[kernel2_ind, *] = 0
;;         mask[*, kernel2_ind] = 0
;;       endif else begin
;;         uselmax[kernel1_ind] = 0
;;         mask[kernel1_ind, *] = 0
;;         mask[*, kernel1_ind] = 0
;;       endelse
;;     endif else begin
;;       mask[kernel1_ind, kernel2_ind] = 0
;;       mask[kernel2_ind, kernel1_ind] = 0
;;     endelse
;;     ind_of_interest = where(mask eq 1, nindofint)
;;   endwhile

; RETURN THE SURVIVING LOCAL MAXIMA  

;;   newlmax = lmaxin[where(uselmax eq 1)]
;;   kernel_ct = strcompress(string(n_elements(newlmax)))
;;   message, /con,  'Number of Kernels after derivative decimation: '+kernel_ct

;;   return, newlmax

;; end                             ; OF DERIV_DECIMATE_KERNELS
