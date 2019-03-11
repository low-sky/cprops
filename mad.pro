function mad, x, window = window, finite = finite, dimension = dimension
;+
; NAME:
;   MAD
; PURPOSE:
;   To calculate the Median Absolute Deviataion of a set of data in
;   order to calucate the RMS of the noise.
;
; CALLING SEQUENCE:
;   sigma = MAD(X)
;
; INPUTS:
;   X -- A data array
;
; KEYWORD PARAMETERS:
;   None
;
; OUTPUTS:
;   Sigma -- The standard deviation of the data.
;
; MODIFICATION HISTORY:
;
;       Thu Sep 15 18:55:37 2016, <erosolow@siglab>
;              Added dimension keyword.
;
;       Mon Oct 4 2004, Adam Leroy <aleroy@astro>
;               Altered MAD to consider only finite values if
;               the finite keyword is on. NB: you may not always want
;               this (only matters for infinities).
;               
;
;       Tue Oct 7 15:59:16 2003, Erik Rosolowsky <eros@cosmic>
;		Added Compatibility for distributions with non-zero
;		mean (oops).
;
;       Mon Oct 6 13:26:11 2003, Erik Rosolowsky <eros@cosmic>
;		Written.
;
;
;-


  mad2gauss = 0.6744897501960817

  if n_elements(dimension) ne 0 then begin
     if n_elements(finite) gt 0 or n_elements(window) gt 0 then $
        message,'Keywords FINITE is not compatible with DIMENSION',/con
     sz = size(x)
     remaining = indgen(sz[0])+1
     remaining = remaining[where(dimension ne remaining)]
     med = median(x, dimension=dimension)
     if sz[0] eq 2 then begin
        med = rebin(med,sz[remaining[0]],sz[dimension])
        sortdim = sort([remaining[0],sz[dimension]])
        med = transpose(med,sortdim)
     endif
     if sz[0] eq 3 then begin
        med = rebin(med,sz[remaining[0]],sz[remaining[1]],sz[dimension])
        sortdim = sort([remaining[0],remaining[1],dimension])
        med = transpose(med,sortdim)
     endif
     mad = median(abs(x-med), dimension=dimension)/mad2gauss
  endif else begin
     if (n_elements(window) eq 0) then begin
        if (keyword_set(finite)) then begin
           ind = where(finite(x) eq 1) 
           mad = median(abs(x[ind]-median(x[ind])))/mad2gauss
        endif else $
           mad = median(abs(x-median(x)))/mad2gauss
     endif else begin
        mad = median(abs(x-median(x, window)),window)/mad2gauss
     endelse 
  endelse


  return, mad
end
