function erf0, x, erftarg = erftarg
; ERF0: Subroutine used by ERRMAP_ROB.pro 
; to find where error function equals a given value.
  return, errorf(x)-erftarg
end

function channelnoise, data, spline = spline
;+
; NAME:
;   channelnoise
; PURPOSE:
;   To estimate the noise across a map in the third dimension using
;   the rms of the planes formed by the first two dimensions.  The
;   algorithm uses a 3.5 sigma rejection of points to calculate the
;   standard deviation.
;
; CALLING SEQUENCE:
;   noise = CHANNELNOISE(data)
;
; INPUTS:
;   DATA -- A datacube with velocity/frequency in the third dimesion.
;
; KEYWORD PARAMETERS:
;   None.
;
; OUTPUTS:
;   NOISE -- Vector with RMS of each plane.
;
; MODIFICATION HISTORY:
;       Trapped NaN values in entire planes.
;       Thu Dec 19 09:24:28 2002, <eros@master>
;
;       Promoted to shared utility.
;       Wed Dec 5 11:06:58 2001, Erik Rosolowsky <eros@cosmic>
;
;		
;
;-

  sz = size(data)
  sig_false = bisection(3., 'erf0', erftarg = (1d0-(5d-1)/(sz[1]*sz[2])), $
                        /double)

  error = fltarr(sz[3])
  for i = 0, sz[3]-1 do begin
    pts = data[*, *, i]
    ind1 = where(pts eq pts, ct)
    if ct eq 0 then begin
      error[i] = !values.f_nan
      continue
    endif
    pts = pts[ind1]

    pts = pts[where(finite(pts))]
    ind = where(pts lt 0, ct)
    if ct lt 10 then begin
      error[i] = !values.f_nan
      continue
    endif
    sigma = stdev([pts[ind], -pts[ind]])
    ind = where(pts lt sig_false*sigma)
    error[i] = stdev(pts[ind])
  endfor

  if keyword_set(spline) then begin
    xax = findgen(sz[3])
    ind = where(error eq error)
    fit = bspline_iterfit(xax[ind], error[ind], nbkpts = sz[3]/30 < 10)
    error = bspline_valu(xax, fit)
  endif 

  return, error
end
