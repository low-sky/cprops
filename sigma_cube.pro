function sigma_cube, data, width = width, emap = emap, spline = spline
;+
; NAME:
;   SIGMA_CUBE
; PURPOSE:
;   Generates a cube containing the estimated errors at the position
;   of each pixel.
;
; CALLING SEQUENCE:
;   rms_cube = SIGMA_CUBE(data [, width=width])
;
; INPUTS:
;   DATA -- a data cube
;
; KEYWORD PARAMETERS:
;   WIDTH -- The width of smoothing to apply to the data cube.
;            Generally the number of pixels in a beam.  Too much
;            smoothing may artificially lower the rms.
;
; OUTPUTS:
;   RMS_CUBE -- A cube of same dimension as DATA containing the RMS at
;               each pixel.
;
; MODIFICATION HISTORY:
;       Written
;       Thu Dec 19 09:11:51 2002, <eros@master>
;-

  sz = size(data)
  if n_elements(width) eq 0 then width = 3
  if n_elements(emap) eq 0 then emap = errmap_rob(data)

;  if width gt 0 then emap = smooth(emap, width, /edge_trun, /nan)
  if width gt 0 then emap = median(emap, width)
  escale = channelnoise(data, spline = spline)
  escale = escale/mean(escale[where(escale eq escale)])
  ecube = fltarr(sz[1], sz[2], sz[3])
  for ii = 0, sz[3]-1 do ecube[*, *, ii] = emap*escale[ii]

  return, ecube
end







