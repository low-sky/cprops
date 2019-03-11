function sigma_cube, data, width = width, emap = emap, $
                     spline = spline, savgol=savgol, $
                     twod = twod
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
;   TWOD  -- Set this to ignore channel variations
;   WIDTH -- The width of smoothing to apply to the data cube.
;            Generally the number of pixels in a beam.  Too much
;            smoothing may artificially lower the rms.
;
; OUTPUTS:
;   RMS_CUBE -- A cube of same dimension as DATA containing the RMS at
;               each pixel.
;
; MODIFICATION HISTORY:
;
;   Thu Jan 24 18:36:32 2019, <erosolow@rick>
;        Borrowed SAVGOL2D functionality from BOLOCAT
;  
;		 Written Thu Dec 19
;09:11:51 2002, <eros@master> -

  sz = size(data)
  if n_elements(width) eq 0 then width = 3
  if n_elements(emap) eq 0 then emap = errmap_rob(data)
;  if width gt 0 then emap = smooth(emap, width, /edge_trun, /nan)
  if width gt 0 then emap = median(emap, width)
  if n_elements(twod) eq 0 then twod = 0b 

  if n_elements(savgol) gt 0 then begin
     filter = savgol2d(nx=savgol, ny=savgol, degree=3)
     bad = emap ne emap
     emaporig = emap
     emap = convol(emap, filter, /edge_trun, invalid = 0, /nan)
     relt = floor(savgol/2 + 1)
     elt = shift(dist(2*relt+1, 2*relt+1), relt, relt) le relt
     bad = 1b-(erode(1b-bad, elt))
     badidx = where(bad, ct)
     if ct gt 0 then emap[badidx] = emaporig[badidx]
  endif

  if 1b-keyword_set(twod) then begin 
     escale = channelnoise(data, spline = spline)
     escale = escale/mean(escale[where(escale eq escale)])
  endif else escale = fltarr(sz[3]) + 1
  ecube = fltarr(sz[1], sz[2], sz[3])
  for ii = 0, sz[3]-1 do ecube[*, *, ii] = emap*escale[ii]


  return, ecube
end







