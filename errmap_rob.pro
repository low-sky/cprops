function erf0, x, erftarg = erftarg
; ERF0: Subroutine used by ERRMAP_ROB.pro 
; to find where error function equals a given value.
  return, errorf(x)-erftarg
end

function errmap_rob, cube
;+
; NAME:
;    errmap_rob
; PURPOSE:
;    To generate a map of the errors pixel-wise in a data cube by
;    taking the RMS along a pixel and then rejecting anything over 
;    a threshold determined by the number of channels in a
;    spectrum, so that there is 25% chance of noise producing a spike
;    above the threshold.
; CALLING SEQUENCE:
;    ERRMAP_ROB, cube
;
; INPUTS:
;    CUBE - A three dimensional data cube
;
; REQUIRES:
;   ERF0.pro
; KEYWORD PARAMETERS:
;    None
;
; OUTPUTS:
;    MAP - A map of the RMS values at each pixel.
;
; MODIFICATION HISTORY:
;       Added MAD values as default for cube.
;       Tue Oct 12 16:28:20 2004, <eros@master>
;		
;       Correctly trapped some
;       bad inputs.  Thu Mar 6 10:51:47 2003, Erik Rosolowsky
;       <eros@cosmic>
;
;       Trapped some bad inputs.
;       Sat Jan 11 09:45:41 2003, <eros@master>
;       
;       Written
;       Tue Oct 23 19:19:09 2001, Erik Rosolowsky <eros@cosmic>
;
;		
;-

sz = size(cube)
if sz[0] ne 3 then begin
  message, 'No Cube not 3-D!', /con
  return, 0
endif
sig_false = bisection(3., 'erf0', erftarg = (1d0-(5d-1)/sz[3]), /double)

; Cube might be really big, so use FOR loops.

map = fltarr(sz[1], sz[2])+!values.f_nan;mad(cube)
for i = 0, sz[1]-1 do begin
  for j = 0, sz[2]-1 do begin
    spec = (cube[i, j, *])
    ind = where(spec eq spec)
    if total(ind) eq -1 then continue    
    ind = where(spec lt 0)
    if n_elements(ind) lt 10 then continue
    sigma = stdev([spec[ind], -spec[ind]])
    ind = where(spec lt sig_false*sigma)
    map[i, j] = stdev(spec[ind])
  endfor 
endfor

return, map
end
