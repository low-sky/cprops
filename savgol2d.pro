;+
; NAME:
;
;      SAVGOL2D()
;
; PURPOSE:
;      Generate the 2D SavGol Smoothing kernels
;
; CALLING SEQUENCE:
;       filter = SAVGOL2D()
;
; KEYWORD PARAMETERS:
;       NX, NY -- Size of the Filter in X and Y.  The filter only
;                 supports square, centered filters (default =101)
;       DEGREE -- The degree of the polynomial (default = 7)
;
; OUTPUTS:
;       A 2D Savitsky - Golay convolution filter kernel
; SIDE EFFECTS:
;       Hair Loss
;
; MODIFICATION HISTORY:
;       Written and documented
;       Fri Apr 24 13:43:30 2009, Erik Rosolowsky <erosolo@A302357>
;-

function savgol2d, nx = nx, ny = ny, degree = degree, order = order


  if n_elements(nx) eq 0 then nx = 101
  if n_elements(ny) eq 0 then ny = 101
  if n_elements(degree) eq 0 then degree = 7

; ORDER > 0 doesn't work... yet
  if n_elements(order) eq 0 then order = 0

  x = (findgen(nx)-nx/2)#replicate(1, ny)
  y = transpose(x)
  x = reform(x, nx*ny)
  y = reform(y, nx*ny)

  n = ((degree+1)*(degree+2))/2
  info = findgen(n, n_elements(x)) 
  
  xpow = indgen(degree+1)#replicate(1, degree+1)
  ypow = transpose(xpow)
  deg = xpow+ypow
  ind = where(deg le degree)
  xpv = xpow[ind]
  ypv = ypow[ind]


; Because nothing says "I thought this through" like nested for loops
  for i = 0, n-1 do begin
      for k = 0, n_elements(x)-1 do begin
        info[i, k] = x[k]^xpv[i]*y[k]^ypv[i]
      endfor
  endfor


  a = transpose(invert(info#transpose(info))#info)
  filter = fltarr(n_elements(x))
  for i = 0, n_elements(x)-1 do begin
    target = fltarr(n_elements(x))
    target[i] = 1
    filter[i] = (a##target)[0]
  endfor
  filter = reform(filter, nx, ny)

  return, filter
end
