function bisection, xinit, fname, itmax = itmax, tol = tol,$
               _extra = extra, radius = radius, double = double
;+
; NAME:
;   BISECTION
; PURPOSE:
;   To find a root by bisecting a region around the given x-value.
;
; CALLING SEQUENCE:
;   root = BISECTION(x0, fname [, itmax=itmax, tol=tol, radius=radius,
;   double=double])
;
; INPUTS:
;   x0 -- initial guess for an x-value
;   FNAME -- function name.
; KEYWORD PARAMETERS:
;   ITMAX -- max number of iterations to search for a root
;   TOL -- tolerance for error in root value
;   RADIUS -- Initial radius of search
;   Extra keywords passed to function FNAME
; OUTPUTS:
;   ROOT -- Hopefully the evil root. (Root of all evil?)
;
; MODIFICATION HISTORY:
;       Written --
;       Wed Oct 24 10:00:07 2001, Erik Rosolowsky <eros@cosmic>
;-

if not keyword_set(itmax) then itmax = 50
if not keyword_set(tol) then tol = 1e-4
if not keyword_set(radius) then radius = 0.5

iter = 0
; First search for x values with y-values on both side of zero
x = double(xinit)
if n_elements(extra) gt 0 then begin
  x1 = x+radius
  x0 = x-radius
  y1 = call_function(fname, x1, _extra = extra)
  y0 = call_function(fname, x0, _extra = extra)
  while y1*y0 gt 0 do begin
    radius = radius*2
    x1 = x+radius
    x0 = x-radius
    y1 = call_function(fname, x1, _extra = extra)
    y0 = call_function(fname, x0, _extra = extra)
    iter = iter+1
    if iter gt itmax then begin
      message, 'Could not locate region with roots on both sides.', /con
      message, 'Restrict radius or initial guess.', /con
      return, 0
    endif
  endwhile
  toler = abs(x0-x1)
  while toler gt tol do begin
    xold = x
    x = (x0+x1)*0.5
    y = call_function(fname, x, _extra = extra)
    x0 = x0+(x-x0)*(y*y0 gt 0) 
    y0 = y0+(call_function(fname, x, _extra = extra)-y0)*(y*y0 gt 0) 
    x1 = x1+(x-x1)*(y*y1 gt 0) 
    y1 = y1+(call_function(fname, x, _extra = extra)-y1)*(y*y1 gt 0) 
    toler = abs(x0-x1)
    iter = iter+1
    if iter gt itmax then begin
      message, 'Could not locate root with number of interations.', /con
      return, 0
    endif
  endwhile
  if not keyword_set(double) then return, float(x)
  return, x
endif else begin
  x1 = x+radius
  x0 = x-radius
  y1 = call_function(fname, x1)
  y0 = call_function(fname, x0)
  while y1*y0 gt 0 do begin
    radius = radius*2
    x1 = x+radius
    x0 = x-radius
    y1 = call_function(fname, x1)
    y0 = call_function(fname, x0)
    iter = iter+1
    if iter gt itmax then begin
      message, 'Could not locate region with roots on both sides.', /con
      message, 'Restrict radius or initial guess.', /con
      return, 0
    endif
  endwhile
  toler = abs(x0-x1)
  while toler gt tol do begin
    xold = x
    x = (x0+x1)*0.5
    y = call_function(fname, x)
    x0 = x0+(x-x0)*(y*y0 gt 0) 
    y0 = y0+(call_function(fname, x)-y0)*(y*y0 gt 0) 
    x1 = x1+(x-x1)*(y*y1 gt 0) 
    y1 = y1+(call_function(fname, x)-y1)*(y*y1 gt 0) 
    toler = abs(x0-x1)
    iter = iter+1
    if iter gt itmax then begin
      message, 'Could not locate root with number of interations.', /con
      return, 0
    endif
  endwhile
  if not keyword_set(double) then return, float(x)
  return, x
endelse


  return, 0
end
