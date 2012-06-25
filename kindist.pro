function kindist, l, b, v, near = near, far = far, rc = rc, $
                  r0 = r0, v0 = v0, rgal = rgal, clemens = clemens
;+
; NAME:
;   KINDIST 
; PURPOSE:
;   To return the distance to an object given l,b,v
;
; CALLING SEQUENCE:
;   dist = KINDIST (L, B, V)
;
; INPUTS:
;   
;
; KEYWORD PARAMETERS:
;
;
; OUTPUTS:
;   DIST -- the kinematic distance in units of R0 (defaults to pc).
;
; MODIFICATION HISTORY:
;
;-


  if n_elements(r0) eq 0 then r0 = 8.5d3
  if n_elements(v0) eq 0 then v0 = 2.2d2


  if keyword_set(clemens) then begin
    






    rvec = (dindgen(190)/100+0.1)*r0
    legcoeff = [663.53792d0, -1652.3148d0, 1881.9831d0, -1455.6453d0, $
                808.84150d0, -315.02388d0, 77.873785d, -9.2156079d0]
    omega = r0/rvec
    for i = 0, n_elements(legcoeff)-1 do omega = omega+$
       legendre(rvec/2/r0, i)*legcoeff[i]
    omega_measure = (1+v/v0/sin(l*!dtor)/cos(b*!dtor))
    rgal = interpol(rvec, omega, omega_measure, /lsquad)
    null = rgal/r0
  endif else begin
    null = 1/(1+v/(v0*sin(l*!dtor)*cos(b*!dtor)))
  endelse

;  The > 0 traps things near the tangent point and sets them to the
;  tangent distance.  So quietly.  Perhaps this should pitch a flag?
  radical = sqrt(((cos(l*!dtor))^2-(1-null^2)) > 0)
  
  fardist = r0*(cos(l*!dtor)+radical)/(cos(b*!dtor))
  
  neardist = r0*(cos(l*!dtor)-radical)/(cos(b*!dtor))
  rgal = null*r0
  ind = where(abs(l-180) lt 90, ct)
  if ct gt 0 then neardist[ind] = !values.f_nan

  if (not keyword_set(near)) then dist = fardist else dist = neardist


  return, abs(dist)
end
