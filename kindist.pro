function kindist, l, b, v, near = near, far = far, rc = rc, $
                  r0 = r0, v0 = v0, rgal = rgal
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

  null = 1/(1+v/(v0*sin(l*!dtor)*cos(b*!dtor)))

  fardist = r0*(cos(l*!dtor)+sqrt((cos(l*!dtor))^2-$
                                  (1-null^2)))/(cos(b*!dtor))

  neardist = r0*(cos(l*!dtor)-sqrt((cos(l*!dtor))^2-$
                                   (1-null^2)))/(cos(b*!dtor))
  rgal = null*r0

  if (not keyword_set(near)) then dist = fardist else dist = neardist


  return, abs(dist)
end
