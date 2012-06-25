function convang,  angle, ra = raflag
;+
; NAME:
;  CONVANG 
; PURPOSE:
;  Converts Sexigesimal to decimal angles and back.
;
; CALLING SEQUENCE:
;  result = CONVANG( angle [, ra = ra])
;
; INPUTS:
;  ANGLE -- Angle to be converted
;
; KEYWORD PARAMETERS:
;  RA -- The angle is RA in hours, minutes, seconds.
;
; OUTPUTS:
;  RESULT -- a scalar if the input was in sexigesimal and a vector if
;            the input was a scalar.
;
; MODIFICATION HISTORY:
;
;       Thu May 1 11:30:10 2003, Erik Rosolowsky <eros@cosmic>
;		Trapped boo-boo for negative dec.
;Documented.  Wed Nov 21 11:37:02 2001, Erik Rosolowsky <eros@cosmic>
;-



if (n_elements(angle) gt 1) then goto, dms2ddeg

if keyword_set(raflag) then begin
 angle= [floor(angle/15), floor(angle*4 mod 60), (angle*240 mod 60)]
endif else begin
 sign = angle lt 0 ? -1 : 1
 angle = abs(angle)
 angle =  [sign*floor(angle),floor(angle*60 mod 60), angle*3600 mod 60]
endelse
return, angle

dms2ddeg:
if keyword_set(raflag) then begin 
  angle = [angle, 0, 0]
  angle = angle[0]*15+angle[1]/4.+angle[2]/240.
endif else begin
  sign = total(angle lt 0) ? -1 : 1
  angle = abs(angle)
  angle = [angle, 0, 0]
  angle = sign*(angle[0]+angle[1]/60.+angle[2]/3600.)
endelse

return, angle
end



