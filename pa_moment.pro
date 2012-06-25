function pa_moment, x, y, t, major = major, minor = minor
;+
; NAME:
;   PA_MOMENT
; PURPOSE:
;   To calculate the position angle of an emission distribution using
;   moment methods (PCA).
;
; CALLING SEQUENCE:
;   posang = PA_MOMENT(x, y, t [, major = major, minor = minor])
;
; INPUTS:
;   X, Y, T -- The X,Y and T_A components of a cloud
;
; KEYWORD PARAMETERS:
;   MAJOR, MINOR -- The major and minor axes of the distribution along
;                   position angles.  
;
; OUTPUTS:
;   The position angle CCW from the X=0 line of the coordinates used.
;
; MODIFICATION HISTORY:
;
;       Fri Jun 18 13:45:39 2004, <eros@master>
;		Written.
;
;-

  major = 0
  minor = 0
  x0 = total(t*x, /nan)/total(t, /nan)
  y0 = total(t*y, /nan)/total(t, /nan)
  wt = t
  mat = 1./(total(wt, /nan))*$
        [[total(wt*(x-x0)^2, /nan), total(wt*(x-x0)*(y-y0), /nan)], $
         [total(wt*(x-x0)*(y-y0), /nan), total(wt*(y-y0)^2, /nan)]]
  if determ(mat, /check) eq 0 then return, !values.f_nan
  evals = hqr(elmhes(mat))
  if (total(transpose(mat) eq mat) eq 4 and mat[1, 0] eq 0) then $
     evec = [[1., 0], [0., 1.]] else $
        evec = eigenvec(mat, evals, residual = res, /double)
  big = max(abs(evals), ind)
  bigvec = float(evec[*, ind])
  posang = atan(bigvec[1], bigvec[0])
  posang = posang+!pi*(posang lt 0)
  major = float(sqrt(evals[ind]))
  minor = float(sqrt(evals[1-ind]))
  posang = float(posang)
  return, posang
end
