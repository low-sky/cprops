function extrap, xin, yin, targett = targett, fast = fast, $
                 square = square, scatter = scatter, weight = weight, $
                 coeffs = coeffs

;+
; NAME:
;  EXTRAP
; PURPOSE:
;  To extrapolate a function T, Y(T) to a set value of T
;
; CALLING SEQUENCE:
;   value = EXTRAP(t, y, [targett = targett, /fast, /square, scatter =
;   scatter, /weight])
;
; INPUTS:
;   T, Y -- The indpendent and dependent variables for extrapolation. 
;   TARGETT -- The target T value for extrapolation
;
; KEYWORD PARAMETERS:
;   /FAST -- Subsample the data to speed up performance
;   /SQUARE -- Do a quadratic extrapolation
;   SCATTER -- The M.A.D. of the extrapolated moment values.
;   /WEIGHT -- Weight the data according to contour level.
;   COEFFS -- Set to a named variable to received the coefficients of
;             the fit.
;
; OUTPUTS:
;   VALUE -- The extrapolated value.
;
; MODIFICATION HISTORY:
;
;       Documented -- Fri Sep 2 16:47:24 2005, Erik Rosolowsky
;                     <erosolow@asgard.cfa.harvard.edu>
;
;-

  scatter = !values.f_nan

; THE TARGET FOR THE EXTRAPOLATION
  if (n_elements(targett) gt 0) then begin
    if ((targett gt min(xin)) AND (targett lt max(xin))) then begin
      dummy = min(abs(xin - targett), mind)
      return, yin[mind]
    endif
  endif

; CAN'T EXTRAPOLATE WITH LESS THAN 5 DATA POINTS
  if (n_elements(xin) lt 5) then begin
    return, !values.f_nan
  endif

; REVERSE THE (SORTED) ARRAYS BEFORE EXTRAPOLATING
  xuse = reverse(xin)
  yuse = reverse(yin)

; CULL OUT BAD VALUES
  goodind = where(finite(yuse), ct)
  if ct gt 0 then begin
    yuse = yuse[goodind]
    xuse = xuse[goodind]
  endif
; TO DO THINGS FAST, WE DECIMATE THE ARRAY DOWN TO 250 ELEMENTS

 if (keyword_set(fast)) and (n_elements(xuse) gt 250) then begin
    useindex = findgen(250)*n_elements(xuse)/250. 
    xuse = xuse[useindex]
    yuse = yuse[useindex]
  endif
  
; KEEP FITTING, INCLUDING PROGRESSIVELY MORE AND MORE OF THE CLOUD FOR
; THE FIT. 
  if (NOT keyword_set(weight)) then begin

    coeffs = dblarr(2, n_elements(xuse))*!values.f_nan
    if (keyword_set(square)) then $
      coeffs = dblarr(3, n_elements(xuse))*!values.f_nan
    
    for i = long(3), n_elements(xuse)-1 do begin
      xfit = xuse[0:i]
      yfit = yuse[0:i]    
      n = n_elements(xfit)

      if (NOT keyword_set(square)) then begin
        M = [[n, total(xfit)], $
             [total(xfit), total(xfit^2)]]
        covar = invert(M)
        coeffs[*, i] = reform(covar##transpose([total(yfit), total(xfit*yfit)]))
      endif else begin
        M = [[n, total(xfit), total(xfit^2)], $
             [total(xfit), total(xfit^2), total(xfit^3)], $
             [total(xfit^2), total(xfit^3), total(xfit^4)]]
        covar = invert(M)
        coeffs[*, i] = reform(covar##transpose([total(yfit), total(xfit*yfit), total(xfit^2*yfit)]))
      endelse    
    endfor

;   DO THE EXTRAPOLATION
    extrap_value = median(coeffs[0, *] + coeffs[1, *]*targett)
    scatter =  mad(coeffs[0, *] + coeffs[1, *]*targett)
    if (keyword_set(square)) then begin
      extrap_value = median(coeffs[0, *] + coeffs[1, *]*targett + $
                            coeffs[2, *]*targett^2)  
      scatter = mad(coeffs[0, *] + coeffs[1, *]*targett + $
                    coeffs[2, *]*targett^2)  

    endif
  endif else begin

    coeffs = dblarr(2)*!values.f_nan
    if (keyword_set(square)) then $
      coeffs = dblarr(3)*!values.f_nan

    xfit = xuse
    yfit = yuse



    wfit = (reverse(findgen(n_elements(yfit))+1.0)); WEIGHT BY CONTOUR NUMBER (BIG = MORE WEIGHT)
    
    if (NOT keyword_set(square)) then begin
      M = [[total(wfit), total(xfit*wfit)], $
           [total(xfit*wfit), total(xfit^2*wfit)]]
      covar = invert(M)
      coeffs = reform(covar##transpose([total(yfit*wfit), total(xfit*yfit*wfit)]))
      yf = coeffs[0]+coeffs[1]*xfit
      chisq = total((yfit-yf)^2)/(n_elements(xfit)-2)
    endif else begin
      M = [[total(wfit), total(xfit*wfit), total(xfit^2*wfit)], $
           [total(xfit*wfit), total(xfit^2*wfit), total(xfit^3*wfit)], $
           [total(xfit^2*wfit), total(xfit^3*wfit), total(xfit^4*wfit)]]
      covar = invert(M)
      coeffs = $
        reform(covar##transpose([total(yfit*wfit), total(xfit*yfit*wfit) $
                                 , total(xfit^2*yfit*wfit)]))
      yf = coeffs[0]+coeffs[1]*xfit+coeffs[2]*xfit^2
      chisq = total((yfit-yf)^2)/(n_elements(xfit)-2)
    endelse

    extrap_value = coeffs[0, *] + coeffs[1, *]*targett
    scatter = sqrt(chisq*covar[1, 1]*targett^2+chisq*covar[0, 0])/extrap_value

    if (keyword_set(square)) then begin
      extrap_value = extrap_value + coeffs[2, *]*targett^2
      scatter = sqrt(chisq*covar[0, 0]+chisq*covar[1, 1]*targett^2+$
                     covar[2, 2]*targett^4)/extrap_value
    endif
  endelse

  return, extrap_value
  
end                             ; OF EXTRAP
