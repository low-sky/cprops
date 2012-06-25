function cloudmom, xin, yin, vin, tin $
                   , mom2x = mom2x, mom2y = mom2y, mom2v = mom2v $
                   , mom0t = mom0t, targett = targett $
                   , forcelin = forcelin, noextrap = noextrap $
                   , find_covariance = covar, covariance = cov, err = err

;+
; NAME:
;  CLOUDMOM
; PURPOSE:
;   To calculate the moments of an emission distribution.
;
; CALLING SEQUENCE: 
;
;  moment = CLOUDMOM(x, y, v, t, [mom2x = mom2x, mom2x = mom2x, mom2y
;   = mom2y, mom2v = mom2v, mom0t = mom0t, targett = targett,
;   /forcelin /noextrap, /find_covariance, covariance = covariance])
;
; INPUTS:
;   X,Y,V,T -- Vectors containing the data (from VECTORIFY.pro)
;
; KEYWORD PARAMETERS: 

;   MOM2X, MOM2Y, MOM2V, MOM0T -- The value of the moments as a function of
;                  decreasing values of input T_A
;   TARGETT -- The target T_A value to extrapolate to (default = 0K)
;   /FORCELIN -- Forces linear extrapolation for the flux (default =
;                quadratic) 
;   /NOEXTRAP -- Force no extrapolation.
;   /FIND_COVARIANCE -- Calculate the extrapolated value of the covariance
;                 between x and y.
;   COVARIANCE -- Named variable to contain the covariance
; OUTPUTS:
;   MOMENT -- A stucture with the following tags:
;     RMSX, RMSY, RMSV, FLUX -- the extapolated moments.
;     ERMSX, ERMSY, ERMSV, EFLUX -- Estimated errors in the moments.
;     COVAR -- The extrapolated covariance.
;   
; MODIFICATION HISTORY:
;
;       Documented -- Fri Sep 2 16:47:39 2005, Erik Rosolowsky
;       <erosolow@asgard.cfa.harvard.edu>
;
;       Thu Dec 9 09:56:04 2004, Erik Rosolowsky <eros@cosmic>
;		Defaulted to TARGETT=0
;
; REQUIRES:
;
;  EXTRAP
;
;-

  if n_elements(targett) eq 0 then targett = 0

  if n_elements(forcelin) eq 0 then square = 1 else square = 0

; COPY THE INPUT ARRAYS SO THAT WE CAN MODIFY THEM WITHOUT CHANGING
; THE ORIGINALS
  x = xin
  y = yin
  v = vin
  t = tin
  
; SORT BY DECREASING ANTENNA TEMPERATURE AND REORDER ALL OF THE
; ARRAYS. THE BRIGHTEST PIXEL IS NOW ELEMENT 0.
  sort_t = reverse(sort(t))
  x = x[sort_t]
  y = y[sort_t]
  v = v[sort_t]
  t = t[sort_t]
  
; CALCULATE THE FIRST MOMENTS (MEAN X,Y,V POSITION). THE CUMULATIVE
; FLAG TO TOTAL RETURNS A CUMULATIVE VALUE FOR THE MOMENT. WE NOW HAVE
; ARRAYS OF MEANX,Y, AND V
  meanx = total(x*t, /cum)/total(t, /cum)
  meany = total(y*t, /cum)/total(t, /cum)
  meanv = total(v*t, /cum)/total(t, /cum)

; ZEROETH MOMENT IS EASY
  mom0t = total(t, /cum)

; SECOND MOMENTS ARE HARDER
;  mom2x = dblarr(n_elements(xin))
;  mom2y = dblarr(n_elements(xin))
;  mom2v = dblarr(n_elements(xin))
;  for i = long(0), n_elements(xin)-1 do begin
;    mom2x[i] = sqrt(total((x[0:i]-meanx[i])^2*t[0:i])/total(t[0:i]))
;    mom2y[i] = sqrt(total((y[0:i]-meany[i])^2*t[0:i])/total(t[0:i]))
;    mom2v[i] = sqrt(total((v[0:i]-meanv[i])^2*t[0:i])/total(t[0:i]))      
;  endfor

  term1x = total(t*double(x)^2., /cum)
  term2x = (total(t*double(x), /cum))^2./mom0t
  mom2x = sqrt((term1x - term2x)/mom0t)
  zeroind = where(abs(term1x - term2x) lt 1.d-10, num)
  if (num gt 0) then $
    mom2x[zeroind] = 0.0
  
  term1y = total(t*double(y)^2., /cum)
  term2y = (total(t*double(y), /cum))^2./mom0t
  mom2y = sqrt((term1y - term2y)/mom0t)
  zeroind = where(abs(term1y - term2y) lt 1.d-10, num)
  if (num gt 0) then $
    mom2y[zeroind] = 0.0  

  term1v = total(t*double(v)^2., /cum)
  term2v = (total(t*double(v), /cum))^2./mom0t
  mom2v = sqrt((term1v - term2v)/mom0t)
  zeroind = where(abs(term1v - term2v) lt 1.d-10, num)
  if (num gt 0) then $
    mom2v[zeroind] = 0.0

  if keyword_set(covar) then begin
    term1cov = total(t*double(x)*double(y), /cum)/mom0t
    term2cov = meanx*meany
    cov = term1cov-term2cov
  endif else cov = !values.d_nan

  xcen = total(x*t)/total(t)
  ycen = total(y*t)/total(t)
  vcen = total(v*t)/total(t)

  ermsx = !values.d_nan
  ermsy = !values.d_nan
  ermsv = !values.d_nan
  eflux = !values.d_nan

  if n_elements(err) eq n_elements(x) then begin
    eflux = double(sqrt(total(err^2)))
 ;    null = wt_moment(x, t, err = err)
;     ermsx = double(null.errsd)
;     null = wt_moment(y, t, err = err)
;     ermsy = double(null.errsd)
;     null = wt_moment(v, t, err = err)
;     ermsv = double(null.errsd)
    count = n_elements(t)
    l = n_elements(t)-1 
    meanx = total(x*t)/total(t)
    meanx_err = total(x^2*err^2)/mom0t[l]^2-$
      2*total(t*x)*total(x*err^2)/mom0t[l]^3+$
      total(err^2)*total(t*x)^2/mom0t[l]^4
    meanx2_err = total(x^4*err^2)/mom0t[l]^2-$
      2*total(t*x^2)*total(x^2*err^2)/mom0t[l]^3+$
      total(err^2)*total(t*x^2)^2/mom0t[l]^4
    xterm = total(x^3*err^2)/mom0t[l]^2-$
      total(x^2*err^2)*total(t*x)/mom0t[l]^3-$
      total(x*err^2)*total(t*x^2)/mom0t[l]^3+$
      total(t*x^2)*total(t*x)*total(err^2)/mom0t[l]^4

    ermsx = sqrt((meanx2_err+4*(meanx^2)*meanx_err-2*meanx*xterm)/mom2x[l]^2/4)
    
    meany = total(y*t)/total(t)
    meany_err = total(y^2*err^2)/mom0t[l]^2-$
      2*total(t*y)*total(y*err^2)/mom0t[l]^3+$
      total(err^2)*total(t*y)^2/mom0t[l]^4
    meany2_err = total(y^4*err^2)/mom0t[l]^2-$
      2*total(t*y^2)*total(y^2*err^2)/mom0t[l]^3+$
      total(err^2)*total(t*y^2)^2/mom0t[l]^4
    yterm = total(y^3*err^2)/mom0t[l]^2-$
      total(y^2*err^2)*total(t*y)/mom0t[l]^3-$
      total(y*err^2)*total(t*y^2)/mom0t[l]^3+$
      total(t*y^2)*total(t*y)*total(err^2)/mom0t[l]^4
    ermsy = sqrt((meany2_err+4*(meany^2)*meany_err-2*meany*yterm)/mom2y[l]^2/4)
    
    meanv = total(v*t)/total(t)
    meanv_err = total(v^2*err^2)/mom0t[l]^2-$
      2*total(t*v)*total(v*err^2)/mom0t[l]^3+$
      total(err^2)*total(t*v)^2/mom0t[l]^4
    meanv2_err = total(v^4*err^2)/mom0t[l]^2-$
      2*total(t*v^2)*total(v^2*err^2)/mom0t[l]^3+$
      total(err^2)*total(t*v^2)^2/mom0t[l]^4
    vterm = total(v^3*err^2)/mom0t[l]^2-$
      total(v^2*err^2)*total(t*v)/mom0t[l]^3-$
      total(v*err^2)*total(t*v^2)/mom0t[l]^3+$
      total(t*v^2)*total(t*v)*total(err^2)/mom0t[l]^4
    ermsv = sqrt((meanv2_err+4*(meanv^2)*meanv_err-2*meanv*vterm)/mom2v[l]^2/4)
  endif  


  if (keyword_set(noextrap)) then $
    return, {rmsx: double(mom2x[n_elements(mom2x)-1]) $
             , rmsy: double(mom2y[n_elements(mom2y)-1]) $
             , rmsv: double(mom2v[n_elements(mom2v)-1]) $
             , flux: double(mom0t[n_elements(mom0t)-1]) $
             , ermsx:ermsx $
             , ermsy:ermsy $
             , ermsv:ermsv $
             , eflux:eflux $
             , covar:cov, number:long(n_elements(t)) $
             , xcen:double(xcen), ycen:double(ycen) $
             , vcen:double(vcen) $
            }

; USE THE "EXTRAP" ROUTINE TO EXTRAPOLATE THE MOMENT VALUES TO OUR
; DESIRED ANTENNA TEMPERATURE CONTOUR, "targett."

  ex_mom2x = extrap(t, mom2x, targett = targett, /fast, scatter = e_mom2x, /weight)
  ex_mom2y = extrap(t, mom2y, targett = targett, /fast, scatter = e_mom2y, /weight)
  ex_mom2v = extrap(t, mom2v, targett = targett, /fast, scatter = e_mom2v, /weight)
  ex_mom0t = extrap(t, mom0t, targett = targett, /fast, square = square, $
                    scatter = e_mom0t, /weight)

;  if keyword_set(log) then begin
;    ex_mom0t = extrap(t, alog(mom0t), targett = targett, $
;                      /fast, square = square, $
;                      scatter = e_mom0t, /weight)
;    ex_mom0t = exp(ex_mom0t)
;  endif

  if keyword_set(covar) then $
    ex_cov = extrap(t, cov, targett = targett, /fast, scatter = e_cov) $
    else ex_cov = !values.d_nan
  if (ex_mom0t lt max(mom0t, /NAN)) then begin
    help, calls = calls
    if n_elements(calls) eq 2 then begin
      print, 'ERROR in EXTRAPOLATION!!! Extrapolated flux is LESS THAN measured flux. That is BAD!'
      print, 'Defaulting to linear extrapolation...'
    endif
    ex_mom0t = extrap(t, mom0t, targett = targett, /fast)
  endif
  return, {rmsx:(double(ex_mom2x))[0], rmsy: (double(ex_mom2y))[0] $
           , rmsv: (double(ex_mom2v))[0], flux: (double(ex_mom0t))[0] $
           , ermsx:(double(e_mom2x/ex_mom2x))[0] $
           , ermsy:(double(e_mom2y/ex_mom2y))[0] $
           , ermsv:(double(e_mom2v/ex_mom2v))[0] $
           , eflux:(double(e_mom0t/ex_mom0t))[0] $
           , covar:(double(ex_cov))[0], number:long(n_elements(t)) $
           , xcen:double(xcen), ycen:double(ycen) $
           , vcen:double(vcen) $
          }
  
end                             ; of cloudmom



