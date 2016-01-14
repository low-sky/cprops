pro deconvolve_gauss $
   , meas_maj = meas_maj $      ; measured major axis
   , meas_min = meas_min $      ; measured minor axis
   , meas_pa = meas_pa $        ; measured position angle
   , beam_maj = beam_maj $      ; beam major axis
   , beam_min = beam_min $      ; beam minor axis
   , beam_pa = beam_pa $        ; beam position angle
   , src_maj = src_maj $        ; source major axis
   , src_min = src_min $        ; source minor axis
   , src_pa = src_pa $          ; source position angle
   , worked = worked $          ; returns 0 if not possible to convolve
   , point = point $            ; returns 1 if close to a point source
   , verbose = verbose          ; set output level

;  WARNING!!! Currently this appears to have issues returning a full
;  range of position angle values. This needs more investigation.

;  ADAPTED FROM gaupar.for in MIRIAD via K. Sandstrom

;
;  Determine the parameters of a gaussian deconvolved with another
;  gaussian.
;
;  Input:
;    bmaj1,bmin1	Major and minor FWHM of the source..
;    bpa1		Position angle of 1st gaussian, in degrees.
;    bmaj2,bmin2	Major and minor FWHM of gaussian to deconvolve with.
;    bpa2		Position angle of 2nd gaussian, in degrees.
;  Output:
;    bmaj,bmin		Major and minor axes of resultant gaussian.
;    bpa		Position angle of the result, in radians.
;    fac		Always 1 (for future use ...).
;    ifail		Success status: 0   All OK.
;					1   Result is pretty close to a
;					    point source.
;					2   Illegal result.

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; DEFAULTS AND DEFINITIONS
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

; FAIL IF MAJOR AXIS OF MEASUREMENT LACKING
  if n_elements(meas_maj) eq 0 then begin
     if keyword_set(verbose) then $
        message, "No measurement supplied. Failing.", /info
     worked = 0B
     return
  endif

; FAIL IF MAJOR AXIS OF BEAM LACKING
  if n_elements(beam_maj) eq 0 then begin
     if keyword_set(verbose) then $
        message, "No beam supplied. Failing.", /info
     worked = 0B
     return
  endif

; NO MINOR AXIS MEASURED FOR MEASUREMENT
  if n_elements(meas_min) eq 0 then begin
     if keyword_set(verbose) then $
        message, "Minor axis not supplied. Assuming round measurement.", /info
     meas_min = meas_maj
  endif

; NO MINOR AXIS BEAM FOR BEAM
  if n_elements(beam_min) eq 0 then begin
     if keyword_set(verbose) then $
        message, "Minor axis not supplied. Assuming round beam.", /info
     beam_min = beam_maj
  endif

; NO POSITION ANGLE - DEFAULT TO 0
  if n_elements(meas_pa) eq 0 then begin
     if keyword_set(verbose) then $
        message, "Position not supplied. Assuming measurement PA = 0.", /info     
     meas_pa = 0.0
  endif

; NO POSITION ANGLE - DEFAULT TO 0
  if n_elements(beam_pa) eq 0 then begin
     if keyword_set(verbose) then $
        message, "Position not supplied. Assuming beam PA = 0.", /info     
     beam_pa = 0.0
  endif

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; CALCULATIONS
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

; CONVERT TO RADIANS
  meas_theta = meas_pa*!dtor
  beam_theta = beam_pa*!dtor
  
; MATH (FROM MIRIAD VIA K. SANDSTROM)
  alpha = $
     (meas_maj*cos(meas_theta))^2d + (meas_min*sin(meas_theta))^2d - $
     (beam_maj*cos(beam_theta))^2d - (beam_min*sin(beam_theta))^2d

  beta = $
     (meas_maj*sin(meas_theta))^2d + (meas_min*cos(meas_theta))^2d - $
     (beam_maj*sin(beam_theta))^2d - (beam_min*cos(beam_theta))^2d
  
  gamma = $
     2d*((meas_min^2-meas_maj^2)*sin(meas_theta)*cos(meas_theta) - $
         (beam_min^2-beam_maj^2)*sin(beam_theta)*cos(beam_theta))

  s = alpha + beta
  t = sqrt((alpha-beta)^2d + gamma^2d)

; FIND THE SMALLEST RESOLUTION
  limit = meas_min < meas_maj < beam_maj < beam_min
  limit = 0.1*limit*limit
  
; TWO CASES:

  if (alpha lt 0 or beta lt 0 or s lt t) then BEGIN

;    ... FAILURE
     src_maj = 0
     src_min = 0
     src_pa = 0

;    ... NOTE THAT IT DIDN'T WORK
     worked = 0B
     if keyword_set(verbose) then $
         message, "Illegal alpha, beta, or s value.", /info     
 
;    ... CLOSE TO A POINT SOURCE
     if (0.5*(s-t) lt limit and $
         alpha gt -1*limit and $
         beta gt -1*limit) then begin
        point = 1B
     endif else begin
;    ... FAILURE BUT NOT A POINT SOURCE
        point = 0B
     endelse
     
  endif else begin

;    ... SUCCESS
     src_maj = sqrt(0.5*(s+t))
     src_min = sqrt(0.5*(s-t))
     
     if (abs(gamma)+abs(alpha-beta) eq 0) then BEGIN
        src_pa = 0
     endif else BEGIN
        src_pa  = (1.0d/!dtor)*0.5*atan(-1*gamma,alpha-beta)
     endelse
     
;    ... NOTE THAT IT WORKED AND IT'S NOT A POINT SOURCE
     worked = 1B
     point = 0B
  endelse

; RETURN!
  return

end
