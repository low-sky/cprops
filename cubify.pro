pro cubify,  x, y, v, t, szin, cube = cube, mask = mask, indcube = indcube $
             , id = id, indvec = indvec, pad = pad, location = location $
             , twod = twod
;+
; NAME:
;   CUBIFY
; PURPOSE:
;   Converts vectorized arrays and IDs into 
;
; CALLING SEQUENCE:
;   CUBIFY, x, y, v, t, sz, CUBE = cube [,ID = id, MASK = mask]
;
; INPUTS:
;   X,Y,V -- Pixel coordinates (zero indexed) of elements in a datacube
;   T -- Brigthness units to fill cube with.
;   SZ -- Size structure characterizing the datacube. 
; KEYWORD PARAMETERS:
;   ID -- Vector of cloud identifications.  
;   TWOD -- Set this flag to force "cubification" to a 2-D matrix.
;           Mostly for compatibility with routines outside the CRPOPS
;           distribution. 
; OUTPUTS:
;   CUBE -- Set this to named variable containing the output data cube
;           restored to full size.
;   MASK -- Set this to named variable containing the output mask
;           array with values equal to ID for each pixel.  If ID is
;           not set, then mask defaults to a BYTE array with a value
;           of 1B where every value of T was set.
;
; MODIFICATION HISTORY:
;
;       Thu Dec 2 12:29:39 2004, Erik Rosolowsky <eros@cosmic>
;		Written.
;
;-

  help, calls = calls
  if n_elements(szin) lt 4 then begin 
    if n_elements(calls) eq 2 then $
    message, 'No size information.  Assuming minimum size.', /con

;   FIGURE OUT THE MINIMUM SIZE OF THE (POSSIBLY PADDED) CUBE
    if (not keyword_set(pad)) then begin
      sz = [-1, max(x)-min(x)+1, max(y)-min(y)+1, max(v)-min(v)+1] 
      minx = min(x)
      miny = min(y)
      minv = min(v)
    endif else begin
      sz = [-1, max(x)-min(x)+1+2*pad $
            , max(y)-min(y)+1+2*pad, max(v)-min(v)+1+2*pad]
      minx = min(x) - pad
      miny = min(y) - pad
      minv = min(v) - pad
    endelse
    if sz[0] eq 2 or keyword_set(twod) then begin
      sz[3] = 1
      minv = 0
    endif
;   NUMBER OF ELEMENTS
    sz[0] = long(sz[1])*sz[2]*sz[3]

;   MAKE EMPTY ARRAYS OF APPROPRIATE TYPE
    cube = make_array(sz[1], sz[2], sz[3], type = size(t, /type))*$
           !values.f_nan
    mask = make_array(sz[1], sz[2], sz[3], type = size(id, /type) > 1)
    indcube = make_array(sz[1], sz[2], sz[3], $
                         type = size(indvec, /type) > 1)-1

;   PLACE THE DATA INTO THE (POSSIBLY PADDED) ARRAY
    cube[x-minx, y-miny, v-minv] = t

;   AND FILL OUT THE MASK (WITH 1s IF NOT SUPPLIED)
    if (n_elements(id) eq n_elements(x)) then $
      mask[x-minx, y-miny, v-minv] = id else $
      mask[x-minx, y-miny, v-minv] = 1b
    
;   AND FILL OUT THE MASK (WITH 1s IF NOT SUPPLIED)
    if (n_elements(indvec) eq n_elements(x)) then $
      indcube[x-minx, y-miny, v-minv] = indvec else $
      indcube[x-minx, y-miny, v-minv] = -1

    indexcube = lindgen(sz[1], sz[2], sz[3])
    location = indexcube[x-minx, y-miny, v-minv]

  endif else begin
    sz = szin
    if sz[0] eq 2 then sz[3] = 1
    if n_elements(minx) eq 0 then minx = 0
    if n_elements(miny) eq 0 then miny = 0
    if n_elements(minv) eq 0 then minv = 0
;   ERROR CHECKING
    if (keyword_set(pad) AND (n_elements(sz) gt 0)) then $
      print, 'Keyword PAD has no effect when SZ is supplied. Ignoring. Dumbass.'

;   MAKE THE EMPTY CUBE
    cube = make_array(sz[1], sz[2], sz[3], type = size(t, /type))*$
           !values.f_nan
    mask = make_array(sz[1], sz[2], sz[3], type = size(id, /type) > 1)
    indcube = make_array(sz[1], sz[2], sz[3], type = size(indvec, /type) > 1)

;   FILL IN THE DATA
    cube[x, y, v] = t

;   AND THE MASK
    if n_elements(id) eq n_elements(x) then $
      mask[x, y, v] = id else $
      mask[x, y, v] = 1B

;   AND THE INDEX CUBE
    if n_elements(indvec) eq n_elements(x) then $
      indcube[x, y, v] = indvec else $
      indcube[x, y, v] = 1B

    indexcube = lindgen(sz[1], sz[2], sz[3])
    location = indexcube[x-minx, y-miny, v-minv]

  endelse

  return
end

