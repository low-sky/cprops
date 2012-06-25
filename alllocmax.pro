function alllocmax, cubein, indcube = indcube, friends = friends, $
                    specfriends = specfriends

;+
;
; NAME:
;   ALLLOCMAX()
; PURPOSE:
;   To establish a candidate set of local maxima within a data cube by
;   searching over a F by F by S box of pixels for a point
;   that's larger than all others in the box.  F and S default to 1
;
; CALLING SEQUENCE:
;   local_maxima = ALLLOCMAX(cube, [friends = friends, specfriends =
;   specfriends, indcube = indcube])
;
; INPUTS:
;   CUBE -- data cube to find local maxima in.
;             
; KEYWORD PARAMETERS:
;   INDCUBE -- (optional) if supplied, then the indices returned are
; drawn from this cube instead of indexing the supplied cube itself.
;   FRIENDS -- Sets the search region in the position dimensions to be
;              2*FRIENDS+1 in size.  
;   SPECFRIENDS -- Sets the search region in the velocity dimensions to be
;              2*SPECFRIENDS+1 in size.  
;
; OUTPUTS:
;   LOCAL_MAXIMA -- indices in CUBE (or in INDCUBE if set) that are
;                   local maxima.
;
; MODIFICATION HISTORY:
;
;       Documented -- Fri Sep 2 15:48:17 2005, Erik Rosolowsky
;                     <erosolow@asgard.cfa.harvard.edu>
;-

; GET THE SIZE AND MAKE AN EMPTY CUBE TO FLAG WHERE LOCAL MAXIMA
; EXIST.
  sz = size(cubein)
  if sz[0] eq 2 then begin
    specfriends = 0
    sz[3] = 1
  endif
; INITIALIZE EVERYTHING TO BE A LOCAL MAX, WE WILL REJECT THEM OVER TIME
  lmaxcube = bytarr(sz[1], sz[2], sz[3]) + 1B

; INITIALIZE THE DEFAULT BOX SIZE TO BE +/- ONE PIXEL (3 x 3 BOX)
  if (n_elements(friends) eq 0) then $
    friends = 1

; SET THE NONSENSICAL INDICES TO NOT-A-NUMBERS
  badind = where(cubein ne cubein, badct)
  cube = cubein
  if (badct gt 0) then begin
    cube[badind] = -!values.f_infinity
    lmaxcube[badind] = 0B
  endif

  if n_elements(specfriends) eq 0 then specfriends = 1 


; A LOCAL MAXIMA IS A POINT WHICH IS GREATER THAN ALL OF THE POINTS
; AROUND IT.

  if specfriends gt 0 then begin
;;     for k = -specfriends, specfriends do $
;;       for j = -friends, friends do $
;;         for i = -friends, friends do $
;;           if NOT ((i eq 0) AND (j eq 0) AND (k eq 0)) then $
;;             lmaxcube = $
;;       lmaxcube*((cube gt shift(cube, i, j, k)) OR $
;;                 (finite(shift(cube, i, j, k)) eq 0)) 
    for k = -specfriends, specfriends do begin
      for j = -friends, friends do begin
        for i = -friends, friends do begin
          if NOT ((i eq 0) AND (j eq 0) AND (k eq 0)) then $
            lmaxcube = $
      lmaxcube*((cube gt shift(cube, i, j, k)) OR $
                (finite(shift(cube, i, j, k)) eq 0)) 
        endfor
      endfor
    endfor
  endif else begin
    for j = -friends, friends do $
      for i = -friends, friends do $
        if NOT ((i eq 0) AND (j eq 0)) then $
          lmaxcube = $
      lmaxcube*((cube gt shift(cube, i, j)) OR $
                (finite(shift(cube, i, j)) eq 0))
  endelse

    lmaxind = where(lmaxcube eq 1B, num)
  if (num eq 0) then begin
    message, 'No true local max found, defaulting to high point in data.', /con
    dummy = max(lmaxcube, lmaxind, /nan)
  endif

; IF THE INDEX CUBE IS SUPPLIED AND THERE ARE LOCAL MAXIMA THEN
; SUBSTITUTE THE INDICES FROM THE CUBE FOR THE ACTUAL INDICES
  if ((n_elements(indcube) gt 0)) then begin
    lmaxind = indcube[lmaxind]
  endif 

  return, lmaxind
end                             ; OF ALLLOCMAX


