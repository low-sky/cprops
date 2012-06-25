pro addnoise, data, sigma, corrnoise = corrnoise, beaminpix = beaminpix

;+
;
; NAME:
;   ADDNOISE
; PURPOSE:
;   Adds noise to a datacube
; CALLING SEQUENCE:
;   ADDNOISE, data, sigma [/corrnoise, beaminpix = beaminpix]
; INPUTS:
;   DATA -- a datacube
;   SIGMA -- The variance of the noise _added_ to the datacube.
; KEYWORD PARAMETERS:
;   /CORRNOISE -- Smooths the noise to be correlated over the size of
;                 the beam.
;   BEAMINPIX -- The length of the beam, measured in pixels.
; OUTPUTS:
;   DATA -- replaces data in situ.
; SOME PROCEDURES THAT AREN'T IDL/GODDARD THAT WE USE:
;
; MODIFICATION HISTORY:
;   Adam's code.
;-

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; I. ADD NOISE TO THE DATA CUBE
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

; IF THE ADDNOISE KEYWORD IS SET JUST ADD AN APPROPRIATE AMOUNT OF
; NOISE IN THE DUMBEST WAY POSSIBLE
  if (n_elements(corrnoise) eq 0) then begin
    data = data + randomn(seed, n_elements(data))*sigma
    return
  endif else begin
; IF THE ADDCORRNOISE KEYWORD IS SET, ADD THE APPROPRIATE LEVEL OF
; *CORRELATED* NOISE (MAKE THE NOISE, CONVOLVE WITH A BEAM AND THEN
; ADD).
    if (n_elements(beaminpix) eq 0) then begin
      return
    endif

;   AN UNCORRELATED CUBE OF NOISE
    corrnoise = finite(data)*0. + randomn(seed, n_elements(data))*sigma

;   MAKE A 2D BEAM. DOESN'T HAVE TO BE PERFECT, JUST CLOSE ENOUGH TO
;   GET THE NOISE CORRELATED ABOUT RIGHT. MAKE IT TWO FWHM ACROSS...
    corrxaxis = findgen(beaminpix*2+1)    
    corryaxis = findgen(beaminpix*2+1)
    corrxaxis = corrxaxis - mean(corrxaxis)
    corryaxis = corryaxis - mean(corryaxis)
    corrxra = (fltarr(n_elements(corrxaxis)) + 1.0) ## corrxaxis
    corryra = corryaxis ## (fltarr(n_elements(corryaxis)) + 1.0)
    corrbeam = exp(-(corrxra^2/2./(beaminpix^2-1)) - $
                   (corryra^2/2./(beaminpix^2-1)))
    corrbeam = corrbeam/total(corrbeam, /nan)
;   CONVOLVE THE CUBE, PLANE BY PLANE WITH THE PSF
    for l = 0, (size(corrnoise))[3]-1 do begin
      corrnoise[*, *, l] = convolve(corrnoise[*, *, l], corrbeam)
    endfor

;   NORMALIZE IT, JUST IN CASE
    corrnoise = corrnoise*sigma/mad(corrnoise)
    data = data + corrnoise
    return
  endelse

  return
end
