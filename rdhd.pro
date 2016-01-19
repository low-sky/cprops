pro rdhd, hd, structure = structure, cmat = cmat, ms = ms, fast = fast, $
          vector = vector, full = full
;+
; NAME:
;    rdhd
; PURPOSE:
;    To extract vital information from the header of a FITS file and
;    return it in a structure.  Designed for Interferometer Data.
;
; CALLING SEQUENCE:
;   READHD, hd, [structure = structure]
;
; INPUTS:
;   HD - A FITS header
;
; KEYWORD PARAMETERS:
;   STRUCTURE - The name of a structure in which the header
;               information is returned.
;   CMAT - A structure containing the matrices of RA and DEC values.
;   MS - Set this keyword to return velocities in m/s.  Defaults to km/s.
;   FAST - Set this keyword to skip calculation of the coordinate
;          matrices, returning only vectors.
;   VECTOR - Set this keyword to indicate a single spectrum is being
;              processed.  RDHD was originally designed for DATA
;              cubes.
;   FULL - Calculate the coordinate matrices with this keyword set.
; OUTPUTS:
;
;
; MODIFICATION HISTORY:
;
;       Wed Jun 2 09:40:27 2004, <eros@master>
;		Added more AIPS compatibility (FREQ0 convention and
;		multiple clean entries)
;
;       Tue May 11 12:29:45 2004, Erik Rosolowsky <eros@cosmic>
;       Added AIPS Beam compatibility.
;
;        Made
;       the /FAST keyword the default Fri Jul 12 13:27:45 2002, Erik
;       Rosolowsky <eros@cosmic>
;       
;       Added compatibility with GSS images.
;       Mon Jun 10 14:04:10 2002, Erik Rosolowsky <eros@cosmic>
;
;       Added /spectrum keyword.
;       Fri Feb 22 11:17:55 2002, Erik Rosolowsky <eros@cosmic>
;
;       Changed velocity indexing to reflect FITS is 1 indexed and not
;       zero indexed.
;       Tue Dec 4 10:52:43 2001, Erik Rosolowsky <eros@cosmic>
;
;      Changed JY/Beam -> Kelvin convervsion factor calculation to a
;      geometric mean from an arithmetic mean.  This is more correct.
;      Wed Aug 8 11:05:41 2001, Erik Rosolowsky <eros@cosmic>
;
;      Implemented FAST keyword.
;      Tue Jul 31 15:13:03 2001, Erik Rosolowsky <eros@cosmic>
;
;      Fixed bug for 3D images to include CDELT information in structure.
;      Tue Jan 23 00:02:39 2001, Erik Rosolowsky <eros@cosmic>
;
;      Fixed bug for 2D images -- Tue Dec 5 01:13:31 2000, Erik
;                                 Rosolowsky <eros@cosmic>
;
;      Added CDELT and CRVAL keywords to output structure -
;      Mon Nov 13 10:57:52 2000, Erik Rosolowsky <eros@cosmic>
;      Written - Mon Oct 30 19:38:50 2000, Erik Rosolowsky
;                <eros@cosmic>
;
;-

  if keyword_set(vector) then begin
    rdhd_sp, hd, structure = structure, ms = ms
    return
  endif 
  
  if keyword_set(full) then fast = 0b else fast = 1b

  dim = sxpar(hd, 'NAXIS')
  naxis1 = sxpar(hd, 'NAXIS1')
  naxis2 = sxpar(hd, 'NAXIS2')
  naxis3 = sxpar(hd, 'NAXIS3')
  bunit = sxpar(hd, 'BUNIT')
  bm_maj = sxpar(hd, 'BMAJ')*3600
  bm_min = sxpar(hd, 'BMIN')*3600
  bpa = sxpar(hd, 'BPA')

; AIPS beam extraction bits.
  if bm_maj eq 0 then begin
    ind = where(strpos(hd, 'CLEAN') gt 0 and $
                strpos(hd, 'BMAJ') gt 0 and $
                strpos(hd, 'AIPS') gt 0,  ct)
    if ct gt 0 then begin
      str = hd[max(ind)]
      str = strsplit(str, /ext)
      bm_maj = float(str[4])*3600
      bm_min = float(str[6])*3600
      bpa = float(str[8])
    endif
  endif
  freq = sxpar(hd, 'RESTFREQ') > sxpar(hd, 'FREQ0')
  date = sxpar(hd, 'DATE-OBS')
  extast, hd, astrom
;nelts = max([naxis1, naxis2])
  lam_mm = 3d11/freq
; Conversion Factors for the map to change from Jy/beam to K.
  jypb2k = 14*(lam_mm/(sqrt(bm_maj*bm_min)))^2
  k2jypb = 1./jypb2k




  strings = astrom.ctype
  if total(strpos(strings, 'GSS') gt 0) then begin
    gsssextast, hd, astrom_gss
    if not keyword_set(fast) then begin
      gsssxyad,  astrom_gss, findgen(naxis1)#replicate(1., naxis2), $
                 replicate(1., naxis1)#findgen(naxis2), ra, dec
      gsssxyad, astrom_gss, 0, indgen(naxis2), null, decvec
      gsssxyad, astrom_gss, indgen(naxis1), 0, ravec, null
    endif else begin
      gsssxyad, astrom_gss, 0, indgen(naxis2), null, decvec
      gsssxyad,  astrom_gss, indgen(naxis1), 0, ravec, null
    endelse
    structure = {naxis1:naxis1, naxis2:naxis2, $
                 ra:ravec, dec:decvec, $
                 ctype:astrom.ctype[0:dim-1]}
  endif else begin
; if requested, calculate coord. matrices, otherwise just generated
; the vector or coordinates through 
    

    if not keyword_set(fast) then begin
      if dim eq 1 then begin
        ravec = sxpar(hd, 'CRVAL2')
        decvec = sxpar(hd, 'CRVAL3')
        ra = ravec
        dec = decvec
      endif else begin
        xy2ad, findgen(naxis1)#replicate(1., naxis2), $
               replicate(1., naxis1)#findgen(naxis2), astrom, ra, dec
        xy2ad, astrom.crpix[0], indgen(naxis2), astrom, null, decvec
        xy2ad, indgen(naxis1), astrom.crpix[1], astrom, ravec, null
      endelse
      cmat = {ra:ra, dec:dec}	
    endif else begin
      if dim eq 1 then begin
        ravec = sxpar(hd, 'CRVAL2')
        decvec = sxpar(hd, 'CRVAL3')
      endif else begin
;        xy2ad, astrom.crpix[0], indgen(naxis2), astrom, null, decvec
;        xy2ad, indgen(naxis1), astrom.crpix[1], astrom, ravec, null
        xy2ad, naxis1/2, indgen(naxis2), astrom, null, decvec
        xy2ad, indgen(naxis1), naxis2/2, astrom, ravec, null
      endelse
    endelse

    if naxis3 gt 1 then begin
; First, regenerate the astrometry structure to include the THIRD
; DIMENSION!!!
      astrom2 = {CD:astrom.cd, $
                 CDELT:[astrom.cdelt, sxpar(hd, 'CDELT3')], $
                 CRPIX:[astrom.crpix, sxpar(hd, 'CRPIX3')], $
                 CRVAL:[astrom.crval, sxpar(hd, 'CRVAL3')], $
                 CTYPE:[astrom.ctype, sxpar(hd, 'CTYPE3')], $
                 LONGPOLE:astrom.longpole, $
                 LATPOLE:astrom.latpole, $
                 PV2:astrom.pv2}
      astrom = astrom2
      velvec = (findgen(naxis3)+1-astrom.crpix[2])*$
               astrom.cdelt[2]+astrom.crval[2]
    endif else begin
      velvec = 0 
      dim = 2
    endelse
    if not keyword_set(ms) then velvec = velvec/1000.

    ppbeam = abs((bm_maj*bm_min/3600.^2)/(astrom.cdelt[0]*astrom.cdelt[1])*$
                 2*!pi/(8*alog(2)))

    structure = {naxis1:naxis1, naxis2:naxis2, $
                 naxis3:naxis3, ra:ravec, dec:decvec, v:velvec, $
                 ctype:astrom.ctype, $
                 crpix:astrom.crpix, $
                 crval:astrom.crval, $
                 cdelt:astrom.cdelt, $
                 k2jypb:k2jypb, jypb2k:jypb2k, freq:freq, $
                 bmaj:bm_maj, bmin:bm_min, date:date, bpa:bpa, $
                 ppbeam:ppbeam}
  endelse
  return
end


