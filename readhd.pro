;...................... START READHD ......................
pro readhd,header
;---------------------
; read in fits header
;---------------------
@header.cb

if (n_params(0) lt 1) then begin 
  print,'PROCEDURE readhd, header'
  print,' header = FITS header' 
  return 
endif  

naxis1=0 & naxis2=0 & naxis3=0
crval1=0.0 & crvax2=0.0 & crval3=0.0
crpix1=0.0 & crpix2=0.0 & crpix3=0.0
cdelt1=0.0 & cdelt2=0.0 & cdelt3=0.0
ctype1='????????' & ctype2='????????' & ctype3='????????'
bmaj=0.0 & bmin=0.0 & bpa=0.0

s=size(header) & nlines=s(1)
for i=0,nlines-1 do begin
  if (strpos(header(i),'NAXIS1') eq 0) then begin
    s=strmid(header(i),10,20) & s=strtrim(s) & naxis1=fix(s)
  endif
  if (strpos(header(i),'NAXIS2') eq 0) then begin
    s=strmid(header(i),10,20) & s=strtrim(s) & naxis2=fix(s)
  endif
  if (strpos(header(i),'NAXIS3') eq 0) then begin
    s=strmid(header(i),10,20) & s=strtrim(s) & naxis3=fix(s)
  endif
  if (strpos(header(i),'CRVAL1') eq 0) then begin
    s=strmid(header(i),10,20) & s=strtrim(s) & crval1=double(s)
  endif
  if (strpos(header(i),'CRVAL2') eq 0) then begin
    s=strmid(header(i),10,20) & s=strtrim(s) & crval2=double(s)
  endif
  if (strpos(header(i),'CRVAL3') eq 0) then begin
    s=strmid(header(i),10,20) & s=strtrim(s) & crval3=double(s)
  endif
  if (strpos(header(i),'CRPIX1') eq 0) then begin
    s=strmid(header(i),10,20) & s=strtrim(s) & crpix1=double(s)
  endif
  if (strpos(header(i),'CRPIX2') eq 0) then begin
    s=strmid(header(i),10,20) & s=strtrim(s) & crpix2=double(s)
  endif
  if (strpos(header(i),'CRPIX3') eq 0) then begin
    s=strmid(header(i),10,20) & s=strtrim(s) & crpix3=double(s)
  endif
  if (strpos(header(i),'CDELT1') eq 0) then begin
    s=strmid(header(i),10,20) & s=strtrim(s) & cdelt1=double(s)
  endif
  if (strpos(header(i),'CDELT2') eq 0) then begin
    s=strmid(header(i),10,20) & s=strtrim(s) & cdelt2=double(s)
  endif
  if (strpos(header(i),'CD1_1') eq 0) then begin
    s=strmid(header(i),10,20) & s=strtrim(s) & cd1_1=double(s)
  endif
  if (strpos(header(i),'CD2_2') eq 0) then begin
    s=strmid(header(i),10,20) & s=strtrim(s) & cd2_2=double(s)
  endif
  if (strpos(header(i),'CDELT3') eq 0) then begin
    s=strmid(header(i),10,20) & s=strtrim(s) & cdelt3=double(s)
  endif
  if (strpos(header(i),'CTYPE1') eq 0) then begin
    s=strmid(header(i),10,20) & ctype1=strtrim(s)
  endif
  if (strpos(header(i),'CTYPE2') eq 0) then begin
    s=strmid(header(i),10,20) & ctype2=strtrim(s)
  endif
  if (strpos(header(i),'CTYPE3') eq 0) then begin
    s=strmid(header(i),10,20) & ctype3=strtrim(s)
  endif
end

; IDL array indices run from 0 to n-1
; FITS array indices run from 1 to n
; hence we must convert ---
crpix1=crpix1-1
crpix2=crpix2-1
crpix3=crpix3-1

return
end
;......................  END READHD  ......................
