function galaxies, name, m31 = m31, m33 = m33
;+
; NAME:
;   GALAXIES
; PURPOSE:
;   Returns structure with information regarding a given galaxy.
;
; CALLING SEQUENCE:
;   gstr = GALAXIES(name, [/M31, /M33])
;
; INPUTS:
;   NAME -- Name of the Galaxy in question.  Currently, I've entered
;           M31, M33, M64, NGC628, NGC1068, NGC2903, NGC3351, NGC3521,
;           NGC3627, NGC4528, NGC4321, NGC4414, NGC4736, NGC5055,
;           NGC5194, NGC6946, NGC7331, NGC5033, NGC5457
; KEYWORD PARAMETERS:
;   M31, M33 -- force values for these galaxies.
;
; OUTPUTS:
;   GSTR -- Structure with Tags:
;             NAME -- Name of galaxy
;             RA_GC -- RA of the Galaxy Center (J2000)
;             DEC_GC -- DEC of Galaxy Center (J2000)
;             RA, DEC -- As RA_GC and DEC_GC
;             DIST -- Distance in parsecs
;             POSANG -- Position angle of the major velocity axis
;             INC -- Inclination angle of galaxy
;             VLSR -- VLSR of the center
; MODIFICATION HISTORY:
;
;       Fri May 9 11:22:24 2003, <eros@master>
;		Documented.
;
;-

  if keyword_set(m31) then name = 'M31'
  if keyword_set(m33) then name = 'M33'
  name = strupcase(name)

  s = {name:'', ra_gc:0.d, dec_gc:0.d, dist:0., posang:0., inc:0., $
       vlsr:0., ra:0.d, dec:0.d}

  if name eq 'M31' then begin 
    s.name = 'M31'
    s.ra_gc = convang(double([00, 42, 44.31]), /ra)
    s.dec_gc = convang(double([41, 16, 09.4]))
    s.ra = s.ra_gc
    s.dec = s.dec_gc
    s.dist = 7.5e5
    s.posang = 34.0
    s.inc = 77.0               
    s.vlsr = -300
  endif  

  if name eq 'M64' then begin
    s.name = 'M64'
    s.ra_gc =  194.18183639
    s.dec_gc = 21.682682
    s.ra = s.ra_gc
    s.dec = s.dec_gc
    s.dist = 4.1e6
    s.posang = -67.6
    s.inc = 58.9             
    s.vlsr = 411.3
  endif

; ; Old M64 Params
;   if name eq 'M64' then begin
;     s.name = 'M64'
; ;   s.ra_gc = convang(double([12, 56, 43.5]), /ra)
; ;   s.dec_gc = convang(double([21, 40, 58.30]))
;     s.ra_gc = 1.94181629536e+02
;     s.dec_gc = 2.16830549517e+01
;     s.ra = s.ra_gc
;     s.dec = s.dec_gc
;     s.dist = 4.1e6
;     s.posang = -67.6
;     s.inc = 55               
;     s.vlsr = 411
;   endif

  if name eq 'M64KIN' then begin
    s.name = 'M64'
    s.ra_gc = convang(double([12, 56, 43.66]), /ra)
    s.dec_gc = convang(double([21, 41, 1]))
    s.ra = s.ra_gc
    s.dec = s.dec_gc
    s.dist = 4.1e6
    s.posang = -66.4
    s.inc = 56.3               
    s.vlsr = 411
  endif

  if name eq 'M33' then begin
    s.name = 'M33'
    s.ra_gc = convang(double([01, 33, 50.8]), /ra)
    s.dec_gc = convang(double([30, 39, 36.7]))
    s.ra = s.ra_gc
    s.dec = s.dec_gc
    s.dist = 8.4e5
    s.posang = 22.0+180         ; Because receding side is in the SW
;    s.inc = 56.0
    s.inc = 52.0     
    s.vlsr = -179.
  endif


  if name eq 'NGC628' then begin
    s.name = 'NGC628'
    s.ra_gc = convang(double([01, 36, 41.7]), /ra)
    s.dec_gc = convang(double([15, 46, 59.4]))
    s.ra = s.ra_gc
    s.dec = s.dec_gc
    s.dist = 7.3e6
    s.posang = 25
    s.inc = 24
    s.vlsr = 657
  endif

  if name eq 'NGC1068' then begin
    s.name = 'NGC1068'
    s.ra_gc = convang(double([02, 42, 40.74]), /ra)
    s.dec_gc = convang(double([00, 00, 47.7]))
    s.ra = s.ra_gc
    s.dec = s.dec_gc
    s.dist = 14.4e6
    s.posang = 13
    s.inc = 33
    s.vlsr = 1136
  endif

  if name eq 'NGC2903' then begin
    s.name = 'NGC2903'
    s.ra_gc = convang(double([09, 32, 10.05]), /ra)
    s.dec_gc = convang(double([21, 30, 0.2]))
    s.ra = s.ra_gc
    s.dec = s.dec_gc
    s.dist = 6.3e6
    s.posang = 17
    s.inc = 61
    s.vlsr = 556
  endif

  if name eq 'NGC3351' then begin
    s.name = 'NGC3351'
    s.ra_gc = convang(double([10, 43, 57.98]), /ra)
    s.dec_gc = convang(double([11, 42, 14.4]))
    s.ra = s.ra_gc
    s.dec = s.dec_gc
    s.dist = 10.1e6
    s.posang = 13
    s.inc = 40
    s.vlsr = 778
  endif

  if name eq 'NGC3521' then begin
    s.name = 'NGC3521'
    s.ra_gc = convang(double([11, 05, 49.26]), /ra)
    s.dec_gc = convang(double([00, 02, 02.3]))
    s.ra = s.ra_gc
    s.dec = s.dec_gc
    s.dist = 7.2e6
    s.posang = 164
    s.inc = 58
    s.vlsr = 805
  endif

  if name eq 'NGC3627' then begin
    s.name = 'NGC3627'
    s.ra_gc = convang(double([11, 20, 15.07]), /ra)
    s.dec_gc = convang(double([12, 59, 21.7]))
    s.ra = s.ra_gc
    s.dec = s.dec_gc
    s.dist = 11.1e6
    s.posang = 176
    s.inc = 63
    s.vlsr = 727
  endif

  if name eq 'NGC4258' then begin
    s.name = 'NGC4258'
    s.ra_gc = convang(double([12, 18, 57.52]), /ra)
    s.dec_gc = convang(double([47, 18, 14.2]))
    s.ra = s.ra_gc
    s.dec = s.dec_gc
    s.dist = 8.1e6
    s.posang = 176
    s.inc = 65
    s.vlsr = 448
  endif

  if name eq 'NGC4321' then begin
    s.name = 'NGC4321'
    s.ra_gc = convang(double([12, 22, 54.84]), /ra)
    s.dec_gc = convang(double([15, 49, 20]))
    s.ra = s.ra_gc
    s.dec = s.dec_gc
    s.dist = 16.7e6
    s.posang = 154
    s.inc = 30
    s.vlsr = 1571
  endif

  if name eq 'NGC4414' then begin
    s.name = 'NGC4414'
    s.ra_gc = convang(double([12, 26, 27.19]), /ra)
    s.dec_gc = convang(double([31, 13, 24]))
    s.ra = s.ra_gc
    s.dec = s.dec_gc
    s.dist = 19.1e6
    s.posang = 159 
    s.inc = 55 
    s.vlsr = 716
  endif

  if name eq 'NGC4736' then begin
    s.name = 'NGC4736'
    s.ra_gc = convang(double([12, 50, 53.06]), /ra)
    s.dec_gc = convang(double([41, 07, 13.6]))
    s.ra = s.ra_gc
    s.dec = s.dec_gc
    s.dist = 4.3e6
    s.posang = 100
    s.inc = 35
    s.vlsr = 308
  endif

  if name eq 'NGC5055' then begin
    s.name = 'NGC5055'
    s.ra_gc = convang(double([13, 15, 49.25]), /ra)
    s.dec_gc = convang(double([42, 01, 49.3]))
    s.ra = s.ra_gc
    s.dec = s.dec_gc
    s.dist = 7.2e6
    s.posang = 81
    s.inc = 56
    s.vlsr = 504
  endif

  if name eq 'NGC5194' then begin
    s.name = 'NGC5194'
    s.ra_gc = convang(double([13, 29, 52.35]), /ra)
    s.dec_gc = convang(double([47, 11, 53.8]))
    s.ra = s.ra_gc
    s.dec = s.dec_gc
    s.dist = 8.4e6
    s.posang = 0
    s.inc = 15
    s.vlsr = 463
  endif

  if name eq 'NGC6946' then begin
    s.name = 'NGC6946'
    s.ra_gc = convang(double([20, 34, 52.33]), /ra)
    s.dec_gc = convang(double([60, 09, 14.2]))
    s.ra = s.ra_gc
    s.dec = s.dec_gc
    s.dist = 5.5e6
    s.posang = 65
    s.inc = 54
    s.vlsr = 48
  endif

  if name eq 'NGC7331' then begin
    s.name = 'NGC7331'
    s.ra_gc = convang(double([22, 37, 4.09]), /ra)
    s.dec_gc = convang(double([34, 24, 56.3]))
    s.ra = s.ra_gc
    s.dec = s.dec_gc
    s.dist = 15.1e6
    s.posang = 172
    s.inc = 62
    s.vlsr = 821
  endif

  if name eq 'NGC5033' then begin
    s.name = 'NGC5033'
    s.ra_gc = convang(double([13, 13, 27.5]), /ra)
    s.dec_gc = convang(double([36, 35, 38]))
    s.ra = s.ra_gc
    s.dec = s.dec_gc
    s.dist = 17.5e6
    s.posang = 353
    s.inc = 68
    s.vlsr = 875
  endif

  if name eq 'NGC5457' then begin
    s.name = 'NGC5457'
    s.ra_gc = convang(double([14, 03, 12.5]), /ra)
    s.dec_gc = convang(double([54, 20, 55]))
    s.ra = s.ra_gc
    s.dec = s.dec_gc
    s.dist = 7.4e6
    s.posang = 42
    s.inc = 21
    s.vlsr = 241
  endif

  if name eq 'M82' then begin
    s.name = 'M82'
    s.ra_gc = convang(double([09, 55, 52.2]), /ra)
    s.dec_gc = convang(double([69, 40, 47]))
    s.ra = s.ra_gc
    s.dec = s.dec_gc
    s.dist = 3.25e6
    s.posang = 0
    s.inc = 0
    s.vlsr = 203
  endif

  if name eq 'M81' then begin
    s.name = 'M81'
    s.ra_gc = convang(double([09, 55, 33.2]), /ra)
    s.dec_gc = convang(double([69, 03, 55]))
    s.ra = s.ra_gc
    s.dec = s.dec_gc
    s.dist = 3.25e6
    s.posang = 0
    s.inc = 0
    s.vlsr = -34
  endif

  if name eq 'LMC' then begin
    s.name = 'LMC'
    s.ra_gc = convang(double([05, 23, 34.5]), /ra)
    s.dec_gc = convang(double([-69, 45, 22.0]))
    s.ra = s.ra_gc
    s.dec = s.dec_gc
    s.dist = 6.5e4
    s.posang = 0 
    s.inc = 0 
    s.vlsr = 278.0 
  endif

;  if name eq '' then begin
;    s.name = ''
;    s.ra_gc = convang(double([]), /ra)
;    s.dec_gc = convang(double([]))
;    s.ra = s.ra_gc
;    s.dec = s.dec_gc
;    s.dist = 
;    s.posang = 
;    s.inc = 
;    s.vlsr = 
;  endif

  return, s
end

