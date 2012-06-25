function cloudalyze, cube, mask, gal = gal, hd = hdin, dist = dist $
                     , xco = xco, inK = inK, noextrap = noextrap $
                     , filename = filename, kindist = usekindist $
                     , bootstrap = bootstrap, ppbeam = ppbeam $
                     , near = near, far = far, $
                     sigma = sigma, ecube = ecube, r0 = r0, v0 = v0
;+
;
; NAME:
;   CLOUDALYZE
; PURPOSE:
;   To assign physical properties to identified clouds using the
;   methods of the PASP article.  See the manual for a detailed
;   explanation of what's going on here.
;
; CALLING SEQUENCE:
;   properties = CLOUDALYZE(data, mask [, gal = gal, hd = hd, dist =
;   dist, xco = xco, /inK, /noextrap, filename = filename,
;   /kindist, bootstrap = bootstrap, ppbeam = ppbeam, /near, /far])
;
; INPUTS:
;   DATA -- The orignal data cube
;   MASK -- An integer assignment mask so that every cloud has a
;           unique integer assignment.
;
; KEYWORD PARAMETERS:
;   /INK -- The orignal data are already in Kelvins
;   /NOEXTRAP -- Do not extrapolate cloud properties to 0 K contour.
;   /KINDIST -- For Milky Way clouds, use the kinematic distance
;   /NEAR, /FAR -- Used for kinematic distances in the Milky Way.
;   GAL -- One of two datatypes.  Either AKL's GAL structure, which
;          contains header information or EWR's GAL strcuture which
;          contains information about the target galaxy.
;   HD -- The FITS header, required with EWR's GAL structure or just
;         the distance , but not AKL's
;   DIST -- the distance, in PARSECS to the target object,
;   PPBEAM -- Number of pixels per beam, defaults to 1..  
;   BOOTSTRAP -- Number of bootstrap iterations to use in estimating
;                cloud uncertainties.
;   XCO -- The CO-to-H2 conversion factor in cm^-2/K km/s (defaults to
;          2 X 10^20.
; OUTPUTS:
;   PROPERTIES -- An array of structures (1 per cloud) that have the
;                 physical properties of the clouds.
;
; MODIFICATION HISTORY:
;       Documented.
;	Mon Sep  5 13:27:34 2005
;   
; REQUIRES:
; 
;  CLOUDMOM (AND EXTRAP)
; 
;  VECTORIFY
;
;  RDHD (maybe)
;
;-

  forward_function kindist

; FOR CONVENIENCE, ASSIGN "not-a-number" TO THE nan VARIABLE
  nan = !values.f_nan
  dnan = !values.d_nan
; CREATE A STRUCTURE THAT WILL HOLD THE MEASURED GMC PROPERTIES
  origcprops = { $
                 npix : long(0), $ ; NUMBER OF PIXELS
                 cloudnum : long(0), $ ; IDENTIFIER -- NOT IN USE
                 xctr_pix : nan, $ ; PIXELS
                 xctr_deg : dnan, $ ; RA DEGREES?
                 yctr_pix : nan, $ ; PIXELS
                 yctr_deg : dnan, $ ; DEC DEGREES?
                 vctr_pix : nan, $ ; PIXELS
                 vctr_kms : dnan, $ ; LSR KM/S?
                 posang : nan, $ ; POSITION ANGLE OF THE MAJOR AXIS
                 momxpix : nan, $ ; PIXELS
                 momxpix_uc : nan, $ ; FRACTIONAL UNCERTAINTY
                 momxpix_noex : nan, $ ; FRACTIONAL UNCERTAINTY
                 momxpix_noex_uc : nan, $ ; FRACTIONAL UNCERTAINTY
                 momypix : nan, $ ; PIXELS
                 momypix_uc : nan, $ ; FRACTIONAL UNCERTAINTY
                 momypix_noex : nan, $ ; FRACTIONAL UNCERTAINTY
                 momypix_noex_uc : nan, $ ; FRACTIONAL UNCERTAINTY
                 mommajpix : nan, $ ; PIXELS
                 mommajpix_uc : nan, $ ; FRACTIONAL UNCERTAINTY
                 mommajpix_noex : nan, $ ; FRACTIONAL UNCERTAINTY
                 mommajpix_noex_uc : nan, $ ; FRACTIONAL UNCERTAINTY   
                 momminpix : nan, $ ; PIXELS
                 momminpix_uc : nan, $ ; FRACTIONAL UNCERTAINTY
                 momminpix_noex : nan, $ ; FRACTIONAL UNCERTAINTY
                 momminpix_noex_uc : nan, $ ; FRACTIONAL UNCERTAINTY   
                 momvpix : nan, $ ; PIXELS
                 momvpix_uc : nan, $ ; FRACTIONAL UNCERTAINTY
                 momvpix_noex : nan, $ ; FRACTIONAL UNCERTAINTY
                 momvpix_noex_uc : nan, $ ; FRACTIONAL UNCERTAINTY
                 rad_pc : nan, $ ; PARSECS
                 rad_uc : nan, $ ; FRACTIONAL UNCERTAINTY
                 rad_nodc : nan, $ ; PARSECS
                 rad_nodc_uc : nan, $ ; FRACTIONAL UNCERTAINTY
                 rad_noex : nan, $ ; PARSECS
                 rad_noex_uc : nan, $ ; FRACTIONAL UNCERTAINTY
                 rad_nodc_noex : nan, $ ; PARSECS
                 rad_nodc_noex_uc : nan, $ ; FRACTIONAL UNCERTAINTY
                 sigv_kms : nan, $ ; KM/S
                 sigv_uc : nan, $ ; FRACTIONAL UNCERTAINTY
                 sigv_nodc : nan, $ ; KM/S
                 sigv_nodc_uc : nan, $ ; FRACTIONAL UNCERTAINTY
                 sigv_noex : nan, $ ; KM/S
                 sigv_noex_uc : nan, $ ; FRACTIONAL UNCERTAINTY
                 sigv_nodc_noex : nan, $ ; KM/S
                 sigv_nodc_noex_uc : nan, $ ; FRACTIONAL UNCERTAINTY
                 flux_kkms_pc2 : nan, $ ; K KM/S
                 flux_uc : nan, $ ; FRACTIONAL UNCERTAINTY
                 flux_noex : nan, $ ; K KM/S
                 flux_noex_uc : nan, $ ; FRACTIONAL UNCERTAINTY
                 mlum_msun : nan, $ ; MSUN                 
                 mlum_uc : nan, $ ; FRACTIONAL UNCERTAINTY
                 mvir_msun : nan, $ ; MSUN
                 mvir_uc : nan, $ ; FRACTIONAL UNCERTAINTY
                 name : '', $   ; STRING OF DATASET NAME
                 beamfwhm_pc : nan, $ ; PARSECS
                 sigchan_kms : nan, $ ; KM/S SIGMA OF ONE CHANNEL
                 galname : '', $ ; NAME OF GALAXY WHERE CLOUD LIVES
                 rmstorad : nan, $ ; RMS -> RADIUS CONVERSION
                 ppbeam : nan, $ ; PIXELS PER BEAM (AREA)
;                 area_pc2 : nan, $ ; AREA in PARSECS^2
;                 area_noex : nan, $
                 s2n: nan, $
                 tmax_K: nan, $
                 distance_pc : nan $ ; DISTANCE TO CLOUD                
               }

; IF AN EMPTY MASK IS PASSED TO CLOUDALYZE PASS BACK A DUMMY
; STRUCTURE.
  if total(mask) eq 0 then return, origcprops

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; DEFAULTS & DEFINITIONS
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

; EXTRACT A VECTOR OF INTERESTING ELEMENTS FROM THE CUBE
  vectorify, cube, mask = mask, x = x, y = y, v = v, t = t, id = id, sz = sz
  if n_elements(ecube) gt 0 then $
      vectorify, ecube, mask = mask, x = x, y = y, $
                 v = v, t = err

  if n_elements(usekindist) eq 0 then usekindist = 0b


; This curious little set of conditionals allows a AKL's GAL structure
; to function to pass header information but EWR's GAL structure only
; passes physical information and is used later.  We assume that EWR
; passes header and a GAL structure whereas AKL only passes a GAL
; structure.

  if (keyword_set(hdin)) then begin
; THINGS WE CAN GET FROM A .fits HEADER
    hd = hdin
    if keyword_set(usekindist) and 1b-stregex(sxpar(hd, 'CTYPE1'), 'GLON', /bool) then begin
        message, /con, 'Warning!  Header not in GLON/GLAT coordinates.  Blindly attempting to convert! '

      heuler, hd, /galactic
    endif
    rdhd, hd, s = hdstruct
    extast, hd, astrom
    beamfwhm = sqrt(hdstruct.bmaj*hdstruct.bmin)
    ppbeam = hdstruct.ppbeam
    raxis = hdstruct.ra
    daxis = hdstruct.dec
    vaxis = hdstruct.v
    bunit = sxpar(hd, 'BUNIT')
    if stregex(strupcase(bunit), 'JY/BEAM', /bool) then begin
      t = t*hdstruct.jypb2k     ; MAKE IT KELVIN!!!
      if n_elements(sigma) gt 0 then  sigma = sigma*hdstruct.jypb2k 
      if n_elements(err) gt 0 then err = err*hdstruct.jypb2k
    endif
    if stregex(bunit, 'K.KM/S', /bool) or $
      stregex(bunit, 'K km/s', /bool) then begin
      t = t/(hdstruct.cdelt[[2]]/1e3)[0] ; MAKE IT KELVIN!!
       if n_elements(sigma) gt 0 then $
         sigma = sigma/(hdstruct.cdelt[[2]]/1e3)[0] ; MAKE IT KELVIN!!
      if n_elements(err) gt 0 then err = $
        err/(hdstruct.cdelt[[2]]/1e3)[0] ; MAKE IT KELVIN!!
    endif
  endif else begin
    if (keyword_set(gal)) then begin
      beamfwhm = gal.beamfwhm   ; ARCSEC
      raxis = gal.raxis         ; DECIMAL DEGREES
      daxis = gal.daxis         ; DECIMAL DEGREES
      vaxis = gal.vaxis         ; KM/S
      if (NOT keyword_set(inK)) then begin
        t = t/abs(gal.vaxis[1] - gal.vaxis[0]) ; GO FROM K KM/S -> K
        if n_elements(err) gt 0 then $
          err = err/abs(gal.vaxis[1] - gal.vaxis[0])
         if n_elements(sigma) gt 0 then sigma = sigma/abs(gal.vaxis[1] - gal.vaxis[0])
      endif
    endif
  endelse

; ERROR CATCHING ON BEAM INFORMATION
  if (n_elements(beamfwhm) eq 0) then begin
    print, "Need some kinda BEAM information, you're killin' us here."
    return, -1
  endif

; THINGS WE CAN'T GET FROM A .fits HEADER
  if (keyword_set(gal)) then begin 
    dist = gal.dist             ; PARSECS
    if total(tag_names(gal) eq 'NAME') then $
      origcprops.galname = gal.name
  endif
; DEFINE THE PHYSICAL AND ASTRONOMICAL CONSTANTS THAT WE WILL NEED
; (AVOIDS THE NEED TO CALL AN EXTERNAL CONSTANTS BATCH FILE)
  mh = 1.673534d-24             ; hydrogen mass CGS
  ms = 1.98900d+33              ; solar mass CGS
  pc = 3.0857d18                ; parsec CGS
  if (n_elements(xco) eq 0) then $
    xco = 2.d20                 ; "Galactic" XCO cm^-2/K km s^-1
  origcprops.rmstorad = 1.91    ; conversion from moment to radius
  gainuc = 0.2                  ; gain uncertainty

  if ppbeam ne ppbeam or ppbeam le 0 then ppbeam = 1d0

  origcprops.ppbeam = ppbeam

; SET THE TARGET CONTOUR FOR THE EXTRAPOLATION. DEFAULTS TO ZERO
; KELVIN CONTOUR (OUR PREFERRED TARGET)
  if (n_elements(targett) eq 0) then $
    targett = 0.0

; GET THE UNIQUE MASK ELEMENTS
  cloudids = id[uniq(id, sort(id))]

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; MEASURE PROPERTIES OF EACH SUPPLIED CLOUD 
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
  
  for i = 0, n_elements(cloudids)-1 do begin

    useind = where(id eq cloudids[i], num)
    xuse = x[useind]
    yuse = y[useind]
    vuse = v[useind]
    tuse = t[useind]

;   CALL THE "cloudmom" ROUTINE TO CALCULATE THE CUMULATIVE MOMENTS
;   AND THEN EXTRAPOLATE THEM TO THE TARGET CONTOUR, "targett"

    mom = cloudmom(xuse, yuse, vuse, tuse, targett = targett)

    mom_noex = cloudmom(xuse, yuse, vuse, tuse, targett = targett, $
                        /noextrap)

;   FIND THE MAJOR AXIS AND ROTATE THEN MEASURE MAJOR/MINOR AXES
    pa = pa_moment(xuse, yuse, tuse)

    xrot = xuse*cos(pa)+yuse*sin(pa)
    yrot = -xuse*sin(pa)+yuse*cos(pa)    
    mom_rot = $ 
      cloudmom(xrot, yrot, vuse, tuse, targett = targett)
    mom_noex_rot = $
      cloudmom(xrot, yrot, vuse, tuse, targett = targett, /noextrap)

    if (n_elements(bootstrap) gt 0) then begin
      npts = n_elements(xuse)
      for j = 0, bootstrap-1 do begin
;       GENERATE A NEW SET OF DATA POINTS FROM RESAMPLING THE DATA
        bootind = fix(randomu(seed, npts)*npts)

;       THE DATA FROM THE BOOTSTRAP
        bootx = xuse[bootind]
        booty = yuse[bootind]
        bootv = vuse[bootind]
        boott = tuse[bootind]        

        bootmom = cloudmom(bootx, booty, bootv, boott, targett = targett)
        bootmom_noex = cloudmom(bootx, booty, bootv, boott $
                                , targett = targett, /noextrap)        

;       BOOTSTRAP THE ROTATED DATA, TOO
        bootx_rot = xrot[bootind]
        booty_rot = yrot[bootind]

        bootmom_rot = cloudmom(bootx_rot, booty_rot $
                               , bootv, boott, targett = targett)
        bootmom_noex_rot = cloudmom(bootx_rot, booty_rot, bootv, boott $
                                    , targett = targett, /noextrap)        
        
        if j eq 0 then begin
          bootmomra = bootmom 
          bootmomra_noex = bootmom_noex
          bootmomra_rot = bootmom_rot
          bootmomra_noex_rot = bootmom_noex_rot
        endif else begin
          bootmomra = [bootmomra, bootmom]        
          bootmomra_noex = [bootmomra_noex, bootmom_noex]
          bootmomra_rot = [bootmomra_rot, bootmom_rot]        
          bootmomra_noex_rot = [bootmomra_noex_rot, bootmom_noex_rot]
        endelse
      endfor
    endif

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; THEN TURN PIXEL-WISE PROPERTIES INTO PHYSICAL PROPERTIES.
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

;   GIVE THE PROPERTIES AN EMPTY STRUCTURE
    cprops = origcprops

;   NUMBER OF PIXELS IN THE CLOUD
    cprops.npix = n_elements(xuse) 

;   LOCATION OF THE DATA
    if n_elements(filename) gt 0 then $
      cprops.name = filename

;   NUMBER OF THIS CLOUD IN THE DATASET
    cprops.cloudnum = (i+1)

;   DISTANCE TO THE CLOUD IN PARSECS
    if (not (usekindist)) then $
      cprops.distance_pc = dist
    
;   INPUT THE KNOWN FWHM BEAM SIZE (IN PC) INTO THE CPROPS STRUCTURE    
    if (not (usekindist)) then $
      cprops.beamfwhm_pc = beamfwhm/3600.*!dtor*dist ; PARSECS

;   FIRST CALCULATE PHYSICAL RESOLUTIONS (X, Y, V PIXEL SIZE IN PC &
;   KM/S)
    if (not (usekindist)) then begin
      deltax_pc = abs(raxis[1]-raxis[0])* $
        cos(!dtor*mean(daxis))*!dtor*dist
      deltay_pc = abs(daxis[1]-daxis[0])*!dtor*dist
      deltamaj_pc = deltay_pc
      deltamin_pc = deltay_pc
    endif
    deltav_kms = abs(vaxis[1]-vaxis[0])
    cprops.sigchan_kms = deltav_kms/sqrt(2*!pi) ; CHECK THIS!!!!!!!

;   STORE THE PIXEL RESULTS FOR THE EXTRAP AND NOEXTRAP MOMENTS
    cprops.momxpix = mom.rmsx
    cprops.momxpix_noex = mom_noex.rmsx
    cprops.momypix = mom.rmsy
    cprops.momypix_noex = mom_noex.rmsy
    cprops.momvpix = mom.rmsv
    cprops.momvpix_noex = mom_noex.rmsv
    cprops.mommajpix = mom_rot.rmsx
    cprops.mommajpix_noex = mom_noex_rot.rmsx
    cprops.momminpix = mom_rot.rmsy
    cprops.momminpix_noex = mom_noex_rot.rmsy
    cprops.posang = pa

;   CALCULATE THE PIXEL-WISE AND ON-THE-SKY CENTERS (NO EXTRAPOLATION
;   IS NEEDED OR APPLIED)
    cprops.xctr_pix = total(tuse*xuse, /NAN)/total(tuse, /NAN)
    cprops.yctr_pix = total(tuse*yuse, /NAN)/total(tuse, /NAN)
    if n_elements(astrom) eq 0 then begin
      cprops.xctr_deg = total(tuse*raxis[xuse], /NAN)/total(tuse, /NAN)
      cprops.yctr_deg = total(tuse*daxis[yuse], /NAN)/total(tuse, /NAN)
    endif else begin
      xy2ad, cprops.xctr_pix, cprops.yctr_pix, astrom, $
             ra, dec
      cprops.xctr_deg = ra
      cprops.yctr_deg = dec
    endelse

    cprops.vctr_pix = total(tuse*vuse, /NAN)/total(tuse, /NAN)
    cprops.vctr_kms = total(tuse*vaxis[vuse], /NAN)/total(tuse, /NAN)

;   IF the KINDIST flag is set, calculate the kinematic distance to
;   each cloud.  This only works for cubes in l,b,v coordinates.
    if keyword_set(usekindist) then begin
      dist = kindist(cprops.xctr_deg, cprops.yctr_deg, $
                     cprops.vctr_kms, near = near, far = far, $
                     r0 = r0, v0 = v0)
      cprops.distance_pc = dist
      cprops.beamfwhm_pc = beamfwhm/3600.*!dtor*dist ; PARSECS
      deltax_pc = abs(raxis[1]-raxis[0])* $
        cos(!dtor*mean(daxis))*!dtor*dist
      deltay_pc = abs(daxis[1]-daxis[0])*!dtor*dist
      deltamaj_pc = deltay_pc
      deltamin_pc = deltay_pc
    endif

;   CONVERT THE MOMENTS TO PHYSICAL UNITS (PC, KM/S, K KM/S PC^2)
    rmsx_pc = cprops.momxpix*deltax_pc
    rmsy_pc = cprops.momypix*deltay_pc
    rmsmaj_pc = cprops.mommajpix*deltamaj_pc
    rmsmin_pc = cprops.momminpix*deltamin_pc
    rmsv_kms = cprops.momvpix*deltav_kms    
    flux_kkms_pc2 = mom.flux*deltav_kms*deltax_pc*deltay_pc

    rmsx_noex = cprops.momxpix_noex*deltax_pc
    rmsy_noex = cprops.momypix_noex*deltay_pc
    rmsmaj_noex = cprops.mommajpix_noex*deltamaj_pc
    rmsmin_noex = cprops.momminpix_noex*deltamin_pc
    rmsv_noex = cprops.momvpix_noex*deltav_kms
    flux_noex = mom_noex.flux*deltav_kms*deltax_pc*deltay_pc

;   CALCULATE QUANTITIES OF INTEREST:

;   ****************************************************************
;   ********************** THE RADIUS ******************************
;   ****************************************************************
;   1A) DECONVOLVED EFFECTIVE RADIUS (PC)

    cprops.rad_pc = sqrt(sqrt(rmsmaj_pc^2-(cprops.beamfwhm_pc/2.354)^2)*$
                         sqrt(rmsmin_pc^2-(cprops.beamfwhm_pc/2.354)^2))*$
      cprops.rmstorad

;    cprops.rad_pc = $
;      cprops.rmstorad*sqrt(rmsx_pc^2 + rmsy_pc^2 $
;                           - 2.*(cprops.beamfwhm_pc/2.354)^2.)*sqrt(0.5)


;   1B) UNEXTRAPOLATED RADIUS
    cprops.rad_noex = sqrt(sqrt(rmsmaj_noex^2-(cprops.beamfwhm_pc/2.354)^2)*$
                           sqrt(rmsmin_noex^2-(cprops.beamfwhm_pc/2.354)^2))*$
      cprops.rmstorad
    
;    cprops.rad_noex = $
;      cprops.rmstorad*sqrt(rmsx_noex^2 + rmsy_noex^2 $
;                           - 2.*(cprops.beamfwhm_pc/2.354)^2.)*sqrt(0.5)

;   1C) NOT-DECONVOLVED EFFECTIVE RADIUS (PC)
    cprops.rad_nodc = $
      cprops.rmstorad*sqrt(rmsmaj_pc*rmsmin_pc)

;    cprops.rad_nodc = $
;      cprops.rmstorad*sqrt(rmsx_pc^2 + rmsy_pc^2)*sqrt(0.5)

;   1D) VIRGINAL EFFECTIVE RADIUS (NO DC, NO EX) (PC)
    cprops.rad_nodc_noex = $
      cprops.rmstorad*sqrt(rmsmaj_noex*rmsmin_noex)

;    cprops.rad_nodc_noex = $
;      cprops.rmstorad*sqrt(rmsx_noex^2 + rmsy_noex^2)*sqrt(0.5)

;   ****************************************************************
;   ********************** THE LINE WIDTH **************************
;   ****************************************************************
;   2A) RMS LINEWIDTH (KM/S)
    cprops.sigv_kms = sqrt(rmsv_kms^2 - cprops.sigchan_kms^2)

;   2B) UNEXTRAPOLATED LINEWIDTH (KM/S)
    cprops.sigv_noex = sqrt(rmsv_noex^2 - cprops.sigchan_kms^2)

;   2C) NOT DECONVOLVED LINEWIDTH (KM/S)
    cprops.sigv_nodc = rmsv_kms

;   2D) VIRGINAL LINEWIDTH (NO DC, NO EX) (KM/S)
    cprops.sigv_nodc_noex = rmsv_noex

;   ****************************************************************
;   ********************** THE LUMINOSITY **************************
;   ****************************************************************
;   3A) FLUX (K KM/S)
    cprops.flux_kkms_pc2 = flux_kkms_pc2

;   3B) UNEXTRAPOLATED FLUX
    cprops.flux_noex = flux_noex

;   ****************************************************************
;   ********************** DERIVED MASSES **************************
;   ****************************************************************
;   4) VIRIAL MASS (MSUN)
    cprops.mvir_msun = $
      1040.*cprops.rad_pc*cprops.sigv_kms^2

;   5) "LUMINOUS MASS" -- XCO INFERRED MASS W/ HELIUM (MSUN)
    cprops.mlum_msun = $
      cprops.flux_kkms_pc2*(xco*(2.*mh)*(1.36)*(pc*pc)/ms)

;   6) MAX and SIGNAL to NOISE
      cprops.tmax_k = max(tuse, /nan)
      if n_elements(err) gt 0 or n_elements(sigma) gt 0 then begin
        if keyword_set(err) then sigma_cld = median(err[useind]) else $
          sigma_cld = sigma
        cprops.s2n = cprops.tmax_k/sigma_cld
      endif


;   6) AREA
;    index_2d = xuse+yuse*sz[1]
;    npix = intarr(n_elements(tuse))
;    sind = reverse(sort(tuse))
;    tvec = tuse[sind]
;     for kk = 0, n_elements(tuse)-1 do npix[kk] = n_elements(uniq(index_2d[sind[0:kk]], sort(index_2d[sind[0:kk]])))
;     cprops.area_pc2 = extrap(tvec,npix, targett = targett)*deltax_pc*deltay_pc
;     cprops.area_noex = npix[kk-1]*deltax_pc*deltay_pc
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
;   IF WE ARE BOOTSTRAPPING, ASSIGN UNCERTAINTIES FOR EACH MOMENT
;   BASED ON THE BOOTSTRAPPING
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

    if (n_elements(bootstrap) gt 0) then begin

;     THE FACTOR TO SCALE THE ERRORS BY TO ACCOUNT FOR CORRELATION IN
;     DATA DUE TO THE OVERSAMPLED BEAM
      indfac = sqrt(cprops.ppbeam)

;     ****************************************************************
;     UNCERTAINTIES IN THE RAW PIXEL VALUES, THESE ARE ALL PERCENTAGES
;     ****************************************************************
      cprops.momxpix_uc = indfac*mad(bootmomra.rmsx)/median(bootmomra.rmsx)
      cprops.momypix_uc = indfac*mad(bootmomra.rmsy)/median(bootmomra.rmsy)
      cprops.momvpix_uc = indfac*mad(bootmomra.rmsv)/median(bootmomra.rmsv)
      cprops.flux_uc = indfac*mad(bootmomra.flux)/median(bootmomra.flux)

      cprops.momxpix_noex_uc = indfac*mad(bootmomra_noex.rmsx) $
        /median(bootmomra_noex.rmsx)
      cprops.momypix_noex_uc = indfac*mad(bootmomra_noex.rmsy) $
        /median(bootmomra_noex.rmsy)
      cprops.momvpix_noex_uc = indfac*mad(bootmomra_noex.rmsv) $
        /median(bootmomra_noex.rmsv)
      cprops.flux_noex_uc = indfac*mad(bootmomra_noex.flux) $
        /median(bootmomra_noex.flux)      

;     ****************************************************************
;     UNCERTAINTIES IN THE RADIUS
;     ****************************************************************

;     THE NOMINAL RADIUS
      bootrad = cprops.rmstorad*sqrt((bootmomra.rmsx*deltax_pc)^2 + $
                                     (bootmomra.rmsy*deltay_pc)^2 $
                                     - 2.*(cprops.beamfwhm_pc/2.354)^2.)*sqrt(0.5)
      cprops.rad_uc = indfac*mad(bootrad)/median(bootrad)

;     THE NOT DECONVOLVED RADIUS
      bootrad_nodc = cprops.rmstorad*sqrt((bootmomra.rmsx*deltax_pc)^2 + $
                                          (bootmomra.rmsy*deltay_pc)^2)*sqrt(0.5)
      cprops.rad_nodc_uc = indfac*mad(bootrad_nodc)/median(bootrad_nodc)

;     THE NOT EXTRAPOLATED RADIUS
      bootrad_noex = cprops.rmstorad*sqrt((bootmomra_noex.rmsx*deltax_pc)^2 + $
                                          (bootmomra_noex.rmsy*deltay_pc)^2 - $
                                          2.*(cprops.beamfwhm_pc/2.354)^2) $
        *sqrt(0.5)
      cprops.rad_noex_uc = indfac*mad(bootrad_noex)/median(bootrad_noex)

;     THE NOT DECONVOLVED NOT EXTRAPOLATED RADIUS
      bootrad_nodc_noex = $
        cprops.rmstorad*sqrt((bootmomra_noex.rmsx*deltax_pc)^2 + $
                             (bootmomra_noex.rmsy*deltay_pc)^2) $
        *sqrt(0.5)
      cprops.rad_nodc_noex_uc = indfac*mad(bootrad_nodc_noex) $
        /median(bootrad_nodc_noex)      

;     ****************************************************************
;     THE VELOCITY LINE WIDTH
;     ****************************************************************

;     THE NOMINAL LINE WIDTH
      bootsigvkms = sqrt((bootmomra.rmsv*deltav_kms)^2 - $
                         (cprops.sigchan_kms)^2)
      cprops.sigv_uc = indfac*mad(bootsigvkms)/median(bootsigvkms)

;     THE NOT EXTRAPOLATED LINEWIDTH      
      bootsigvkms_noex = sqrt((bootmomra_noex.rmsv*deltav_kms)^2 - $
                              (cprops.sigchan_kms)^2)
      cprops.sigv_noex_uc = indfac*mad(bootsigvkms_noex)/median(bootsigvkms_noex)

;     THE NOT DECONVOLVED LINEWIDTH
      bootsigvkms_nodc = bootmomra.rmsv*deltav_kms
      cprops.sigv_nodc_uc = indfac*mad(bootsigvkms_nodc)/median(bootsigvkms_nodc)

;     THE NOT DECONVOLVED NOT EXTRAPOLATED LINEWIDTH
      bootsigvkms_nodc_noex = bootmomra_noex.rmsv*deltav_kms
      cprops.sigv_nodc_noex_uc = $
        indfac*mad(bootsigvkms_nodc_noex)/median(bootsigvkms_nodc_noex)

;     ****************************************************************
;     THE DERIVED MASSES
;     ****************************************************************

;     LUMINOUS MASS UNCERTAINTY
      cprops.mlum_uc = cprops.flux_uc

;     VIRIAL MASS UNCERTAINTY
      bootvmass = 1040.*bootrad*bootsigvkms^2.
      cprops.mvir_uc = indfac*mad(bootvmass)/median(bootvmass)

    endif
    
;   ADD THE CLOUD PROPERTIES TO THE ARRAY OF CLOUD PROPERTIES FOR THIS CUBE
    if (n_elements(cpropsra) eq 0) then $
      cpropsra = cprops else $
      cpropsra = [cpropsra, cprops]

  endfor

; RETURN THE ARRAY OF CLOUD PROPERTIES
  return, cpropsra

end                             ; OF CLOUDALYZE
