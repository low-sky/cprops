pro decomp_wrap, x, y, v, t, assignment, kernels = kernels $
                 , newkernels = newkernels $
                 , decomp = decomp, subcloud = subcloud $
                 , minpix = minpix $
                 , sigdiscont = sigdiscont $
                 , nredun = nredun $
                 , nodecomp = nodecomp, clfind = clfind $
                 , sigma = sigma, eclump = eclump, fscale = fscale $
                 , noextrap = noextrap, friends = friends $
                 , specfriends = specfriends, delta = delta, ppbeam = ppbeam

;+
; NAME:
;
;    DECOMP_WRAP
;
; PURPOSE:
;
;    THE WRAPPER FOR THE CLOUD DECOMPOSITION ROUTINES.
;
; CALLING SEQUENCE:
;    DECOMP_WRAP, x, y, v, t, assignment, subcloud =
;    subcloud, newkernels = newkernels
;     ... plus keyword listed below.
; INPUTS:
;   X,Y,V,T -- vectors containing the signal of the cube.
;   ASSIGNMENT -- vector containing the cloud (island) assignment of
;                 the data.  
; KEYWORD PARAMETERS:
;   /CLFIND -- Set this flag to run J. Williams' clumpfind on the
;              data.
;   /ECL -- Set this keyword to run E. Rosolowsky's modified clumpfind
;           on the data.
;   PPBEAM -- Input number of pixels per beam.
;   DECOMP -- Binary array with same number of elt. as X,Y,V,T that
;             has 1 if cloud is to be decomposed else 0.

;   DELTA, MINPIX -- passed to decimate_kernels()
;   NOEXTRAP -- Perform all moment analysis w/o extrapolation to 0 K.
;   
;   FRIENDS, SPECFRIENDS -- passed to alllocmax
;   SIGDISCONT, NREDUN, FSCALE -- Passed to deriv_decimate_kernels
;   KERNELS -- Defunct.  don't believe this keyword.
; OUTPUTS:
;   SUBCLOUD -- The assignment WITHIN a cloud that the pixel has been
;                assigned to.
;   NEWKERNELS -- Keyword named variable containing the new kernels of
;                 the data cube.
;   
; MODIFICATION HISTORY:
;       Wrote a blithe header.
;	Sat Sep  3 20:11:03 2005, Erik 
;
;-

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; DEFAULTS, DEFINITIONS, AND INPUTS
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

; FORWARD FUNCTION DECLARATION. IDL IS FLAKY.
  forward_function propdecomp

; DEFAULT TO DECOMPOSING EVERYTHING IF THE USER DOESN'T SUPPLY FLAGS.
  if n_elements(decomp) eq 0 then $
    decomp = bytarr(n_elements(x))+1b

; If CLUMPFIND is desired, check that we have an RMS estimate.
  if n_elements(sigma) eq 0 and keyword_set(clfind) then begin
    message, 'CLUMPFIND requires an RMS estimate.', /con
    return
  endif

; Create a vector of uniq assignment values, in essence a list of
; clouds.
  avals = assignment[uniq(assignment, sort(assignment))]

; Now make a vector containing clump (i.e. sub-cloud) assignments. All
; CLUMPs start out as assigned to 0.
  clump_asgn = lonarr(n_elements(x)) ;+1L

; Default to a property-based decomposition. Do this if no flags are
; set or more than one flag is set.
  if (keyword_set(nodecomp) + keyword_set(eclump) + $
      keyword_set(clfind)) ne 1 then $
        propdecomp = 1
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; LOOP OVER CLOUDS AND DECOMPOSE
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

; The main loop. Tackle each of the clouds in succession.
  for k = 0, n_elements(avals)-1 do begin
;   First extract the cloud from the cube (get X,Y,V,T)
    clindex = where(assignment eq avals[k])
    xcld = x[clindex]
    ycld = y[clindex]
    vcld = v[clindex]
    cloud = t[clindex]

;   Check whether the 'decomposition flag' is set, indicating that we
;   should attempt to decompose this cloud.
    decompflag = total(decomp[clindex] eq 1b) gt 0

    pad = 1 > sqrt(ppbeam/!pi) > friends

;     First make a cube out of the indices. It's a small cube, just
;     big enough to hold the data plus a little convenient padding.
        cubify, xcld, ycld, vcld, cloud, cube = minicube, $
                pad = pad, $
                indcube = indcube, indvec = clindex

    if decompflag then begin
;     -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
;     NOW DECOMPOSE THE CLOUD, IF DESIRED
;     -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

;     -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
;     1. FLAG SET TO FORCE NO DECOMPOSITION
;     -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
      if (keyword_set(nodecomp)) then begin
        print, 'No decomposition.'
      endif else if keyword_set(clfind) then begin
;     -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
;     2. WILLIAMS CLUMPFIND
;     -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+    
        print, 'Williams Clumpfind.'
        writefits, 'dummy_file.fits', minicube
        if n_elements(delta) gt 0 then sigma = delta/2
        clfind, file = 'dummy_file', low = 0.0, inc = 2*sigma
        dummy = readfits('dummy_file.fits.clf', hd)
        spawn, 'rm dummy_file.fits.clf', /sh
        spawn, 'rm dummy_file.fits', /sh
        flag = total(dummy ge 2, /nan) gt 0
        if flag then dummy = dummy+1*(dummy gt 0)
        indices = indcube[where(indcube ne -1)]
        clump_asgn[indices] = dummy[where(indcube ne -1)]
      endif else begin
;     
;     The remaining options for decomposition require a set of kernels
;     which are a subset of the local maxima.  So, let's get
;     started making some Local Maxima!


;     Determine a master set of contour values with the contour_values
;     function.
        levels = contour_values(cloud, nlevels = (n_elements(xcld)/50 > 250) < 1000)
;  Search for local maxima by comparing each datum with its neighbors.
        lmax = alllocmax(minicube, friends = friends, $
                         specfriends = specfriends)
        message, 'Number of Local Maxima Identified:'+$
                 strcompress(string(n_elements(lmax))), /con

; Only continue if there are 2+ maxima.
        if n_elements(lmax) gt 1 then begin
;  Pare down this kernels according to area and contrast restrictions.

          decimkern = decimate_kernels(lmax, minicube, $
                                     all_neighbors = all_neighbors $
                                     , delta = delta, sigma = sigma $
                                     , minpix = minpix);, levels = levels)
          message, 'Number of kernels after area and contrast decimation:'+$
                   strcompress(string(n_elements(decimkern))), /con

; Only continue if there are 2+ valid kernels. 
; AL --- SEEMS LIKE A BUG TO ME, I'M RENAMING KERNELS TO DECIMKERN
          if n_elements(decimkern) gt 1 then begin
; Calculate the merger matrix
            merger = mergefind(minicube, decimkern, levels = levels)
            kernels2 = deriv_decimate_kernel(minicube, merger, $
                                             decimkern $
                                             , levels = levels $
                                             , sigdiscont = sigdiscont $
                                             , nredun = nredun $ 
                                             , fscale = fscale)
            message, 'Number of kernels after derivative decimation:'+$
              strcompress(string(n_elements(kernels2))), /con

; Armed with a good set of kernels, let's get to decomposing.
            if n_elements(kernels2) gt 1 then begin
              if keyword_set(eclump) then begin
;     -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
;     3. ROSOLOWSKY CLUMPFIND
;     -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
                print, 'Rosolowsky Clumpfind.'
                print, 'Kernel Indices: ', kernels2
                cube_partition, minicube $
                                , level = 0.5*sigma $
                                , terminate = 0, astr = astr $
                                , kernels = kernels2, ppbeam = ppbeam
                sz = size(minicube)
                dummy = intarr(sz[1], sz[2], sz[3])
                dummy[astr.index] = astr.assignment+2
                indices = indcube[where(indcube ne -1)]
                clump_asgn[indices] = dummy[where(indcube ne -1)]
              endif else if keyword_set(propdecomp) then begin
;     -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
;     4. PROPERTY-BASED DECOMPOSITION
;     -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
                print, 'Property based decomposition.'
; Now go through and find all well-defined clouds with smooth behavior
; in the clouds.
                dummy = where(minicube eq minicube)
                vector_kernels2 = lonarr(n_elements(kernels2)) 
                for kk = 0, n_elements(vector_kernels2)-1 do $
                  vector_kernels2[kk] = where(dummy eq kernels2[kk])
                clump_asgn[clindex] = $
                  propdecomp(xcld, ycld, vcld, cloud, vector_kernels2 $
                             , merger = merger $
                             , levels = levels, all_neighbors = all_neighbors)
              endif             ; End of PROPDECOMP check
            endif               ; Check number after DERIV_DECIMATE
          endif                 ; Check number of KERNELS 
        endif                   ; Check number of LMAX 
      endelse                   ; Rosolowsky CLUMPFIND or PROPDECOMP
    endif                       ; if decompflag


;     IF *ONLY* WATERSHED IS RETURNED THEN GO AHEAD AND ZERO
;     EVERYTHING (I.E. TREAT IT AS ONE CLOUD).
    if (total(clump_asgn[clindex] ne 1) eq 0) then begin
      print, 'All watershed. Defaulting to one cloud.'
      clump_asgn[clindex] = 0
    endif
  endfor                        ; of loop over clouds

  if n_elements(decimkern) eq 0 and n_elements(kernels) gt 0 then $
    decimkern = kernels
  if n_elements(kernels2) gt 0 then $
    newkernels = kernels2 else if n_elements(decimkern) gt 0 then $
    newkernels = decimkern

  subcloud = clump_asgn

  return
end

