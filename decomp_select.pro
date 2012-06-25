pro decomp_select, cube, mask, all_neighbors = all_neighbors, $
                   minpix = minpix, $
                   x = x, y = y, v = v, t = t, index = index $
                   , decompose = decomp, assignment = asgn 
;+
;
; NAME:
;    DECOMP_SELECT
; PURPOSE:
;    Select candidate clouds for decomposition into constituent GMCs.  
;
; CALLING SEQUENCE:
;    DECOMP_SELECT, cube, mask
;
; INPUTS:
;    CUBE -- Original data cube
;    MASK -- Binary mask of emission in the data set.
;
; KEYWORD PARAMETERS:
;   ALL_NEIGHBORS -- Set this keyword to make pixels have 26 neighbors
;                    and not 6.
;   MINPIX -- The minimum MEAN number of pixels per kernel in the
;              decomposition.  Defaults to zero, i.e. no minimum.
;
; OUTPUTS:  Set keywords to variable names that will contain:
;
;   X,Y,V,T -- The position and data value for each cloud
;
;   ASSIGNMENT -- Integer assignment to of pixel to cloud N
;
;   DECOMPOSE -- A binary flag for each pixel: 1 to decompose, 0 to
;                not.
;
;   KERNELS -- Indices of X,Y,V,T to indicate the coordinates of the
;              local maxima
;
;   INDEX -- The 1D index in the original array CUBE of the clouds in
;            question.
;
; REQUIRES:
;    ALLLOCMAX, MERGEFIND, DECIMATE_KERNELS
;
; MODIFICATION HISTORY:
;
;       Mon Dec 6 12:53:15 2004, Erik Rosolowsky <eros@cosmic>
;		Written.
;
;-

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; HANDLE INPUTS AND DEFAULTS
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

; DEFAULT TO NO PIXEL MINIMUM
  if n_elements(minpix) eq 0 then $
    minpix = 0.0

; SET DEFAULT MINIMUM DYNAMIC RANGE TO BE A FACTOR OF TWO SIGMA
  if n_elements(delta) eq 0 then $
    delta = 2.0
  
; SIZE OF THE DATA CUBE
  sz = size(cube)

; GET THE MEAN ABSOLUTE DEVIATION OF DATA IN THE CUBE
  sigma = mad(cube)


; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; LOCATE LOCAL MAXIMA
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

; 1. LABEL SIMPLY CONNECTED REGIONS
  regions = label_region(mask, all_neighbors = all_neighbors)

; 2. HISTOGRAM REGION-LABELS TO GET A COUNT OF ELEMENTS IN EACH CLOUD
  h = histogram(regions, binsize = 1, omax = omax, $
                omin = omin, reverse = ri)

; 3. TOTAL NUMBER OF PIXELS INCLUDED IN REGION LABELS  
  nelts = total(h[1:*])

; 4. AN ARRAY OF INDICES CORRESPONDING TO LABELLED PIXELS IN THE DATA
; CUBE
  index = where(regions)

; 5. THE X, Y, V, T, AND ASSIGNMENT VALUES ASSOCIATED WITH EACH
; LABELLED PIXEL
  coordinates = transpose(array_indices(cube, index))
  x =  coordinates[*, 0]
  y =  coordinates[*, 1]
  v =  coordinates[*, 2]
  t =  cube[index]
  asgn = regions[index]

; 6. A BINARY FLAG INDICATING WHETHER THE CLOUD SHOULD BE DECOMPOSED
  decomp = bytarr(n_elements(x))+1b

  print, 'Identified ', omax, ' island(s).'

; 7. LOOP OVER ALL CLOUDS AND FLAG THEM TO BE EITHER DECOMPOSED OR NOT
; DEPENDING ON CLOUD SIZE AND CONTRAST BETWEEN PEAK AND SHARED
; CONTOURS

  for k = 1, omax do begin
;   GET THE PIXELS ASSOCIATED WITH CLOUD k
    mindex = where(asgn eq k, count)
    if count lt 2*minpix then  decomp[mindex] = 0b
  endfor
  return
end





