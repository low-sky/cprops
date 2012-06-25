function reject_region, cube, oldmask $
                        , minarea = minarea, mindeltav = mindeltav $
                        , minpksig = minpksig, edge = edge

;+
;
; NAME:
;   REJECT_REGION
;
; PURPOSE:
;
;   Removes regions from a mask that are not interesting and returns a
;   pared-down mask. Masks a certain amount from the edge and then
;   rejects regions with insufficient area, velocity spanning, and
;   peak signal.
;
; CALLING SEQUENCE:
;
;   newmask = reject_region(cube, oldmask, MINAREA = minarea $,
;              MINDELTAV = mindeltav, MINPKSIG = minpksig, EDGE = edge)
;
; INPUTS:
;
;   CUBE      - input dataset
;   OLDMASK   - original mask (output from signal extraction algorithm)
;   MINAREA   - minimum area of an interesting region in pixels.
;   MINDELTAV - minimum velocity spanning of an interesting region
;               (pixels).
;   MINPKSIG  - minimum peak signal in an interesting region (data
;               values).
;   EDGE      - pixels from the edge to mask (should be >= 1 to make
;               LABEL_REGION work correctly, defaults to 1).
;
; KEYWORD PARAMETERS:
;
;
; OUTPUTS:
;
;   NEWMASK   - the new, pared-down mask.
;
; MODIFICATION HISTORY:
;
;-

; COPY THE OLD MASK TO THE NEW MASK!!!
  mask = oldmask

; GET THE SIZE VECTOR FOR THE CUBE
  sz = size(cube)

; SET THE REJECTION CRITERIA
  if (n_elements(minarea) eq 0) then $
    minarea = 1L

  if (n_elements(mindeltav) eq 0) then $
    mindeltav = 1L

  if (n_elements(minpksig) eq 0) then $
    minpksig = 0.0D

  if (n_elements(edge) eq 0) then $
    edge = 1L

; ZERO THE EDGE REGIONS
  if (edge gt 0) then begin
    mask[0:(edge-1), *, *] = 0B
    mask[*, 0:(edge-1), *] = 0B
    mask[*, *, 0:(edge-1)] = 0B
    mask[(sz[1]-edge):*, *, *] = 0B
    mask[*, (sz[2]-edge):*, *] = 0B
    mask[*, *, (sz[3]-edge):*] = 0B
  endif

; IDENTIFY REGIONS OF THE MASK USING THE assign2 FUNCTION ...
  astr = assign2(mask)

; CHECK THAT WE DO INDEED HAVE SOME CLOUDS DEFINED. IF NOT THEN RETURN
  if size(astr, /type) ne 8 then $
    return, mask


; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; CALCULATE PROPERTIES OF INTEREST FOR EACH CLOUD.
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

; CALCULATE THE AREA OF EACH CLOUD
  xyind = astr.x + astr.y*sz[1]
  area = lonarr(max(astr.assignment)+1)
  for i = 0L, max(astr.assignment) do begin
    ind = where(astr.assignment eq i)
    area[i] = n_elements(uniq(xyind[ind], sort(xyind[ind])))
  endfor

; CALCULATE THE VELOCITY SPAN OF EACH CLOUD
  deltav = lonarr(max(astr.assignment)+1)
  for i = 0L, max(astr.assignment) do begin
    ind = where(astr.assignment eq i)
    deltav[i] = n_elements(uniq(astr.v[ind], sort(astr.v[ind])))
  endfor

; CALCULATE THE PEAK SIGNAL FOR EACH CLOUD
  pksig = dblarr(max(astr.assignment)+1)
  for i = 0L, max(astr.assignment) do begin
    ind = where(astr.assignment eq i)
    pksig[i] = max(cube[astr.index[ind]], /NAN)
  endfor

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; PARE REGIONS BASED ON A SERIES OF CRITERIA, KEEPING ONLY REGIONS
; FOR WHICH WE MIGHT WANT TO MEASURE CLOUD PROPERTIES.
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

; LOOP OVER ALL CLOUDS
  for i = 0L, max(astr.assignment) do begin
    
; 1) REJECT ON AREA -- SPATIALLY UNRESOLVED REGIONS ARE UNINTERESTING
    if (area[i] lt minarea) then begin
      ind = where(astr.assignment eq i, num)
      mask[astr.x[ind], astr.y[ind], astr.v[ind]] = 0B   
    endif  

; 2) REJECT ON VELOCITY SPAN -- REGIONS NOT VELOCITY-RESOLVED ARE
;                               UNINTERESTING
    if (deltav[i] lt mindeltav) then begin
      ind = where(astr.assignment eq i, num)
      mask[astr.x[ind], astr.y[ind], astr.v[ind]] = 0B   
    endif  
    
; 3) REJECT ON PEAK SIGNAL -- INSIGNIFICANT PEAKS STINK.
    if (pksig[i] lt minpksig) then begin
      ind = where(astr.assignment eq i, num)
      mask[astr.x[ind], astr.y[ind], astr.v[ind]] = 0B   
    endif

  endfor 
  return, mask

end                             ; OF REJECT_REGION
