function assign2, mask, all_neighbors = all_neighbors, success = success
;+
; NAME:
;   ASSIGN2
; PURPOSE: 
;
;  USE THE LABEL REGIONS COMMAND TO SLICKLY PARTITION A MASK
;  INTO INDIVIDUAL REGIONS. RETURN THESE REGIONS IN AN "ASSIGNMENT"
;  STRUCTURE.
;
; CALLING SEQUENCE:
;   asgn = ASSIGN2(mask [,/all_neighbors, success = success])
;
; INPUTS:
;   MASK -- An integer mask where each distinct cloud has its own
;           integer label.
;
; KEYWORD PARAMETERS:
;   /ALL_NEIGHBORS -- Passed to label_region()
;   SUCCESS -- Named variable containing a 1 if successful or 0 otherwise
; OUTPUTS:
;   ASGN -- An assignment structure with the following tags:
;     X,Y,V,T -- Vectors of the X,Y,V positions and T amplitude of the
;                data
;     ASSIGNMENT -- Unique integer assignment (ZERO Indexed) for each
;                   cloud.
;     INDEX -- 1D index of position in original data cube.
;     XSZ, YSZ, VSZ -- Dimensions of original datacube.
; MODIFICATION HISTORY:
;       Documented.
;	Mon Sep  5 13:03:39 2005, Erik 
;
;-

  sz = size(mask)
  success = 1B

; ERROR CHECKING -- THE MASK MUST CONTAIN EMISSION (LABELLED AS 1's)
  If total(mask) eq 0 then begin
    message, 'Mask contains no valid clouds.', /con
    success = 0B
    return, {index:0, x:0, y:0, v:0, counts:0 $
             , assignment:0 $
             , xsz:sz[1], ysz:sz[2], vsz:sz[3]}
  endif

; USE "LABEL REGIONS" TO IDENTIFY SIMPLY CONNECTED
; REGIONS. ALL-NEIGHBORS ALLOWS DIAGONAL NEIGHBORING.
  assignment = label_region(mask, all_neighbors = all_neighbors, /ulong)
  
; ERROR CHECKING -- EVEN IF WE HAVE A NONZERO MASK WE MAY NOT HAVE ANY
;                   REGIONS B/C OF EDGE EFFECTS
  if (total(assignment) eq 0) then begin
    message, 'Mask contains no valid clouds.', /con
    success = 0B
    return, {index:0, x:0, y:0, v:0, counts:0 $
             , assignment:0 $
             , xsz:sz[1], ysz:sz[2], vsz:sz[3]}    
  endif
  
; BY HISTOGRAMING THE ASSIGNMENT, WE GET THE NUMBER OF PIXELS ASSIGNED
; TO EACH REGION.
  counts = histogram(assignment, min = 1) ;, reverse_indices = ri)

; TAKE ALL PIXELS WHICH HAVE BEEN ASSIGNED TO A CLOUD...
  index = where(assignment gt 0)
; ... GET THE X PIXEL FROM THE INDEX,
  x = index mod sz[1]
; AND THE Y PIXEL FROM THE INDEX
  y = (index mod (sz[1]*sz[2]))/sz[1]
; AND THE V PIXEL FROM THE INDEX
  v = index/(sz[1]*sz[2])
; AND FINALLY GET THE CLOUD ASSIGNMENT FROM THE ASSIGNMENT CUBE
; (DROPPING DOWN ONE SO THAT THE CLOUDS BEGIN AT CLOUD 0)
  assignment = assignment[index]-1

; ORGANIZE THE DATA INTO A STRUCTURE AND HAND IT BACK, NOW EACH PIXEL
; OF EMISSION HAS AN ENTRY WITH ITS X, Y, AND V LOCATION, AND HOW MANY
; COUNTS THAT CLOUD HAS, ALONG WITH THE ONE DIMENSIONAL INDEX.
  return, {index:index, x:x, y:y, v:v $
           , assignment:assignment, counts:counts $
           , xsz:sz[1], ysz:sz[2], vsz:sz[3]}

end
