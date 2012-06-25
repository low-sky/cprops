function dilator3d,  image, kernel, constraint = constraint, $
                     ppbeam = ppbeam, vwt = vwt, $
                     distance = distance
;+
; NAME:
;  DILATOR3D
;
; PURPOSE:
;  This routine dilates an assignment mask out from kernels through an
;  entire data cube, assigning conflicts based on geometric proximity
;  to the original kernel.
;
; CALLING SEQUENCE:
;   assignment = DILATOR3D( cube, kernels, [ constraint = constraint,
;   ppbeam = ppbeam, vwt = vwt, distance = distance])
;
; INPUTS:
;   CUBE -- A data cube
;   KERNELS -- A vector of indices within image that are the beginning
;              of the dilation.
; KEYWORD PARAMETERS:
;   CONSTRAINT -- a mask of emission to consider.
;   PPBEAM -- The number of pixels per beam
;   VWT -- The relative weight of the velocity distance in breaking
;          conflicted ownership of a pixel.
;   DISTANCE -- Named variable that returns the distance (along a
;               path) to the kernel that owns a pixel.
;
; OUTPUTS:
;   ASSIGNMENT -- An assignment cube with unique integer labels for
;                 each unique cloud.
;
; MODIFICATION HISTORY:
;       Gutted & Documented
;	Sat Sep  3 15:55:53 2005, Erik
;
;-

  
  if n_elements(vwt) eq 0 then vwt = 1.0
; Set this keyword to loop until saturated.
  if keyword_set(loop) then looper = 1 else looper = 0

; Input image with pixel assignments and indices of valid maxima
; Output: dilation with conflicts assigned to nearest kernels.

; Error checking...
  sz = size(image)
  if sz[0] ne 3 then begin
    message, 'This is for 3D data.  Try DILATOR.PRO.', /con
    return, 0
  endif
  if n_elements(distance) eq 0 then begin
    distance = intarr(sz[1], sz[2], sz[3])+32767
    ind = where(image, ct)
    if ct gt 0 then distance[where(image)] = 1
    distance[kernel] = 0
  endif

  if n_elements(ppbeam) gt 0 then begin
    radbeam = floor(sqrt(ppbeam/!pi)) > 1
    structure = bytarr(2*radbeam+1, 2*radbeam+1, 3)
    structure[radbeam, radbeam, 0] = 1b
    structure[radbeam, radbeam, 2] = 1b
    structure[*, *, 1] = shift(dist(2*radbeam+1, 2*radbeam+1), $
                              radbeam, radbeam) le radbeam
  endif else begin
; Default to basic cross (dilation is same in all directions)
    structure = bytarr(3, 3, 3)
    structure[1, 1, 0] = 1b
    structure[1, 1, 2] = 1b
    structure[0:2, 1, 1] = 1b
    structure[1, 0:2, 1] = 1b
  endelse 


  if n_elements(kernel) gt 0 then begin
    x_kernel = kernel mod sz[1]
    y_kernel = (kernel mod (sz[1]*sz[2]))/sz[1]
    v_kernel = kernel/(sz[1]*sz[2])
  endif

  if n_elements(constraint) eq 0 then  imout = image else $
    imout = image*constraint

  iter = 0
  finishing = 0
  repeat begin
    iter = iter+1
    if n_elements(new_m) eq 0 then m = (image gt 0) else m = imout gt 0
; Dilate the Structure
    new_m = dilate(m, structure)
; Mask to current region
    if n_elements(constraint) gt 0 then new_m = new_m*constraint
; Find new points in the mask 
; Sometimes m[x]=1 but new_m[x]=0.  Why?  I dunno...
    additions = new_m-m
    ind_adds = where(additions eq 1, break_ct)
; Start assigning if there are new points.
    if break_ct gt 0 then begin
; If there are new pixels, we aren't finishing.
      finishing = 0 
; Get the indices of neighbors to the new pixels
      nbrs_adds = replicate(1, 18)#ind_adds+$
                  [sz[1], -sz[1], sz[1]+1, sz[1]-1, $
                   +1, -1, sz[1]*sz[2]+1, $
                   sz[1]*sz[2]-1, -sz[1]*sz[2]+1, -sz[1]*sz[2]-1, $
                   sz[1]*sz[2]+sz[1], sz[1]*sz[2]-sz[1], $
                   sz[1]*sz[2], -sz[1]*sz[2], $
                   -sz[1]*sz[2]+sz[1], -sz[1]*sz[2]-sz[1], $
                   -sz[1]-1, -sz[1]+1]#$
                  replicate(1, n_elements(ind_adds))
; Get the assignments of the neighbors
      vals = imout[nbrs_adds]
      dmat = distance[nbrs_adds]
      nbval = lonarr(1, n_elements(ind_adds)) ; Value of NBR
      nnbrs = intarr(1, n_elements(ind_adds)) ; Num. of NBRs
      nbdist = intarr(1, n_elements(ind_adds))+32767 

      for k = 0, 17 do begin
        nnbrs = nnbrs+(vals[k, *] ne nbval and vals[k, *] ne 0)
        nbval = nbval > vals[k, *]
        nbdist = nbdist < dmat[k, *]
      endfor
; Note that NNBRS is not the number of neighbors (precisely) but is
; bigger than 1 if there is more than one neighbor.
      if n_elements(kernel) gt 0 then begin
; If there are kernel points to use for partitioning use them.
; Otherwise just let the local max be the new assignment.
        ind = where(nnbrs gt 1, ct)
        if ct gt 0 then begin
          xpos = ind_adds[ind] mod sz[1]
          ypos = (ind_adds[ind] mod (sz[1]*sz[2]))/sz[1]
          vpos = ind_adds[ind]/(sz[1]*sz[2])
          for k = 0L, ct-1 do begin
; This finds all possible kernels that are path-connected to the pixel
; in question.
            possible = vals[*, ind[k]]
            i1 = where(possible gt 0)
            possible = possible[i1]
            i2 = uniq(possible, sort(possible))
            possible = possible[i2]

; Establish the distance metric to those kernels.
            darr = sqrt((xpos[k]-x_kernel[possible])^2+$
                        (ypos[k]-y_kernel[possible])^2+$
                        vwt*(vpos[k]-v_kernel[possible])^2)
            null = min(darr, winner)
; Assign to the "closest"
            imout[ind_adds[ind[k]]] = possible[winner]
            i3 = where(vals[*, ind[k]] eq possible[winner])
            distance[ind_adds[ind[k]]] = min(dmat[i3, ind[k]])+1
          endfor
        endif
; If there's only 1 nbr, then give the new pixel the neighbor's assgn.
        ind = where(nnbrs eq 1, ct)
        if ct gt 0 then begin 
          imout[ind_adds[ind]] = nbval[ind]
          distance[ind_adds[ind]] = nbdist[ind]+1
        endif
; Ignore kernels with 0 nbrs.
        ind = where(nnbrs eq 0, ct)
        break_ct = break_ct-ct
      endif else begin
; If no kernels just give the max neighboring value.
        imout[ind_adds] = nbval
      endelse
    endif


  endrep until (break_ct eq 0)

  return, imout
end
