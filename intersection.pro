function intersection, a_in, b_in, count, disjoint = disjoint, $
  presort = presort, index_tag = index_tag, a_index = a_index, $
  b_index = b_index
;+
; NAME:
;    INTERSECTION
; PURPOSE:
;    To return the values common to two arrays.
;
; CALLING SEQUENCE:
;    vals = INTERSECTION(A, B [,count, disjoint=disjoint])
;
; INPUTS:
;    A,B -- Two arrays for comparison.
;
; KEYWORD PARAMETERS:
;    DISJOINT -- Values that do not appear in both sets.
;    PRESORT -- The input arrays are unique, sorted values already.
;    INDEX_TAG -- Passes where the elements appear in the original
;                 arrays into the variables, A_INDEX and B_INDEX
;
; OUTPUTS:
;    VALS -- The VALUES common to both arrays.
;    COUNT -- The number of common values.
;    a_index, b_index -- Indices of where the COMMON elements appear
;                        in arrays A and B respectively.
; NOTES: 
;    For arrays of length P and Q respectively, the operation performs
;    2(P+Q) comparisons and P+Q subtractions.  A naive comparison
;    will take P*Q/2 operations.  There probably remains a better way
;    of doing this, but this is a start.
;
;    
; MODIFICATION HISTORY:
;
;       Wed Dec 8 10:55:11 2004, Erik Rosolowsky <eros@cosmic>
;		Worked around the problems with SORT by tracking which
;		element came from which vector.
;
;       Mon Apr 12 10:34:02 2004, <eros@master>
;		Added the TAG_INDEX functionality
;
;       Tue Nov 4 15:47:33 2003, <eros@master>
;		Added compatibility for 1 element arrays.
;       
;       Included string compatility and disjoint keyword.
;       Fri Jun 20 16:28:42 2003, Erik Rosolowsky <eros@cosmic>
;		
;       Added
;       empty set check.  Mon Dec 16 13:08:09 2002, <eros@master>
;
;       Written --
;       Fri Dec 13 19:54:38 2002, <eros@master>
;-

  if keyword_set(index_tag) then begin
    tagalong1 = lindgen(n_elements(a_in))
    tagalong2 = lindgen(n_elements(b_in))
  endif 

  if n_elements(a_in) eq 0 or n_elements(b_in) eq 0 then begin
    count = 0
    return, !values.f_nan
  endif

  if not keyword_set(presort) then begin
    a = a_in[uniq(a_in, sort(a_in))]
    b = b_in[uniq(b_in, sort(b_in))]
    if keyword_set(index_tag) then begin
      tagalong1 = tagalong1[uniq(a_in, sort(a_in))]
      tagalong2 = tagalong2[uniq(b_in, sort(b_in))]
    endif
  endif else begin
    a = a_in
    b = b_in
  endelse

  merge = [a, b]
  if n_elements(tagalong1) gt 0 then begin
    tagmerge = [tagalong1, tagalong2]
    from_a = bytarr(n_elements(tagmerge))
    from_a[0:n_elements(tagalong1)-1] = 1b
  endif
  if n_elements(a) eq 1 and n_elements(b) eq 1 then begin
    if a eq b then begin
      ncount = 0
      count = 1
      return, a
    endif else begin
      count = 0
      ncount = 1
      return, !values.f_nan
    endelse
  endif
  sind = sort(merge)
  merge = merge[sind]
  if keyword_set(index_tag) then begin
    tagmerge = tagmerge[sind]
    from_a = from_a[sind]
  endif
  ind = where(merge eq shift(merge, 1), count)
  ind2 = where(merge ne shift(merge, 1) and merge ne shift(merge, -1), ncount)

  if keyword_set(index_tag) then begin
    if count eq 0 then begin
      a_index = -1
      b_index = -1
    endif else begin
      b_index = tagmerge[(ind-1)*(from_a[ind])+ind*(1b-from_a[ind])]
      a_index = tagmerge[ind*(from_a[ind])+(ind-1)*(1b-from_a[ind])]
    endelse
  endif
  if ncount gt 0 then disjoint = merge[ind2] else disjoint = !values.f_nan
  if count eq 0 then return, !values.f_nan else return, merge[ind]
end
