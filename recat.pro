function recat, assignment
;+
; NAME:
;    RECAT
; PURPOSE:
;    Reindices an assignment vector to start at 1 and increment by 1.
;
; CALLING SEQUENCE:
;   new = RECAT( old )
;
; INPUTS:
;   OLD -- An assignment vector.  Generally a vector of assignments
;          corresponding to an array of indicies.
;
; KEYWORD PARAMETERS:
;   NONE
;
; OUTPUTS:
;   NEW -- A recataloged assignment vector that starts with 1.
;
; MODIFICATION HISTORY:
;
;       Mon Jun 7 11:03:21 2004, <eros@master>
;		Written.
;
;-

  uniq_vals = assignment[uniq(assignment, sort(assignment))]
  new = assignment
  new[*] = 0
  for k = 0, n_elements(uniq_vals)-1 do begin
    index = where(assignment eq uniq_vals[k])
    new[index] = k+1
  endfor

  return, new
end
