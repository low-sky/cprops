;+
; NAME:
;          COUNTER
; PURPOSE:
;          Print a progress status to the screen on a single line.
;
; CALLING SEQUENCE:
;          COUNTER, NUMBER, OUTOF  [,INFOSTRING, /PERCENT, WAIT_TIME=variable]
;
; INPUTS:
;          NUMBER:  The current number. Usually the loop index variable
;          OUTOF:   The total number of iterations. The last loop index
;
; OPTIONAL INPUTS:
;
;          INFOSTRING: A string telling the user what is being
;                      counted e.g. 'Flat '
;
; KEYWORD PARAMETERS:
;         
;          PERCENT: Set to output update in percent completed 
;                   percent = rount(number/outof) * 100
;
;          WAIT_TIME:  Used for test purposes for fast loops. Don't use
;                     this if time is of the essense.
; OUTPUTS:
;          Status is printed to the screen and updated on a single line.
;
; SIDE EFFECTS:
;          This program doesn't take much longer than a simple print
;          statement. But keep in mind that this amount of time is not
;          0 (zero). So be careful where you use this. If your loop
;          consists of only two relatively quick commands, updating
;          the status with this program could take up a significant
;          portion of the loop time!
;
; PROCEDURE:
;          Put counter statement inside your loop, preferably at the end.
;
; PROCEDURES CALLED:
;            
; EXAMPLE:
;          Try this to see how it works:
;
;          IDL> for i = 0,4 do counter,i,4,'test ',wait=.5
;
;
; MODIFICATION HISTORY:
;
;       Thu Aug 14 13:19:59 2003, Erik Rosolowsky <eros@cosmic>
;		Hacked to stand like a man!
;       Written by JohnJohn, Berkeley 06 January 2003 
;-
pro counter, num, outof, infostring $
             , wait_time = waittime, percent = percent

  on_error, 2
  clearline = string("15b)
  if n_elements(infostring) eq 0 then $
    st = 'Number '              ;set default info string
  
  lenst = strcompress(strlen(infostring), /rem)
  leni = strcompress(strlen(strcompress(num, /rem)), /rem)
  leno = strcompress(strlen(strcompress(outof, /rem)), /rem)

  if keyword_set(percent) then begin
    per = strcompress(round(float(num)*100./outof), /rem)
    lenp = str(strlen(str(per)))
    form = "($,a"+lenp+",' % Completed',a)"
    print, form = form, per, clearline
  endif else begin
    form = "($,a"+lenst+",i"+leni+",' out of ',i"+leno+",a)"
    print, form = form, infostring, num, outof, clearline
  endelse
  
  if n_elements(waittime) gt 0 then wait, waittime
end
