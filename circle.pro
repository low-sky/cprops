pro circle, radius = radius, color = color, fill = fill, thick = thick
;+
; NAME:
;   CIRCLE 
; PURPOSE:
;   To make the plotting symbol a circle using USERSYM (psym = 8)
;
; CALLING SEQUENCE:
;   CIRCLE [, radius = radius, color = color, fill = fill, thick = thick]
;
; INPUTS:
;   None
;
; KEYWORD PARAMETERS:
;   RADIUS -- Radius of the circle in character sizes.
;   COLOR -- Color to draw the symbol in (or to FILL with)
;   FILL -- Fill the circle
;   THICK -- Thickness of the border to draw.
; OUTPUTS:
;
;
; MODIFICATION HISTORY:
;       Documented --
;       Wed Nov 21 11:27:36 2001, Erik Rosolowsky <eros@cosmic>
;-


if n_elements(thick) eq 0 then thick = !p.thick
if not keyword_set(radius) then radius = 1
phi = findgen(41)/40*!pi*2
usersym, radius*cos(phi), radius*sin(phi), fill = fill, thick = thick, $
         color = color

  return
end
