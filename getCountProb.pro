;; use a results fits file; output from countAll or countRegion

function getCountProb, counts, values, $
                       inverse = inverse, $
                       print = print

  cts = counts[sort(counts)]
  ncts = n_elements(cts)
  y = findgen(ncts)/(ncts-1)

  if NOT keyword_set(INVERSE) then $
     probs = y[value_locate(cts, values)] $
  else $
     probs = cts[value_locate(y, values)]

  if keyword_set(PRINT) then begin
     print, f = '(%"Median count value: %i")', cts[ceil(ncts/2)-1]
     if NOT keyword_set(INVERSE) then begin
        for ii = 0, n_elements(values) - 1 do $
           print, f = '(%"Prob < %i: %4.2f")', values[ii], probs[ii]
     endif else begin
        for ii = 0, n_elements(values) - 1 do $
           print, f = '(%"Value @ p = %4.2f: %i")', values[ii], probs[ii]
     endelse
  endif
  
  RETURN, probs
end
