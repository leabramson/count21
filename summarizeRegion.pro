function summarizeRegion, resStruct

  if n_elements(resStruct) gt 1 then begin
     cts = resStruct.COUNTS
     cts = total(total(cts, 3),1)
  endif else begin
     cts = total(resStruct.COUNTS, 1)
  endelse
  cts = cts[sort(cts)]

  pctles = [0.05,0.16,0.25,0.50,0.75,0.84,0.95]
  
  RETURN, cts[ceil(pctles*n_elements(cts))-1]
end
