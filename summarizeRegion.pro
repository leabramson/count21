function summarizeRegion, resStruct

  cts = resStruct.COUNTS
  cts = total(total(cts, 3),1)
  cts = cts[sort(cts)]

  pctles = [0.05,0.16,0.25,0.50,0.75,0.84,0.95]
  
  RETURN, cts[ceil(pctles*n_elements(cts))-1]
end
