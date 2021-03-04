function trans2020official, struct, $
                            wts = wts

  if not keyword_set(wts) then $
     wts = [1.51,1.77,1.42,1.48,1.68]               ;; 2020 Spa4 interp into CD13
  oth = (struct.TOTAL - struct.PERSONS) / mean(wts) ;; approximation here
  tot = struct.PERSONS + oth

  output = {TRACT: '0000.00', $
            TOT_OBJ: 0., $
            TOT_PPL: 0., $
            TOT_IND: 0.}
  output = replicate(output, n_elements(Struct))  
  
  for ii = 0, n_elements(struct) - 1 do begin
     output[ii].TRACT   = struct[ii].TRACT
     output[ii].TOT_OBJ = tot[ii]
     output[ii].TOT_PPL = struct[ii].TOTAL
     output[ii].TOT_IND = struct[ii].PERSONS
  endfor

  return, output
  
end

pro fitsCSV

  readcol, '2020impliedCounts.csv', tract, raw, tot, ind, $
           f = 'A,F,F,F'
  tract = strmid(tract,0,4)+'.'+strmid(tract,1,/rev)

  savedata = {TRACT: '0', $
              PERSONS: 0., $
              TOTAL: 0.}
  savedata = replicate(saveData, n_elements(tract))
  for ii = 0, n_elements(tract) - 1 do begin
     savedata[ii].TRACT   = tract[ii]
     savedata[ii].TOTAL   = tot[ii]
     savedata[ii].PERSONS = ind[ii]
  endfor
  mwrfits, savedata, 'official2020occupancies.fits', /create
     
end
