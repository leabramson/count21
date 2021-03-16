function trans2020official, struct, $
                            wts = wts

  if not keyword_set(wts) then $
     wts = [1.51,1.77,1.42,1.48,1.68]               ;; 2020 Spa4 interp into CD13

  if total(strcompress(tag_names(struct),/rem) eq 'C') then begin
     oth = transpose([[struct.C/wts[0]], $
                      [struct.V/wts[1]], $
                      [struct.R/wts[2]], $
                      [struct.T/wts[3]], $
                      [struct.M/wts[4]]])
  
     tot = struct.PERSONS + total(oth, 1)
  endif else begin
     oth = (struct.TOTAL - struct.PERSONS) / mean(wts) ;; approximation here
     tot = struct.PERSONS + oth
  endelse
     
  output = {TRACT: '0000.00', $
            TOT_PPL: 0., $
            TOT_OBJ: 0., $
            TOT_IND: 0., $
            C: 0., $
            V: 0., $
            R: 0., $
            T: 0., $
            M: 0.}
  output = replicate(output, n_elements(Struct))  
  
  for ii = 0, n_elements(struct) - 1 do begin
     output[ii].TRACT   = struct[ii].TRACT
     output[ii].TOT_PPL = struct[ii].TOTAL
     output[ii].TOT_OBJ = tot[ii]
     output[ii].TOT_IND = struct[ii].PERSONS
     if total(strcompress(tag_names(struct),/rem) eq 'C') then begin
        output[ii].C       = oth[0,ii]
        output[ii].V       = oth[1,ii]
        output[ii].R       = oth[2,ii]
        output[ii].T       = oth[3,ii]
        output[ii].M       = oth[4,ii]
     endif
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

pro fitsCSVactual

  readcol, '2020actualCounts.csv', tract,pop,p,c,v,r,t,m,cts, $
           f = 'A,F,F,F,F,F,F,F,F'
  tract = strmid(tract,0,4)+'.'+strmid(tract,1,/rev)

  savedata = {TRACT: '0', $
              PERSONS: 0., $
              C: 0., $
              V: 0., $
              R: 0., $
              T: 0., $
              M: 0., $
              TOTAL: 0.}
  savedata = replicate(saveData, n_elements(tract))
  for ii = 0, n_elements(tract) - 1 do begin
     savedata[ii].TRACT   = tract[ii]
     savedata[ii].TOTAL   = pop[ii]
     savedata[ii].PERSONS = p[ii]
     savedata[ii].C = c[ii]
     savedata[ii].V = v[ii]
     savedata[ii].R = r[ii]
     savedata[ii].T = t[ii]
     savedata[ii].M = m[ii]
  endfor
  mwrfits, savedata, 'official2020CompleteOccupancies.fits', /create
     
end

pro fitsCSVredux

  readcol, '2020actualCountsW191902.csv', tract,pop,p,c,v,r,t,m,cts, $
           f = 'A,F,F,F,F,F,F,F,F'
  tract = strmid(tract,0,4)+'.'+strmid(tract,1,/rev)

  savedata = {TRACT: '0', $
              PERSONS: 0., $
              C: 0., $
              V: 0., $
              R: 0., $
              T: 0., $
              M: 0., $
              TOTAL: 0.}
  savedata = replicate(saveData, n_elements(tract))
  for ii = 0, n_elements(tract) - 1 do begin
     savedata[ii].TRACT   = tract[ii]
     savedata[ii].TOTAL   = pop[ii]
     savedata[ii].PERSONS = p[ii]
     savedata[ii].C = c[ii]
     savedata[ii].V = v[ii]
     savedata[ii].R = r[ii]
     savedata[ii].T = t[ii]
     savedata[ii].M = m[ii]
  endfor
  mwrfits, savedata, 'official2020CompleteOccupanciesW191902.fits', /create
     
end
