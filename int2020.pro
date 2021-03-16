pro int2020

  data = mrdfits('official2020completeOccupanciesW191902.fits', 1)
  d = trans2020official(data)

  eh = fltarr(n_elements(data))
  for ii = 0, n_elements(data) - 1 do $
     eh[ii] = eHoLookup(data[ii].TRACT)
  eho = where(eh, compl = ho)
  
  ;; weighted
  pop = transpose([[data.PERSONS],$
                   [data.C],[data.V],[data.R], $
                   [data.T],[data.M]])
  tepop = total(pop[*,eho], 2)
  thpop = total(pop[*,ho], 2)
  
  ;; unweighted
  cts = transpose([[d.TOT_IND], $
                   [d.C],[d.V],[d.R], $
                   [d.T],[d.M]])
  tects = total(cts[*,eho], 2)
  thcts = total(cts[*,ho], 2)

  ;; official community summaries
  oeppl =  656. * [0.25,0.067,0.157,0.024,0.212,0.29]
  ohppl = 1058. * [0.388,0.078,0.08,0.043,0.31,0.101]

  oects = [656.*0.25  ,29,58,11,94,113]
  ohcts = [1058.*0.388,55,48,32,222,64]

  print, 'HOLLYWOOD'
  print, 'ppl', thpop - ohppl
  print, 'cts', thcts - ohcts
;  print, '%cts', (thcts - ohcts) / total(ohcts)
  print, ''
  print, 'EAST HOLLYWOOD'
  print, 'ppl', tepop - oeppl
  print, 'cts', tects - oects
;  print, '%cts', (tects - oects) / total(oects)
  
end
