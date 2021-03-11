pro tallyTentPersonSwap

  data = mrdfits('countHollywoodResults2021.fits', 1)
  persons = total(data.RAWCOUNTS[0:2], 1)
  tents   = total(data.RAWCOUNTS[6:7], 1)
  vehicles = total(data.RAWCOUNTS[3:5], 1)
  dwellings = tents + vehicles

  ;; get 2020
  d2 = mrdfits('official2020occupancies.fits', 1)
  d2 = trans2020official(d2, wts = data[0].WTS)

  p2020 = d2.TOT_IND
  d2020 = d2.TOT_OBJ - p2020

  pdelta = persons - p2020
  ddelta = dwellings - d2020

  q1 = where(ddelta / d2020 gt 0.1, nq1)
  q2 = where(ddelta / d2020 gt 0.5, nq2)
  q3 = where(ddelta / d2020 gt 1.0, nq3)
  print, data[q1].TRACT
  print, 'F tracts', nq1 / float(n_elements(ddelta))
  print, 'F T+M   ', total(data[q1].RAWCOUNTS[6:7])/total(data.RAWCOUNTS)
  print, 'F Dwel  ', total(data[q1].RAWCOUNTS[3:7])/total(data.RAWCOUNTS)
  print, 'F Cts   ', total(data[q1].RAWCOUNTS)/total(data.RAWCOUNTS)
  print, data[q2].TRACT
  print, 'F tracts', nq2 / float(n_elements(ddelta))
  print, 'F T+M   ', total(data[q2].RAWCOUNTS[6:7])/total(data.RAWCOUNTS)
  print, 'F Dwel  ', total(data[q2].RAWCOUNTS[3:7])/total(data.RAWCOUNTS)
  print, 'F Cts   ', total(data[q2].RAWCOUNTS)/total(data.RAWCOUNTS)
  print, data[q3].TRACT
  print, 'F tracts', nq3 / float(n_elements(ddelta))
  print, 'F T+M   ', total(data[q3].RAWCOUNTS[6:7])/total(data.RAWCOUNTS)
  print, 'F Dwel  ', total(data[q3].RAWCOUNTS[3:7])/total(data.RAWCOUNTS)
  print, 'F Cts   ', total(data[q3].RAWCOUNTS)/total(data.RAWCOUNTS)
  
  swap = where(pdelta lt 0 AND ddelta gt 0, nswap)
  
  d = data[swap]

  s = sort(ddelta(swap))
  for ii = 0, nswap - 1 do $
     print, f = '(%"Tract: %7.2f | Tent increase: %i (%4.2f)")', $
            data[swap[s[ii]]].TRACT, ddelta[swap[s[ii]]], $
            ddelta[swap[s[ii]]] / d2020[swap[s[ii]]]
  
  print, float(nswap) / n_elements(p2020)
  print, total(d.RAWCOUNTS)/total(data.RAWCOUNTS)
  print, total(mean(d.COUNTS, dim = 2)) $
         / total(mean(data.COUNTS, dim = 2))
 
end
