pro dumpYearYear

  d2020 = mrdfits('official2020CompleteOccupanciesW191902.fits',1)
  d2020 = trans2020Official(d2020)
  d2021 = mrdfits('countHollywoodResults2021w191902.fits', 1)

  tc2021 = total(d2021.RAWCOUNTS,1)
  tp2021 = median(total(d2021.COUNTS, 1), dim = 1)

  tc2020 = d2020.TOT_OBJ
  tp2020 = d2020.TOT_PPL

  print, transpose([[d2020.TRACT],[d2021.TRACT]])
  
  close, 1
  openw, 1, 'GreaterHollywoodUnshelteredData.csv', width = 256
  printf, 1, '#TRACT, 2021_TOTAL_PPL, 2020_TOTAL_PPL, 2021_TOTAL_COUNTS, 2020_TOTAL_COUNTS'
  for ii = 0, n_elements(tc2020) - 1 do $
     printf, 1, f = '(%"%7.2f, %5.1f, %5.1f, %3i, %3i")', $
             d2021[ii].TRACT, tp2021[ii], tp2020[ii], tc2021[ii], tc2020[ii]
  close, 1

  print, total(tp2021), total(tc2021), total(tp2020), total(tc2020)
  
end
