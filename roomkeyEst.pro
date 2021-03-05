pro roomkeyEst, results

  ly = 1058. + 656.
  
  if N_ELEMENTS(results) eq 0 then $
     results = 'countHollywoodResults2021.fits'

  cd13Seniors = 324.  ;; from 2020 LAHSA sheet
  laCoSeniots = 4939. ;; from 2020 LAHSA sheet

  prkRooms = 1608. ;; https://projectroomkeytracker.com/
  
  cd13Frac = cd13Seniors/laCoSeniots

  cd13Prk = cd13Frac * prkRooms

  d = mrdfits(results, 1)
  cts = median(total(total(d.COUNTS, 3), 1))

  delta = cts - ly

  print, cd13Prk, delta, cd13Prk/delta
  
end
