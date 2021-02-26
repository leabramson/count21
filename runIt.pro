pro runit, csv

  fitsCount, csv, 'countHollywood2021.fits'
  data = mrdfits('countHollywood2021.fits', 1)
  hoTractLookup, data.TRACT
  countAll, 'countHollywood2021.fits', 1d4, output = 'countHollywoodResults2021.fits'
  summarize, 'countHollywoodResults2021.fits'
  makeplots, 'countHollywoodResults2021.fits'
  data = mrdfits('countHollywoodResults2021.fits', 1)
  plotBarNewOld, data[where(~data.EASTFLAG)], /hwood, $
                   output = 'Hwood2120Bars.eps'
  plotBarNewOld, data[where(data.EASTFLAG)], /eho, $
                   output = 'Eho2120Bars.eps'
  lastYear = 1058.
  td = data[where(~data.EASTFLAG)]
  mwrfits, td, 'hwood2021results.fits', /create
  findNullWeights, 'hwood2021Results.fits', lastYear
  cts = total(total(td.COUNTS, 3), 1)
  print, getCountProb(cts, lastYear)

end
