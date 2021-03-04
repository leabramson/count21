pro charVolOnly

  data = mrdfits('countHollywoodResults2021.fits', 1)
  h    = data[where(~data.EASTFLAG)]
  e    = data[where(data.EASTFLAG)]
  hpro = idProTracts(h.TRACT)
  epro = idProTracts(e.TRACT)

  vh = h[where(~hpro)]
  ve = e[where(~epro)]

  print, 1-(total(vh.RAWCOUNTS)/total(h.RAWCOUNTS)), $
         1-(total(ve.RAWCOUNTS)/total(e.RAWCOUNTS))
  
  d2020 = mrdfits('official2020occupancies.fits',1)
  h2020 = d2020[where(~data.EASTFLAG)]
  e2020 = d2020[where(data.EASTFLAG)]

  vh2020 = h2020[where(~hpro)]
  ve2020 = e2020[where(~epro)]
  
  hcts = [total(total(h.RAWCOUNTS[0:2], 1)), total(total(h.RAWCOUNTS[3:7], 1))]
  vhcts = [total(total(vh.RAWCOUNTS[0:2], 1)), total(total(vh.RAWCOUNTS[3:7], 1))]
  ects = [total(total(e.RAWCOUNTS[0:2], 1)), total(total(e.RAWCOUNTS[3:7], 1))]
  vects = [total(total(ve.RAWCOUNTS[0:2], 1)), total(total(ve.RAWCOUNTS[3:7], 1))]

  base2020 = [total(h2020.PERSONS), total(h2020.TOTAL - h2020.PERSONS), $
              total(e2020.PERSONS), total(e2020.TOTAL - e2020.PERSONS)]
  plot2020 = [total(vh2020.PERSONS), total(vh2020.TOTAL - vh2020.PERSONS), $
              total(ve2020.PERSONS), total(ve2020.TOTAL - ve2020.PERSONS)]

  base2021 = [hcts,ects]
  plot2021 = [vhcts,vects]

  cgbarplot, base2020, col = 'cccccc'x, $
             barname = ['HI', 'HD', 'EI', 'ED'], barcoord = bx
  cgbarplot, base2021, col = 'ff5500'x, /over
  cgbarplot, plot2020, /over, col = '777777'x
  cgbarplot, plot2021, /over, col = 'ffa500'x
  for ii = 0, 3 do $
     oploterror, bx[ii], plot2020[ii], sqrt(plot2020[ii]), errcol = 0, errthick = 2, psym = 3
  for ii = 0, 3 do $
     oploterror, bx[ii], plot2021[ii], sqrt(plot2021[ii]), errcol = 'ff0000'x, errthick = 2, psym = 3
  plotsym, 8, /fill
  legend, /top, /right, box = 0, $
          ['2020', '2021'], psym = 8, col = ['777777'x,'ffa500'x]
  
  stop
  
end
