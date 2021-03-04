pro charVolOnly

  data = mrdfits('countHollywoodResults2021.fits', 1)
  h    = data[where(~data.EASTFLAG)]
  e    = data[where(data.EASTFLAG)]
  hpro = idProTracts(h.TRACT)
  epro = idProTracts(e.TRACT)

  vh = h[where(~hpro)]
  ve = e[where(~epro)]
  ph = h[where(hpro)]
  pe = e[where(epro)]
  
  d2020 = mrdfits('official2020occupancies.fits',1)
  d2020 = trans2020official(d2020, wts = data[0].WTS[3:7])
  h2020 = d2020[where(~data.EASTFLAG)]
  e2020 = d2020[where(data.EASTFLAG)]
  vh2020 = h2020[where(~hpro)]
  ve2020 = e2020[where(~epro)]
  ph2020 = h2020[where(hpro)]
  pe2020 = e2020[where(epro)]
  
  print, total(hpro), total(epro)
  print, total(hpro)/n_elements(h), total(epro)/n_elements(e)
  
  print, '2021'
  print, total(ph.RAWCOUNTS)/total(h.RAWCOUNTS), $
         total(pe.RAWCOUNTS)/total(e.RAWCOUNTS)
  print, '2020'
  print, total(ph2020.TOT_OBJ)/total(h2020.TOT_OBJ), $
         total(pe2020.TOT_OBJ)/total(e2020.TOT_OBJ)
  
  hcts = [total(total(h.RAWCOUNTS[0:2], 1)), total(total(h.RAWCOUNTS[3:7], 1))]
  vhcts = [total(total(vh.RAWCOUNTS[0:2], 1)), total(total(vh.RAWCOUNTS[3:7], 1))]
  phcts = [total(total(ph.RAWCOUNTS[0:2], 1)), total(total(ph.RAWCOUNTS[3:7], 1))]
  ects = [total(total(e.RAWCOUNTS[0:2], 1)), total(total(e.RAWCOUNTS[3:7], 1))]
  vects = [total(total(ve.RAWCOUNTS[0:2], 1)), total(total(ve.RAWCOUNTS[3:7], 1))]
  pects = [total(total(pe.RAWCOUNTS[0:2], 1)), total(total(pe.RAWCOUNTS[3:7], 1))]
  
  base2020 = [total(h2020.TOT_IND), total(h2020.TOT_OBJ - h2020.TOT_IND), $
              total(e2020.TOT_IND), total(e2020.TOT_OBJ - e2020.TOT_IND)]
  plot2020 = [total(vh2020.TOT_IND), total(vh2020.TOT_OBJ - vh2020.TOT_IND), $
              total(ve2020.TOT_IND), total(ve2020.TOT_OBJ - ve2020.TOT_IND)]
  pro2020 = [total(ph2020.TOT_IND), total(ph2020.TOT_OBJ - ph2020.TOT_IND), $
              total(pe2020.TOT_IND), total(pe2020.TOT_OBJ - pe2020.TOT_IND)]

  base2021 = [ hcts, ects]
  plot2021 = [vhcts,vects]
  pro2021  = [phcts,pects]

  if 0 then begin
  cgbarplot, base2020, col = 'cccccc'x, barcoord = bx, $
             pos = [0.1,0.25,0.9,0.9]
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

  cgplot, bx, base2021/base2020, /noEr, $
          pos = [0.1,0.1,!X.WINDOW[1],!Y.WINDOW[0]], $
          xtickname = [' ','HI', 'HD', 'EI', 'ED', ' '], yr = [0.2,1.2], xticks = 4, $
          xr = !X.CRANGE, xtickint = 1, ytickint = 0.25, xtickv = bx
  oploterror, bx, plot2021/plot2020, plot2021/plot2020 * sqrt(1./plot2020 + 1./plot2021), col = 'ffa500'x
  oploterror, bx, pro2021/pro2020, pro2021/plot2020 * sqrt(1./pro2020 + 1./pro2021), col = '00a5ff'x
  endif

  x = [0,1,2,3]
  plot, x, base2021/base2020, $;
        xtickname = [' ','HI', 'HD', 'EI', 'ED', ' '], yr = [0.2,1.2], xticks = 4, $
        xr = [-0.5,3.5], xtickint = 1, ytickint = 0.25, xtickv = bx, /ys, $
        ytitle = '% decline vs 2020'
  oploterror, x, base2021/base2020, base2021/base2020 * sqrt(1./base2020 + 1./base2021), $
              col = 'ffffff'x, errcol = 'ffffff'x
  oploterror, x, plot2021/plot2020, plot2021/plot2020 * sqrt(1./plot2020 + 1./plot2021), $
              col = 'ffa500'x, errcol = 'ffa500'x
  oploterror, x, pro2021/pro2020, pro2021/plot2020 * sqrt(1./pro2020 + 1./pro2021), $
              col = '00a5ff'x, errcol = '00a5ff'x
  oplot, !X.CRANGE, [1,1], col = 255
  legend, /bottom, /left, $
          ['pro', 'vol', 'all'], $
          col = ['00a5ff'x,'ffa500'x,'ffffff'x], $
          pspacing = 0.5, linesty = 0, box = 0
  
  stop
  
end
