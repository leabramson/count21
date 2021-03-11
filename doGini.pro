pro doGini

  data = mrdfits('countHollywoodResults2021.fits', 1)
  readcol, 'hwood2020usCensus.csv', tract, pop, f = 'A,F', delim = ','
  p = pop/total(pop)
  s = sort(p)
  gGeneral = total(p[s], /cum)
     
  struct     = total(data.RAWCOUNTS[3:7], 1)
  cts        = total(data.RAWCOUNTS, 1)
  tot_struct = total(data.RAWCOUNTS[3:7])
  tot_cts    = total(data.RAWCOUNTS)
  
  scounts = sort(cts / tot_cts)
  b = total(cts[scounts], /cum) / tot_cts
  x = findgen(n_elements(b))/(n_elements(b) - 1)
  dx = x[1]-x[0]

  allPpl = total(data.COUNTS, 1)
  cap = fltarr(n_elements(data),n_elements(allPpl[*,0]))
  for ii = 0, 1d4-1 do begin
     tp = reform(allppl[ii,*])
     s = sort(tp / total(tp))
     cap[*,ii] = total(tp[s], /cum) / total(tp)
  endfor
  ag = 1. - 2*total(cap*dx, 1)
  b = mean(cap, dim = 2)
  
  gini = mean(ag); 1.-2*total(b*dx)
  giniGeneral = 1.-2*total(gGeneral*dx)
  
  print, 'Gini_unsh = ', gini
  print, 'Gini_tot  = ', giniGeneral

  set_plot, 'PS'
  device, filename = 'gini.eps', $
          /col, /encap, /decomp, bits_per_pix = 8
  plot, x, b, psym = 10, $
        xtitle = 'fraction of tracts', $
        ytitle = 'fraction of population', $
        /iso, xthick = 4, ythick = 4, charthick = 5, charsize = 1.5, $
        thick = 6, /nodat, xr = [0,1], yr = [0,1]
  polyfill, [x,reverse(x)], [x, reverse(b)], col = 'cccccc'x
  polyfill, [x,reverse(x)], [x, reverse(gGeneral)], col = 'ffa500'x
  oplot, x, gGeneral, thick = 6, col = 'ff5500'x
  oplot, x, b, thick = 6
  one_one, thick = 6
  oplot, x[value_locate(b, 0.5)] * [1,1], [0,b[value_locate(b, 0.5)]], $
         thick = 6, linesty = 1
  oplot, x[value_locate(gGeneral, 0.5)] * [1,1], [0,gGeneral[value_locate(gGeneral, 0.5)]], $
         thick = 6, linesty = 1, col = 'ff5500'x
  plotsym, 8, /fill
  legend, /top, /left, box = 0, $
          ['Equitable', 'Unsheltered, Gini='+string(gini, f = '(F4.2)')+$
           textoidl('\pm')+string(stddev(ag),f='(F4.2)'), $
           'All (2020), Gini='+string(giniGeneral, f = '(F4.2)'), $
           'Median'], $
          linesty = [0,0,0,1], psym = [0,0,0,0], $
          col = [255, 0, 'ffa500'x, 0], pspacing = 0.5, $
          charsize = 1.25, charthick = 4, thick = 6
  device, /close
  set_plot, 'X'
  spawn, 'open gini.eps &'
end
