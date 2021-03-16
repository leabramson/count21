pro doGini2020

  data = mrdfits('official2020completeOccupancies.fits', 1)
  d = trans2020official(data)
  readcol, 'hwood2020usCensus.csv', tract, pop, income, f = 'A,F,F', delim = ','
  p = pop/total(pop)
  s = sort(p)
  gGeneral = total(p[s], /cum)

  ;; Sort stuff
  rawcounts = transpose([[d.TOT_IND],[d.C],[d.V],[d.R],[d.T],[d.M]])
  upops     = transpose([[data.PERSONS],[data.C],[data.V],[data.R],[data.T],[data.M]])
  
  cts        = total(RAWCOUNTS, 1)
  tot_cts    = total(RAWCOUNTS)
  ppl        = total(UPOPS, 1)
  tot_ppl    = total(UPOPS)
  dwel       = total(UPOPS[4:*,*], 1)
  tot_dwel   = total(UPOPS[4:*,*])
  
  scounts = sort(cts / tot_cts)
  bc = total(cts[scounts], /cum) / tot_cts
  scounts = sort(ppl / tot_ppl)
  bp = total(ppl[scounts], /cum) / tot_ppl
  scounts = sort(dwel / tot_dwel)
  bd = totaL(dwel[scounts], /cum) / tot_dwel
  
  x = findgen(n_elements(bc))/(n_elements(bc) - 1)
  dx = x[1]-x[0]
  
  gini        = 1.-2*total(bp*dx)
  giniDwel    = 1.-2*total(bd*dx)
  giniGeneral = 1.-2*total(gGeneral*dx)
  
  print, 'Gini_unsh = ', gini
  print, 'Gini_tent = ', giniDwel
  print, 'Gini_tot  = ', giniGeneral

  cgloadct, 27, /brewer, ncol = 9
  
  tfill = cgcolor('2');'00a5ff'x; 
  ufill = cgcolor('4');'5c4cea'x; long('0055ff'x)
  
  set_plot, 'PS'
  device, filename = 'gini2020.eps', $
          /col, /encap, /decomp, bits_per_pix = 8
  plot, x, bp, psym = 10, $
        xtitle = 'fraction of tracts', $
        ytitle = 'fraction of population', $
        /iso, xthick = 4, ythick = 4, charthick = 5, charsize = 1.5, $
        thick = 6, /nodat, xr = [0,1], yr = [0,1]
  polyfill, [x,reverse(x)], [x, reverse(bp)], col = ufill;'cccccc'x
  polyfill, [x,reverse(x)], [x, reverse(gGeneral)], col = tfill;'ffa500'x
  oplot, x, gGeneral, thick = 8;, col = 'ff5500'x
  oplot, x, bp, thick = 8;, col = cgcolor('3')
  one_one, thick = 10, col = cgcolor('7')
  oplot, x[value_locate(bp, 0.5)] * [1,1], [0,bp[value_locate(bp, 0.5)]], $
         thick = 6, linesty = 2, col = ufill
  oplot, x[value_locate(gGeneral, 0.5)] * [1,1], [0,gGeneral[value_locate(gGeneral, 0.5)]], $
         thick = 6, linesty = 2, col = tfill
  plotsym, 8, /fill
  legend, pos = [!X.WINDOW[0],!Y.WINDOW[1]-0.025], /norm, box = 0, $
          ['Equitable', 'Unshelt., !18c!X!DGini!N='+string(gini, f = '(F4.2)'), $
           'All, !18c!X!DGini!N='+string(giniGeneral, f = '(F4.2)'), $
           'Median'], $
          linesty = [0,0,0,1], psym = [0,8,8,0], $
          col = [cgcolor('7'), ufill, tfill, 0], pspacing = 0.5, $
          charsize = 1.25, charthick = 4, thick = 6, spacing = 1.75
  device, /close
  set_plot, 'X'
  spawn, 'open gini2020.eps &'

;  stop
  
end
