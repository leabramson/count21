pro show1907

  ;; 2021
  data = mrdfits('countHollywoodResults2021.fits', 1)
  d = data[where(data.TRACT eq 1907.00)]

  ppl = total(d.COUNTS[0:2,*], 1)
  dwellings = total(d.COUNTS[3:7,*], 1)
  tents     = total(d.COUNTS[6:7,*], 1)
  
  outs = fltarr(3,5)
  pctles = [0.05,0.16,0.5,0.84,0.95]
  outs[0,*] = getCountProb(ppl, pctles, /inv)
  outs[1,*] = getCountProb(dwellings, pctles, /inv)
  outs[2,*] = getCountProb(tents, pctles, /inv)

  print, outs
  
  ;; 2020 from
  ;; https://www.lahsa.org/data?id=45-2020-homeless-count-by-community-city
  p2020 = [73., 104.-73.]
  p2020Err = sqrt(p2020)
  foo  = mrdfits('official2020CompleteOccupancies.fits', 1)
  foo = trans2020official(foo, wts = d.WTS[3:7])
  hit = where(foo.TRACT eq 1907.00)
  t2020 = total(transpose([[foo[hit].T],[foo[hit].M]]))
  t2020Err = sqrt(t2020) * d.WTS[6]
  t2020 *= d.WTS[6]
  
;  print, t2020, t2020 * d.WTS[6], p2020[1]
  
  cgloadct, 27, /brewer, ncol = 9
  
  pcol = cgcolor('3');'aa00aa'x
  dcol = cgcolor('7');'00aa00'x
  
  plotsym, 0, /fill
  x = [0.1,0.9]

  set_plot, 'PS'
  device, filename = 'tract1907comp.eps', $
          /col, /encap, /decomp, bits_per_pix = 8
  !P.CHARSIZE = 1.5
  !P.CHARTHICK = 5
  !X.THICK = 4
  !Y.THICK = 4
  
  plot, [-0.25,1.25], [0,150], /nodat, xr = [-0.25,1.25], /xsty, $
        xtitle = 'year', ytitle = 'total unsheltered people', $
        xtickname = ['2020', '2021'], $
        xminor = 1, xtickv = x, xticks = 1, yr = [0,170], $
        title = 'Tract 1907.00 (Central Hollywood)', /ys
  oplot, x, [total(p2020), total(outs[0:1,2])], col = 'aaaaaa'x, thick = 10
  oplot, x, [total(p2020), total(outs[0:1,2])], psym = 8, symsize = 2
  oplot, x, [total(p2020), total(outs[0:1,2])], psym = 8, col = 'aaaaaa'x, symsize = 1.5
  oplot, x, [p2020[0], outs[0,2]], thick = 8, col = pcol
  oploterror, x, [p2020[0], outs[0,2]], [2*p2020Err[0], outs[0,4]-outs[0,2]], $
              /hibar, /nohat, errthick = 6, psym = 3
  oploterror, x, [p2020[0], outs[0,2]], [2*p2020Err[0], outs[0,2]-outs[0,0]], $
              /lobar, /nohat, errthick = 6, psym = 3
  oploterror, x, [p2020[0], outs[0,2]], [p2020Err[0], outs[0,3]-outs[0,2]], $
              /hibar, /nohat, errthick = 10, psym = 3
  oploterror, x, [p2020[0], outs[0,2]], [p2020Err[0], outs[0,2]-outs[0,1]], $
              /lobar, /nohat, errthick = 10, psym = 3
  oplot, x, [p2020[0], outs[0,2]], psym = 8, symsize = 2
  oplot, x, [p2020[0], outs[0,2]], psym = 8, symsize = 1.5, col = pcol
  oplot, x, [p2020[1], outs[1,2]], col = dcol, thick = 8, linesty = 2
  oplot, x, [t2020, outs[2,2]], col = dcol, thick = 8
  oploterror, x, [t2020, outs[2,2]], [2*t2020Err, outs[2,4]-outs[2,2]], $
              /hibar, /nohat, errthick = 6, psym = 3
  oploterror, x, [t2020, outs[2,2]], [2*t2020Err, outs[2,2]-outs[2,0]], $
              /lobar, /nohat, errthick = 6, psym = 3
  oploterror, x, [t2020, outs[2,2]], [t2020Err, outs[2,3]-outs[2,2]], $
              /hibar, /nohat, errthick = 10, psym = 3
  oploterror, x, [t2020, outs[2,2]], [t2020Err, outs[2,2]-outs[2,1]], $
              /lobar, /nohat, errthick = 10, psym = 3
  oplot, x, [t2020, outs[2,2]], psym = 8, symsize = 2
  oplot, x, [t2020, outs[2,2]], psym = 8, symsize = 1.5, col = dcol
  legend, pos = [!X.WINDOW[0]+0.025, !Y.WINDOW[1]-0.025], /norm, box = 0, $
          ['Total', 'Persons', 'Tents+Mkshfts', 'Dwellings'], $
          linesty = [0,0,0,2], col = ['aaaaaa'x, pcol, dcol, dcol], $
          thick = 6, pspacing = 0.75
  
  device, /close
  set_plot, 'X'
  spawn, 'open tract1907comp.eps &'
  
  
end
