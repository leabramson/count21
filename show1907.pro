pro show1907

  ;; 2021
  data = mrdfits('countHollywoodResults2021.fits', 1)
  d = data[where(data.TRACT eq 1907.00)]

  ppl = total(d.COUNTS[0:2,*], 1)
  dwellings = total(d.COUNTS[3:7,*], 1)
  
  outs = fltarr(2,5)
  pctles = [0.05,0.16,0.5,0.84,0.95]
  outs[0,*] = getCountProb(ppl, pctles, /inv)
  outs[1,*] = getCountProb(dwellings, pctles, /inv)

  print, outs
  
  ;; 2020 from
  ;; https://www.lahsa.org/data?id=45-2020-homeless-count-by-community-city
  p2020 = [73., 104.-73.]
  p2020Err = sqrt(p2020)

  pcol = 'aa00aa'x
  dcol = '00aa00'x
  
  plotsym, 0, /fill
  x = [0,1]

  set_plot, 'PS'
  device, filename = 'tract1907comp.eps', $
          /col, /encap, /decomp, bits_per_pix = 8
  !P.CHARSIZE = 1.25
  !P.CHARTHICK = 4
  !X.THICK = 4
  !Y.THICK = 4
  
  plot, [-1,2], [0,150], /nodat, xr = [-0.25,1.25], /xsty, $
        xtitle = 'year', ytitle = 'total unsheltered people', $
        xtickint = 1, xtickname = ['2020', '2021'], $
        xminor = 1, xtickv = [0,1], xticks = 1, yr = [0,150], $
        title = 'Tract 1907.00 (Central Hollywood)'
  oplot, x, [total(p2020), total(outs[*,2])], col = 'aaaaaa'x, thick = 10
  oplot, x, [total(p2020), total(outs[*,2])], psym = 8, symsize = 2
  oplot, x, [total(p2020), total(outs[*,2])], psym = 8, col = 'aaaaaa'x, symsize = 1.5
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
  oplot, x, [p2020[1], outs[1,2]], col = dcol, thick = 8
  oploterror, x, [p2020[1], outs[1,2]], [2*p2020Err[1], outs[1,4]-outs[1,2]], $
              /hibar, /nohat, errthick = 6, psym = 3
  oploterror, x, [p2020[1], outs[1,2]], [2*p2020Err[1], outs[1,2]-outs[1,0]], $
              /lobar, /nohat, errthick = 6, psym = 3
  oploterror, x, [p2020[1], outs[1,2]], [p2020Err[1], outs[1,3]-outs[1,2]], $
              /hibar, /nohat, errthick = 10, psym = 3
  oploterror, x, [p2020[1], outs[1,2]], [p2020Err[1], outs[1,2]-outs[1,1]], $
              /lobar, /nohat, errthick = 10, psym = 3
  oplot, x, [p2020[1], outs[1,2]], psym = 8, symsize = 2
  oplot, x, [p2020[1], outs[1,2]], psym = 8, symsize = 1.5, col = dcol
  legend, /top, /left, box = 0, $
          ['Total', 'Persons', 'Dwellings'], $
          linesty = 0, col = ['aaaaaa'x, pcol, dcol], $
          thick = 6, pspacing = 0.75
  
  device, /close
  set_plot, 'X'
  spawn, 'open tract1907comp.eps &'
  
  
end
