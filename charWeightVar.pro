pro charWeightVar

  ;; SPA4/CD13

  wts = fltarr(5,2,3)

  wts[*,0,0] = [1.58, 1.52, 1.77, 1.49, 1.56]
  wts[*,0,1] = [1.43, 1.70, 1.96, 1.33, 1.87]
  wts[*,0,2] = [1.51, 1.77, 1.42, 1.48, 1.68]

  wts[*,1,0] = [0.20, 0.31, 0.34, 0.19, 0.15]
  wts[*,1,1] = [0.24, 0.28, 0.30, 0.16, 0.32]
  wts[*,1,2] = [0.25, 0.42, 0.28, 0.11, 0.31]

  cgloadct, 27, /brewer, ncol = 5

  !p.multi = [0,2,0]
  x = [0,1,2]

  offs = randomn(seed, 5) * 0.1

  window, 0, xsize = 800, ysize = 500
  
  plot, x, wts[0,0,*], xr = [-0.5,2.5], $
        yr = [1,2.5], xtickname = ['2018', '2019', '2020'], $
        xticks = 2, xtickv = x, /nodat, $
        ytitle = 'weight'
  for ii = 0, 4 do $
     oplot, x+offs[ii], wts[ii,0,*], $;wts[ii,1,*], $
                 col = cgcolor(string(fix(ii)));, errcol = cgcolor(string(fix(ii)))
  legend, /top, /left, box = 0, $
          ['C','V','R','T','M'], $
          col = [cgcolor('0'), cgcolor('1'), cgcolor('2'), $
                  cgcolor('3'), cgcolor('4')], $
          pspacing = 0.5, linesty = 0

  plot, x, wts[0,0,*], xr = [-0.5,2.5], $
        yr = [0.5,1.5], xtickname = ['2018', '2019', '2020'], $
        xticks = 2, xtickv = x, /nodat, ytitle = 'weight / 2020'
  for ii = 0, 4 do $
     oplot, x+offs[ii], wts[ii,0,*]/wts[ii,0,0], $;wts[ii,1,*], $
                 col = cgcolor(string(fix(ii)));, errcol = cgcolor(string(fix(ii)))
  
  stop
  

end  
