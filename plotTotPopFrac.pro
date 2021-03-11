pro plotTotPopFrac

  readcol, 'hwood2020usCensus.csv', tract, pop, f = 'A,F', delim = ','
  data = mrdfits('countHollywoodResults2021.fits', 1)


  up = total(data.COUNTS, 1)
  up = arrstats(transpose(up))

  pf = (up.P50 / pop) > 0

  s = sort(pf)

  data = data[s]
  tract = tract[s]
  pop = pop[s]
  eh = where(data.EASTFLAG, compl = h)  
  pf = pf[s]
  up = up[s]
  
  mask = fltarr(n_elements(pf))
  emask = mask
  emask[eh] = 1
  hmask = mask
  hmask[h] = 1
  
  epf = pf * emask
  hpf = pf * hmask

  meanOcc = total(up.P50)/total(pop)
  
  set_plot, 'PS'
  device, filename = 'totPopFrac.eps', $
          /col, /encap, /decomp, bits_per_pix = 8
  !P.CHARTHICK = 4 & !P.CHARSIZE = 1.25 & !X.THICK = 4 & !Y.THICK = 4

  plotsym, 8, /fill
  cgbarplot, pf, barcoord = bx, col = 'aaaaaa'x, $
             ytitle = 'unsheltered fraction of total tract population', $
             yminor = 2, yr = [0,0.05];, xtitle = 'tract'
  OPLOT, !X.CRANGE, meanOcc * [1,1], linesty = 5, thick = 8
  cgbarplot, epf, /over, col = '00a5ff'x
  cgbarplot, hpf, /over, col = 'ffa500'x
  for ii = 0, total(emask) - 1 do begin
     oploterror, bx[eh[ii]], epf[eh[ii]], ((up.P95-up.P50)/pop)[eh[ii]], psym = 3, /hibar, $
                 errcol = long('0000ff'x), /nohat, errthick = 6
     oploterror, bx[eh[ii]], epf[eh[ii]], ((up.P50-up.P05)/pop)[eh[ii]], psym = 3, /lobar, $
                 errcol = long('0000ff'x), /nohat, errthick = 6
  endfor
  for ii = 0, total(hmask) - 1 do begin
     oploterror, bx[h[ii]], hpf[h[ii]], ((up.P95-up.P50)/pop)[h[ii]], psym = 3, /hibar, $
                 errcol = long('ff0000'x), /nohat, errthick = 6
     oploterror, bx[h[ii]], hpf[h[ii]], ((up.P50-up.P05)/pop)[h[ii]], psym = 3, /lobar, $
                 errcol = long('ff0000'x), /nohat, errthick = 6
  endfor
  for ii = 0, n_elements(bx) - 1 do $
     if up[ii].P95/pop[ii] lt meanOcc then $
        cgtext, bx[ii], meanOcc + 0.0015, /data, tract[ii], charsize = 0.8, col = 'aaaaaa'x, $
                orien = 75 $
     else $
        cgtext, bx[ii]+0.3*(bx[1]-bx[0]), meanOcc - 0.0025, /data, tract[ii], charsize = 0.8, col = 'ffffff'x, $
                orien = 90, align = 1       
  legend, /top, /left, box = 0, $
          ['Hollywood', 'East Hollywood', '90% confidence', 'mean'], $
          psym = [8,8,0,0], linesty = [0,0,0,2], $
          col = ['ffa500'x, '00a5ff'x, 0, 0], $
          pspacing = 1, thick = 6
  
  device, /close
  set_plot, 'X'
  spawn, 'open totPopFrac.eps &'
  
end
