;; Tract by tract comparison between years

pro plotTractTract, newData, oldData, $
                    test = test

  oldd = mrdfits(oldData, 1)
  newd = mrdfits(newData, 1)

  ;; sort, align, cull

  keep    = []
  newInds = []
  for ii = 0, n_elements(oldD) - 1 do begin
     hit = where(newd.TRACT eq oldd[ii].TRACT, nhit)
     if nhit eq 1 then begin
        keep    = [keep, ii]
        newInds = [newInds, hit]
     endif
  endfor
  oldD = oldD[keep]
  newD = newD[hit]
  
  oldRaw    = total(oldD.RAWCOUNTS, 1)
  oldRawErr = sqrt(oldRaw / oldD.NCOUNTERS)
  oldTot    = mean(total(oldD.COUNTS, 1), dim = 1)
  oldTotErr = stddev(total(oldD.COUNTS, 1), dim = 1)

  newRaw    = total(newD.RAWCOUNTS, 1)
  newTot    = mean(total(newD.COUNTS, 1), dim = 1)
  newRawErr = sqrt(newRaw / newD.NCOUNTERS)
  newTotErr = stddev(total(newD.COUNTS, 1), dim = 1)
  if keyword_set(test) then begin
     newRaw += randomn(seed, n_elements(oldRaw)) * oldRaw/2
     newTot = 1.4 * newRaw
  endif
  
  delRaw    = newRaw - oldRaw
  delRawErr = sqrt(newRawErr^2 + oldRawErr^2)
  delTot    = newTot - oldTot
  delTotErr = sqrt(newTotErr^2 + oldTotErr^2)

;  window, 0, xsize = 800, ysize = 800
  !p.multi = [0,2,0]
  
  plot, oldRaw, newRaw, psym = 1, /iso, $
        xran = minmax(oldRaw), yran = minmax(newRaw), $
        xtitle = 'old counts [obj]', ytitle = 'new counts [obj]', /nodat
  one_one
  oploterror, oldRaw, newRaw, oldRawErr, newRawErr, $
              psym = 1, /nohat

  plot, oldTot, newTot, psym = 1, /iso, $
        xran = minmax(oldTot), yran = minmax(newTot), $
        xtitle = 'old pop [ppl]', ytitle = 'new pop [ppl]', /nodat
  one_one
  oploterror, oldTot, newTot, oldTotErr, newTotErr, $
              psym = 1, /nohat

  !p.multi = 0
  s = sort(delRaw)

  nUpRaw    = total(delRaw gt 0)
  nUpRawSig = total(delRaw gt delRawErr)
  nUpTot    = total(delTot gt 0)
  nUpTotSig = total(delTot gt delTotErr)

  netRaw    = total(newRaw - oldRaw)
  netRawErr = sqrt(total(newRawErr^2 + oldRawErr^2))
  netTot    = total(newTot - oldTot)
  netTotErr = sqrt(total(newTotErr^2 + oldTotErr^2))
  
  window, 0, xsize = 1000, ysize = 800
  cgbarplot, delRaw[s], ytitle = greek('Delta')+' [counts or ppl]', $
             barcoord = bx, baroffset = 1, barwidth = 0.25, barspace = 0.75, col = '777777'x, yr = minmax([delTot,delRaw])
  cgbarplot, delTot[s], /over, $
             barcoord = bx2, baroffset = 2, barwidth = 0.25, barspace = 0.75, col = 'ffa500'x
  oploterror, bx, delRaw[s], delRawErr[s], psym = 3, /nohat, errcol = 0
  oploterror, bx2, delTot[s], delTotErr[s], psym = 3, /nohat, errcol = 'ff0000'x
  oplot, !X.CRANGE, [0,0], thick = 4, col = 0
  for ii = 0, n_elements(bx) - 1 do $
     cgtext, 0.5 * (bx + bx2)[ii], !Y.CRANGE[0] - 0.1 * (!Y.CRANGE[1]-!Y.CRANGE[0]), $
             string(oldD[s[ii]].TRACT, f = '(F7.2)'), align = 0.5, orien = 45, /data, col = 0
  plotsym, 8, /fill
  legend, /top, /left, $
          ['counts', 'people', $
           '!18N!X!Draw, up!N: '+string(nupRaw, f = '(I0)')+' ('+string(nupRawSig, f = '(I0)')+')', $
           '!18N!X!Dttot, up!N: '+string(nupTot, f = '(I0)')+' ('+string(nupTotSig, f = '(I0)')+')', $
           'Net raw: '+string(netRaw, f = '(I0)')+texToIdl('\pm')+string(netRawErr, f = '(I0)'), $
           'Net tot: '+string(netTot, f = '(I0)')+texToIdl('\pm')+string(netTotErr, f = '(I0)')], $
          col = ['777777'x, 'ffa500'x, replicate(0, 4)], psym = [8,8,replicate(0,4)], linesty = replicate(0, 6), $
          pspacing = 0.5, textcol = 0
  
  
end
