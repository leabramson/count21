;; Tract by tract comparison between years

pro plotTractTract, newData, oldData, test = test

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
  newD = newD[newInds]
  
  oldRaw    = total(oldD.RAWCOUNTS, 1)
  oldRawErr = sqrt(oldRaw / oldD.NCOUNTERS)
  oldTot    = mean(total(oldD.COUNTS, 1), dim = 1)
  oldTotErr = stddev(total(oldD.COUNTS, 1), dim = 1)

  oldTypes    = total(oldD.RAWCOUNTS, 2)
  oldTypesErr = sqrt(oldTypes / oldD.NCOUNTERS)
  newTypes    = total(newD.RAWCOUNTS, 2)
  newTypesErr = sqrt(newTypes / newD.NCOUNTERS)
  
  newRaw    = total(newD.RAWCOUNTS, 1)
  newTot    = mean(total(newD.COUNTS, 1), dim = 1)
  newRawErr = sqrt(newRaw / newD.NCOUNTERS)
  newTotErr = stddev(total(newD.COUNTS, 1), dim = 1)
  
  delRaw    = newRaw - oldRaw
  delRawErr = sqrt(newRawErr^2 + oldRawErr^2)
  delTot    = newTot - oldTot
  delTotErr = sqrt(newTotErr^2 + oldTotErr^2)

;  window, 0, xsize = 800, ysize = 800
  !p.multi = [0,2,0]
  plot, oldRaw, newRaw, psym = 1, /iso, $
        xran = minmax(oldRaw)>0, yran = minmax(newRaw)>0, $
        xtitle = 'old counts [obj]', ytitle = 'new counts [obj]', /nodat
  one_one
  oploterror, oldRaw, newRaw, oldRawErr, newRawErr, $
              psym = 1, /nohat
  oplot, oldRaw, newRaw, psym = 1, col = 'ffa500'x

  plot, oldTot, newTot, psym = 1, /iso, $
        xran = minmax(oldTot)>0, yran = minmax(newTot)>0, $
        xtitle = 'old pop [ppl]', ytitle = 'new pop [ppl]', /nodat
  one_one
  oploterror, oldTot, newTot, oldTotErr, newTotErr, $
              psym = 1, /nohat, col = 'ffa500'x
  oplot, oldTot, newTot, psym = 1, col = 'ffa500'x
  
  print, f= '(%"Old total obj: %i")', total(oldRaw)
  print, f= '(%"New total obj: %i")', total(newRaw)
  print, f= '(%"Old total ppl: %i")', total(oldTot)
  print, f= '(%"New total pp;: %i")', total(newTot)

  stop

  !p.multi = 0
  cgbarplot, oldTypes, barwidth = 0.25, baroffset = 1, $
             col = '777777'x, barspace = 0.75, barcoord = bx      
  cgbarplot, newTypes, barwidth = 0.25, baroffset = 2, $
             col = 'ffa500'x, barspace = 0.75, /over, barcoord = bx2
  oploterror, bx, oldTypes, oldTypesErr, psym = 3, barcol = 0
  oploterror, bx2, newTypes, newTypesErr, psym = 3, barcol = 'ff0000'x
  for ii = 0, n_elements(bx) - 1 do $
     cgtext, bx[ii], -20, oldD[0].TAGS[ii], align = 0.5, orien = 15           
  legend, /top, /right, $
          ['OLD counts', 'NEW counts'], $
          col = ['777777'x, 'ffa500'x], psym = 8, $
          pspacing = 0.5, textcol = 0
  stop

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
