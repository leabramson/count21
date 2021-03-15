;; Tract by tract comparison between years

pro plotTractTractOfficial, newData;, $
;                            region = region

;  foo  = mrdfits('official2020Occupancies.fits', 1)
  foo  = mrdfits('official2020CompleteOccupancies.fits', 1)
  newd = mrdfits(newData, 1)
  oldD = trans2020official(foo, wts = newD[0].WTS[3:7]) ;; ensure same weights used year on year at least
  
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
  
  oldRaw    = oldD.TOT_OBJ
  oldRawErr = sqrt(oldRaw)
  oldTot    = oldD.TOT_PPL
;  oldTotErr = oldRawErr * mean(newD[0].WTS[3:7])

  wts = newD[0].WTS[3:7]
  oldTotErr = sqrt(total([[oldD.C*wts[0]^2], [oldD.V*wts[1]^2], $
                          [oldD.R*wts[2]^2], [oldD.T*wts[3]^2], [oldD.M*wts[4]^2]], 2) $
                   + oldD.TOT_IND)
  
  
;  if strupcase(region) eq 'HWOOD' then $
;     oldtypes = [410,0,0,55,48,32,222,64,0] $
;  else if strupcase(REGION) eq 'EHO' then $
;     oldtypes = [164,0,0,29,58,11,94,113,0]

;  oldtypes = [410,0,0,55,48,32,222,64,0] + [164,0,0,29,58,11,94,113,0]
  oldTypes = [total(oldD.TOT_IND), 0, 0, $
              total(oldD.C), total(oldD.V), total(oldD.R), total(oldD.T), total(oldD.M), $
              0]
  
  oldTypesErr = sqrt(oldTypes)
  newTypes    = total(newD.RAWCOUNTS, 2)
  newTypesErr = sqrt(newTypes / newD.NCOUNTERS)
  newTypes[0] += total(newTypes[[1,2,8]])
  newTypes[[1,2,8]] = 0
  newTypesErr[0] += sqrt(total(newTypesErr[[1,2,8]]^2))
  newTypesErr[[1,2,8]] = 0
  
  newRaw    = total(newD.RAWCOUNTS, 1)
  newTot    = median(total(newD.COUNTS, 1), dim = 1)
  newRawErr = sqrt(newRaw / newD.NCOUNTERS)
  newTotErr = stddev(total(newD.COUNTS, 1), dim = 1)
  
  delRaw    = newRaw - oldRaw
  delRawErr = sqrt(newRawErr^2 + oldRawErr^2)
  delTot    = newTot - oldTot
  delTotErr = sqrt(newTotErr^2 + oldTotErr^2)

;  window, 0, xsize = 800, ysize = 800
  !p.multi = [0,2,0]
  plot, oldRaw, newRaw, psym = 1, /iso, $
        xran = [0,max([oldRaw,newRaw])], yran = [0,max([oldRaw,newRaw])], $
        xtitle = 'old counts [obj]', ytitle = 'new counts [obj]', /nodat
  one_one
  oploterror, oldRaw, newRaw, oldRawErr, newRawErr, $
              psym = 1, /nohat
  oplot, oldRaw, newRaw, psym = 1, col = 'ffa500'x

  plot, oldTot, newTot, psym = 1, /iso, $
        xran = [0,max([oldToT,newTot])], yran = [0,max([oldTot,newTot])], $
        xtitle = 'old pop [ppl]', ytitle = 'new pop [ppl]', /nodat
  one_one
  oploterror, oldTot, newTot, oldTotErr, newTotErr, $
              psym = 1, /nohat, col = 'ffa500'x
  oplot, oldTot, newTot, psym = 1, col = 'ffa500'x
  
  print, f= '(%"Old total obj: %i")', total(oldRaw)
  print, f= '(%"New total obj: %i")', total(newRaw)
  print, f= '(%"Old total ppl: %i")', total(oldTot)
  print, f= '(%"New total ppl: %i")', total(newTot)

;  stop

  tags = ['Persons', newD[0].TAGS[3:7]]

  inds = [0,3,4,5,6,7]
  
  !p.multi = 0
  cgbarplot, oldTypes[inds], barwidth = 0.25, baroffset = 1, $
             col = '777777'x, barspace = 0.75, barcoord = bx      
  cgbarplot, newTypes[inds], barwidth = 0.25, baroffset = 2, $
             col = 'ffa500'x, barspace = 0.75, /over, barcoord = bx2
  oploterror, bx, oldTypes[inds], oldTypesErr[inds], psym = 3, barcol = 0
  oploterror, bx2, newTypes[inds], newTypesErr[inds], psym = 3, barcol = 'ff0000'x
  for ii = 0, n_elements(bx) - 1 do $
     cgtext, 0.5*(bx+bx2)[ii], -50, strcompress(tags[ii],/rem), align = 0.5;, orien = 15           
  plotsym, 8, /fill
  legend, /top, /right, $
          ['2020 raw counts', '2021 raw counts'], $
          col = ['777777'x, 'ffa500'x], psym = 8, $
          pspacing = 0.5, textcol = 0
;  stop

  s = sort(delRaw)

  nUpRaw    = total(delRaw gt 0)
  nUpRawSig = total(delRaw gt delRawErr)
  nUpTot    = total(delTot gt 0)
  nUpTotSig = total(delTot gt delTotErr)

  nDnRaw    = total(delRaw lt 0)
  nDnRawSig = total(delRaw lt -delRawErr)
  nDnTot    = total(delTot lt 0)
  nDnTotSig = total(delTot lt -delTotErr)
  
  netRaw    = total(newRaw - oldRaw)
  netRawErr = sqrt(total(newRawErr^2 + oldRawErr^2))
  netTot    = total(newTot - oldTot)
  netTotErr = sqrt(total(newTotErr^2 + oldTotErr^2))

  eho = where(newD[s].EASTFLAG)

  foo = idProTracts(newD[s].Tract)
  pros = where(foo)
  
  set_plot, 'PS'
  device, filename = 'tractsYrYr.eps', $
          /col, /encap, /decomp, bits_per_pix = 8, $
          xsize = 8, /in
  !X.THICK = 4
  !Y.THICK = 4
  !P.CHARTHICK = 4
  !P.CHARSIZE = 1.25

  cgbarplot, delRaw[s]>(-50), $
             ytitle = 'change from 2020 [counts or ppl]', $
             barcoord = bx, baroffset = 1, $
             barwidth = 0.3, barspace = 0.75, $
             col = '777777'x, yr = [-50,50], /ys, $
             title = 'Tract-by-tract Comparison'
  cgbarplot, delTot[s]>(-50), /over, $
             barcoord = bx2, baroffset = 2, $
             barwidth = 0.3, barspace = 0.75, col = 'ff00ff'x

;  bx = findgen(n_elements(s))
;  bx2 = bx
;  plot, bx, delRaw[s], thick = 6, $
;        xr = [-1,n_elements(s)+1], $
;        xtickname = replicate(' ', 60), $
;        ytitle = 'change from 2020 [counts or ppl]', $
;        yr = [-40,40], /xs
;  oplot, bx, delTot[s], thick = 6, col = 'ff00ff'ax
  oploterror, bx, delRaw[s], delRawErr[s], $
              psym = 3, /nohat, errcol = 0, errthick = 4
  oploterror, bx2, delTot[s], delTotErr[s], $
              psym = 3, /nohat, errcol = 'aa00aa'x, errthick = 4
  plotsym, 0, /fill
  oplot, bx[pros], delRaw[s[pros]], psym = 8, symsize = 1
  oplot, bx2[pros], delTot[s[pros]], psym = 8, symsize = 1, col = 'aa00aa'x
  oplot, !X.CRANGE, [0,0], thick = 4, col = 0
;  for ii = 0, total(foo) - 1 do $
;     cgtext, 0.5*(bx+bx2)[pros[ii]], delRaw[s[pros[ii]]]+5, /data, "pro", charsize = 1, charthick = 2, col = 0, align = 0.5
  for ii = 0, n_elements(bx) - 1 do $
     cgtext, 0.5 * (bx + bx2)[ii], !Y.CRANGE[0] - 0.05 * (!Y.CRANGE[1]-!Y.CRANGE[0]), $
             string(oldD[s[ii]].TRACT, f = '(F7.2)'), align = 0.75, orien = 45, /data, col = 0, $
             charsize = 0.8
  legend, /top, /left, box = 0, $
          ['counts', 'people', 'pro counters'], $
          col = [0, 'ff00ff'x, 0], psym = [0,0,8], linesty = [0,0,0], thick = 6, pspacing = 0.5
  legend, /bottom, /right, box = 0, $
          ['!18N!X!Dcounts, up!N: '+string(nupRaw, f = '(I0)')+' ('+string(nupRawSig, f = '(I0)')+')', $
           '!18N!X!Dcounts, dn!N: '+string(ndnRaw, f = '(I0)')+' ('+string(ndnRawSig, f = '(I0)')+')', $
           '!18N!X!Dpeople, up!N: '+string(nupTot, f = '(I0)')+' ('+string(nupTotSig, f = '(I0)')+')', $
           '!18N!X!DPeople, dn!N: '+string(ndnTot, f = '(I0)')+' ('+string(ndnTotSig, f = '(I0)')+')', $
           'Net counts: '+string(netRaw, f = '(I0)')+texToIdl('\pm')+string(netRawErr, f = '(I0)'), $
           'Net people: '+string(netTot, f = '(I0)')+texToIdl('\pm')+string(netTotErr, f = '(I0)')], $
          linesty = -1, pspacing = 0.5
  device, /close
  spawn, 'open tractsYrYr.eps &'
  
  stop
  
end
