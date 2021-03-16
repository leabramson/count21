pro charVolOnly, $
   cut1927 = cut1927

;  data = mrdfits('countHollywoodResults2021.fits', 1)
  data = mrdfits('countHollywoodResults2021w191902.fits', 1)
  if keyword_set(cut1927) then $
     ;; chuck 1927 and see what happens
     data = data[where(data.TRACT ne 1927.00)]

  h    = data[where(~data.EASTFLAG)]
  e    = data[where(data.EASTFLAG)]
  hpro = idProTracts(h.TRACT)
  epro = idProTracts(e.TRACT)

  vh = h[where(~hpro)]
  ve = e[where(~epro)]
  ph = h[where(hpro)]
  pe = e[where(epro)]
  
;  d2020 = mrdfits('official2020occupancies.fits',1)
;  d2020 = mrdfits('official2020completeOccupancies.fits',1)
  d2020 = mrdfits('official2020completeOccupanciesW191902.fits',1)
  d2020 = trans2020official(d2020, wts = data[0].WTS[3:7])
  if keyword_set(cut1927) then $
     ;; chuck 1927 and see what happens
     d2020 = d2020[where(d2020.TRACT ne 1927.00)]

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
  print, total(mean(ph.COUNTS,dim=2))/total(mean(h.COUNTS,dim=2)), $
         total(mean(pe.COUNTS,dim=2))/total(mean(e.COUNTS,dim=2))
  print, '2020'
  print, total(ph2020.TOT_OBJ)/total(h2020.TOT_OBJ), $
         total(pe2020.TOT_OBJ)/total(e2020.TOT_OBJ)
  print, total(ph2020.TOT_PPL)/total(h2020.TOT_PPL), $
         total(pe2020.TOT_PPL)/total(e2020.TOT_PPL)
  
  hcts = [total(total(h.RAWCOUNTS[0:2], 1)), total(total(h.RAWCOUNTS[3:7], 1))]
  vhcts = [total(total(vh.RAWCOUNTS[0:2], 1)), total(total(vh.RAWCOUNTS[3:7], 1))]
  phcts = [total(total(ph.RAWCOUNTS[0:2], 1)), total(total(ph.RAWCOUNTS[3:7], 1))]
  ects = [total(total(e.RAWCOUNTS[0:2], 1)), total(total(e.RAWCOUNTS[3:7], 1))]
  vects = [total(total(ve.RAWCOUNTS[0:2], 1)), total(total(ve.RAWCOUNTS[3:7], 1))]
  pects = [total(total(pe.RAWCOUNTS[0:2], 1)), total(total(pe.RAWCOUNTS[3:7], 1))]
  
  base2020 = [total(h2020.TOT_IND), total(h2020.TOT_OBJ - h2020.TOT_IND), $
              total(e2020.TOT_IND), total(e2020.TOT_OBJ - e2020.TOT_IND)]
  vol2020 = [total(vh2020.TOT_IND), total(vh2020.TOT_OBJ - vh2020.TOT_IND), $
             total(ve2020.TOT_IND), total(ve2020.TOT_OBJ - ve2020.TOT_IND)]
  pro2020 = [total(ph2020.TOT_IND), total(ph2020.TOT_OBJ - ph2020.TOT_IND), $
              total(pe2020.TOT_IND), total(pe2020.TOT_OBJ - pe2020.TOT_IND)]

  base2021 = [ hcts, ects]
  vol2021  = [vhcts,vects]
  pro2021  = [phcts,pects]

  if keyword_set(cut1927) then begin
     savedata = {base2020: base2020[2:3], $
                 vol2020: vol2020[2:3], $
                 pro2020: pro2020[2:3], $
                 base2021: base2021[2:3], $
                 vol2021: vol2021[2:3], $
                 pro2021: pro2021[2:3]}
     mwrfits, savedata, 'volCompNo1927.fits', /create
  endif

  set_plot, 'PS'
  device, filename = 'volProfComp.eps', $
          /col, /encap, /decomp, bits_per_pix = 8
  !p.charsize = 1.5
  !p.charthick = 5
  !X.thick = 4
  !y.thick = 4

  cgloadct, 27, /brewer, ncol = 9
  
  vc = cgcolor('2');'00a5ff'x; 
  pc = cgcolor('4');'5c4cea'x; long('0055ff'x)
  
;  vc = '00bb00'x;;'898f2c'x; 
;  pc = 'aa00cc'x;'5c4cea'x; long('0055ff'x)
  
  x = [-0.1,1,2,3.1]
  plotsym, 0, /fill
  plot, x, base2021/base2020, $
        xtickname = ["Hollywood!CPersons", "Hollywood!CDwellings", $
                     'E. Ho.!CPersons', 'E. Ho.!CDwellings'], $
        yr = [0.3,1.3], xticks = 3, $
        xr = [-0.5,3.5], ytickint = 0.25, xtickv = x, /ys, $
        ytitle = 'fraction of 2020 raw counts', xminor = 1, /nodat
  oplot, !X.CRANGE, [1,1], col = 'aaaaaa'x, thick = 20
  oplot, mean(!X.CRANGE)*[1,1], !Y.CRANGE, thick = 10;, linesty = 2
  td = mrdfits('volCompNo1927.fits', 1)
  oploterror, x[2:3], td.BASE2021/td.BASE2020, td.BASE2021/td.BASE2020 * sqrt(1./td.BASE2021+1./td.BASE2020), $
              linesty = 2, thick = 6, errthick = 4
  oploterror, x, base2021/base2020, base2021/base2020 * sqrt(1./base2020 + 1./base2021), $
              thick = 8, errthick = 8, psym = 8, symsize = 1.5, /nohat
  oplot, x[0:1], (base2021/base2020)[0:1], thick = 8
  oplot, x[2:3], (base2021/base2020)[2:3], thick = 8
  oplot, x[0:1]-0.1, (vol2021/vol2020)[0:1], col = vc, thick = 10
  oplot, x[2:3]-0.1, (vol2021/vol2020)[2:3], col = vc, thick = 10
  oplot, x[0:1]+0.1, (pro2021/pro2020)[0:1], col = pc, thick = 10
  oplot, x[2:3]+0.1, (pro2021/pro2020)[2:3], col = pc, thick = 10
  oploterror, x-0.1, vol2021/vol2020, vol2021/vol2020 * sqrt(1./vol2020 + 1./vol2021), $
              thick = 10, errthick = 8, psym = 8, symsize = 1.5, /nohat ;, col = vc, errcol = vc
  oploterror, x[2:3]+0.1, td.PRO2021/td.PRO2020, td.PRO2021/td.PRO2020 * sqrt(1./td.PRO2021+1./td.PRO2020), $
              linesty = 2, thick = 8, errthick = 6, col = pc, errcol = pc
  oploterror, x+0.1, pro2021/pro2020, pro2021/pro2020 * sqrt(1./pro2020 + 1./pro2021), $
              thick = 10, errthick = 8, psym = 8, symsize = 1.5, /nohat;, col = pc, errcol = pc
  oplot, x, base2021/base2020, psym = 8, symsize = 1.2
  oplot, x-0.1, vol2021/vol2020, psym = 8, symsize = 1.2, col = vc
  oplot, x+0.1, pro2021/pro2020, psym = 8, symsize = 1.2, col = pc
  legend, /bottom, /left, $
          ['Pro tracts', 'Vol tracts', 'All tracts'], $
          col = [pc,vc,0], $
          pspacing = 1, linesty = 0, box = 0, thick = 10, charsize = 1.25
  legend, /bottom, /right, $
          'w/o 1927.00', $
          pspacing = 1, linesty = 2, box = 0, thick = 8, charsize = 1.25
  device, /close
  spawn, 'open volProfComp.eps &'
  set_plot, 'X'
  
end

;;
;;
;;

pro foo
     
  cgbarplot, base2020, col = 'cccccc'x, barcoord = bx, $
             pos = [0.1,0.25,0.9,0.9]
  cgbarplot, base2021, col = 'ff5500'x, /over
  cgbarplot, vol2020, /over, col = '777777'x
  cgbarplot, vol2021, /over, col = 'ffa500'x
  for ii = 0, 3 do $
     oploterror, bx[ii], vol2020[ii], sqrt(vol2020[ii]), errcol = 0, errthick = 2, psym = 3
  for ii = 0, 3 do $
     oploterror, bx[ii], vol2021[ii], sqrt(vol2021[ii]), errcol = 'ff0000'x, errthick = 2, psym = 3
  plotsym, 8, /fill
  legend, /top, /right, box = 0, $
          ['2020', '2021'], psym = 8, col = ['777777'x,'ffa500'x]

  cgplot, bx, base2021/base2020, /noEr, $
          pos = [0.1,0.1,!X.WINDOW[1],!Y.WINDOW[0]], $
          xtickname = [' ','HI', 'HD', 'EI', 'ED', ' '], yr = [0.2,1.2], xticks = 4, $
          xr = !X.CRANGE, xtickint = 1, ytickint = 0.25, xtickv = bx
  oploterror, bx, vol2021/vol2020, vol2021/vol2020 * sqrt(1./vol2020 + 1./vol2021), col = 'ffa500'x
  oploterror, bx, pro2021/pro2020, pro2021/pro2020 * sqrt(1./pro2020 + 1./pro2021), col = '00a5ff'x

end
