pro plotEPstuff

  ;; Baseline data from 2020 CVRTM LAHSA website
  d2020f = [0.15,0,0,0.092,0.060,0.038,0.143,0.015,0]
  d2020t = [15., 0, 0, 41., 27., 17., 64.,7.,0]
  
  ntrials = 3
  ef = 11./30. + [-1,0,1] * 1.96/2 * sqrt(11.*19./30^2)/sqrt(30)

  spawn, 'ls epSandboxResults*.fits > test.list'
  readcol, 'test.list', files, f = 'A'
  nfiles = n_elements(files)

  out = fltarr(5,nfiles)
  pctles = [0.05,0.25,0.5,0.75,0.95]
  lastYear = 174.
  overP = fltarr(nfiles)
  for ii = 0, nfiles - 1 do begin
     d = mrdfits(files[ii], 1)
     if ii eq 0 then $
        d2020c = d2020t / [1,1,1,d.WTS,1]
     cts = total(d.COUNTS,1)
     cts = cts[sort(cts)]
     out[*,ii] = cts[ceil(n_elements(cts) * pctles) - 1]
     overp[ii] = 1-getcountprob(cts, lastyear)
     print, files[ii], overP[ii]
  endfor

  set_plot, 'PS'
  !X.THICK = 4
  !Y.THICK = 4
  !P.CHARSIZE = 1.4
  !P.CHARTHICK = 5

  ;; Bar plot
  ;; A T U C V R T M F  
  
  device, filename = 'epSandbox/EPsummaryBars.eps', $
          /col, /encap, /decomp, bits_per_pix = 8
  title = 'Tract '+string(d.TRACT, f = '(F7.2)')+' (EPL+) Data'
  cgbarplot, d2020c, barwidth = 0.4, barspace = 0.65, $
             baroffset=0.5, barcoord = bx, col = '777777'x, $
             yr = [0,max(d.RAWCOUNTS)*1.1], ytitle = 'People or dwellings counted', $
             title = title, xticklen = 1d-6
  oploterror, bx, d2020c, sqrt(d2020c) * 2, psym = 3, errthick = 6, /nohat
  cgbarplot, d.RAWCOUNTS, barwidth = 0.4, barspace = 0.65, $
             baroffset = 1.5, col = 'ff5500'x, /over, barcoord = bx2, xticklen = 1d-6
  oploterror, bx2, d.RAWCOUNTS, sqrt(d.RAWCOUNTS/d.NCOUNTERS) * 2, psym = 3, $
              errcol = 'ffa500'x, errthick = 6, /nohat
  tl = 0.5*(bx + bx2)
  for ii = 0, n_elements(d2020c) - 1 do $
     cgtext, tl[ii], -7, orien = 15, align = 0.5, strcompress(d.TAGS[ii])
  plotsym, 8, /fill
  legend, /top, /left, box=0, $
          ['LAHSA 2020 PIT', 'SELAH fall 2020', '95% CI'], $
          psym = [8,8,0], linesty = [0,0,0], col = ['777777'x, 'ff5500'x,0], $
          thick = [1,1,6], pspacing = 0.5
  cgtext, tl[0], 67, 'LAHSA PIT total . '+string(total(d2020c), f='(I0)')+textoidl('\pm')+$
          string(sqrt(total(d2020c)) * 2, f = '(I0)')
  cgtext, tl[0], 60, 'SELAH ................ '+string(total(d.RAWCOUNTS), f='(I0)')+textoidl('\pm')+$
          string(sqrt(total(d.RAWCOUNTS/d.NCOUNTERS)) * 2, f = '(I0)')
  device, /close
  spawn, 'open epSandbox/EPsummaryBars.eps &'
  
  ;; Summary plot  
  device, filename = 'epSandbox/EPwEmpties.eps', $
          /col, /encap, /decomp, bits_per_pix = 8
  title = 'Tract '+string(d.TRACT, f = '(F7.2)')+' (EPL+) Population Estimates'
  plot, [0,nfiles], minmax(out), /nodat, $
        xtickname = replicate(' ', 60), yr = [110,270], $
        xtickint = 1, yminor = 5, ytitle = 'Unsheltered persons', $
        xminor = 1, title = title, ysty = 8+1, $
        pos = [0.125,0.15,0.875,0.9]
  axis, yaxis = 1, yran = !Y.CRANGE - lastYear, /ysty, ytitle = 'Increase from 2020'
  oplot, nfiles-1.5*[1,1], !Y.CRANGE, thick = 20, col = '555555'x
  oplot, !X.CRANGE, replicate(lastYear,2), thick = 8, linesty = 5         
  xxx = !X.CRANGE[[0,0,1,1]]
  yyy = out[*,0]
  polyfill, xxx, yyy[[0,4,4,0]], col = 'ffff00'x, $
            /line_fill, orien = 45, thick = 1, spacing = 0.025
  polyfill, xxx, yyy[[1,3,3,1]], col = 'ffa500'x, $
            /line_fill, orien = -45, thick = 1, spacing = 0.025
  oplot, !X.CRANGE, yyy[2] * [1,1], col = 'ff5500'x, thick = 8
  cgtext, 0.1, yyy[4]-15, /data, "SELAH fall 2020 w/ SPA4 !18CVRTM!X", col = 'ff0000'x
  plotsym, 0, /fill
  out = out[*,1:*]
  for ii = 0, nfiles - 2 do begin
     x = ii+1
     col = cgcolor(string(fix(ii)))
     oploterror, x, out[2,ii], out[4,ii]-out[2,ii], /hibar, psym = 8, $
                 errthick = 6, symsize = 1.5
     oploterror, x, out[2,ii], out[2,ii]-out[0,ii], /lobar, $
                 errthick = 6, symsize = 1.5
     if ii le 2 then $
        col = '00a5ff'x $
;        oploterror, x, out[2,ii], symsize = 1.2
     else $
        col = 'ff00ff'x
     oplot, [x], [out[2,ii]], psym = 8, symsize = 1.2, col = col

     if ii eq 0 then $
        cgtext, x, out[0,ii]-10, /data, $
                string(overp[ii+1]*100,f='(I2)')+'% chance!Cincrease!Cover 2020', $
                align = 0.5, charsize = 1.1, col = '777777'x $
     else $
        cgtext, x, out[0,ii]-10, /data, $
                string(overp[ii+1]*100,f='(I2)')+'%', $
                align = 0.5, charsize = 1.1, col = '777777'x 
     
     if ii le 2 then $
        cgtext, x, !Y.CRANGE[0]-10, /data, string(ef[ii]*100,f='(I2)')+'%', align = 0.5 $
     else begin
        cgtext, x, !Y.CRANGE[0]-10, /data, "!18T!X=1.1", align = 0.5
        cgtext, x, !Y.CRANGE[0]-25, /data, 'Hwood tent wt', align = 0.5
     endelse
     if ii eq 1 then $
        cgtext, x, !Y.CRANGE[0]-25, /data, 'Empty tent fraction', align = 0.5
  endfor
  cgtext, 0.1,lastyear+4,/data, 'LAHSA', col = '555555'x, align =0
  cgtext, 0.1,lastyear-10,/data, '2020 PIT', col = '555555'x, align =0
  device, /close
  spawn, 'open epSandbox/EPwEmpties.eps &'
  
end
