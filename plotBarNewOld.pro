pro plotBarNewOld, struct, lastYearRaw, $
                   output = output, $
                   eho = eho, $
                   hwood = hwood, $
                   mo = mo

  cgloadct, 25, /brewer, ncol = 8
  col = cgcolor('1')
  errcol = cgcolor('2')
    
  if keyword_set(EHO) then begin
     lastYearRaw = [164.,29.,58.,11.,94.,113.,0] ;; from the LAHSA community sheets
     title = 'East Hollywood CoC'
     col = '00a5ff'x
     errcol = long('0055ff'x)
     tcol = long('0000ff'x)
  endif else if keyword_set(HWOOD) then begin
     lastYearRaw = [411.,55.,48.,32.,222.,64.,0]
;     lastYearRaw -= [10, 2, 2, 3, 2, 2, .0]
     title = 'Hollywood CoC'
     col = 'ff5500'x
     errcol = 'ffa500'x
     tcol = 'ffff00'x
  endif
  
  ;; A C V R T M F
    
  thisYearRaw = total(struct.RAWCOUNTS, 2)
  thisYearRawErr = struct.RAWCOUNTS $
                   / (struct.NCOUNTERS ## replicate(1,n_elements(thisYearRaw))) ;; variance per tract
  thisYearRawErr = sqrt(total(thisYearRawErr,2))                           ;; total RMS

  tyR       = fltarr(7)
  tyR[0]    = total(thisYearRaw[0:2])
  tyr[1:*]  = thisYearRaw[3:*]
  tyRe      = fltarr(7)
  tyRe[0]   = sqrt(total(thisYearRawErr[0:2]^2))
  tyRe[1:*] = thisYearRawErr[3:*]
  ;; Collapse the new stuff into individuals

  set_plot, 'PS'
  !X.THICK = 4
  !Y.THICK = 4
  !P.CHARSIZE = 1.4
  !P.CHARTHICK = 5

  ;; Bar plot
  ;; A+TAY+UM C V R T M F  

  tags = ['Persons', 'Cars', 'Vans', 'RVs', 'Tents', 'Makeshifts', 'Families']
  
  device, filename = output, $
          /col, /encap, /decomp, bits_per_pix = 8
  cgbarplot, lastYearRaw, barwidth = 0.4, barspace = 0.65, $
             baroffset=0.5, barcoord = bx, col = '777777'x, $
             yr = [0,max([tyR,lastYearRaw])*1.1], ytitle = 'People or dwellings counted (raw)', $
             xticklen = 1d-6, title = title
  oploterror, bx, lastYearRaw, sqrt(lastYearRaw) * 2, psym = 3, errthick = 6, /nohat
  cgbarplot, tyR, barwidth = 0.4, barspace = 0.65, $
             baroffset = 1.5, col = col, /over, barcoord = bx2, xticklen = 1d-6
  oploterror, bx2, tyR, tyRe * 2, psym = 3, $
              errcol = errcol, errthick = 6, /nohat
  tl = 0.5*(bx + bx2)
  for ii = 0, n_elements(lastYearRaw) - 1 do begin
     if lastYearRaw[ii] lt 20 then $
        cgtext, bx[ii], lastYearRaw[ii]+20, string(lastYearRaw[ii], f = '(I0)'), $
                align = 0.5, charsize = 1 $
     else $
        cgtext, bx[ii], lastYearRaw[ii]-20, string(lastYearRaw[ii], f = '(I0)'), $
                align = 0.5, col = 'ffffff'x, charsize = 1
     if tYR[ii] lt 20 then $
        cgtext, bx2[ii], tYR[ii]+20, string(tYR[ii], f = '(I0)'), $
                align = 0.5, charsize = 1, col = errcol $
     else $
        cgtext, bx2[ii], tYR[ii]-20, string(tYR[ii], f = '(I0)'), $
                align = 0.5, col = tcol, charsize = 1
     cgtext, tl[ii], -0.1*!Y.CRANGE[1], orien = 15, align = 0.5, strcompress(tags[ii])
  endfor
  plotsym, 8, /fill
  if not keyword_set(MO) then mo = 'Feb'
  if keyword_set(HWOOD) then begin
     legend, /top, /right, box=0, $
             ['LAHSA 2020 PIT', 'Feb 2021', '95% CI'], $
             psym = [8,8,0], linesty = [0,0,0], col = ['777777'x, col,0], $
             thick = [1,1,6], pspacing = 0.5
     cgtext, !X.WINDOW[1]-0.025, 0.65, 'LAHSA PIT total . '+string(total(lastYearRaw), f='(I0)')+textoidl('\pm')+$
             string(sqrt(total(lastYearRaw)) * 2, f = '(I0)'), /norm, align = 1
     cgtext, !X.WiNDOW[1]-0.025, 0.6, 'Volunteers ......... '+string(total(tyR), f='(I0)')+textoidl('\pm')+$
             string(sqrt(total(tyRe^2)) * 2, f = '(I0)'), /norm, align = 1
  endif else if keyword_set(EHO) then begin
     legend, pos = [bx[1]-0.2*(bx[1]-bx[0]),170], /data, box=0, $
             ['LAHSA 2020 PIT', 'Feb 2021', '95% CI'], $
             psym = [8,8,0], linesty = [0,0,0], col = ['777777'x, col,0], $
             thick = [1,1,6], pspacing = 0.5
     cgtext, bx[1]*1.05, 125, 'LAHSA PIT total . '+string(total(lastYearRaw), f='(I0)')+textoidl('\pm')+$
             string(sqrt(total(lastYearRaw)) * 2, f = '(I0)'), /data, align = 0
     cgtext, bx[1]*1.05, 115, 'Volunteers ......... '+string(total(tyR), f='(I0)')+textoidl('\pm')+$
             string(sqrt(total(tyRe^2)) * 2, f = '(I0)'), /data, align = 0
  endif else begin
     legend, /top, /right, box=0, $
             ['LAHSA 2020 PIT', mo+' 2021', '95% CI'], $
             psym = [8,8,0], linesty = [0,0,0], col = ['777777'x, col,0], $
             thick = [1,1,6], pspacing = 0.5
     cgtext, !X.WINDOW[1]-0.025, 0.65, 'LAHSA PIT total . '+string(total(lastYearRaw), f='(I0)')+textoidl('\pm')+$
             string(sqrt(total(lastYearRaw)) * 2, f = '(I0)'), /norm, align = 1
     cgtext, !X.WiNDOW[1]-0.025, 0.6, 'Volunteers ......... '+string(total(tyR), f='(I0)')+textoidl('\pm')+$
             string(sqrt(total(tyRe^2)) * 2, f = '(I0)'), /norm, align = 1
  endelse
  device, /close
  spawn, 'open '+output+' &'
  
end
