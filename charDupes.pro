pro charDupes, dataFits, $
               outdir = outdir, $
               excl = excl

  if NOT keyword_set(outDir) then outdir = './'

  data = mrdfits(dataFits, 1)
  data = data[where(~data.FLAG)]
  nlines = n_elements(data)     
  
  tract = data.TRACT
  s = sort(tract)
  stract = tract[s]
  utract = stract[uniq(stract)]
  nutract = n_elements(utract)

  ;; find ncounters
  ncount = fltarr(nutract)
  for ii = 0, nutract - 1 do $
     ncount[ii] = total(tract eq utract[ii])

  forComp = where(ncount gt 1, nutract)
  utract = utract[forComp]
  ncount = ncount[forComp]

  ;; find total number of pairs
  ;; one pair for all tracts (t1-t0), then 2 additional pairs for the
  ;; triplets (t2-t0, t2-t1)
  ntrips = total(ncount eq 3)
  np = nutract + 2 * ntrips
  
  ;; make it individuals
  input = transpose([[data.ADULT+data.TAY+data.MINOR],$
                     [data.CAR],[data.VAN],[data.RV],$
                     [data.TENT],[data.MAKESHIFT]]);, [data.FAMILY]])
  ncat = n_elements(input[*,0])

  ;; outputs
  team1  = fltarr(ncat, nutract)
  team2  = fltarr(ncat, nutract)
  team3  = fltarr(ncat, nutract)
  rcomp  = fltarr(ncat, np)
  rerr   = fltarr(ncat, np)
  rtcomp = fltarr(np)
  rterr  = fltarr(np)

  meanCounts = fltarr(ncat, np)
  meanTots = fltarr(np)
  ptracts = fltarr(np)

  jj = 0
  ii = 0
  while ii lt nutract do begin     

     hit = where(tract eq utract[ii], nhit)
     
     if nhit eq 2 then ptracts[jj] = utract[ii]
     if nhit eq 3 then ptracts[jj:jj+nhit-1] = utract[ii]
     
     master = hit[0]
     slave  = hit[1]
     if nhit eq 3 then $
        trip = hit[2]
     
     team1[*,ii] = input[*,master]
     team2[*,ii] = input[*,slave]
     rcomp[*,jj] = input[*,slave] - input[*,master]
     rerr[*,jj]  = sqrt(input[*,master] + input[*,slave])
     tmp0        = total(input[*,master])
     tmp1        = total(input[*,slave])
     rtcomp[jj]  = (tmp1-tmp0)
     rterr[jj]   = sqrt(tmp1+tmp0)

     meanCounts[*,jj] = (input[*,master] + input[*,slave])/2.
     meanTots[jj]   = (tmp0+tmp1)/2.
     
     if n_elements(hit) eq 3 then begin

        tmp2 = total(input[*,trip])

        team3[*,ii]   = input[*,trip]
        rcomp[*,jj+1] = team3[*,ii] - team1[*,ii]
        rcomp[*,jj+2] = team3[*,ii] - team2[*,ii]
        rerr[*,jj+1]  = sqrt(team3[*,ii] + team1[*,ii])
        rerr[*,jj+2]  = sqrt(team3[*,ii] + team2[*,ii])
        rtcomp[jj+1]  = tmp2 - tmp0
        rtcomp[jj+2]  = tmp2 - tmp1
        rterr[jj+1]   = sqrt(tmp2+tmp0)
        rterr[jj+2]   = sqrt(tmp2+tmp1)
        
        meanCounts[*,jj+1] = (input[*,master]+input[*,trip])/2.
        meanCounts[*,jj+2] = (input[*,slave]+input[*,trip])/2.
        meanTots[jj+1]   = (tmp0+tmp2)/2.
        meanTots[jj+2]   = (tmp1+tmp2)/2.
     endif
     
     ii++
     if nhit eq 3 then jj+=nhit else jj++

  endwhile

  ;; normalized spread in counter differences, total
  ;; this should be equivalent to root(meanTots)
  nsp = sqrt(rtcomp^2) / sqrt(2)
  
  if keyword_set(excl) then $
     norm = where(ptracts ne excl, compl = oddTract) $;; known shite tract ;;1901.00 for Hwood
  else $
     norm = indgen(n_elements(ptracts))
            
  ;; cut the sample into nths and do the exercise w/ and w/o hot tract
  nsplit = 4
  xs = fltarr(nsplit,2)
  ys = fltarr(nsplit,2,2)
  for ii = 0, 1 do begin
     case ii of
        0: begin
           tmt = meanTots
           trt = rtcomp
           f0 = fn_stats(tmt, sqrt(trt^2/2), n=10, /meanloc)
        end
        1: begin
           tmt = meanTots[norm]
           trt = rtcomp[norm]
           f1 = fn_stats(tmt, sqrt(trt^2/2), n=10, /meanloc)
        end
     endcase

     nt = n_elements(tmt)
     s = sort(tmt)

     tmt = tmt[s]
     trt = trt[s]

     print, mean(sqrt(trt^2/2)/sqrt(tmt), /nan)
     
     for jj = 0, nsplit - 1 do begin
        stt = (ceil(jj*nt/nsplit)-1)>0
        stp = ceil((jj+1)*nt/nsplit)-1
        nq = (stp-stt)+1
        xs[jj,ii] = mean(tmt[stt:stp])
        ys[jj,*,ii] = [mean(sqrt(trt[stt:stp]^2/2)), stddev(sqrt(trt[stt:stp]^2/2), /nan)/sqrt(nq)]
     endfor
  endfor

  set_plot, 'PS'
  device, filename = outdir+'intDupeChar.eps', $
          /col, /encap, /decomp, bits_per_pix = 8
  !P.CHARTHICK = 5
  !P.CHARSIZE = 1.5
  !X.THICK = 4
  !Y.THICK = 4

  cgloadct, 18, /brewer, ncol = 2, clip = [47,220], /rev
  c1 = 'da8428'x                ;cgcolor(string(1));          
  if keyword_set(excl) then $
     c2 = '69b498'x $           ;cgcolor(string(0));'8de0b9'x;
  else $
     c2 = c1
  nu = greek('nu')
  
  plotsym, 0, /fill
  plot, meanTots, nsp, psym = 8, $
        xtitle = '(!18'+nu+'!X!D1!N+!18'+nu+'!X!D2!N)!18/!X2 [avg. total people + dwellings]', $
        ytitle = 'sqrt[(!18'+nu+'!X!D1!N-!18'+nu+'!X!D2!N)!E2!N!18/!X2]  [people or dwellings]', $
        /nodat, yr = [0,2*sqrt(max(meanTots))], xr = [1,120], /xlog, /xsty
  x = findgen((10.^!X.CRANGE[1]-10.^!X.CRANGE[0])/0.25+1)*0.25 + !X.CRANGE[0]
;  oplot, xs[*,1], ys[*,0,1], psym = 8, symsize = 2.5, col = 'aaaaaa'x
;  oplot, xs[*,0], ys[*,0,0], psym = 8, symsize = 2.5, col = 'aaaaaa'x
;  oploterror, xs[*,1], ys[*,0,1], ys[*,1,1], psym = 8, $
;              errcol = c1, col = c1, symsize = 1.5, errthick = 6
;  oploterror, xs[*,0], ys[*,0,0], ys[*,1,0], psym = 8, $
;              errcol = c2, col = c2, symsize = 2, errthick = 8
;  oplot, xs[*,1], ys[*,0,1], psym = 8, symsize = 1.5, col = 'aa00aa'x
  oplot, x, sqrt(x), col = '0055ff'x, thick = 10, linesty = 4
  oploterror, meanTots, nsp, rterr/sqrt(2), errthick = 3, errcol = 'cccccc'x, /nohat, psym = 3
  oplot, meanTots, nsp, psym = 8, symsize = 1.7, col = 'cccccc'x
  oplot, meanTots, nsp, psym = 8, symsize = 1.1
  oploterror, f0.LOC, f0.MEAN, f0.SIGMA/sqrt(f0.COUNT), $
              psym = 8, errcol = c1, col = c1, symsize = 2.7, errthick = 6
  oploterror, f1[-1].LOC, f1[-1].MEAN, f1[-1].SIGMA/sqrt(f1[-1].COUNT), $
              psym = 8, errcol = c2, col = c2, symsize = 2.7, errthick = 6
  if keyword_set(excl) then $
     legend, /top, /left, box = 0, $
             ['tract-level', 'mean', 'mean excl '+excl, 'Poisson'], $
             psym = [8,8,8,0], linesty = [0,0,0,4], $
             col = [0,c1,c2,'0055ff'x], $
             symsize = [1,2,2,1], pspacing = 1, $
             thick = [1,1,1,6], spacing = 1.7 $
  else $
     legend, /top, /left, box = 0, $
             ['tract-level', 'mean', 'Poisson'], $
             psym = [8,8,0], linesty = [0,0,4], $
             col = [0,c1,'0055ff'x], $
             symsize = [1,2,1], pspacing = 1, $
             thick = [1,1,6], spacing = 1.7

  device, /close
;  spawn, 'open intDupeChar.eps &'
  
  ;'P C V R T M'
  ;; actual mean-square difference from the counters and its uncertainty
  sqrd = mean(sqrt(rcomp^2/2),dim=2)
  esqrd = stddev(sqrt(rcomp^2/2),dim=2)/sqrt(np)
  
  ;; expected poisson noise from avg. occupancy per tract
  sigd = sqrt(mean(meanCounts, dim = 2, /nan))
;  foo = mrdfits('countHollywoodResults2021.fits', 1)
;  sigd2 = sqrt(mean(foo.COUNTS, dim = 2))
  
  x = findgen(ncat)
  device, filename = outDir+'catDupeChar.eps', $
          /col, /encap, /decomp, bits_per_pix = 8

  plot, x, sigd, xr = [-1,6], /xsty, $
        xtickname = [' ', '!18P!X', '!18C!X', '!18V!X', $
                     '!18R!X', '!18T!X', '!18M!X', ' '], $
        xticks = 7, xtickint = 1, $
        ytitle = 'sqrt[(!18n!X!D1!N-!18n!X!D2!N)!E2!N!18/!X2] [people or dwellings]', $
        xminor = 1, /nodat, yr = [0,5], yminor = 2 
;  polyfill, [x,reverse(x)], [sqrd+esqrd, reverse(sqrd-esqrd)], $
;            col = 'ffa500'x
;  oplot, x, sqrd, col = 'ff0000'x
  oplot, x, sigd, col = '0055ff'x, thick = 10, psym = 10, linesty = 4
  oploterror, x, sqrd, esqrd, psym = 8, $
              symsize = 2.7, errthick = 6, errcol = c1, col = c1
  legend, /top, /right, box = 0, $
          ['mean, all pairs, w/ std. err', 'Poisson'], $
          psym = [8,0], linesty = [0,4], $
          col = [c1,'0055ff'x], $
          symsize = [1.4,0], pspacing = 1, thick = [1,6], $
          spacing = 1.7
  device, /close
;  spawn, 'open catDupeChar.eps &'
  set_plot, 'X'

  !X.thick = 1
  !y.thick = 1
  !p.charthick = 1

  s = sort(meanTots)
  qui = s[where(ptracts[s] eq excl, compl = use)]
  use = s[use]
  plotsym, 0, /fill
  plot, meanTots[s], sqrt(rtcomp[s]^2)/rterr[s], psym = 1, $
        xtitle = '!18<n>!X', $
        ytitle = 'sqrt[(!18n!X!D1!N-!18n!X!D2!N)!E2!N!18/(n!X!D1!N+!18n!D!X2!N)]', $
        /xlog, xr = [1,100]
  f = fn_Stats(meanTots[s], sqrt(rtcomp[s]^2)/rterr[s], n=5, /meanLoc)
  g = fn_Stats(meanTots[use], sqrt(rtcomp[use]^2)/rterr[s], n=5, /meanLoc)
  oplot, 10.^!X.CRANGE, [1,1], col = 255
  oplot, g.LOC, g.MEAN, psy = 8, col = '00a500'x
  oplot, f.LOC, f.MEAN, psy = 8, col = 'ffa500'x
  oploterror, g.loc, g.mean, g.sigma/sqrt(g.count), psym = 3, errcol = '00a500'x, /nohat
  oploterror, f.loc, f.mean, f.sigma/sqrt(f.count), psym = 3, errcol = 'ffa500'x, /nohat
  print, mean(sqrt(rtcomp[s]^2)/rterr[s], /nan), median(sqrt(rtcomp[s]^2)/rterr[s]), n_elements(s)

  outlier = where(sqrt(rtcomp[s]^2)/rterr[s] gt 3)
  if n_elements(outlier) gt 0 then $
     print, ptracts[s[outlier]]
  
;  stop
  
  if 0 then begin
  pctles = [0.05,0.16,0.50,0.84,0.95]
  stats = fltarr(ncat,n_elements(pctles))
  for ii = 0, ncat - 1 do begin
     x = comp[ii,where(finite(comp[ii,*]))]
     x = x[sort(x)]
     stats[ii,*] = x[ceil(pctles*nutract)-1]
  endfor

  slp = fltarr(100)
  int = fltarr(100)
  for ii = 0, 99 do begin
     dx = randomn(seed, n_elements(bb[*,0])) * sqrt(bb[*,0])
     dy = randomn(seed, n_elements(bb[*,1])) * sqrt(bb[*,1])
     foo = poly_fit((bb[*,0]+dx)>0, (bb[*,1]+dy) > 0, 1)
     slp[ii] = foo[1]
     int[ii] = foo[0]
  endfor
  
  trips = where(ncount eq 3, ntrips)
  plot, bb[*,0], bb[*,1], psym = 1, xr = [0,120], yr = [0,120], /iso, $
        xtitle = 'count red', ytitle = 'count blue'
  for ii = 0, 99 do $
     oplot, bb[*,0], int[ii] + slp[ii] * bb[*,0], col = '00ff00'x, thick = 1
  one_one
  oploterror, bb[*,0], bb[*,1], sqrt(bb[*,0]), sqrt(bb[*,1]), psym = 1, /nohat
  oplot, bb[trips,0], bbb, psym = 1, col = 255
  oploterror, bb[trips,0], bbb, sqrt(bb[trips,0]), sqrt(bbb), psym = 1, /nohat, $
              errcol = 'ff5500'x
  
  d = findgen(101)-50
  del = bb[*,1] - bb[*,0] ;; duplicates
  del = [del, bbb - bb[trips,0]];, bbb-bb[trips,1]] ;; with triplicates to 0th and 1st measurement
  ebar = sqrt(bb[*,0] + bb[*,1])
  ebar = [ebar, sqrt(bbb + bb[trips,0])];, sqrt(bbb + bb[trips,1])]

  pdf = fltarr(n_elements(d))
;  npdf = fltarr(n_elements(d))
  for ii = 0, n_elements(ebar) - 1 do $
     pdf += (1./sqrt(2*!pi*ebar[ii]^2) * exp(-0.5 * (d-del[ii])^2/ebar[ii]^2))>0
  endif
  
end
;chardupes, 'countHollywood2021w191902.fits'
;chardupes, 'mcw/countMidCity2021.fits', outdir = 'mcw/'
