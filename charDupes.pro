pro charDupes, dataFits

  data = mrdfits(dataFits, 1)
  data = data[where(~data.FLAG)]
  nlines = n_elements(data)

  tract = data.TRACT
  s = sort(tract)
  stract = tract[s]
  utract = tract[uniq(stract)]
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
        rcomp[*,jj+1] = team3[*,jj] - team1[*,jj]
        rcomp[*,jj+2] = team3[*,jj] - team2[*,jj]
        rerr[*,jj+1]  = sqrt(team3[*,jj] + team1[*,jj])
        rerr[*,jj+2]  = sqrt(team3[*,jj] + team2[*,jj])
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

  norm = where(ptracts ne 1901.00, compl = oddTract) ;; known shite tract
  
  ;; cut the sample into nths and do the exercise w/ and w/o hot tract
  nsplit = 4
  xs = fltarr(nsplit,2)
  ys = fltarr(nsplit,2,2)
  for ii = 0, 1 do begin
     case ii of
        0: begin
           tmt = meanTots
           trt = rtcomp
           f0 = fn_stats(tmt, sqrt(trt^2/2), n=10)
        end
        1: begin
           tmt = meanTots[norm]
           trt = rtcomp[norm]
           f1 = fn_stats(tmt, sqrt(trt^2/2), n=10)
        end
     endcase

     nt = n_elements(tmt)
     s = sort(tmt)

     tmt = tmt[s]
     trt = trt[s]
     
     for jj = 0, nsplit - 1 do begin
        stt = (ceil(jj*nt/nsplit)-1)>0
        stp = ceil((jj+1)*nt/nsplit)-1
        nq = (stp-stt)+1
        xs[jj,ii] = mean(tmt[stt:stp])
        ys[jj,*,ii] = [mean(sqrt(trt[stt:stp]^2/2)), stddev(sqrt(trt[stt:stp]^2/2), /nan)/sqrt(nq)]
     endfor
  endfor

  set_plot, 'PS'
  device, filename = 'intDupeChar.eps', $
          /col, /encap, /decomp, bits_per_pix = 8
  !P.CHARTHICK = 4
  !P.CHARSIZE = 1.25
  !X.THICK = 4
  !Y.THICK = 4

  cgloadct, 18, /brewer, ncol = 2, clip = [47,220], /rev
  c1 = cgcolor(string(0));'69b498'x;'8de0b9'x;
  c2 = cgcolor(string(1));'da8428'x;          
  
  plotsym, 0, /fill
  plot, meanTots, nsp, psym = 8, $
        xtitle = '(!18n!X!D1!N+!18n!X!D2!N)!18/!X2 [avg. total people + dwellings]', $
        ytitle = 'sqrt[(!18n!X!D1!N-!18n!X!D2!N)!E2!N!18/!X2]', $
        /nodat, yr = [0,2*sqrt(max(meanTots))], xr = [0,60]
  x = findgen((!X.CRANGE[1]-!X.CRANGE[0])/0.25+1)*0.25 + !X.CRANGE[0]
;  oplot, xs[*,1], ys[*,0,1], psym = 8, symsize = 2.5, col = 'aaaaaa'x
;  oplot, xs[*,0], ys[*,0,0], psym = 8, symsize = 2.5, col = 'aaaaaa'x
;  oploterror, xs[*,1], ys[*,0,1], ys[*,1,1], psym = 8, $
;              errcol = c1, col = c1, symsize = 1.5, errthick = 6
;  oploterror, xs[*,0], ys[*,0,0], ys[*,1,0], psym = 8, $
;              errcol = c2, col = c2, symsize = 2, errthick = 8
;  oplot, xs[*,1], ys[*,0,1], psym = 8, symsize = 1.5, col = 'aa00aa'x
  oplot, meanTots, nsp, psym = 8, symsize = 1.25, col = 'cccccc'x
  oplot, meanTots, nsp, psym = 8, symsize = 0.75
  oploterror, f0.LOC, f0.MEAN, f0.SIGMA/sqrt(f0.COUNT), $
              psym = 8, errcol = c1, col = c1, symsize = 2, errthick = 6
  oploterror, f1[-1].LOC, f1[-1].MEAN, f1[-1].SIGMA/sqrt(f1[-1].COUNT), $
              psym = 8, errcol = c2, col = c2, symsize = 2, errthick = 6
  oplot, x, sqrt(x), col = 255, thick = 6, linesty = 4
  legend, /top, /left, box = 0, $
          ['tract-level', 'mean', 'mean excl 1901.00', 'Poisson'], $
          psym = [8,8,8,0], linesty = [0,0,0,4], $
          col = [0,c1,c2,'0000ff'x], $
          symsize = [0.7,1.5,1.5,1], pspacing = 1, $
          charsize = 1.25, thick = [1,1,1,6]

  device, /close
;  spawn, 'open intDupeChar.eps &'
  
  ;'P C V R T M'
  ;; actual mean-square difference from the counters and its uncertainty
  sqrd = mean(sqrt(rcomp^2/2),dim=2)
  esqrd = stddev(sqrt(rcomp^2/2),dim=2)/sqrt(np)
  
  ;; expected poisson noise from avg. occupancy per tract
  sigd = sqrt(mean(meanCounts, dim = 2, /nan))

  x = findgen(ncat)
  device, filename = 'catDupeChar.eps', $
          /col, /encap, /decomp, bits_per_pix = 8

  plot, x, sigd, xr = [-1,6], /xsty, $
        xtickname = [' ', '!18P!X', '!18C!X', '!18V!X', $
                     '!18R!X', '!18T!X', '!18M!X', ' '], $
        xticks = 7, xtickint = 1, $
        ytitle = 'sqrt[(!18n!X!D1!N-!18n!X!D0!N)!E2!N!18/!X2] (persons or dwellings)', $
        xminor = 1, /nodat, yr = [0,5], yminor = 2 
;  polyfill, [x,reverse(x)], [sqrd+esqrd, reverse(sqrd-esqrd)], $
;            col = 'ffa500'x
;  oplot, x, sqrd, col = 'ff0000'x
  oplot, x, sigd, col = '0000ff'x, thick = 10, psym = 10, linesty = 4
  oploterror, x, sqrd, esqrd, psym = 8, $
              symsize = 2, errthick = 6, errcol = c1, col = c1
  legend, /top, /right, box = 0, $
          ['mean, all pairs, w/ std. err', 'Poisson'], $
          psym = [8,0], linesty = [0,4], $
          col = [c1,'0000ff'x], $
          symsize = [1.4,0], pspacing = 1, thick = [1,6], $
          charsize = 1.25, spacing = 1.5
  device, /close
  spawn, 'open catDupeChar.eps &'
  
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
;chardupes, 'countHollywood2021.fits'
