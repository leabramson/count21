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

  ;; make it individuals
  
  input = transpose([[data.ADULT+data.TAY+data.MINOR],$
                     [data.CAR],[data.VAN],[data.RV],$
                     [data.TENT],[data.MAKESHIFT]]);, [data.FAMILY]])
  ncat = n_elements(input[*,0])
  
  team1  = fltarr(ncat, nutract)
  team2  = fltarr(ncat, nutract)
  team3  = fltarr(ncat, nutract)
  rcomp  = fltarr(ncat, nutract,3)
  rerr   = fltarr(ncat, nutract,3)
  rtcomp = fltarr(nutract,3)
  rterr  = fltarr(nutract,3)

  ntrips = total(ncount eq 3)
  
  bb = fltarr(nutract,2)
  bbb = fltarr(ntrips)
  for ii = 0, nutract - 1 do begin     

     hit = where(tract eq utract[counter], nhit)
     
     master = 0
     slave = 1
     
     team1[*,ii] = input[*,hit[master]]
     team2[*,ii] = input[*,hit[slave]]
     rcomp[*,ii,0] = input[*,hit[slave]] - input[*,hit[master]]
     rerr[*,ii,0]  = sqrt(input[*,hit[master]] + input[*,hit[slave]])
     
     tmp0 = total(input[*,hit[master]])
     tmp1 = total(input[*,hit[slave]])

     rtcomp[ii] = (tmp1-tmp0)
     rterr[ii]  = sqrt(tmp1+tmp0)

     bb[ii,*] = [tmp0,tmp1]     
     if n_elements(hit) eq 3 then begin
        tmp2 = total(input[*,hit[2]])
        
        team3[*,ii]   = input[*,hit[2]]
        rcomp[*,ii,1] = input[*,hit[2]] - input[*,hit[master]]
        rcomp[*,ii,2] = input[*,hit[2]] - input[*,hit[slave]]
        rerr[*,ii,1]  = sqrt(input[*,hit[2]] + input[*,hit[master]])
        rerr[*,ii,2]  = sqrt(input[*,hit[2]] + input[*,hit[slave]])
        rtcomp[ii,1]  = tmp2 - tmp0
        rtcomp[ii,2]  = tmp2 - tmp1
        rterr[ii,1]   = sqrt(tmp2+tmp0)
        rterr[ii,2]   = sqrt(tmp2+tmp1)
        
        bbb[jj] = tmp2
        jj++
     endif
     
  endfor
  meanCounts = 0.5 * (team1 + team2)
  meanTots   = 0.5 * (total(team1,1) + total(team2,1))

  ;; normalized spread in counter differences, total
  ;; this should be equivalent to root(meanTots)
  nsp = sqrt(rtcomp^2) / sqrt(2)

  norm = where(utract ne 1901.00, compl = oddTract) ;; known shite tract
  
  ;; cut the sample into nths and do the exercise w/ and w/o hot tract
  nsplit = 2
  xs = fltarr(nsplit,2)
  ys = fltarr(nsplit,2,2)
  for ii = 0, 1 do begin
     case ii of
        0: begin
           tmt = meanTots
           trt = rtcomp
        end
        1: begin
           tmt = meanTots[norm]
           trt = rtcomp[norm]
        end
     endcase

     ;; add the triplicates in, which thankfully do not cover the hot
     ;; tract
     
     
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

  plotsym, 0, /fill
  plot, meanTots, nsp, psym = 8, $
        xtitle = '(!18n!X!D1!N+!18n!X!D2!N)!18/!X2 [avg. total people + dwellings]', $
        ytitle = 'sqrt(!18n!X!D1!N-!18n!X!D2!N)!E2!N!18/!Xsqrt(2)', $
        /nodat, yr = [0,2*sqrt(max(meanTots))], xr = [0,60]
  x = findgen(!X.CRANGE[1]-!X.CRANGE[0]) + !X.CRANGE[0]
;  oplot, meanTots[norm], nsp[norm], psym = 8, symsize = 1, col = '777777'x
  oplot, meanTots, nsp, psym = 8, symsize = 0.7
  oplot, x, sqrt(x), col = 255
  oploterror, xs[*,0], ys[*,0,0], ys[*,1,0], psym = 8, $
              errcol = '00a5ff'x, col = '00a5ff'x, symsize = 2
  oploterror, xs[*,1], ys[*,0,1], ys[*,1,1], psym = 8, $
              errcol = 'ffa500'x, col = 'ffa500'x, symsize = 2
  legend, /top, /left, box = 0, $
          ['tract-level', 'mean', 'mean excl 1901.00', 'Poisson'], $
          psym = [8,8,8,0], linesty = [0,0,0,0], $
          col = ['ffffff'x,'00a5ff'x,'ffa500'x,'0000ff'x], $
          symsize = [0.7,1.5,1.5,1], pspacing = 0.5, $
          charsize = 1.25

  stop
  
;  print, 'P C V R T M'
;  print, sqrt(mean(rcomp^2, dim = 2))/sqrt(2)
;  print, stddev(rcomp, dim = 2, /nan)

  ;; square of the differences per category
  sqd = sqrt(mean(rcomp^2, dim = 2, /nan))/sqrt(2)

  ;; expected poisson noise from avg. occupancy per tract
  sigd = sqrt(mean(means, dim = 2, /nan))
  
  plot, findgen(ncat), sqd, xr = [-1,6], $
        xtickname = [' ', 'P', 'C', 'V', 'R', 'T', 'M', ' '], $
        xticks = 7, xtickint = 1, $
        ytitle = 'sqrt[(n!D1!N-n!D0!N)!E2!N]/sqrt(2)', xminor = 1
  oplot, findgen(ncat), sigd, linesty = 2

  stop
  
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

  stop
  
  d = findgen(101)-50
  del = bb[*,1] - bb[*,0] ;; duplicates
  del = [del, bbb - bb[trips,0]];, bbb-bb[trips,1]] ;; with triplicates to 0th and 1st measurement
  ebar = sqrt(bb[*,0] + bb[*,1])
  ebar = [ebar, sqrt(bbb + bb[trips,0])];, sqrt(bbb + bb[trips,1])]

  pdf = fltarr(n_elements(d))
;  npdf = fltarr(n_elements(d))
  for ii = 0, n_elements(ebar) - 1 do $
     pdf += (1./sqrt(2*!pi*ebar[ii]^2) * exp(-0.5 * (d-del[ii])^2/ebar[ii]^2))>0
;  for ii = 0, n_elements(ebar) - 1 do $
;     npdf += (1./sqrt(2*!pi*del[ii]^2) * exp(-0.5 * (d-del[ii])^2/del[ii]^2))>0
  
;  plot, d/mean(ebar), pdf/max(pdf), $
;        xtitle = greek('Delta')+'!18/<!X'+greek('sigma')+'!18>!X', $
;        ytitle = 'prob', xr = [-3,3]
;  null = exp(-0.5 * (d*0.1)^2)
;  oplot, d*(0.1), null/max(null), col = 255
  
;  stop
  
end
