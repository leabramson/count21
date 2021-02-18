function drawDist, x, pdf, n, $
                   normalize = normalize

  dx = x[1] - x[0]
  cdf = total(pdf * dx, /cum)
  
  if keyword_set(normalize) then cdf /= max(cdf)

  out = x[value_locate(CDF, randomu(seed, n))]

  RETURN, out
end

;;
;;
;;

function readIn, table

  readcol, table, $
           tract, adult, tay, minor, $
           car, van, rv, tent, makeshift, $
           f = "F,F,F,F,F,F,F,F,F", /preserve_null, $
           delim = ','

  output = {ADULT: 0., $
            TAY: 0., $
            MINOR: 0., $
            CAR: 0., $
            VAN: 0., $
            RV: 0., $
            TENT: 0., $
            MAKESHIFT: 0.}
  
  for ii = 0, 7 do begin
     case ii of
        0: data = adult
        1: data = tay
        2: data = minor
        3: data = car
        4: data = van
        5: data = rv
        6: data = tent
        7: data = makeshift
     endcase

     output.(ii) = total(data[0:6])
     
  endfor

  RETURN, output
end

;;
;;
;;

function genFunc, x, peak, sigma, cutoff

  func = 1./sqrt(sigma^2 * 2 * !pi) $
         * exp(-0.5 * (x-peak)^2/sigma^2)

  func[where(x lt cutoff)] = 0

  func /= total(func * (x[1]-x[0]))
  
  RETURN, func
end

;;
;;
;;

function findSig, x, peak, sigi, sigf, $
                  nsteps = nsteps, $
                  ntile = ntile, $
                  cutoff = cutoff, $
                  target = target

  if NOT keyword_set(ntile) then ntile = 0.95
  
  dsig = (sigf - sigi)/(nsteps-1)
  sigran = findgen(nsteps) * dsig + sigi
  nout = fltarr(nsteps)
  for ii = 0, nsteps - 1 do begin
     foo = genFunc(x, peak, sigran[ii], 1)
     cdf = total(foo, /cum)
     nout[ii] = x[value_locate(cdf, ntile * cdf[-1])]
  endfor

  sigma = (sigran[value_locate(nout, target)] > min(sigran)) < max(sigran)
  
  RETURN, sigma
end

;;
;;
;;

pro monteCount, datafile, n, $
                wtTargets = wtTargets
  

;  peaks = [1.39,1.66,1.91,1.26,1.80] ;; 2019
  peaks = [1.38,1.68,1.32,1.45,1.64] ;; 2020

  if NOT keyword_set(wtTargets) then $
     wtTargets = peaks + 2 * [0.1, 0.2, 0.15, 0.1, 0.16]; [1.5 * peaks
  
  data = readIn(dataFile)
  
  pplPer = findgen(1001)/200.

  carSig = findSig(pplPer, peaks[0], 0.1, 1.0, $
                   nsteps = 20, $
                   cutoff = 1, $
                   target = wtTargets[0])
  vanSig = findSig(pplPer, peaks[1], 0.1, 1.0, $
                   nsteps = 20, $
                   cutoff = 1, $
                   target = wtTargets[1])
  rvSig = findSig(pplPer, peaks[2], 0.1, 1.0, $
                  nsteps = 20, $
                  cutoff = 1, $
                  target = wtTargets[2])
  tentSig = findSig(pplPer, peaks[3], 0.1, 1.0, $
                    nsteps = 20, $
                    cutoff = 1, $
                    target = wtTargets[3])
  makeSig = findSig(pplPer,  peaks[4], 0.1, 1.0, $
                    nsteps = 20, $
                    cutoff = 1, $
                    target = wtTargets[4])

  print, carSig, vanSig, rvSig, tentSig, makeSig
;  stop
  
  carPDF  = genFunc(pplPer, peaks[0], carSig, 1) ;; 2020 wts = maxLike
  vanPDF  = genFunc(pplPer, peaks[1], vanSig, 1)
  rvPDF   = genFunc(pplPer, peaks[2], rvSig, 1)
  tentPDF = genFunc(pplPer, peaks[3], tentSig, 1)
  makePDF = genFunc(pplPer, peaks[4], makeSig, 1)

  cgloadct, 33, ncol = 5, clip = [10,240]
  plot, pplPer, carPdf / max(carpdf), /nodat, $
        xran = [0.8,3], yr = [0,1.1], /xsty, /ysty, $
        xtitle = 'Avg. people per structure', $
        ytitle = 'relative probability'
  oplot, pplPer, carPdf / max(carpdf), col = cgcolor('0')
  oplot, pplPer, vanPdf / max(vanpdf), col = cgcolor('1')
  oplot, pplPer, rvPdf / max(rvpdf), col = cgcolor('2')
  oplot, pplPer, tentPdf / max(tentpdf), col = cgcolor('3')
  oplot, pplPer, makePdf / max(makePdf), col = cgcolor('4')
  plotsym, 0, /fill
  for ii = 0, n_elements(peaks) - 1 do begin
     oplot, [peaks[ii]], [0], psym = 8, symsize = 2, $
            col = cgcolor(string(ii,f='(I1)'))
     oplot, [peaks[ii]], [!Y.CRANGE[1]], psym = 8, symsize = 2, $
            col = cgcolor(string(ii,f='(I1)'))
  endfor
  key = ['Car', 'Van', 'RV', 'Tent', 'Makeshift', 'LAHSA SPA4']
  col = [cgcolor('0'), cgcolor('1'), cgcolor('2'), $
         cgcolor('3'), cgcolor('4'), 'ffffff'x]
  lsty = replicate(0, 6)
  psym = [replicate(0, 5), 8]
  legend, /top, /right, $
          key, pspacing = 0.5, $
          col = col, linesty = lsty, psym = psym
  
  carWts  = drawDist(pplPer, carPDF , n, /norm)
  vanWts  = drawDist(pplPer, vanPDF , n, /norm)
  rvWts   = drawDist(pplPer, rvPDF  , n, /norm)
  tentWts = drawDist(pplPer, tentPDF, n, /norm)
  makeWts = drawDist(pplPer, makePDF, n, /norm)

  pplWts = replicate(1,n)

  wts = transpose([[pplWts],[pplWts],[pplWts],$
                   [carWts],[vanWts],[RVwts], $
                   [tentWts],[makeWts]])

  input = transpose([[data.ADULT],[data.TAY],[data.MINOR],$
                     [data.CAR],[data.VAN],[data.RV],$
                     [data.TENT],[data.MAKESHIFT]])
  input = input # replicate(1,n)
  
  output = (input + randomn(seed, [8,n]) * (sqrt(input)>1)) $
           * wts

  totals = total(output, 1)
  totals = totals[sort(totals)]

  bs = 5
  h = histogram(totals, bins = bs, loc = bins)
  mode = bins[where(h eq max(h))]
  
  set_plot, 'PS'
  device, filename = 'CHNCrecountDist.eps', $
          /col, /encap, /decomp, bits_per_pix = 8
  plot, bins, h / total(h), psym = 10, $
        ytitle = 'probability', $
        xtitle = 'inferred unsheltered people', $
        ysty = 8, $
        pos = [0.15,0.15,0.9,0.9], $
        xr = [400,600], $;totals[[0.05,0.95]*n-1] + [-20,20], $
        ymin = 2, /nodat, xthick = 4, ythick = 4, $
        charthick = 4, charsize = 1.25
  hf = histofill(bins, h/total(h), /pad)
  qui = where(hf.BINS ge totals[0.05*n-1] $
              and hf.BINS le totals[0.95*n-1])
  qui2 = where(hf.BINS ge totals[0.16*n-1] $
               and hf.BINS le totals[0.84*n-1])
  ty = findgen(n)/(n-1) * !Y.CRANGE[1]
;  dx = 1
;  tx = findgen((!X.CRANGE[1] - !X.CRANGE[0])/dx+1) * dx + !X.CRANGE[0]
;  ty = interpol(h/total(h), bins, tx)
;  qui = where(tx ge totals[0.05*n-1] $
;              and tx le totals[0.95*n-1])
;  qui2 = where(tx ge totals[0.16*n-1] $
;               and tx le totals[0.84*n-1])
  oplot, mode[0] * [1,1], !Y.CRANGE, col = '0000ff'x, thick = 8
  oplot, replicate(totals[n/2-1], 2), !Y.CRANGE, $
         col = '555555'x, thick = 8
  oplot, replicate(totals[0.16*n-1], 2), !Y.CRANGE, $
         col = '555555'x, linesty = 5, thick = 8
  oplot, replicate(totals[0.84*n-1], 2), !Y.CRANGE, $
         col = '555555'x, linesty = 5, thick = 8
  oplot, replicate(totals[0.05*n-1], 2), !Y.CRANGE, $
         col = '555555'x, linesty = 2, thick = 6
  oplot, replicate(totals[0.95*n-1], 2), !Y.CRANGE, $
         col = '555555'x, linesty = 2, thick = 6
  legend, pos = [!X.CRANGE[1]-5,0.048], /right, /clear, box = 0, $
          ['mode ('+string(mode[0],f='(I3)')+')', $
           'median ('+string(totals[n/2-1],f='(I3)')+')', $
           '16/84 pctle', '5/95 pctle', 'cumulative'], $
          linesty = [0,0,5,2,4], $
          col = ['0000ff'x, replicate('555555'x, 3), '00a5ff'x], $
          pspacing = 2, charsize = 1.2, charthick = 4, thick = 6, $
          spacing = 1.7
  polyfill, [hf.bins[qui[0]], hf.bins[qui], hf.bins[qui[-1]]], $
            [0,hf.hist[qui],0], col = 'ffa500'x, $
            /line_fill, spacing = 0.025, thick = 1, orien = 45
  polyfill, [hf.bins[qui2[0]], hf.bins[qui2], hf.bins[qui2[-1]]], $
            [0,hf.hist[qui2],0], col = 'ff0000'x, $
            /line_fill, spacing = 0.025, thick = 1, orien = -45
  oplot, bins, h / total(h), psym = 10, thick = 8
  axis, yaxis = 1, yr = [0,1], /ysty, col = '00a5ff'x, $
        ytitle = '!18P!X(<!18X!X)', ythick = 4, $
        charthick = 4, charsize = 1.25
  oplot, totals, ty, col = '00a5ff'x, thick = 6, linesty = 4
  device, /close
  set_plot, 'X'
  spawn, 'open CHNCrecountDist.eps &'

  tags = tag_Names(data)
  
  print, f = '(%"\nMedian: %i")', totals[n/2-1]
  print, f = '(%"Mode  : %i")', mode
  print, f = '(%"68\% CI: %i--%i")', totals[0.16*n-1], totals[0.84*n-1]
  print, f = '(%"90\% CI: %i--%i")', totals[0.05*n-1], totals[0.95*n-1]
  print, f = '(%"90\% CL: +/-%i (%i\%)")', $
         (totals[0.05*n-1] - totals[0.95*n-1])/(-2), $
         (totals[0.05*n-1] - totals[0.95*n-1])/(-2)/totals[n/2-1]*100
  print, f = '(%"\n >>> DATA BY CLASS <<<\n")'
  for ii = 0, n_elements(tags) - 1 do $
     print, f = '(%"%9s: %3i +/- %2i")', tags[ii], data.(ii), $
            sqrt(data.(ii)) > 1
  print, f = '(%"\n >>> CONTRIBUTION BY CLASS (95\% CL) <<<\n")'
  for ii = 0, n_elements(tags) - 1 do $
     print, f = '(%"%9s: %2i\% +/- %2i\%")', tags[ii], $
            mean(output[ii,*]/total(output, 1)) * 100, $
            2 * stddev(output[ii,*]/total(output, 1)) * 100

;  pctles = [0.05,0.16,0.25,0.50,0.75,0.84,0.95]
;  RETURN, [totals[ceil(pctles * n) - 1], bins[where(h eq max(h))]]
  
end
; monteCount, 'chncRecount.csv', 10000, wtTarg = [2,2,2.5,2,2.5]

;;

pro do10

  answers = fltarr(8,4)
  answers[*,0] = monteCount('chncRecount.csv', 10000)
  answers[*,1] = monteCount('chncRecount.csv', 10000, $
                            wtTargets = [2,2,2.5,2,2.5])
  answers[*,2] = monteCount('chncRecount.csv', 10000, $
                            wtTargets = replicate(2, 5))
  answers[*,3] = monteCount('chncRecount.csv', 10000, $
                            wtTargets = replicate(2.5,5))
  
end
