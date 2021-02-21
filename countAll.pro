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

pro countAll, datafile, niter, $
              output = output, $
              wtTargets = wtTargets, $ ;; DEPRECATED
              Peaks = Peaks

  ;; SPA 4 CVRTM weights & Stats
  ;; https://www.lahsa.org/documents?id=4693-2020-greater-los-angeles-homeless-count-cvrtm-conversion-factors
  if NOT keyword_set(PEAKS) then $
     peaks = [1.38,1.68,1.32,1.45,1.64] ;; 2020 SPA 4
;  peaks = [1.10,1.16,1.74,1.40,1.22] ;; 2020 CD13
  nobj  = [521., 694., 735., 2126., 1235.]
  nppl  = [674., 1117., 930., 3050., 1990.]
  se    = [50., 141., 102., 122., 190.]
  
  ;; Standard errors on the means
  ;; They seem to find the weight by N_tot/N_obj, and then use the
  ;; "standard error" on N_tot to define the error on the weight.
  ;; The error on N_tot, however, is like 2-4x sqrt(N_tot),
  ;; so one could assume the 95% limit is "1" standard error on the mean...
  ;; I dunno; I'll default 2 for now.
  stderrs = [0.11, 0.22, 0.15, 0.06, 0.16] ;; 2020 SPA4
  ;; [0.47,1.08,1.26,0.84,0.93] CD 13 2020
  if NOT keyword_set(wtTargets) then $
     wtTargets = peaks + 2 * stderrs
                  
  ;; Read in count data and find unique tracts and Ncounters
  data = mrdfits(dataFile, 1)
  nlines = n_elements(data)
  tracts = data.TRACT
  utracts = tracts[UNIQ(tracts)]
  nutracts = n_elements(utracts)
  ncounters = fltarr(nutracts)
  for ii = 0, nutracts - 1 do $
     ncounters[ii] = total(tracts eq utracts[ii])
  
  ;; Generate the PDFs on the weights
  pplPer = findgen(1001)/200.

  carSig  = stderrs[0]
  vanSig  = stderrs[1]
  rvSig   = stderrs[2]
  tentSig = stderrs[3]
  makeSig = stderrs[4]
  
;  carSig = findSig(pplPer, peaks[0], 0.05, 1.0, $
;                   nsteps = 200, $
;                   cutoff = 1, $
;                   target = wtTargets[0])
;  vanSig = findSig(pplPer, peaks[1], 0.05, 1.0, $
;                   nsteps = 200, $
;                   cutoff = 1, $
;                   target = wtTargets[1])
;  rvSig = findSig(pplPer, peaks[2], 0.05, 1.0, $
;                  nsteps = 200, $
;                  cutoff = 1, $
;                  target = wtTargets[2])
;  tentSig = findSig(pplPer, peaks[3], 0.01, 1.0, $
;                    nsteps = 2000, $
;                    cutoff = 1, $
;                    target = wtTargets[3])
;  makeSig = findSig(pplPer,  peaks[4], 0.05, 1.0, $
;                    nsteps = 200, $
;                    cutoff = 1, $
;                    target = wtTargets[4])
;
;  print, carSig, vanSig, rvSig, tentSig, makeSig
  
  ;; Draw niter realizations for the weights from the PDFs
  carPDF  = genFunc(pplPer, peaks[0], carSig, 1) ;; 2020 wts = maxLike
  vanPDF  = genFunc(pplPer, peaks[1], vanSig, 1)
  rvPDF   = genFunc(pplPer, peaks[2], rvSig, 1)
  tentPDF = genFunc(pplPer, peaks[3], tentSig, 1)
  makePDF = genFunc(pplPer, peaks[4], makeSig, 1)

  ;; Plot to see if they're sensible distributions
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

  carWts  = drawDist(pplPer, carPDF , niter, /norm)
  vanWts  = drawDist(pplPer, vanPDF , niter, /norm)
  rvWts   = drawDist(pplPer, rvPDF  , niter, /norm)
  tentWts = drawDist(pplPer, tentPDF, niter, /norm)
  makeWts = drawDist(pplPer, makePDF, niter, /norm)

  ;; Compile everything into master arrays
  pplWts = replicate(1,niter)
  wts = transpose([[pplWts],[pplWts],[pplWts],$
                   [carWts],[vanWts],[RVwts], $
                   [tentWts],[makeWts],[pplWts]])
  input = transpose([[data.ADULT],[data.TAY],[data.MINOR],$
                     [data.CAR],[data.VAN],[data.RV],$
                     [data.TENT],[data.MAKESHIFT], [data.FAMILY]])
  nStruct = n_elements(input[*,0])

  ;; Establish a base error rate by distributing all raw counts evenly
  ;; over the entire region, but ceil it to 1 since we're
  ;; pretty good at identifying structures as opposed to ppl.
  foo = fltarr(nstruct, nutracts)
  for ii = 0, nutracts - 1 do begin
     hit = where(tracts eq utracts[ii], nhit)
     if nhit gt 1 then $
        foo[*,ii] = mean(input[*,hit], dim = 2) $
     else $
        foo[*,ii] = input[*,hit]
  endfor
  bkgRates = (sqrt(total(foo, 2) / n_elements(foo[0,*]))) ;< 1
  bkgRates[where(bkgRates eq 0)] = min(bkgRates[where(bkgRates gt 0)])

  ;; Count, looping over unique tracts and averaging where
  ;; mutiple counters.
  baseCounts  = fltarr(nStruct, nutracts)
  finalCounts = fltarr(nStruct, niter, nutracts)
  for ii = 0, nutracts - 1 do begin
     hit = where(tracts eq utracts[ii], nhit)
     if nhit gt 1 then $
        baseCounts[*,ii] = mean(input[*,hit], dim = 2) $
     else $
        baseCounts[*,ii] = input[*,hit]

     baseErr = sqrt(baseCounts[*,ii] / nhit) ;; poison err on mean
     baseErr[where(baseErr eq 0)] = bkgRates[where(baseErr eq 0)] ;; fill nulls w/ bkg
     base    = baseCounts[*,ii] # replicate(1,niter)     
     bErr    = (baseErr # replicate(1,niter)) * randomn(seed, [nStruct,niter]) ;; poisson noise, accounting for multiple counters, with the bkg Rate
;     stop
;     baseErr             = (sqrt(base/nhit)>1) $
;                           * randomn(seed, [nStruct,niter])
     finalCounts[*,*,ii] = ((base + bErr) * wts)>0 ;; a draw from the underlying real distribution boosted by a draw from the weight PDF
  endfor

  ;; Produce a summary file
  savedata = {TRACT: '0', $
              NCOUNTERS: 0, $
              EASTFLAG: 0, $
              COUNTS: fltarr(nStruct,niter), $
              RAWCOUNTS: fltarr(nStruct), $
              TAGS: ['Adults','TAY','Minors',$
                     'Cars', 'Vans', 'RVs', $
                     'Tents', 'Makeshifts', 'Families'], $
              WTS: pks, $
              WTERRS: [0,0,0,$
                       carSig,vanSig,rvSig,$
                       tentSig,makeSig,0]}
  savedata = replicate(savedata, nUtracts)
  for ii = 0, nutracts - 1 do begin
     savedata[ii].TRACT     = utracts[ii]
     savedata[ii].NCOUNTERS = ncounters[ii]
     savedata[ii].EASTFLAG  = eHoLookup(utracts[ii])
     savedata[ii].COUNTS    = finalCounts[*,*,ii]
     savedata[ii].RAWCOUNTS = baseCounts[*,ii]
  endfor
  mwrfits, savedata, output, /create
  
  print, ''
  print, f = '(%">>> Results output to %s <<<\n")', output
   
end
;; countAll, 'test.fits', 1000, output = 'testResults.fits'

;;
;;
;;

pro summarize, resfile

  data     = mrdfits(resfile, 1)
  eHoDat   = data[where(data.EASTFLAG, compl = Hwood, nEho)]
  hWoodDat = data[Hwood]
  nHwood   = n_elements(hWoodDat)

  ;; Print summary stats
  
  for ii = 0, 2 do begin
     case ii of
        0: begin $
           d = data
           region = 'GREATER HOLLYWOOD'
        end
        1: begin
           d = eHoDat
           region = 'EAST HOLLYWOOD'
        end
        2: begin
           d = hWoodDat
           region = 'HOLLYWOOD'
        end
     endcase

     means = mean(d.COUNTS, dim = 2)
     errs  = stddev(d.COUNTS, dim = 2)

     tot    = total(means)
     totErr = 2*sqrt(total(errs^2))

     s = size(means, /dim)
     if n_elements(s) gt 1 then begin
        totStruct     = total(means, 2)
        totStructErrs = 2*sqrt(total(errs^2, 2))
     endif else begin
        totStruct     = means
        totStructErrs = 2*errs
     endelse
     
     tstring = "*** SUMMARY FOR "+region+" (95% conf.) ***"
     print, ''
     print, tstring
     print, replicate('-', strlen(tstring)/2)
     print, f = '(%"Total: %i+/-%i")', round(tot), round(totErr)
     for jj = 0, n_elements(d[0].TAGS) - 1 do $
        print, f = '(%" > %s: %i+/-%i (%4.1f\%)")', d[0].TAGS[jj], $
               round(totStruct[jj]), round(totStructErrs[jj]), $
               round(totStruct[jj]/tot*100)
     
  endfor

  
;  print, f = '(%"\nMedian: %i")', totals[n/2-1]
;  print, f = '(%"Mode  : %i")', mode
;  print, f = '(%"68\% CI: %i--%i")', totals[0.16*n-1], totals[0.84*n-1]
;  print, f = '(%"90\% CI: %i--%i")', totals[0.05*n-1], totals[0.95*n-1]
;  print, f = '(%"90\% CL: +/-%i (%i\%)")', $
;         (totals[0.05*n-1] - totals[0.95*n-1])/(-2), $
;         (totals[0.05*n-1] - totals[0.95*n-1])/(-2)/totals[n/2-1]*100
;  print, f = '(%"\n >>> DATA BY CLASS <<<\n")'
;  for ii = 0, n_elements(tags) - 1 do $
;     print, f = '(%"%9s: %3i +/- %2i")', tags[ii], data.(ii), $
;            sqrt(data.(ii)) > 1
;  print, f = '(%"\n >>> CONTRIBUTION BY CLASS (95\% CL) <<<\n")'
;  for ii = 0, n_elements(tags) - 1 do $
;     print, f = '(%"%9s: %2i\% +/- %2i\%")', tags[ii], $
;            mean(output[ii,*]/total(output, 1)) * 100, $
;            2 * stddev(output[ii,*]/total(output, 1)) * 100
  

end

;;
;;
;;

pro makePlots, resFile

  data     = mrdfits(resfile, 1)
  eHoDat   = data[where(data.EASTFLAG, compl = Hwood, nEho)]
  hWoodDat = data[Hwood]
  nHwood   = n_elements(hWoodDat)

  structs = strcompress(data[0].TAGS,/rem)
  nstruct = n_elements(structs)
  outBoxes = fltarr(nstruct*4,5)
  
  bs = 2
  
  ;; Do the big bar plot first
  counts  = total(data.COUNTS, 3)
  tc      = total(counts, 1) ## replicate(1,nstruct)
  csums   = arrstats(counts/tc)
  hCounts = total(hWoodDat.COUNTS, 3)
  htc     = total(hcounts, 1) ## replicate(1,nstruct)
  hCsums  = arrstats(hcounts/htc)
  eCounts = total(eHoDat.COUNTS, 3)
  etc     = total(ecounts, 1) ## replicate(1,nstruct)
  eCsums  = arrstats(eCounts/etc)

  ;; Make text summaries
  dumpPremade, data, 'allTracts/breakdown_summary.dat', region = 'Greater Hollywood'
  dumpPremade, eHoDat, 'eHo/breakdown_summary.dat', region = 'East Hollywood'
  dumpPremade, hWoodDat, 'hWood/breakdown_summary.dat', region = 'Hollywood'       

  ;; Plot
  inds = findgen(nstruct)
  outboxes[4*inds,*]   = [csums.P05, csums.P25, csums.P50, $
                          csums.P75, csums.P95]
  outboxes[4*inds+1,*] = [hcsums.P05, hcsums.P25, hcsums.P50, $
                          hcsums.P75, hcsums.P95]
  outboxes[4*inds+2,*] = [ecsums.P05, ecsums.P25, ecsums.P50, $
                          ecsums.P75, ecsums.P95]

  set_plot, 'PS'
  device, filename = 'allTracts/allBreakdown.eps', $
          /col, /encap, /decomp, bits_per_pix = 8

  !x.thick = 4
  !y.thick = 4
  !p.charsize = 1.25
  !p.charthick = 4

  allFillCol = 'aaaaaa'x
  allCol     = 0
  hoFillCol  = 'ffa500'x
  hoCol      = 'ff0000'x
  eHoFillCol = '00a5ff'x
  eHoCol     = long('0000ff'x)
  
  plot, [-3,4*nstruct], [0,1.2*max(outboxes)], $
        /nodat, xtickname = replicate(' ', 60), $
        /xsty, /ysty, $
        xminor = 1, yminor = 5, $
        xtickv = 4*inds+1, xticks = nstruct-1, $
        ytitle = 'unsheltered fraction', $
        title = 'Unsheltered by Dwelling'
  for ii = 0, nstruct - 1 do $
     cgtext, 4*inds[ii]+1, !Y.CRANGE[0]-0.03, structs[ii], $
             align = 0.7, /data, orien = 20
  width = ((!X.CRange[1] - !X.Crange[0]) / (40))
  for ii = 0, 2 do begin
     case ii of
        0: begin
           use = 4*inds
           bcol = allFillCol
           ocol = allCol
        end
        1: begin
           use++
           bcol = hoFillCol
           ocol = hoCol
        end
        2: begin
           use++
           bcol = eHoFillCol
           ocol = eHoCol
        end
     endcase
     cgboxplot, outboxes[use, *], $
                xloc = use, /over, width = width, $
                /fillbox, boxCol = bcol
     cgboxplot, outboxes[use, *], $
                xloc = use, /over, width = width, $
                col = ocol, boxthick = 6
  endfor
  plotsym, 8, /fill
  legend, /top, /left, $
          ['Greater Hollywood', 'Hollywood', 'E. Hollywood'], $
          psym = 8, col = [allFillCol,hoFillCol,eHoFillCol], $
          pspacing = 0.5
  device, /close
;  spawn, 'open allTracts/allBreakdown.eps &'

  device, filename = 'allTracts/allBreakdownBar.eps', $
          /col, /encap, /decomp, bits_per_pix = 8
  cgbarplot, csums.P50, $
             barwidth = 0.25, $
             barspace = 0.75, $
             baroffset = 1, $
             col = allFillcol, $
             yr = [0, max(outboxes)], /ysty, $
             ymin = 5, barcoord = allx, $
             title = 'Unsheltered by Dwelling'
  oploterror, allx, csums.P50, csums.P95 - csums.P50, $
              /hibar, thick = 6, /nohat, errcol = allCol, psym = 3
  oploterror, allx, csums.P50, csums.P50 - csums.P05, $
              /lobar, thick = 6, /nohat, errcol = allCol, psym = 3
  cgbarplot, hcsums.P50, $
             barwidth = 0.25, $
             barspace = 0.75, $
             baroffset = 2, $
             col = hoFillcol, /over, barcoord = hox
  oploterror, hox, hcsums.P50, hcsums.P95 - hcsums.P50, $
              /hibar, thick = 6, /nohat, errcol = hoCol, psym = 3
  oploterror, hox, hcsums.P50, hcsums.P50 - hcsums.P05, $
              /lobar, thick = 6, /nohat, errcol = hoCol, psym = 3
  cgbarplot, ecsums.P50, $
             barwidth = 0.25, $
             barspace = 0.75, $
             baroffset = 3, $
             col = ehoFillcol, /over, barcoord = ehox
  oploterror, ehox, ecsums.P50, ecsums.P95 - ecsums.P50, $
              /hibar, thick = 6, /nohat, errcol = ehoCol, psym = 3
  oploterror, ehox, ecsums.P50, ecsums.P50 - ecsums.P05, $
              /lobar, thick = 6, /nohat, errcol = ehoCol, psym = 3
  for ii = 0, nstruct - 1 do $
     cgtext, hox[ii], !Y.CRANGE[0]-0.03, structs[ii], $
             align = 0.7, /data, orien = 20
  legend, /top, /left, box = 0, $
          ['Greater Hollywood', 'Hollywood', 'E. Hollywood', $
           '90% confidence interval'], $
          psym = [8,8,8,0], col = [allFillCol,hoFillCol,eHoFillCol,0], $
          pspacing = 0.5, linesty = [0,0,0,0], thick = 6
  device, /close
;  spawn, 'open allTracts/allBreakdownBar.eps &'

  ;; Now do the regional and tract-level subcounts.
  for ii = 0, 2 do begin
     case ii of
        0: begin $
           d = data
           region = 'Greater Hollywood'
           outdir = 'allTracts'
        end
        1: begin
           d = eHoDat
           region = 'East Hollywood'
           outdir = 'eHo'
        end
        2: begin
           d = hWoodDat
           region = 'Hollywood'
           outdir = 'hWood'
        end
     endcase

     ;; sum over tracts
     counts = total(d.COUNTS, 3)
     csums = arrstats(counts)

     ;; sum over dwellings
     tCounts = total(counts, 1)
     tCsum   = arrstats(tCounts)
     niter = n_elements(tCounts)
     
     h = float(histogram(tCounts, bins = bs, loc = bins))
     mode = bins[where(h eq max(h))]

     ;; plot an aggregate sum
     histName = outdir+'/'+strcompress(region, /rem)+'Dist.eps'     
     
     device, filename = histName, $
             /col, /encap, /decomp, bits_per_pix = 8
     plot, bins, h / total(h), psym = 10, $
           ytitle = 'probability', $
           xtitle = 'inferred unsheltered people', $
           ysty = 8, $
           pos = [0.15,0.15,0.9,0.9], $
           title = region, $
;        xr = [400,600], $;allTotal[[0.05,0.95]*n-1] + [-20,20], $
           ymin = 2, /nodat, xthick = 4, ythick = 4, $
           charthick = 4, charsize = 1.25
     hf = histofill(bins, h/total(h), /pad)
     qui = where(hf.BINS ge tcSum.P05 and hf.BINS le tcSum.P95)
     qui2 = where(hf.BINS ge tcSum.P25 and hf.BINS le tcSum.P75)
     ty = findgen(niter)/(niter-1) * !Y.CRANGE[1]
     oplot, mode[0] * [1,1], !Y.CRANGE, col = '0000ff'x, thick = 8
     oplot, replicate(tcSUm.P50, 2), !Y.CRANGE, $
            col = '555555'x, thick = 8
     oplot, replicate(tcSum.P25, 2), !Y.CRANGE, $
            col = '555555'x, linesty = 5, thick = 8
     oplot, replicate(tcSum.P75, 2), !Y.CRANGE, $
            col = '555555'x, linesty = 5, thick = 8
     oplot, replicate(tcSum.P05, 2), !Y.CRANGE, $
            col = '555555'x, linesty = 2, thick = 6
     oplot, replicate(tcSum.P95, 2), !Y.CRANGE, $
            col = '555555'x, linesty = 2, thick = 6
     legend, /top, /left, /clear, box = 0, $
             ['mode ('+string(mode[0],f='(I0)')+')', $
              'median ('+string(tcSum.P50,f='(I0)')+')', $
              '25/75 pctle ('+$
              string(tcSum.P25,f='(I0)')+'/'+$
              string(tcSum.P75,f='(I0)')+')', $
              '5/95 pctle ('+$
              string(tcSum.P05,f='(I0)')+'/'+$
              string(tcSum.P95,f='(I0)')+')', 'cumulative'], $
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
     oplot, tCounts[sort(tCounts)], ty, col = '00a5ff'x, $
            thick = 6, linesty = 4
     device, /close
     
     ;; Now tracts
     ntracts = n_elements(d)

     for jj = 0, ntracts - 1 do begin
        tract  = d[jj].TRACT
        tCounts = total(d[jj].COUNTS, 1)
        tCsum   = arrstats(tCounts)

        ;; Text summary
        dumpText, d[jj], outdir+'/'+strcompress(tract, /rem)+'_summary.dat'
        
        h = float(histogram(tCounts, bins = bs, loc = bins))
        mode = bins[where(h eq max(h))]

        histName = outdir+'/'+strcompress(tract, /rem)+'_Dist.eps'
        device, filename = histName, $
                /col, /encap, /decomp, bits_per_pix = 8
        plot, bins, h / total(h), psym = 10, $
              ytitle = 'probability', $
              xtitle = 'inferred unsheltered people', $
              ysty = 8, $
              pos = [0.15,0.15,0.9,0.9], $
              title = 'Tract '+strcompress(tract,/rem), $
              xr = [0,200], $
;        xr = [400,600], $;allTotal[[0.05,0.95]*n-1] + [-20,20], $
              ymin = 2, /nodat, xthick = 4, ythick = 4, $
              charthick = 4, charsize = 1.25
        hf = histofill(bins, h/total(h), /pad)
        qui = where(hf.BINS ge tcSum.P05 and hf.BINS le tcSum.P95)
        qui2 = where(hf.BINS ge tcSum.P25 and hf.BINS le tcSum.P75)
        ty = findgen(niter)/(niter-1) * !Y.CRANGE[1]
        oplot, mode[0] * [1,1], !Y.CRANGE, col = '0000ff'x, thick = 8
        oplot, replicate(tcSUm.P50, 2), !Y.CRANGE, $
               col = '555555'x, thick = 8
        oplot, replicate(tcSum.P25, 2), !Y.CRANGE, $
               col = '555555'x, linesty = 5, thick = 8
        oplot, replicate(tcSum.P75, 2), !Y.CRANGE, $
               col = '555555'x, linesty = 5, thick = 8
        oplot, replicate(tcSum.P05, 2), !Y.CRANGE, $
               col = '555555'x, linesty = 2, thick = 6
        oplot, replicate(tcSum.P95, 2), !Y.CRANGE, $
               col = '555555'x, linesty = 2, thick = 6
        legend, /top, /right, /clear, box = 0, $
                ['mode ('+string(mode[0],f='(I0)')+')', $
                 'median ('+string(tcSum.P50,f='(I0)')+')', $
                 '25/75 pctle ('+$
                 string(tcSum.P25,f='(I0)')+'/'+$
                 string(tcSum.P75,f='(I0)')+')', $
                 '5/95 pctle ('+$
                 string(tcSum.P05,f='(I0)')+'/'+$
                 string(tcSum.P95,f='(I0)')+')', 'cumulative'], $
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
        oplot, tCounts[sort(tCounts)], ty, col = '00a5ff'x, $
               thick = 6, linesty = 4
        device, /close
     endfor
     
  endfor
  set_plot, 'X'

  
end

;;
;;
;;

pro runit, csv
  fitsCount, csv, 'test.fits'
  data = mrdfits('test.fits', 1)
  hoTractLookup, data.TRACT
  countAll, 'test.fits', 1d4, output = 'testResults.fits'
  summarize, 'testResults.fits'
  makeplots, 'testResults.fits'
  findNullWeights, 'test.fits', lastYear
end

;;
;;
;;

pro runitRetry, csv
  fitsCount, csv, '2020sandbox/retry2020_hwoodOnly.fits'
  data = mrdfits('2020sandbox/retry2020_hwoodOnly.fits', 1)
  hoTractLookup, data.TRACT
  countAll, '2020sandbox/retry2020_hwoodOnly.fits', 1d4, $
            output = 'retryHwood2020Results.fits';, $
;            peaks = [1.38,1.68,1.32,1.12,1.64]
  summarize, 'retryHwood2020Results.fits'
  makeplots, 'retryHwood2020Results.fits'
  lastyear = 0.9 * 1067
  findNullWeights, '2020sandbox/retry2020_hwoodOnly.fits', lastYear
end
;runitRetry, '2020sandbox/retry2020_hwoodOnly.csv'

;;
;;
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
