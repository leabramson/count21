pro countRegion, datafile, niter, $
                 output = output, $
                 wtTargets = wtTargets, $
                 RegName = regName, $
                 peaks = peaks, $
                 stdErrs = stderrs

  ;; SPA 4 CVRTM weights & Stats
  ;; https://www.lahsa.org/documents?id=4693-2020-greater-los-angeles-homeless-count-cvrtm-conversion-factors
  if not keyword_Set(PEAKS) then $
     peaks = [1.38,1.68,1.32,1.45,1.64] ;; 2020
  nobj  = [521., 694., 735., 2126., 1235.]
  nppl  = [674., 1117., 930., 3050., 1990.]
  se    = [50., 141., 102., 122., 190.]

  if NOT keyword_set(STDERRS) then $
     stderrs = [0.11, 0.22, 0.15, 0.06, 0.16] ;; 2020 SPA4
  ;; [0.47,1.08,1.26,0.84,0.93] CD 13 2020
  if NOT keyword_set(wtTargets) then $
     wtTargets = peaks + 2 * stderrs
  
  ;; Read in count data and find unique tracts and Ncounters
  data = mrdfits(dataFile, 1)

    ;; cull bad data
  if total(tag_names(data) eq "FLAG") gt 0 then begin
     bad = where(data.FLAG, compl = good, nbad)
     if nbad gt 0 then $
        for ii = 0, nbad-1 do $
           print, f = '(%"Bad data detected at tract %7.2f by counter %s")', $
                  data[bad[ii]].TRACT, data[bad[ii]].COUNTER
     data = data[good]
  endif

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
  bkgRates = (sqrt(total(input, 2) / n_elements(input[0,*]))); < 1
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
     base                = baseCounts[*,ii] # replicate(1,niter)
     bkg                 = bkgRates # replicate(1, niter) 
     baseErr             = (sqrt(base/nhit)>bkg) $
                           * randomn(seed, [nStruct,niter])
;     baseErr             = (sqrt(base/nhit)>1) $
;                           * randomn(seed, [nStruct,niter])
     finalCounts[*,*,ii] = ((base + baseErr) * wts);>0
  endfor

;  stop
  
  ;; Produce a summary file
  savedata = {TRACT: '0', $
              NCOUNTERS: 0, $
              REGION: regName, $
              COUNTS: fltarr(nStruct,niter), $
              RAWCOUNTS: fltarr(nStruct), $
              TAGS: ['Adults','TAY','Minors',$
                     'Cars', 'Vans', 'RVs', $
                     'Tents', 'Makeshifts', 'Families'], $
              WTS: [1,1,1,peaks,1], $
              WTERRS: [0,0,0,$
                       carSig,vanSig,rvSig,$
                       tentSig,makeSig,0]}
  savedata = replicate(savedata, nUtracts)
  for ii = 0, nutracts - 1 do begin
     savedata[ii].TRACT     = utracts[ii]
     savedata[ii].NCOUNTERS = ncounters[ii]
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

  d  = mrdfits(resfile, 1)
  nlines = n_elements(d)

  ;; Print summary stats
  region = d[0].REGION

  means = mean(d.COUNTS, dim = 2)
  errs  = stddev(d.COUNTS, dim = 2)
  
  tot    = total(means)
  totErr = 2*sqrt(total(errs^2))

  if n_elements(struct) gt 1 then begin
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
    
end

;;
;;
;;

pro makePlots, resFile, $
               outdir = outdir
    

  data     = mrdfits(resfile, 1)
  nLines   = n_elements(data)
  if NOT keyword_set(OUTDIR) then begin
     outdir = data[0].REGION
     spawn, 'mkdir '+outdir
  endif
     
  structs = strcompress(data[0].TAGS,/rem)
  nstruct = n_elements(structs)
  
  bs = 2
  
  ;; Do the big bar plot first
  if n_elements(data) gt 1 then $
     counts  = total(data.COUNTS, 3) $
  else $
     counts  = data.COUNTS
  tc      = total(counts, 1) ## replicate(1,nstruct)
  csums   = arrstats(counts/tc)


  ;; Make text summaries
  dumpPremade, data, outdir+'/breakdown_summary.dat', region = data[0].REGION

  ;; Plot
  outboxes  = [csums.P05, csums.P25, csums.P50, csums.P75, csums.P95]

  set_plot, 'PS'

  !x.thick = 4
  !y.thick = 4
  !p.charsize = 1.25
  !p.charthick = 4

  allFillCol = 'aaaaaa'x
  allCol     = 0
 
  device, filename = outDir+'/allBreakdownBar.eps', $
          /col, /encap, /decomp, bits_per_pix = 8
  cgbarplot, csums.P50, $
             col = allFillcol, $
             yr = [0, max(outboxes)], /ysty, $
             ymin = 5, barcoord = allx, $
             title = 'Unsheltered by Dwelling', $
             ytitle = 'people'
  oploterror, allx, csums.P50, csums.P95 - csums.P50, $
              /hibar, thick = 6, /nohat, errcol = allCol, psym = 3
  oploterror, allx, csums.P50, csums.P50 - csums.P05, $
              /lobar, thick = 6, /nohat, errcol = allCol, psym = 3
  for ii = 0, nstruct - 1 do $
     cgtext, allx[ii], !Y.CRANGE[0]-0.03, structs[ii], $
             align = 0.7, /data, orien = 20
  legend, /top, /right, box = 0, $
          [data[0].Region, $
           '90% confidence interval'], $
          psym = [8,0], col = [allFillCol, allcol], $
          pspacing = 0.5, linesty = [0,0], thick = 6
  device, /close
  spawn, 'open '+outDir+'/allBreakdownBar.eps &'

  ;; Now do the regional and tract-level subcounts.
  d = data
  region = d[0].REGION
  outdir = outDir

  ;; sum over tracts
  if n_elements(d) gt 1 then $
     counts = total(d.COUNTS, 3) $
  else $
     counts = d.COUNTS
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
          pspacing = 2, charsize = 1.1, charthick = 4, thick = 6, $
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
;           xr = [0,200], $
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
  
  set_plot, 'X'

end

;;
;;
;;

pro runmcw, csv

  cd5Errs = [0.24,0.32,0.26,0.11,0.23]
  
  fitscount, csv, 'mcw/unmerged_countMidCity2021.fits'
  mergesplits, 'mcw/unmerged_countMidCity2021.fits', 'mcw/countMidCity2021.fits'
  countRegion, 'mcw/countMidCity2021.fits', 1d4, $
               output = 'mcw/countMidCityResults2021.fits', $
               regName = 'Mid City West', $
               peaks = getMCWwts(), $
               stderrs = cd5errs

  data = mrdfits('mcw/countMidCityResults2021.fits', 1)
  d = mrdfits('mcw/2020actualCounts.fits', 1)
  d = trans2020Official(d, wts = getMCWwts())  
  
  ly = total(d.TOT_PPL)
  printNeedToKnow, data, lastYear = ly, region = 'Mid City West'
  plotHist, data, out = 'mcw/mcw2021Hist.eps', compval = ly
  plotTractTractOfficial, 'mcw/countMidCityResults2021.fits', $
                          oldData = 'mcw/2020actualCounts.fits', outdir = 'mcw/'

  r = transpose([[d.TOT_IND], [d.C], [d.V], [d.R], [d.T], [d.M]])
  plotBarNewOld, data, [total(r, 2),0], output = 'mcw/bars.eps', mo = 'Mar'

end

;;
;;
;;

pro runit, csv
  fitsCount, csv, 'epSandbox/epSandbox.fits'
  countRegion, 'epSandbox/epSandbox.fits', 1d4, $
               output = 'epSandboxResults.fits', regName = 'Echo Pk Lake Area Fall 2020', $
               peaks = [1.38,1.68,1.32,1.45,1.64] ;;
  summarize, 'epSandboxResultsEmptyT.fits'
  
;  makeplots, 'epSandboxResultsT1.fits', outdir = 'EpSandbox'
;  findNullWeights, 'test.fits', lastYear
end

;runit, 'epSandbox/epSandbox.csv'

;;
;;
;;

pro epTentExp

  ntrials = 3.
  emptyfrac = [0,0,0,11./30,0] ;; From EP census 21 feb
  emptyTentErr = 1.96 * sqrt(11.*19./30.^2)/sqrt(30)/2 ;; 1-sigma
  csv = 'epSandbox/epSandbox.csv'
  ef = fltarr(ntrials+1)
  sign = [-1,0,1]
  spawn, 'rm epSandboxResults??.fits'
  spawn, 'rm epSandbox/epSandbox??.fits'
  for ii = 0, ntrials - 1 do begin
     emptyf = emptyfrac
;     emptyf[3] += randomn(Seed, 1) * emptyTentErr/2.
     emptyf[3] += sign[ii] * emptyTentErr
     ef[ii] = emptyf[3]
     fitsCount, csv, 'epSandbox/epSandbox'+string(ii,f='(I02)')+'.fits', $
                emptyFrac = emptyF
     countRegion, 'epSandbox/epSandbox'+string(ii,f='(I02)')+'.fits', 1d4, $
                  output = 'epSandboxResults'+string(ii,f='(I02)')+'.fits', $
                  regName = 'Echo Pk Lake Area Fall 2020', $
                  peaks = [1.38,1.68,1.32,1.45,1.64]

;     stop
     
  endfor
  countRegion, 'epSandbox/epSandbox.fits', 1d4, $
               output = 'epSandboxResults99.fits', $
               regName = 'Echo Pk Lake Area Fall 2020', $
               peaks = [1.38,1.68,1.32,1.1,1.64] ;; HwoodWts
  
  spawn, 'ls epSandboxResults*.fits > test.list'
  readcol, 'test.list', files, f = 'A'
  nfiles = n_elements(files)

  out = fltarr(5,nfiles)
  pctles = [0.05,0.25,0.5,0.75,0.95]
  lastYear = 174.
  for ii = 0, nfiles - 1 do begin
     d = mrdfits(files[ii], 1)
     cts = total(d.COUNTS,1)
     cts = cts[sort(cts)]
     out[*,ii] = cts[ceil(n_elements(cts) * pctles) - 1]
     print, files[ii], getcountprob(cts, lastyear)
  endfor

  plot, [0,nfiles], minmax(out), /nodat, $
        xtickname = replicate(' ', 60), yr = [100,300], $
        xtickint = 1, yminor = 5, ytitle = 'Unsheltered persons', $
        xminor = 1
  oplot, !X.CRANGE, replicate(lastYear,2), thick = 4, linesty = 5         
  xxx = !X.CRANGE[[0,0,1,1]]
  yyy = out[*,0]
  polyfill, xxx, yyy[[0,4,4,0]], col = 'ffa500'x
  polyfill, xxx, yyy[[1,3,3,1]], col = 'ff5500'x
  oplot, !X.CRANGE, yyy[2] * [1,1], col = 'ff0000'x
  cgtext, 0.025, yyy[4]-10, /data, "SPA-4 2020 CVRTM", col = 'ff0000'x, $
          charthick = 2
  plotsym, 0, /fill
  out = out[*,1:*]
  cgloadct, 33, ncol = nfiles, clip = [10,240], /rev
  for ii = 0, nfiles - 2 do begin
     x = ii+1
     col = cgcolor(string(fix(ii)))
     oploterror, x, out[2,ii], out[4,ii]-out[2,ii], /hibar, psym = 8, $
                 col = col, errcol = col
     oploterror, x, out[2,ii], out[2,ii]-out[0,ii], /lobar, $
                 col = col, errcol = col
     if ii le 2 then $
        cgtext, x, !Y.CRANGE[0]-10, /data, string(ef[ii],f='(F4.2)'), align = 0.5 $
     else $
        cgtext, x, !Y.CRANGE[0]-10, /data, "T=1.1 (Hwood)", align = 0.5
  endfor
  cgtext, nfiles*0.975,lastyear-10,/data, "last year's estimate", align = 1

  stop
  

end

;;
;;
;;

