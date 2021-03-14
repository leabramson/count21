pro countAll, datafile, niter, $
              output = output, $
              wtTargets = wtTargets, $ ;; DEPRECATED
              Peaks = Peaks, $
              stdErrs = stderrs

  ;; SPA 4 CVRTM weights & Stats
  ;; https://www.lahsa.org/documents?id=4693-2020-greater-los-angeles-homeless-count-cvrtm-conversion-factors
  if NOT keyword_set(PEAKS) then $
     peaks = [1.38,1.68,1.32,1.45,1.64] ;; 2020 SPA 4
;  peaks = [1.10,1.16,1.74,1.40,1.22] ;; 2020 CD13
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

     baseErr = sqrt(baseCounts[*,ii] / nhit)                      ;; poison err on mean
     baseErr[where(baseErr eq 0)] = bkgRates[where(baseErr eq 0)] ;; fill nulls w/ bkg
     base    = baseCounts[*,ii] # replicate(1,niter)     
     bErr    = (baseErr # replicate(1,niter)) * randomn(seed, [nStruct,niter]) ;; poisson noise, accounting for multiple counters, with the bkg Rate
;     stop
;     baseErr             = (sqrt(base/nhit)>1) $
;                           * randomn(seed, [nStruct,niter])
     finalCounts[*,*,ii] = ((base + bErr) * wts);>0 ;; a draw from the underlying real distribution boosted by a draw from the weight PDF
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
              WTS: [1,1,1,peaks,1], $
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

pro makePlots, resFile

  data     = mrdfits(resfile, 1)
  eHoDat   = data[where(data.EASTFLAG, compl = Hwood, nEho)]
  hWoodDat = data[Hwood]
  nHwood   = n_elements(hWoodDat)

;  ;; debias the categories from the tract flooring
;  honorms  = scaleToMean(hWoodDat)
;  ehonorms = scaleToMean(eHoDat)
;  norms    = scaleToMean(data)
;  data.COUNTS *= norms.GLOBAL
;  ehodat.COUNTS *= ehoNorms.GLOBAL
;  hWoodDat.C.OUNTS *= hoNorms.GLOBAL
  
  structs = strcompress(data[0].TAGS,/rem)
  nstruct = n_elements(structs)
  outBoxes = fltarr(nstruct*4,5)
  
  bs = 2
  
  ;; Do the big bar plot first
  
  counts  = total(data.COUNTS, 3); * norms.GLOBAL
  tc      = total(counts, 1) ## replicate(1,nstruct)
  csums   = arrstats(counts/tc)
  hCounts = total(hWoodDat.COUNTS, 3); * honorms.GLOBAL
  htc     = total(hcounts, 1) ## replicate(1,nstruct)
  hCsums  = arrstats(hcounts/htc)
  eCounts = total(eHoDat.COUNTS, 3); * ehonorms.GLOBAL
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
             yr = [0, max(outboxes)], ysty=1, $
             ymin = 5, barcoord = allx, $
             title = 'Unsheltered by Dwelling', ytitle = 'People'
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
;           norms = norms.GLOBAL
        end
        1: begin
           d = eHoDat
           region = 'East Hollywood'
           outdir = 'eHo'
;           norms = ehoNorms.GLOBAL
        end
        2: begin
           d = hWoodDat
           region = 'Hollywood'
           outdir = 'hWood'
;           norms = hoNorms.GLOBAL
        end
     endcase

     ;; sum over tracts
     counts = total(d.COUNTS, 3)
     csums = arrstats(counts)

     ;; sum over dwellings
     tCounts = total(counts, 1); * norms
     tCsum   = arrstats(tCounts)
     niter = n_elements(tCounts)

     ;; deBias
;     norms = scaleToMean(d)
     
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
             pspacing = 1, charsize = 1, charthick = 4, thick = 6, $
             spacing = 1.
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

;        norms = scaleToMean(d[jj])

        ;; Text summary
        ;; THESE DO NOT GET DE-BIASED B/C WE DON'T WANT TO
        ;; QUOTE NEGATIVE OCCUPANCIES! YOU JUST HAVE TO NOTE THAT THE
        ;; TRACT- AND GLOBAL-GEOGRAPHY ESTIMATES ARE NOT CONSISTENT!

        dumpText, d[jj], outdir+'/'+strcompress(tract, /rem)+'_summary.dat'

        ;; Histogram
        
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

  cd13spa4 = [1.51,1.77,1.42,1.48,1.68]  ;; SPA4-derived! -- this should be the default
  cd13spa4errs = [0.25,0.42,0.28,0.11,0.31] ;; these might be 2 sigma...

  fitsCount, csv, 'countHollywood2021.fits'
  data = mrdfits('countHollywood2021.fits', 1)
  hoTractLookup, data.TRACT
  countAll, 'countHollywood2021.fits', 1d4, $
            output = 'countHollywoodResults2021.fits', $
            peaks = cd13Spa4, stderrs = cd13spa4errs
  summarize, 'countHollywoodResults2021.fits'
  makeplots, 'countHollywoodResults2021.fits'
  data = mrdfits('countHollywoodResults2021.fits', 1)
  plotBarNewOld, data[where(~data.EASTFLAG)], /hwood, $
                   output = 'Hwood2021Bars.eps'
  plotBarNewOld, data[where(data.EASTFLAG)], /eho, $
                   output = 'Eho2021Bars.eps'
  lastYear = 1058.
  td = data[where(~data.EASTFLAG)]
  mwrfits, td, 'hwood2021results.fits', /create
;  findNullWeights, 'hwood2021Results.fits', lastYear
  cts = total(total(td.COUNTS, 3), 1)
  print, getCountProb(cts, lastYear)

  lastYear = 656.
  td = data[where(data.EASTFLAG)]
  mwrfits, td, 'eho2021results.fits', /create
  cts = total(total(td.COUNTS, 3), 1)
  print, getCountProb(cts, lastYear)
  
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
  data = mrdfits('retryHwood2020Results.fits', 1)
  lastyear = 917.
  findNullWeights, data[where(~data.EASTFLAG)], lastYear
end
;runitRetry, '2020sandbox/retry2020_hwoodOnly.csv'

;;
;;
;;

pro multiWts

  tw = gen2021t() ;; from 2/28 data by brian and me; 38 tents, 47 total w/ unkns
  
  p2020 = [1.381, 1.682, 1.323, 1.453, 1.643] ;; SPA4 in toto
  p2020errs = [0.11, 0.22, 0.15, 0.06, 0.16]
  
  cd13spa4 = [1.51,1.77,1.42,1.48,1.68]  ;; SPA4-derived! -- this should be the default
  cd13spa4errs = [0.25,0.42,0.28,0.11,0.31]

  p2021 = cd13spa4
  p2021errs = cd13spa4errs
  p2021[3] = tw[0]
  p2021errs[3] = tw[1]
  
  p2021m = cd13spa4
  p2021m[3] = tw[2]
  p2021merrs = cd13spa4errs
  p2021merrs[3] = tw[3]

  
  countAll, 'countHollywood2021.fits', 1d4, $
            output = 'hwoodCVRTMtests/base2020.fits', $
            peaks = cd13spa4, stderrs = cd13spa4errs
  countAll, 'countHollywood2021.fits', 1d4, $
            output = 'hwoodCVRTMtests/spa42020.fits', $
            peaks = p2020, stderrs = p2020Errs
  countAll, 'countHollywood2021.fits', 1d4, $
            output = 'hwoodCVRTMtests/tent2021.fits', $
            peaks = p2021, stderrs = p2021errs
  countAll, 'countHollywood2021.fits', 1d4, $
            output = 'hwoodCVRTMtests/tmod2021.fits', $
            peaks = p2021m, stderrs = p2021merrs

    countAll, 'countHollywood2021.fits', 1d4, $
            output = 'countHollywoodResults2021_revisedWtErr.fits', $
            peaks = cd13Spa4, stderrs = cd13spa4errs/1.96
  
;  cd4Wts = [0.92,2.10,1.77,1.67,2.06]
;  cd4errs = [1.52,2.45,1.99,2.40,2.56]/1.96
;  hwoodWts = (cd13wts/cd13errs^2 + cd4wts/cd4errs^2) / (1/cd13errs^2 + 1/cd4errs^2)
;  hwoodErrs = sqrt((cd13Errs^2 + cd4Errs^2)/4)
;  cityWts = [1.58,1.90,1.64,1.45,1.47]
;  cityErrs = [0.21,0.29,0.27,0.09,0.21]/1.96
;  cd13wts = [1.10,1.16,1.74,1.40,1.22]         ;; LOCAL!
;  cd13errs = [0.47,1.08,1.26,0.84,0.93]/1.96

  spawn, 'rm hwoodCVRTMtests/*.fits'
  
  plotsym, 0, /fill
  cgbarplot, cd13spa4, barcoord = bx, $
             ytitle = 'CVRTM weights', barname = ['C','V','R','T','M'], $
             yran = [0,2.5], /ysty, col = 'ffffff'x
  oploterror, bx, cd13Spa4, cd13Spa4errs, symsize = 2, psym = 8
  oploterror, bx, p2021, p2021Errs, symsize = 1.2, $
              col = 'ffa500'x, errcol = 'ffa500'x, psym = 8
  oploterror, bx, p2021m, p2021merrs, symsize = 1.2, $
              col = '00a5ff'x, errcol = '00a5ff'x, psym = 8
;  oploterror, bx + replicate(randomn(seed, 1) * (bx[1]-bx[0])/5.,5), $
;              cityWts, cityErrs, symsize = 1.2, $
;              col = '00a5ff'x, errcol = '00a5ff'x, psym = 8
;  oploterror, bx + replicate(randomn(seed, 1) * (bx[1]-bx[0])/5.,5), $
;              cd13Wts, cd13Errs, symsize = 1.2, $
;              col = long('0000ff'x), errcol = long('0000ff'x), psym = 8
;  oploterror, bx + randomn(seed, 1) * (bx[1]-bx[0])/5., cd4Wts, cd4Errs, symsize = 1.2, $
;              col = 'bbbbbb'x, errcol = 'bbbbbb'x, psym = 8
;  oploterror, bx + randomn(seed, 1) * (bx[1]-bx[0])/5., hwoodWts, hwoodErrs, symsize = 1.2, $
;              col = 'ff0000'x, errcol = 'ff0000'x, psym = 8
  legend, /top, /left, box = 0, $
          ['CD13/SPA4', 'w/ New T', 'w/ Modeled T'], $;, 'CD4', 'Mean(CD4,13)'], $
          psym = 8, col = [0,'ffa500'x,'00a5ff'x], $;,'bbbbbb'x,'ff0000'x], $
          charsize = 1, textcol = 0

;  countAll, 'countHollywood2021.fits', 1d4, $
;            output = 'hwoodCVRTMtests/cd13.fits', $
;            peaks = cd13Wts, stderrs = cd13Errs
;  countAll, 'countHollywood2021.fits', 1d4, $
;            output = 'hwoodCVRTMtests/cd13spa4.fits', $
;            peaks = cd13spa4, stderrs = cd13spa4Errs
;  cd13spa4t = cd13spa4
;  cd13Spa4t[3] = 1.1
;  cd13spa4terrs = cd13spa4errs
;  cd13spa4terrs[3] = 0.07
;  countAll, 'countHollywood2021.fits', 1d4, $
;            output = 'hwoodCVRTMtests/13s4t.fits', $
;            peaks = cd13spa4t, stderrs = cd13spa4tErrs

end

;;
;;
;;

pro plotMultiWts, lastyear, $
                  region = region

  spawn, 'ls hwoodCVRTMtests/*.fits > test.list'
  readcol, 'test.list', files, f = 'A'
  nfiles = n_elements(files)

  if strupcase(region) eq 'HWOOD' then begin
     output = 'hwoodFinal.eps'
     lycol = 'ff5500'x
     stt = 0.75
  endif else if strupcase(region) eq 'EHO' then begin
     output = 'ehoFinal.eps'
     lycol = long('0055ff'x)
     stt = 0.65
  endif
  
  fn = []
  for ii = 0, nfiles - 1 do $
     fn = [fn, (strsplit(files[ii], '/',/extr))[1]]
  
  wtnames = strmid(fn, 0 ,4)

  out = fltarr(5,nfiles)
  means = fltarr(nfiles)
  pctles = [0.05,0.16,0.5,0.84,0.95]
  dcv = fltarr(nfiles)
  ps = fltarr(nfiles,2)
  for ii = 0, nfiles - 1 do begin

     d = mrdfits(files[ii], 1)
;     print, d[0].WTS
     
     if strupcase(region) eq 'HWOOD' then begin        
        d = d[where(d.EASTFLAG eq 0)]
        lastYear = 1058.
        if ii eq 0 then begin
           ;; boost the C & V numbers to 2020
           dc = 55 - total(d.RAWCOUNTS[3])
           dv = 48 - total(d.RAWCOUNTS[4])
        endif
     endif else if strupcase(region) eq 'EHO' then begin
        d = d[where(d.EASTFLAG)]
        lastYear = 656.
        if ii eq 0 then begin
           ;; boost the C & V numbers to 2020
           dc = 29 - total(d.RAWCOUNTS[3])
           dv = 58 - total(d.RAWCOUNTS[4])
        endif
     endif

     dcv[ii] = total([dc,dv] * d[0].WTS[[3,4]])
     
     cts = total(total(d.COUNTS,3),1)
     
     p1 = getCountProb(cts, lastYear)
     p2 = getCountProb(cts + dcv[ii],lastYear)

     ps[ii,*] = [p1,p2]
     
     cts = cts[sort(cts)]
     out[*,ii] = cts[ceil(n_elements(cts) * pctles) - 1]
     print, wtnames[ii], [out[*,ii], p1]
     print, wtnames[ii], [out[*,ii] + dcv[ii], p2]
     means[ii] = total(total(d.RAWCOUNTS, 2) * d[0].WTS)
  endfor
  null = where(wtnames eq '2020')
;  nullOut = out[*,null]
;  nullName = 'SPA4 2020'
  case strupcase(region) of
     'HWOOD': title = 'Hollywood CoC'
     'EHO': title = 'East Hollywood CoC'
  endcase
  
  print, [4,6] - [8,5.5]
  print, [0,4] - [0.5,3.5]

  lye = 0.7*lastYear/2
  tye = 0.5 * (out[3,*]-out[1,*])
  sig = sqrt(lye^2 + tye^2)

  cchange = gauss_pdf((lastYear - out[2,*])/sig)
  
  print, out[2,0] / lastYear, 0.5 * (out[4,0]-out[0,0]) / lastYear
  print, cchange
  
  key = ['2020!CSPA4-CD13', '2020!CSPA4', $
         '2021 !18w!DT!X!N', '2021 !18w!DT!X!N!Cnon-resp model']

  set_plot, 'PS'
  device, filename = output, $
          /col, /encap, /decomp, bits_per_pix = 8
  !X.THICK = 4
  !Y.THICK = 4
  !P.CHARSIZE = 1.25
  !P.CHARTHICK = 4
  !X.TICKLEN = 1d-6
  
  plotsym, 0, /fill
  plot, [-1,nfiles], minmax(out), /nodat, $
        xtickname = replicate(' ', 60), yr = lastyear * [stt,1.1], $
        xtickint = 1, yminor = 5, ytitle = 'Unsheltered people', $
        ysty = 8+1, title = title, pos = [0.13,0.2,0.85,0.9], $
        xr = [-0.5,nfiles-0.5], /xs, xminor = 1
  oplot, !X.CRANGE, lastYear * [1,1], thick = 10, col = lycol
  oplot, !X.CRANGE, out[2,0] * [1,1], linesty = 2, thick = 6, col = '777777'x
  axis, yaxis = 1, yr = (!y.CRANGE / lastYear - 1)*100, /ysty, $
        ytitle = '% change from 2020', yminor = 5
  cgloadct, 19, ncol = nfiles, clip = [47,200], /brewer, /rev 
  for ii = 0, nfiles - 1 do begin
     col = cgcolor(string(fix(ii)))
     x = ii
     oploterror, x, out[2,ii], out[4,ii]-out[2,ii], /hibar, psym = 8, $
                 symsize = 2.5, /nohat, thick = 5;, col = '555555'x, errcol = '555555'x
     oploterror, x, out[2,ii], out[2,ii]-out[0,ii], /lobar, $
                 symsize = 2, /nohat, thick = 5;, col = '555555'x, errcol = '555555'x
     oploterror, x, out[2,ii], out[3,ii]-out[2,ii], /hibar, psym = 8, $
                 symsize = 2.5, /nohat, thick = 10;, col = '555555'x, errcol = '555555'x
     oploterror, x, out[2,ii], out[2,ii]-out[1,ii], /lobar, $
                 symsize = 2, /nohat, thick = 10;, col = '555555'x, errcol = '555555'x
     oplot, [x], [out[2,ii]], psym = 8, symsize = 2, col = col
     cgtext, x, !Y.CRANGE[0]-22, /data, key[ii], align = 0.5
     if ii ne 0 then $
        cgtext, x, out[0,ii]-15, string(ps[ii,0]*100, f='(I0)')+'%', $
                charsize = 1, col = '777777'x, align = 0.5 $
     else $
        cgtext, x, out[0,ii]-20, string(ps[ii,0]*100, f='(I0)')+'%!Cchance of!Cdecrease', $
                charsize = 1, col = '777777'x, align = 0.5
  endfor
  plotsym, 0, thick = 6
  oplot, [0], [out[2,0]], psym = 8, symsize = 3.5, thick = 6
  cgtext, !X.CRANGE[0]+0.1,lastyear+15,/data, "LAHSA 2020 estimate", align = 0, col = lycol
  cgtext, 0+0.15, out[2,0]-22,/data, "Volunteer!C2021!Cbaseline", align = 0
  cgtext, 0+0.1, out[2,0]+20, /data, string((out[2,0]/lastYear-1)*100, f = '(I0)')+'%'+texToIdl('\pm')+$
          string(0.5*(out[4,0]-out[0,0])/lastYear*100, f='(I0)')+'%'
  cgtext, mean(!X.WINDOW), 0.05, /norm, '!18CVRTM!X weight choice', align = 0.5
  device, /close
  spawn, 'open '+output+' &'
  set_plot, 'X'
  
;  c = c[sort(c)]
;  h = histogram(c, min = 850, max = 1100, bins = 5, loc = bins)
;  plot, bins, h, psym = 10, xtitle = 'unsheltered people', ytitle = 'probability', xsty = 8+1
;  axis, xaxis = 1, xr = !Y.CRANGE/1058., /xsty
;  plot, bins, h, psym = 10, xtitle = 'unsheltered people', ytitle = 'probability', xsty = 8+1
;  axis, xaxis = 1, xr = !X.CRANGE/1058., /xsty                                               
;  oplot, 1058.*[1,1], !Y.CRANGE, col = 255
;  oplot, median(c) * [1,1], !Y.CRANGE, col = 'ff0000'x
;  foo = getCountProb(c, [0.05,0.25,0.75,0.95], /inv)
;  oplot, foo[-1] * [1,1], !Y.CRANGE, col = 'ffa500'x, linesty = 2
;  oplot, foo[-2] * [1,1], !Y.CRANGE, col = 'ffa500'x, linesty = 5
;  oplot, foo[1] * [1,1], !Y.CRANGE, col = 'ffa500'x, linesty = 5 
;  oplot, foo[0] * [1,1], !Y.CRANGE, col = 'ffa500'x, linesty = 2
  
end

;;
;;
;;

pro plotRange, region

  if strupcase(region) eq 'HWOOD' then $
     output = 'hwoodFinal.eps' $
  else if strupcase(region) eq 'EHO' then $
     output = 'ehoFinal.eps'

  d0 = mrdfits('hwoodCVRTMtests/base2020.fits', 1)
  d1 = mrdfits('hwoodCVRTMtests/spa42020.fits', 1)
  d2 = mrdfits('hwoodCVRTMtests/tent2021.fits', 1)
  d3 = mrdfits('hwoodCVRTMtests/tmod2021.fits', 1)

  case strupcase(region) of
     'HWOOD': begin
        lastYear = 1058.
        use = where(~d0.EASTFLAG)
        ;; C & V numbers from 2020
        c20 = 55
        v20 = 48
        title = 'Hollywood CoC'
        tscol = 'ffa500'x
        scol  = 'ff5500'x
        mcol  = 'ff0000'x
        lcol = 'ff0000'x
        ll = 0.7 & hl = 1.1
     end
     'EHO': begin
        lastYear = 656.
        use = where(d0.EASTFLAG)
        ;; C & V numbers from 2020
        c20 = 29
        v20 = 58
        title = 'East Hollywood CoC'
        tscol = 'ffa500'x
        scol  = long('0055ff'x)
        mcol  = long('0000ff'x)
        lcol = '0000ff'x
        ll = 0.6 & hl = 1.2
     end
  endcase

  d0 = d0[use]
  d1 = d1[use]
  d2 = d2[use]
  d3 = d3[use]
  
  ;; boost the C & V numbers to 2020
  dc = c20 - total(d0.RAWCOUNTS[3])
  dv = v20 - total(d0.RAWCOUNTS[4]) ;; raw counts are the same in d1 and d0

   ;; FOR LOOP GOES HERE
  
  dcv = [total([dc,dv] * d0.WTS[[3,4]]), $
         total([dc,dv] * d1.WTS[[3,4]]), $
         total([dc,dv] * d2.WTS[[3,4]]), $
         total([dc,dv] * d3.WTS[[3,4]]) $
        ]

  ;; total counts
  cts0 = total(total(d0.COUNTS,3),1)
  cts1 = total(total(d1.COUNTS,3),1)
  cts2 = total(total(d2.COUNTS,3),1)
  cts3 = total(total(d3.COUNTS,3),1)
  
  p0   = getCountProb(cts0, lastYear)
  p0cv = getCountProb(cts0 + dcv[0], lastYear)
  p1   = getCountProb(cts1, lastYear)
  p1cv = getCountProb(cts1 + dcv[1], lastYear)
  p2   = getCountProb(cts2, lastYear)
  p2cv = getCountProb(cts2 + dcv[0], lastYear)
  p3   = getCountProb(cts3, lastYear)
  p3cv = getCountProb(cts3 + dcv[1], lastYear)

  pctles = [0.05,0.16,0.50,0.84,0.95]
  
  cts0 = cts0[sort(cts0)]
  cts1 = cts1[sort(cts1)]
  cts2 = cts2[sort(cts2)]
  cts3 = cts3[sort(cts3)]
  
  out = fltarr(n_elements(pctles),2)
  out[*,0] = cts0[ceil(n_elements(cts0) * pctles) - 1]
  out[*,1] = cts1[ceil(n_elements(cts1) * pctles) - 1]

  print, 'CD13/SPA4 Baseline                        ', [out[[0,2,4],0], p0]
  print, 'CD13/SPA4 Baseline with 2020 CV           ', [out[[0,2,4],0]+dcv[0], p0CV]
  print, 'CD13/SPA4 Baseline with 2020 CV and 2021 T', [out[[0,2,4],1]+dcv[1], p1cv]
  print, 'CD13/SPA4 Baseline with 2021 T            ', [out[[0,2,4],1], p1]

  ps = [p0,p0cv,p1cv,p1]
  
  medians = [out[2,0], out[2,0]+dcv[0], out[2,1]+dcv[1], out[2,1]]
  maxs    = [out[4,0], out[4,0]+dcv[0], out[4,1]+dcv[1], out[4,1]]
  mins    = [out[0,0], out[0,0]+dcv[0], out[0,1]+dcv[1], out[0,1]]
  ih      = [out[3,0], out[3,0]+dcv[0], out[3,1]+dcv[1], out[3,1]]
  il      = [out[1,0], out[1,0]+dcv[0], out[1,1]+dcv[1], out[1,1]]

  set_plot, 'PS'
  device, filename = output, $
          /col, /encap, /decomp, bits_per_pix = 8
  !p.charthick = 4
  !p.charsize = 1.25
  !X.Thick = 4
  !Y.Thick = 4

  lye = sqrt(lastYear) * mean(d0.WTS)

  ;; find chance of being the same as 2020
;  run = findgen(1200)
;  pe = fltarr(n_elements(medians))
  err = 0.5 * (ih - il)
  print, 1-gauss_pdf(abs((medians - lastYear) / sqrt(err^2 + lye^2))) ;; chance the things are the same

  plotsym, 0, /fill
  plot, findgen(6), medians, yr = lastyear * [ll,hl], xr = [0.3,4.5], /xs, $
        xtickname = replicate(' ', 60), $
        xminor = 1, xtickint = 1, /nodat, ytitle = 'Unsheltered people', $
        yminor = 5, ysty = 8+1, title = title, pos = [0.13,0.15,0.85,0.9]
  axis, yaxis = 1, yr = (!y.CRANGE / lastYear - 1)*100, /ysty, $
        ytitle = '% change from 2020', yminor = 5
;  polyfill, !X.CRANGE[[0,1,1,0]], lastYear + lye * 2 * [-1,-1,1,1], $
;            col = 'aaaaaa'x, /line_fill, orien = 45, thick = 1, spacing = 0.025
;  polyfill, !X.CRANGE[[0,1,1,0]], lastYear + lye * [-1,-1,1,1], $
;            col = '777777'x, /line_fill, orien = -45, thick = 1, spacing = 0.025
  oplot, !X.CRANGE, replicate(lastYear, 2), thick = 10, col = lcol
  x = [1,2,3,4]
;  polyfill, [x,reverse(x)], [mins, reverse(maxs)], col = tscol, $
;            /line_fill, orien = -60, thick = 1, spacing = 0.025
;  polyfill, [x,reverse(x)], [ih, reverse(il)], col = scol, $
;            /line_fill, orien = 60, thick = 1, spacing = 0.025
;  oplot, x, medians, linesty = 2, thick = 6, col = mcol
  cols = [0,'0055ff'x,'0055ff'x,'ffa500'x]
  ecols = [0,'777777'x,'777777'x,0]
  ss   = [2,1.25,1.25,2]
  key = ['2020!C!18CVRTM!X', 'match!Ccar/van', 'match!Ccar/van!Cmod !18T!X', '2020!C!18CVRTM!X!Cmod !18T!X']
  for ii = 0, n_elements(medians) - 1 do begin
     oploterror, x[ii], medians[ii], (ih-medians)[ii], $
                 /hibar, psym = 3, errthick = 10, /nohat, errcol = ecols[ii]
     oploterror, x[ii], medians[ii], (maxs-medians)[ii], $
                 /hibar, psym = 3, errthick = 5, /nohat, errcol = ecols[ii]
     oploterror, x[ii], medians[ii], (medians-il)[ii], $
                 /lobar, psym = 3, errthick = 10, /nohat, errcol = ecols[ii]
     oploterror, x[ii], medians[ii], (medians-mins)[ii], $
                 /lobar, psym = 3, errthick = 5, /nohat, errcol = ecols[ii]

     oplot, [x[ii]], [medians[ii]], psym = 8, symsize = ss[ii], col = cols[ii]
     if ii eq 2 then $
        oplot, [x[ii]], [medians[ii]], psym = 8, symsize = 0.75, col = cols[ii+1]
     if ii ne 0 then $
        cgtext, x[ii], mins[ii]-15, string((1-ps[ii])*100, f='(I0)')+'%', $
                charsize = 1, col = '777777'x, align = 0.5 $
     else $
        cgtext, x[ii], mins[ii]-20, 'chance increase!Cover 2020: '+string((1-ps[ii])*100, f='(I0)')+'%', $
                charsize = 1, col = '777777'x, align = 0.5
     cgtext, x[ii], !y.CRANGE[0]-20, key[ii], align = 0.5, charsize = 1.25, charthick = 4
  endfor  
  legend, /top, /left, box = 0, /clear, $
          ['2020 LAHSA weights', 'modified weights', 'modified data'], $
          psym = 8, col = [0,'ffa500'x,'0055ff'x], $
          pspacing = 0.5, charsize = 1.1
  plotsym, 0, thick = 6
  cgtext, !X.CRANGE[1]-0.05, lastyear+25, /data, 'LAHSA 2020', align = 1
  oplot, x[[0,3]], medians[[0,3]], psym = 8, symsize = 3
  cgtext, x[0]+0.35, medians[0]-0.6*(medians-il)[0], 'Baseline', align = 0.5
  cgtext, x[3]-0.35, medians[3]-0.6*(medians-il)[3], 'Best!Cestimate', align = 0.5
  
;  oploterror, x, medians, maxs-medians, /hibar, psym = 8, errthick = 6, symsize = 2
;  oploterror, x, medians, medians-mins, /lobar, psym = 8, errthick = 6, symsize = 2

  device, /close
  spawn, 'open '+output+' &'
  set_plot, 'X'
  
  
end
