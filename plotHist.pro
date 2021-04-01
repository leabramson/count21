pro plotHist, struct, $
              region = region, $
              compval = compval, $
              output = output, $
              binsize = binsize

  if NOT keyword_set(binsize) then binsize = 5

  cts = total(total(struct.COUNTS, 3), 1)
  cts = cts[sort(cts)]
  h = histogram(cts, bins = binsize, loc = bins)
  ints = getCountProb(cts, [0.05,0.16,0.5,0.84,0.95], /inv)
  niter = n_elements(cts)

  cgloadct, 1, /brewer, ncol = 12

  s1col = cgcolor('10')
  s2col = cgcolor('3')
    
  set_plot, 'PS'
  !p.CHARTHICK = 5
  !p.CHARSIZE = 1.5
  !X.THICK = 4
  !Y.THICK = 4

  device, filename = output, $
          /col, /encap, /decomp, bits_per_pix = 8
  plot, bins, h / total(h), psym = 10, $
        ytitle = 'probability', $
        xtitle = 'inferred unsheltered people', $
        ysty = 8, $
        pos = [0.15,0.15,0.9,0.9], $
        title = region, $
;        xr = [400,600], $;allTotal[[0.05,0.95]*n-1] + [-20,20], $
        ymin = 2, /nodat, xthick = 4, ythick = 4
  hf = histofill(bins, h/total(h), /pad)
  qui = where(hf.BINS ge ints[0] and hf.BINS le ints[4])
  qui2 = where(hf.BINS ge ints[1] and hf.BINS le ints[3])
  ty = findgen(niter)/(niter-1) * !Y.CRANGE[1]
  oplot, replicate(ints[2], 2), !Y.CRANGE, $
         col = '555555'x, thick = 8
  oplot, replicate(ints[1], 2), !Y.CRANGE, $
         col = '555555'x, linesty = 5, thick = 8
  oplot, replicate(ints[3], 2), !Y.CRANGE, $
         col = '555555'x, linesty = 5, thick = 8
  oplot, replicate(ints[0], 2), !Y.CRANGE, $
         col = '555555'x, linesty = 2, thick = 6
  oplot, replicate(ints[4], 2), !Y.CRANGE, $
         col = '555555'x, linesty = 2, thick = 6
  if keyword_set(COMPVAL) then $
     oplot, compval * [1,1], !Y.CRANGE, col = 255, thick = 10
  if NOT keyword_set(compval) then $
     legend, /top, /left, /clear, box = 0, $
             ['median ('+string(ints[2],f='(I0)')+')', $
              '25/75 pctle ('+$
              string(ints[1],f='(I0)')+'/'+$
              string(ints[3],f='(I0)')+')', $
              '5/95 pctle ('+$
              string(ints[0],f='(I0)')+'/'+$
              string(ints[4],f='(I0)')+')', 'cumulative'], $
             linesty = [0,5,2,4], $
             col = [replicate('555555'x, 3), '00a5ff'x], $
             pspacing = 1, charsize = 1, charthick = 4, thick = 6, $
             spacing = 1. $
  else $
     legend, /top, /left, /clear, box = 0, $
             ['last year ('+string(compval,f='(I0)')+')', $, $
              'median ('+string(ints[2],f='(I0)')+')', $
              '16/84 pctle ('+$
              string(ints[1],f='(I0)')+'/'+$
              string(ints[3],f='(I0)')+')', $
              '5/95 pctle ('+$
              string(ints[0],f='(I0)')+'/'+$
              string(ints[4],f='(I0)')+')', 'cumulative'], $
             linesty = [0, 0,5,2,4], $
             col = [255, replicate('555555'x, 3), '00a5ff'x], $
             pspacing = 1, charsize = 1, charthick = 4, thick = 6, $
             spacing = 1.
  polyfill, [hf.bins[qui[0]], hf.bins[qui], hf.bins[qui[-1]]], $
            [0,hf.hist[qui],0], col = s2col, $
            /line_fill, spacing = 0.025, thick = 1, orien = 45
  polyfill, [hf.bins[qui2[0]], hf.bins[qui2], hf.bins[qui2[-1]]], $
            [0,hf.hist[qui2],0], col = s1col, $
            /line_fill, spacing = 0.025, thick = 1, orien = -45
  oplot, bins, h / total(h), psym = 10, thick = 8
  axis, yaxis = 1, yr = [0,1], /ysty, col = '00a5ff'x, $
        ytitle = '!18P!X(<!18X!X)', ythick = 4
  oplot, cts, ty, col = '00a5ff'x, thick = 6, linesty = 4
  device, /close
  set_plot, 'X'
  
end

pro plotStuff

  data = mrdfits('countHollywoodResults2021w191902.fits', 1)
  plotHist, data[where(~data.EASTFLAG)], $
            region = 'Hollywood', output = 'hwoodHist.eps', $
            compval = 1058.
  plotHist, data[where(data.EASTFLAG)], $
            region = 'East Hollywood', output = 'ehoHist.eps', $
            compval = 656.
  plotHist, data, $
            region = 'Greater Hollywood', output = 'allHoHist.eps', $
            compval = 656.+1058., bins = 10
  
  data = mrdfits('countHollywoodResults2021_revisedWtErr.fits', 1)
  plotHist, data[where(~data.EASTFLAG)], $
            region = 'Hollywood', output = 'hwoodHist_revErr.eps', $
            compval = 1058.
  plotHist, data[where(data.EASTFLAG)], $
            region = 'East Hollywood', output = 'ehoHist_revErr.eps', $
            compval = 656.

  
end
