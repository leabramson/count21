;; To de-bias the aggregate totals while allowing individual tracts to
;; not be characterized by negative numbers

function scaleToMean, struct

  rawCts = struct.RAWCOUNTS
  wts    = [1,1,1,struct[0].WTS,1]
  catRaw = total(rawCts, 2)
  raw    = total(catRaw * wts)

  catTot = total(struct.COUNTS, 3)
  cts    = total(catTot, 1)
  tot    = median(cts)
  catTot = median(catTot, dim = 2)
  
  norm  = tot / raw ;; median_tot / mean_tot
  cnorm = catTot / catRaw
  
  output = {GLOBAL: 1./norm, $
            CATEGORY: 1./cnorm}
  
;  print, norm
  
  return, output
    
end
