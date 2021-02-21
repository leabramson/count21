pro dumpPremade, results, output, $
                 REGION = region

  if n_elements(results) gt 1 then begin
     raw = total(results.RAWCOUNTS, 2)
     sum = total(results.COUNTS, 3)
  endif else begin
     raw = results.RAWCOUNTS
     sum = results.COUNTS
  endelse
  ans = mean(sum, dim = 2)
;  errs = stddev(sum, dim = 2)
  ferrs = fltarr(n_elements(ans))
  for ii = 0, n_elements(ans) - 1 do $
     ferrs[ii] = 1.96 * sqrt(ans[ii]*(total(ans)-ans[ii])/total(ans))/total(ans) ;; binomial
  errs = ferrs * total(ans)
  frac = ans / total(ans)

  outs = []
  for ii = 0, n_elements(ans) - 1 do $
     outs = [outs, ans[ii], errs[ii]]
  fouts = []
  for ii = 0, n_elements(ans) - 1 do $
     fouts = [fouts, frac[ii], ferrs[ii]]
  
  close, 1
  openw, 1, output, width = 1024
  printf, 1, f = '(%"#COUNT RESULTS FOR %s ")', region
  printf, 1, f = '(%"#%s %s %s %s %s %s %s %s %s")', results[0].TAGS
  printf, 1, f = '(%"%i %i %i %i %i %i %i %i %i")', raw
  printf, 1, f = '(%"%i(%i) %i(%i) %i(%i) %i(%i) %i(%i) %i(%i) %i(%i) %i(%i) %i(%i) %i(%i)")', round(outs)
  printf, 1, f = '(%"%4.2f(%4.2f) %4.2f(%4.2f) %4.2f(%4.2f) %4.2f(%4.2f) %4.2f(%4.2f) %4.2f(%4.2f) %4.2f(%4.2f) %4.2f(%4.2f) %4.2f(%4.2f)")', fouts
  close, 1
  

end
