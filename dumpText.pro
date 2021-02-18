pro dumpText, structure, output

  tract = structure.TRACT

  nitem = n_elements(structure.TAGS)  
  cts = structure.COUNTS

  ans   = mean(cts, dim = 2)
;  errs = stddev(cts, dim = 2)
  fracs = ans / total(ans)
  ferrs = fltarr(nitem)
  for ii = 0, nitem - 1 do $
     ferrs[ii] = 1.96 * sqrt(ans[ii]*(total(ans)-ans[ii])/total(ans))/total(ans) ;; binomial
  errs = ferrs * total(ans)
  
  outs = []
  for ii = 0, n_elements(ans) - 1 do $
     outs = [outs, ans[ii], errs[ii]]
  fouts = []
  for ii = 0, n_elements(ans) - 1 do $
     fouts = [fouts, fracs[ii], ferrs[ii]]
  
  close, 1
  openw, 1, output, width = 1024
  printf, 1, f = '(%"#COUNT RESULTS FOR TRACT %7.2f ")', tract
  printf, 1, f = '(%"#%s %s %s %s %s %s %s %s %s")', structure.TAGS
  printf, 1, f = '(%"%i %i %i %i %i %i %i %i %i")', structure.RAWCOUNTS
  printf, 1, f = '(%"%i(%i) %i(%i) %i(%i) %i(%i) %i(%i) %i(%i) %i(%i) %i(%i) %i(%i) %i(%i)")', round(outs)
  printf, 1, f = '(%"%4.2f(%4.2f) %4.2f(%4.2f) %4.2f(%4.2f) %4.2f(%4.2f) %4.2f(%4.2f) %4.2f(%4.2f) %4.2f(%4.2f) %4.2f(%4.2f) %4.2f(%4.2f)")', fouts
  close, 1
  

end
