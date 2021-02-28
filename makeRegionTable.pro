;; Sum over tracts

pro makeRegionTable, struct

  close, 2
  openw, 2, 'header.txt', width = 1024
  printf, 2, "\begin{table*}[]"
  printf, 2, f = '(%"\\caption{Unsheltered Data}")'
  printf, 2, "\resizebox{\linewidth}{!}{%"
  printf, 2, "\begin{tabular}{lcccccccccc}"
  printf, 2, "\toprule"
  printf, 2, $
          " & Adult & TAY & Unacc Minor & Car & Van & RV & Tent & Makeshift & Family & {\bf Total} \\ \cmidrule{1-11}"
  close, 2

  nitem = n_elements(struct[0].TAGS)  
  cts = total(struct.COUNTS, 3)

  ans   = mean(cts, dim = 2)
  fracs = ans / total(ans)
  ferrs = fltarr(nitem)
  for ii = 0, nitem - 1 do $
     ferrs[ii] = 1.96 * sqrt(ans[ii]*(total(ans)-ans[ii])/total(ans))/total(ans) ;; binomial
  errs = ferrs * total(ans)

  tot = arrstats(total(struct.COUNTS,3))
  out = tot.P50
  ci  = 0.5 * (tot.P95 - tot.P05) ;; 90%
  
  to = getCountProb(total(cts,1), [0.05,0.5,0.95], /inv)
  errs[where(tot.P50 le ci)] = tot[where(tot.P50 le ci)].P95
  
  outs = []
  for ii = 0, n_elements(ans) - 1 do $
     outs = [outs, ans[ii], errs[ii]]
  fouts = []
  for ii = 0, n_elements(ans) - 1 do $
     fouts = [fouts, fracs[ii], ferrs[ii]]

  
  close, 1
  openw, 1, 'tmp.tab', width = 1024
  printf, 1, f = '(%"Counts & %i & %i & %i & %i & %i & %i & %i & %i & %i & {\\bf %i} \\")', $
          total(struct.RAWCOUNTS,2), total(struct.RAWCOUNTS)
  printf, 1, f = '(%"Inhabitants & %i (%i) & %i (%i) & %i (%i) & %i (%i) & %i (%i) & %i (%i) & %i (%i) & %i (%i) & %i (%i) & {\\bf %i (%i)} \\")', $
          round(outs), round(to[1]), round(0.5*(to[2]-to[0]))
  printf, 1, f = '(%"Category share & %4.2f (%4.2f) & %4.2f (%4.2f) & %4.2f (%4.2f) & %4.2f (%4.2f) & %4.2f (%4.2f) & %4.2f (%4.2f) & %4.2f (%4.2f) & %4.2f (%4.2f) & %4.2f (%4.2f) & - ")', fouts
  close, 1
  
  openw, 2, 'footer.txt', width = 1024
  printf, 2, "\\ \bottomrule"
  printf, 2, "\end{tabular}"
  printf, 2, "}"
  printf, 2, "\caption*{Quantities in parentheses denote 95\% uncertainties (binomial in the case of the categories). Uncertainties larger than estimates imply that only upper limits can be stated confidently.}"
  printf, 2, "\label{tbl:}"
  printf, 2, "\end{table*}"
  close, 2

  spawn, "cat tmp.tab >> header.txt"
  spawn, "cat footer.txt >> header.txt"
 

end
