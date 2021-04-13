pro printHeader  

  close, 2
  openw, 2, 'header.txt', width = 1024
  printf, 2, "\begin{table}[]"
  printf, 2, "\caption{Census Tract-level Unsheltered Data}"
  printf, 2, "\resizebox{\textwidth}{!}{%"
  printf, 2, "\begin{tabular}{cccccc}"
  printf, 2, "\toprule"
  printf, 2, "Tract & Community & Counter & Passes & Median Est. & 90\% CI \\ \cmidrule{1-6}"
  close, 2

  openw, 2, 'footer.txt', width = 1024
  printf, 2, "\\ \bottomrule"
  printf, 2, "\end{tabular}"
  printf, 2, "}"
  printf, 2, "\caption*{}"
  printf, 2, "\label{tbl:}"
  printf, 2, "\end{table}"
  close, 2
  
end

;;
;;
;;

pro makeTable, resFile, output

  data      = mrdfits(resfile, 1)
  ntracts   = n_elements(data)
  ncounters = data.NCOUNTERS
  eastflag  = data.EASTFLAG
  tract     = data.TRACT

  pts = idProTracts(data.TRACT)
  pf = replicate('V', ntracts)
  pf[where(pts)] = 'P'

  comm = replicate('H', ntracts)
  comm[where(eastFlag)] = 'E'

  cts = reform(total(data.COUNTS, 1))
  sum = getCountProb(total(cts,2), [0.05,0.5,0.95], /inv)

  printHeader
  
  close, 1
  openw, 1, output, width = 1024
  for ii = 0, ntracts - 1 do begin
     trs = getCountProb(cts[*,ii], [0.05,0.5,0.95], /inv)
     printf, 1, f = '(%"%7.2f & %s & %s & %i & %i & %i--%i \\")', $
             tract[ii], comm[ii], pf[ii], $
             ncounters[ii], round(trs[1]), trs[0] > 0, trs[-1]
  endfor
  printf, 1, f = '(%" & & & %i & %i & %i--%i")', $
          total(ncounters), sum[1], sum[0], sum[2]
  close, 1

  spawn, "cat "+output+" >> header.txt"
  spawn, "cat footer.txt >> header.txt"
  
end

;;
;;
;;

pro makeBigTables, resfile

  data      = mrdfits(resfile, 1)
  ntracts   = n_elements(data)
  eastflag  = data.EASTFLAG
  tract     = data.TRACT

  pts = idProTracts(data.TRACT)
  pf = replicate('V', ntracts)
  pf[where(pts)] = 'P'

  comm = replicate('H', ntracts)
  comm[where(eastFlag)] = 'E'

  cts = data.COUNTS
  rcts = data.RAWCOUNTS

  close, 2
  openw, 2, 'header.txt', width = 1024
  printf, 2, "\begin{table}[]"
  printf, 2, "\caption{Census Tract-level Unsheltered Population Inferences}"
  printf, 2, "\resizebox{\textwidth}{!}{%"
  printf, 2, "\begin{tabular}{ccccccccccc}"
  printf, 2, "\toprule"
  printf, 2, "Tract & Community & Counter & $A$ & {\it TAY} & $C$ & $V$ & $R$ & $T$ & $M$ & {\bf Total}\\ \cmidrule{1-11}"
  close, 2
  
  close, 1
  openw, 1, 'ppl.txt', width = 1024
  for ii = 0, ntracts - 1 do begin
     c = cts[*,*,ii]
     trs = arrstats(c)
     err = (trs.P95 - trs.P05) / 2
     mid = trs.P50

     ct = total(c, 1)
     ct = ct[sort(ct)]

     mid[where(mid lt err)] = 0;trs[where(mid lt err)].P95
     printf, 1, f = '(%"%7.2f & %s & %s & %4.1f (%4.1f) & %4.1f (%4.1f) & %4.1f (%4.1f) & %4.1f (%4.1f) & % 4.1f (%4.1f) & %4.1f (%4.1f) & %4.1f (%4.1f) & %5.1f (%5.1f) \\")', $
             tract[ii], comm[ii], pf[ii], $
             mid[0], err[0], mid[1], err[1], $
             mid[3], err[3], mid[4], err[4], mid[5], err[5], mid[6], err[6], mid[7], err[7], $
             ct[0.5*n_elements(ct)-1]>0, (ct[0.95*n_elements(ct)-1] - ct[0.05*n_elements(ct)-1]) / 2
  endfor
  close, 1

  spawn, "cat ppl.txt >> header.txt"
  spawn, "cat footer.txt >> header.txt"
  spawn, "mv header.txt ppl.txt"

  close, 2
  openw, 2, 'header.txt', width = 1024
  printf, 2, "\begin{table}[]"
  printf, 2, "\caption{Census Tract-level Unsheltered Counts}"
  printf, 2, "\resizebox{\textwidth}{!}{%"
  printf, 2, "\begin{tabular}{ccccccccccc}"
  printf, 2, "\toprule"
  printf, 2, "Tract & Community & Counter & $A$ & {\it TAY} & $C$ & $V$ & $R$ & $T$ & $M$ & {\bf Total} \\ \cmidrule{1-11}"
  close, 2
  
  close, 1
  openw, 1, 'cts.txt', width = 1024
  for ii = 0, ntracts - 1 do begin
     mid = rcts[*,ii]
     printf, 1, f = '(%"%7.2f & %s & %s & %4.1f & %4.1f & %4.1f & %4.1f & %4.1f & %4.1f & %4.1f & %5.1f \\")', $
             tract[ii], comm[ii], pf[ii], $
             mid[[0,1,3,4,5,6,7]], total(mid[[0,1,3,4,5,6,7]])
  endfor
  close, 1

  spawn, "cat cts.txt >> header.txt"
  spawn, "cat footer.txt >> header.txt"
  spawn, "mv header.txt cts.txt"
  
  
end
