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
  
  printHeader
  
  close, 1
  openw, 1, output, width = 1024
  for ii = 0, ntracts - 1 do begin
     trs = getCountProb(cts[*,ii], [0.05,0.5,0.95], /inv)
     printf, 1, f = '(%"%7.2f & %s & %s & %i & %i & %i--%i \\")', $
             tract[ii], comm[ii], pf[ii], $
             ncounters[ii], trs[1], trs[0] > 0, trs[-1]
  endfor
  close, 1

  spawn, "cat "+output+" >> header.txt"
  spawn, "cat footer.txt >> header.txt"
  
end

