pro printHeader  

  close, 2
  openw, 2, 'header.txt', width = 1024
  printf, 2, "\begin{table}[]"
  printf, 2, "\caption{Census Tract-level Unsheltered Data}"
  printf, 2, "\resizebox{\textwidth}{!}{%"
  printf, 2, "\begin{tabular}{ccccc}"
  printf, 2, "\toprule"
  printf, 2, "Tract & Team & $n_{\rm teams}$ & Median Est. & 90\% CI \\ \cmidrule{1-5}"
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

  printHeader
  
  close, 1
  openw, 1, output, width = 1024
  for ii = 0, ntracts - 1 do begin
     trs = summarizeRegion(data[ii])
     printf, 1, f = '(%"%7.2f & %s & %i & %i & %i--%i \\")', $
             tract[ii], 'V', ncounters[ii], trs[3], trs[0], trs[-1]
     print, tract[ii], eastflag[ii]
  endfor
  close, 1

  spawn, "cat "+output+" >> header.txt"
  spawn, "cat footer.txt >> header.txt"
  
end

