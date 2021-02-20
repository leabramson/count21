pro makeBreakoutTable, struct

  close, 2
  openw, 2, 'header.txt', width = 1024
  printf, 2, "\begin{table}[]"
  printf, 2, f = '(%"\caption{Tract %7.2f Unsheltered Data}")', $
          struct.TRACT
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
