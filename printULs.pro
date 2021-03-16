pro printULs

  data = mrdfits('countHollywoodResults2021w191902.fits', 1)
  e = where(data.EASTFLAG, compl = h)

  hc = total(data[h].COUNTS, 3)
  ec = total(data[e].COUNTS, 3)

  hs = arrstats(hc)
  es = arrstats(ec)

  print, 'HOLLYWOOD'
  for ii = 0, n_elements(data[0].TAGS) - 1 do $
     print, f = '(%"%s: <%f")', data[0].TAGS[ii], hs[ii].P95
  print, 'EAST HOLLYWOOD'
  for ii = 0, n_elements(data[0].TAGS) - 1 do $
     print, f = '(%"%s: <%f")', data[0].TAGS[ii], es[ii].P95
  
end
