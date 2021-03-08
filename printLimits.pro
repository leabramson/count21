pro printLimits, struct

  cts = struct.COUNTS

  if n_elements(struct) gt 1 then $
     cts = total(cts, 3)
  lims = arrstats(cts)

  print, struct[0].TAGS
  print, lims.P95
  print, lims.P05
  

end
