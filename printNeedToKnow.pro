pro printNeedToKnow, struct, $
                     region = region, $
                     lastYear = lastYear

  region = strupcase(region)
  if region eq 'HWOOD' then $
     lastYear = 1058. $
  else if region eq 'EHO' then $
     lastYear = 656.
  
  raw       = struct.RAWCOUNTS
  rawTypes  = total(raw, 2)
  rawTracts = total(raw, 1)
  rawTot    = total(rawTracts)

  nTypes = n_elements(rawTypes)
  nTracts = n_elements(rawTracts)
  
  ppl       = struct.COUNTS
  pplTypes  = total(ppl, 3)
  pplTracts = total(ppl, 1)
  pplTot    = total(pplTracts, 2)

  pctles = [0.05,0.25,0.50,0.75,0.95]
  pplTotSummary = getCountProb(pplTot, pctles, /inv)
  pplTotProb    = getCountProb(pplTot, lastYear)
  pplTotFracs   = [pplTotSummary[2]/lastYear, $
                   0.5*(pplTotSummary[4]-pplTotSummary[0])/lastYear]

  typeSummary  = fltarr(nTypes ,5)
  tractSummary = fltarr(nTracts,5)
  for ii = 0, ntypes - 1 do $
     typeSummary[ii,*] = getCountProb(pplTypes[ii,*], pctles, /inv)
  for ii = 0, ntracts - 1 do $
     tractSummary[ii,*] = getCountProb(pplTracts[*,ii], pctles, /inv)

  ;; print the stuff
  print, ''
  print, f = '(%" ------ SUMMARY for %s ------ ")', region
  print, ''
  print, f = '(%"Total People (90\%CI) . %5i+/-%i")', pplTotFracs * lastYear
  print, f = '(%"Fraction vs. last yr . %6.2f+/-%4.2f")', pplTotFracs
  print, f = '(%"Total counts.......... %5i")', rawTot
  for ii = 0, ntypes - 1 do $
     print, f = '(%" > %s: %3i (%2i\%), %3i, (%2i\%)")', $
            struct[0].TAGS[ii], rawTypes[ii], 100*rawTypes[ii]/rawTot, $
            typeSummary[ii,2], 100*typeSummary[ii,2]/pplTotSummary[2] 
  print, f = '(%"Min/max ppl/tract ..... %i, %i")', minmax(tractSummary[*,2])
  print, f = '(%"Min/max counts/tract .. %i, %i")', minmax(rawTracts)
  
  
end
