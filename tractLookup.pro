pro tractLookUp, inTracts, tractfile

  readcol, tractfile, tracts, f = 'F', /silent
  tracts = tracts[sort(tracts)]

  ntracts = n_elements(inTracts)
  
  totHit = intarr(n_elements(tracts))
  check  = intarr(ntracts)
  for ii = 0, nTracts - 1 do begin
     t = float(inTracts[ii])
     hit = where(tracts eq t, nhit)
     if nhit gt 0 then $
        totHit[hit] += 1 $
     else $
        check[ii] += 1
  endfor

  qui = where(totHit gt 1, nqui)  
  ext = where(check, nExt)

  uintracts = inTracts[sort(inTracts)]
  uintracts = uintracts[UNIQ(uintracts)]
  
  print, ''
  print, f = '(%"Total tracts counted .............. %i")', $
         n_elements(uintracts)
  print, f = '(%"Total in-domain tracts counted .... %i of %i")', $
         total(totHit gt 0), n_elements(tracts)
  print, f = '(%" > Tracts with multiple counters .. %i")', $
         total(totHit gt 1)
  if nqui gt 0 then for ii = 0, nqui - 1 do $
     print, f = '(%"      %7.2f")', tracts[qui[ii]]
  if nExt gt 0 then begin
     print, f = '(%"Out-domain tracts ................. %i")', $
            nExt
     for ii = 0, total(check) - 1 do $
        print, f = '(%"      %7.2f")', tracts[ext[ii]]
  endif
  
end
