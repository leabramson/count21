pro hoTractLookUp, inTracts

  ntracts = n_elements(inTracts)
  
  eHoTracts = [1916.10, $
               1916.20, $
               1905.20, $
               1911.10, $
               1911.20, $
               1925.10, $
               1925.20, $
               1926.20, $
               1926.10, $
               1915.00, $
               1912.04, $
               1912.03, $
               1912.01, $
               1913.02, $
               1913.01, $
               1914.10, $
               1914.20, $
               1927.00]

  hoTracts = [1905.10, $
              1898.00, $
              1899.02, $
              1899.03, $
              1899.04, $
              1899.05, $
              1901.00, $
              1902.01, $
              1902.02, $
              1903.01, $
              1907.00, $
              1908.01, $
              1908.02, $
              1909.01, $
              1909.02, $
              1910.00, $
              1917.10, $
              1917.20, $
              1918.10, $
              1918.20, $
              1919.01]

  tracts = [eHoTracts, hoTracts]
  
  eHoHit = intarr(n_elements(eHoTracts))
  hoHit  = intarr(n_elements(HoTracts))
  totHit = intarr(n_elements(tracts))
  for ii = 0, nTracts - 1 do begin
     t = float(inTracts[ii])

     hit = where(tracts eq t, nhit)
     if nhit gt 0 then $
        totHit[hit] += 1

     if nhit gt 0 then begin
        hit = where(hoTracts eq t, nHoHit)
        if nHoHit gt 0 then $
           hoHit[hit] += 1 $
        else $
           eHoHit[where(eHoTracts eq t)] +=1
     endif
     
  endfor

  qui = where(totHit gt 1, nqui)

  print, ''
  print, f = '(%"Total tracts counted .............. %i of %i")', $
         total(totHit gt 0), n_elements(tracts)
  print, f = '(%" > Tracts with multiple counters .. %i")', $
         total(totHit gt 1)
  if nqui gt 0 then for ii = 0, nqui - 1 do $
     print, f = '(%"      %7.2f")', tracts[qui[ii]]
  print, f = '(%"Total E. Hollywood tracts counted . %i of %i")', $
         total(eHoHit gt 0), n_elements(eHotracts)
  print, f = '(%"Total Hollywood tracts counted .... %i of %i")', $
         total(HoHit gt 0), n_elements(Hotracts)
  
end

