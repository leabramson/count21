pro compCHNC

  chncTracts = [1901.00, 1907.00, 1908.01, 1908.02, 1918.10, 1918.20, 1919.01]
  ntracts = n_elements(chncTracts)

  data = mrdfits('countHollywoodResults2021.fits', 1)
  raw = fltarr(9,ntracts)
  for ii = 0, ntracts - 1 do begin
     hit = where(data.TRACT eq chncTracts[ii])
     raw[*,ii] = data[hit].RAWCOUNTS
  endfor

  chncRaw = transpose([[43,59,34,13,9,5,29], $
                       [0,15,0,1,0,0,7], $
                       [0,0,0,0,0,0,0], $
                       [5,0,0,2,1,0,4], $
                       [2,1,1,5,5,3,1], $
                       [2,0,2,3,3,1,6], $
                       [5,22,15,18,12,1,11], $
                       [4,5,8,10,4,1,4], $
                       [0,0,0,0,0,0,0]])

  plot, chncRaw, raw, psym = 1, $
        xr = [0,80], yr = [0,80], /iso, $
        xtitle = 'July 2020 RAW', $
        ytitle = 'Feb 2021 RAW'
  one_one
  oploterror, chncRaw, raw, sqrt(chncRaw), sqrt(raw), $
              psym = 1, /nohat

  print, f = '(%"Abramson 7/2020 ... %i+/-%i")', total(chncRaw), sqrt(total(chncRaw))
  print, f = '(%"Volunteers 2/2021 . %i+/-%i")', total(raw), sqrt(total(Raw))

  print, raw - chncRaw
  
end
