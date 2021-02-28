pro compBID

  chncTracts = [1901.00, 1902.02, 1902.01, 1907.00, 1908.01, 1908.02, 1909.02, 1910.00]
;  chncTracts = [1907.00, 1908.02]
  ntracts = n_elements(chncTracts)

  data = mrdfits('countHollywoodResults2021.fits', 1)
  raw = fltarr(9,ntracts)
  for ii = 0, ntracts - 1 do begin
     hit = where(data.TRACT eq chncTracts[ii])
     raw[*,ii] = data[hit].RAWCOUNTS
  endfor

  print, total(raw[3,*] + raw[4,*] + raw[5,*]), total(raw[0,*]+raw[1,*]), total(raw[6,*]+raw[7,*])
  
  
;  print, f = '(%"Abramson 7/2020 ... %i+/-%i")', total(chncRaw), sqrt(total(chncRaw))
;  print, f = '(%"Volunteers 2/2021 . %i+/-%i")', total(raw), sqrt(total(Raw))

;  print, raw - chncRaw
  
end
