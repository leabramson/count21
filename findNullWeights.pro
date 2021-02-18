pro findNullWeights, dataFile, target

  data = mrdfits(dataFile, 1)
  peaks = [1.38,1.68,1.32,1.45,1.64] ;; 2020
  bestWeights = [1,1,1, peaks, 1]
  input = transpose([[data.ADULT],[data.TAY],[data.MINOR],$
                     [data.CAR],[data.VAN],[data.RV],$
                     [data.TENT],[data.MAKESHIFT], [data.FAMILY]])

  ;; Best simple guess at this year's count
  rawTot = total(input, 2)
  newTot = total(rawTot * bestWeights)

  print, ''
  print, f = '(%"Current estimate . %i")', newTot
  print, f = '(%" > %3i cars")', rawTot[3]
  print, f = '(%" > %3i vans")', rawTot[4]
  print, f = '(%" > %3i RVs")', rawTot[5]
  print, f = '(%" > %3i tents")', rawTot[6]
  print, f = '(%" > %3i makeshifts")', rawTot[7]
  print, f = '(%"Past estimate .... %i\n")', target
  
  avgRuns = findgen((5.-0.)/0.05+1)*0.05
  nRuns = n_elements(avgRuns)

  ;; Find CVRTM weights needed to get this year to equal last year
  nullOuts = fltarr(5)
  for ii = 0, 4 do begin
     newWeights = bestWeights
     tmpTot = fltarr(nRuns)
     for jj = 0, nRuns - 1 do begin
        newWeights[3+ii] = avgRuns[jj]
        tmpTot[jj] = total(rawTot * newWeights)        
     endfor
;     if ii eq 0 then stop
     ind = value_locate(tmpTot, target)
     if ind ge 0 then $
        nullOuts[ii] = avgRuns[ind] $
     else $
        nullOuts[ii] = -1
  endfor
  print, f = '(%"Zero delta requires CVRTM mean occupancies of:")'
  tags = ['car', 'van', 'RV', 'tent', 'mkshft']
  for ii = 0, n_elements(peaks) - 1 do $
     print, f = '(%" > %4.2f ppl/%s (from %4.2f)")', $
            nullOuts[ii], tags[ii], peaks[ii]
  
end
