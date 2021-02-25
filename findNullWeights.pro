pro findNullWeights, data, target

;  data = mrdfits(resFile, 1)
  peaks = [1.38,1.68,1.32,1.45,1.64] ;; 2020
  bestWeights = [1,1,1, peaks, 1]
;  input = transpose([[data.ADULT],[data.TAY],[data.MINOR],$
;                     [data.CAR],[data.VAN],[data.RV],$
;                     [data.TENT],[data.MAKESHIFT], [data.FAMILY]])
  input = data.RAWCOUNTS
  
  ;; Best simple guess at this year's count
  if n_elements(data) gt 1 then $
     rawTot = total(input, 2) $
  else $
     rawtot = input
  newTot = total(rawTot * bestWeights)

  print, ' ** The following are generated from raw counts + mean weights **'
;  print, ' **                  No de-biasing necessary                   ** '
  print, ''
  print, f = '(%"Current estimate . %i")', newTot
  print, f = '(%" > %3i cars")', rawTot[3]
  print, f = '(%" > %3i vans")', rawTot[4]
  print, f = '(%" > %3i RVs")', rawTot[5]
  print, f = '(%" > %3i tents")', rawTot[6]
  print, f = '(%" > %3i makeshifts")', rawTot[7]
  print, f = '(%"Past estimate .... %i\n")', target

  ;; Find CVRTM weights needed to get this year to equal last year
  nullOuts = fltarr(5)
  nullFracs = fltarr(5)
  avgRuns = findgen((5.-0.)/0.025+1)*0.025
  nRuns = n_elements(avgRuns)
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

     ;; Find empty fractions
     nonTent = newTot - (rawTot * bestWeights)[3+ii] ;; the non-tent contribution
     noTentReach = target - nonTent                  ;; last year's total accounted for w/o this year's tents
     nTargTents  = noTentReach / bestweights[3+ii] ;; the number of tents needed at this year's weight to reach last year's total count
     nEmpty = rawTot[3+ii] - nTargTents            ;; number of this year's tents that have to be empty
     nullFracs[ii] = nEmpty / rawTot[3+ii]
  endfor
  tags = ['car', 'van', 'RV', 'tent', 'mkshft']
  print, f = '(%"Zero delta requires CVRTM mean occupancies of:")'
  for ii = 0, n_elements(peaks) - 1 do $
     print, f = '(%" > %4.2f ppl/%s (from %4.2f)")', $
            nullOuts[ii], tags[ii], peaks[ii]
  print, f = '(%"\nZero delta requires CVRTM unoccupied fractions of:")'
  for ii = 0, n_elements(peaks) - 1 do $
     print, f = '(%" > %4.2f (at %4.2f ppl/inhabited %s)")', $
            nullFracs[ii], peaks[ii], tags[ii]
  
end
