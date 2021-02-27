function findSig, x, peak, sigi, sigf, $
                  nsteps = nsteps, $
                  ntile = ntile, $
                  cutoff = cutoff, $
                  target = target

  if NOT keyword_set(ntile) then ntile = 0.95
  
  dsig = (sigf - sigi)/(nsteps-1)
  sigran = findgen(nsteps) * dsig + sigi
  nout = fltarr(nsteps)
  for ii = 0, nsteps - 1 do begin
     foo = genFunc(x, peak, sigran[ii], 1)
     cdf = total(foo, /cum)
     nout[ii] = x[value_locate(cdf, ntile * cdf[-1])]
  endfor

  sigma = (sigran[value_locate(nout, target)] > min(sigran)) < max(sigran)
  
  RETURN, sigma
end

