function gen2021t, maxOcc

  if n_elements(maxOcc) eq 0 then $
     maxOcc = 4.

  cts = [1,1,1,2,1,4,2,1,0,1,1,2,1,1, $
         1,1,3,2,2,0,1,2,1,4,2,1,1,1, $
         1,2,2,0,2,1,1,1,1,1]

  nnull = 9.
  nct = n_elements(cts)
  ntot = nct + nnull

  ;; Proceed in 2 ways
  ;; 1. use the known quantities
  ;; 2. model the distribution of the unknown tents

  ;; 1 -- use known stuff
  t = mean(cts)
  terr = stddev(cts)/sqrt(nct)

  ;; 2 -- model unknown stuff
  niter = 1000
  ms = fltarr(niter)
  for ii = 0, niter - 1 do begin
     occ = randomu(seed,9)*maxOcc
     ms[ii] = mean([cts+randomn(seed,nct)*cts,occ])
  endfor
  es = stddev(ms)
  
  RETURN, [t,terr,mean(ms),es]
end
