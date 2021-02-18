pro wtErrsCalc

  car     = [1.84,1.48,1.25,1.38,1.37,1.26,1.20,1.21]
  ecar    = [0.35,0.12,0.17,0.11,0.13,0.14,0.12,0.17]
  
  van     = [1.38,1.76,1.53,1.68,1.66,1.59,1.67,1.80]
  evan    = [0.68,0.29,0.41,0.22,0.23,0.12,0.32,0.33]

  rv      = [2.12,1.92,2.30,1.32,1.51,1.51,1.96,1.56]
  erv     = [0.13,0.13,0.52,0.15,0.20,0.10,0.57,0.17]

  tent    = [1.66,1.57,1.50,1.45,1.47,1.51,1.39,1.44]
  etent   = [0.23,0.11,0.11,0.06,0.07,0.09,0.12,0.23]

  mkshift  = [2.19,1.38,2.18,1.64,1.65,1.46,1.78,1.26]
  emkshift = [0.38,0.12,0.75,0.16,0.17,0.11,0.35,0.10]

  wts = transpose([[car],[van],[rv],[tent],[mkshift]])
  ewts = transpose([[ecar],[evan],[erv],[etent],[emkshift]])
  
  ncar     = [282,402,225,521,305,529,349,324]
  nvan     = [128,346,181,694,379,638,287,331]
  nrv      = [540,947,269,735,371,1439,427,475]
  ntent    = [216,449,153,2126,396,417,62,201]
  nmkshift = [312,434,153,1235,383,562,219,319]

  types = transpose([[ncar],[nvan],[nrv],[ntent],[nmkshift]])
  
  print, ecar/car
  print, evan/van
  print, erv/rv
  print, etent/tent
  print, emkshift/mkshift

  terr = total(ewts/wts * types, 1) / total(types, 1)
  for ii = 0, n_elements(ncar) - 1 do $
     print, f = '(%"SPA%i base weighting error: %5.2f")', $
            ii+1, terr[ii]
  print, f = '(%"Total base weighting error: %5.2f")', $
         total(total(types, 1) / total(types) * terr) ;; each spa's structure/vehicle contribution weighted by its overall error
  
end
