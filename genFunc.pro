function genFunc, x, peak, sigma, cutoff

  func = 1./sqrt(sigma^2 * 2 * !pi) $
         * exp(-0.5 * (x-peak)^2/sigma^2)

  if N_ELEMENTS(cutoff) eq 0 then $
     cutoff = min(x)
  
  func[where(x lt cutoff)] = 0

  func /= total(func * (x[1]-x[0]))
  
  RETURN, func
end
