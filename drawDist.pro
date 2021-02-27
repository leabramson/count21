function drawDist, x, pdf, n, $
                   normalize = normalize

  dx = x[1] - x[0]
  cdf = total(pdf * dx, /cum)
  
  if keyword_set(normalize) then cdf /= max(cdf)

  out = x[value_locate(CDF, randomu(seed, n))]

  RETURN, out
end
