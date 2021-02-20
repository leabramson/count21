function idProTracts, tracts

  ptlist = [1903.01, $
            1905.10, $
            1908.02, $
            1909.01, $
            1910.00, $
            1911.20, $
            1916.10, $
            1916.20, $
            1927.00]
  
  output = bytarr(n_elements(tracts))
  for ii = 0, n_elements(tracts) - 1 do $
     if total(ptlist eq tracts[ii]) gt 0 then $
        output[ii] = 1

  return, output

end
