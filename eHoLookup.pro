function eHoLookup, tract

  ; Assigning all of 1905.10 to Hollywood
  
  tracts = [1916.10, $
            1916.20, $
            1905.20, $
            1911.10, $
            1911.20, $
            1925.10, $
            1925.20, $
            1926.20, $
            1926.10, $
            1915.00, $
            1912.04, $
            1912.03, $
            1912.01, $
            1913.02, $
            1913.01, $
            1914.10, $
            1914.20, $
            1927.00]

  result = total(tracts eq float(tract))
  
  RETURN, result
end

