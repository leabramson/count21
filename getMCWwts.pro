function getMCWwts

  ;; totUnshelt * [p,c,v,r,t,m]
;  midCityPop = 224 * [0.223,0.098,0.098,0.165,0.277,0.139]
;  midCityCts = [224. * 0.223, 15., 13., 25., 42., 19.]
;  midCityWts = midCityPop / midCityCts

  bevGrovePop = 147 * [0.517,0.055,0.061,0.02,0.15,0.197]
  bevGroveCts = [147 * 0.517, 5., 5., 2., 15., 17.]
  bevGroveWts = bevGrovePop / bevGroveCts
  
  fairfaxPop = 87 * [0.563,0.023,0.126,0.184,0.069,0.035]
  fairfaxCts = [87 * 0.563, 1., 6., 11., 4., 2.]
  fairfaxWts = fairfaxPop / fairfaxCts

  midWilshirePop = 100 * [0.38,0.12,0.18,0.09,0.10,0.13]
  midWilshireCts = [100 * 0.38, 8., 10., 6., 7., 8.]
  midWilshireWts = midWilshirePop / midWilshireCts
 
;  print, bevGroveWts;midCityWts
;  print, fairfaxWts
;  print, midWilshireWts

  ;; Areal fraction weights from BK.
  wts = 0.3 * fairfaxWts + $
        0.3 * midWilshireWts + $
        0.4 * bevGroveWts

  RETURN, wts[1:*] ;; kill 1 for ppl

 end
