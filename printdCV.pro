pro printdCV, lastYear, $
              region = region

  spawn, 'ls hwoodCVRTMtests/*.fits > test.list'
  readcol, 'test.list', files, f = 'A'
  nfiles = n_elements(files)

  if strupcase(region) eq 'HWOOD' then begin
     output = 'hwoodFinal.eps'
     lycol = 'ff5500'x
     stt = 0.75
  endif else if strupcase(region) eq 'EHO' then begin
     output = 'ehoFinal.eps'
     lycol = long('0055ff'x)
     stt = 0.65
  endif
  
  fn = []
  for ii = 0, nfiles - 1 do $
     fn = [fn, (strsplit(files[ii], '/',/extr))[1]]
  
  wtnames = strmid(fn, 0 ,4)

  out = fltarr(5,nfiles)
  means = fltarr(nfiles)
  pctles = [0.05,0.16,0.5,0.84,0.95]
  dcv = fltarr(nfiles)
  ps = fltarr(nfiles,2)
  for ii = 0, nfiles - 1 do begin

     d = mrdfits(files[ii], 1)
;     print, d[0].WTS
     
     if strupcase(region) eq 'HWOOD' then begin        
        d = d[where(d.EASTFLAG eq 0)]
        lastYear = 1058.
        if ii eq 0 then begin
           ;; boost the C & V numbers to 2020
           dc = 55 - total(d.RAWCOUNTS[3])
           dv = 48 - total(d.RAWCOUNTS[4])
        endif
     endif else if strupcase(region) eq 'EHO' then begin
        d = d[where(d.EASTFLAG)]
        lastYear = 656.
        if ii eq 0 then begin
           ;; boost the C & V numbers to 2020
           dc = 29 - total(d.RAWCOUNTS[3])
           dv = 58 - total(d.RAWCOUNTS[4])
        endif
     endif

     dcv[ii] = total([dc,dv] * d[0].WTS[[3,4]])
     
     cts = total(total(d.COUNTS,3),1)
     
     p1 = getCountProb(cts, lastYear)
     p2 = getCountProb(cts + dcv[ii],lastYear)

     ps[ii,*] = [p1,p2]
     
     cts = cts[sort(cts)]
     out[*,ii] = cts[ceil(n_elements(cts) * pctles) - 1]
     print, wtnames[ii], [out[*,ii], p1]
     print, wtnames[ii], [out[*,ii] + dcv[ii], p2]
     means[ii] = total(total(d.RAWCOUNTS, 2) * d[0].WTS)
  endfor
  null = where(wtnames eq '2020')
;  nullOut = out[*,null]
;  nullName = 'SPA4 2020'
  case strupcase(region) of
     'HWOOD': title = 'Hollywood CoC'
     'EHO': title = 'East Hollywood CoC'
  endcase

end
