pro cleancount

  readcol, 'HCdataSandbox.tsv', $
           tract, ad, tay, min, car, van, rv, tent, makeshift, $
           f = 'A,F,F,F,F,F,F,F,F', comment = '#'

  trct = strmid(tract, 0, 6)
  s = sort(trct)
  trct = trct[s]
  data = transpose([[ad[s]],[tay[s]],[min[s]],[car[s]],$
                    [van[s]],[rv[s]],[tent[s]],[makeshift[s]]])

  t = trct[uniq(trct)]
  d = fltarr(8, n_elements(t))
  
  for ii = 0, n_elements(t) - 1 do begin
     hit = where(trct eq t[ii], nhit)
     if nhit gt 1 then $
        d[*,ii] = total(data[*,hit],2) $
     else $
        d[*,ii] = data[*,hit]
  endfor

  close, 1
  openw, 1, 'retry2020.csv', width = 1024
  printf, 1, 'header'
  for ii = 0, n_elements(t) - 1 do $
     printf, 1, f = $
             '(%"time,email,%s,%i,%i,%i,%i,%i,%i,%i,%i,0,file,file")', $
             strmid(string(t[ii],f='(I0)'),0,4)+'.'+$
             strmid(string(t[ii],f='(I0)'), 1, /rev), $
                    d[0,ii], d[1,ii], d[2,ii], $
                    d[3,ii], d[4,ii], d[5,ii], $
                    d[6,ii], d[7,ii]
  close, 1
end
