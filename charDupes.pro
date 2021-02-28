pro charDupes, dataFits

  data = mrdfits(dataFits, 1)
  data = data[where(~data.FLAG)]
  nlines = n_elements(data)

  tract = data.TRACT
  s = sort(tract)
  stract = tract[s]
  utract = tract[uniq(stract)]
  nutract = n_elements(utract)

  ;; find ncounters
  ncount = fltarr(nutract)
  for ii = 0, nutract - 1 do $
     ncount[ii] = total(tract eq utract[ii])

  forComp = where(ncount gt 1, nutract)
  utract = utract[forComp]

  ;; make it individuals
  
  input = transpose([[data.ADULT+data.TAY+data.MINOR],$
                     [data.CAR],[data.VAN],[data.RV],$
                     [data.TENT],[data.MAKESHIFT]]);, [data.FAMILY]])
  ncat = n_elements(input[*,0])

  comp = fltarr(ncat, nutract)
  tcomp = fltarr(nutract)
  bb = fltarr(nutract,2)
  bbb = fltarr(total(ncount eq 3))
  jj = 0
  for ii = 0, nutract - 1 do begin
     hit = where(tract eq utract[ii])
     comp[*,ii] = (input[*,hit[1]] - input[*,hit[0]]) $
                  / sqrt(input[*,hit[0]] + input[*,hit[1]])
     tmp0 = total(input[*,hit[1]])
     tmp1 = total(input[*,hit[0]])
     bb[ii,*] = [tmp0,tmp1]
     tcomp[ii] = (tmp1-tmp0) / sqrt(tmp0 + tmp1)

     if n_elements(hit) eq 3 then begin
        tmp2 = total(input[*,hit[2]])
        bbb[jj] = tmp2
        jj++
     endif
  endfor

  stop

;  qui = where(data.TRACT ne 1901.00)
  
  plot, findgen(ncat), abs(mean(comp, dim = 2, /nan)), yr = [0,2], xr = [-1,6], $
        xtickname = [' ', 'P', 'C', 'V', 'R', 'T', 'M', ' '], $
        xticks = 7, xtickint = 1, $
        ytitle = '|n!D1!N-n!D0!N|/sqrt(n!D1!N+n!D0!N)', xminor = 1
  oplot, findgen(ncat), stddev(comp, dim = 2, /nan), linesty = 2
  oplot, findgen(ncat), stddev(comp, dim = 2, /nan)/sqrt(nutract), col = 255  

  stop
  
  pctles = [0.05,0.16,0.50,0.84,0.95]
  stats = fltarr(ncat,n_elements(pctles))
  for ii = 0, ncat - 1 do begin
     x = comp[ii,where(finite(comp[ii,*]))]
     x = x[sort(x)]
     stats[ii,*] = x[ceil(pctles*nutract)-1]
  endfor

  trips = where(ncount eq 3, ntrips)
  plot, bb[*,0], bb[*,1], psym = 1, xr = [0,120], yr = [0,120], /iso, $
        xtitle = 'team 1', ytitle = 'team 2'
  one_one
  oploterror, bb[*,0], bb[*,1], sqrt(bb[*,0]), sqrt(bb[*,1]), psym = 1, /nohat
  oplot, bb[trips,0], bbb, psym = 1, col = 255
  oploterror, bb[trips,0], bbb, sqrt(bb[trips,0]), sqrt(bbb), psym = 1, /nohat, $
              errcol = 'ff5500'x

  d = findgen(101)-50
  del = bb[*,1] - bb[*,0] ;; duplicates
  del = [del, bbb - bb[trips,0]];, bbb-bb[trips,1]] ;; with triplicates to 0th and 1st measurement
  ebar = sqrt(bb[*,0] + bb[*,1])
  ebar = [ebar, sqrt(bbb + bb[trips,0])];, sqrt(bbb + bb[trips,1])]

  pdf = fltarr(n_elements(d))
;  npdf = fltarr(n_elements(d))
  for ii = 0, n_elements(ebar) - 1 do $
     pdf += (1./sqrt(2*!pi*ebar[ii]^2) * exp(-0.5 * (d-del[ii])^2/ebar[ii]^2))>0
;  for ii = 0, n_elements(ebar) - 1 do $
;     npdf += (1./sqrt(2*!pi*del[ii]^2) * exp(-0.5 * (d-del[ii])^2/del[ii]^2))>0
  
;  plot, d/mean(ebar), pdf/max(pdf), $
;        xtitle = greek('Delta')+'!18/<!X'+greek('sigma')+'!18>!X', $
;        ytitle = 'prob', xr = [-3,3]
;  null = exp(-0.5 * (d*0.1)^2)
;  oplot, d*(0.1), null/max(null), col = 255
  
;  stop
  
end
