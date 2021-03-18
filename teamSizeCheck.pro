pro teamSizeCheck

  readcol, 'teamSizes.csv', id, nc, f = 'A,F', $
           delim = ','
  nc = strcompress(nc, /rem)
  
  data = mrdfits('countHollywood2021w191902.fits',1)
  data = data[where(~data.FLAG)]
  nlines = n_elements(data)

  teamSize = intarr(nlines)
  for ii = 0, nlines - 1 do begin
     hit = where(id eq strcompress(data[ii].COUNTER, /rem), nhit)
     if nhit gt 0 then teamSize[ii] = nc[hit[0]]
  endfor

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
  ncount = ncount[forComp]
  
  ;; find total number of pairs
  ;; one pair for all tracts (t1-t0), then 2 additional pairs for the
  ;; triplets (t2-t0, t2-t1)
  ntrips = total(ncount eq 3)
  np = nutract + 2 * ntrips
  
  ;; make it individuals
  input = transpose([[data.ADULT+data.TAY+data.MINOR],$
                     [data.CAR],[data.VAN],[data.RV],$
                     [data.TENT],[data.MAKESHIFT]]);, [data.FAMILY]])
  ncat = n_elements(input[*,0])

  ;; outputs
  team1  = fltarr(ncat, nutract,2)
  team2  = fltarr(ncat, nutract,2)
  team3  = fltarr(ncat, nutract,2)
  rcomp  = fltarr(ncat, np)
  rerr   = fltarr(ncat, np)
  rtcomp = fltarr(np)
  rterr  = fltarr(np)

  meanCounts = fltarr(ncat, np)
  meanTots = fltarr(np)
  ptracts = fltarr(np)

  smallTeam = fltarr(np)
  bigTeam   = fltarr(np)
  eqTeam    = bytarr(np)

  nteams = intarr(np)
  
  jj = 0
  ii = 0
  while ii lt nutract do begin     

     hit = where(tract eq utract[ii], nhit)

     if nhit eq 2 then ptracts[jj] = utract[ii]
     if nhit eq 3 then ptracts[jj:jj+nhit-1] = utract[ii]
     
     master = hit[0]
     slave  = hit[1]
     if nhit eq 3 then $
        trip = hit[2]
     
     team1[*,ii,0] = input[*,master]
     team2[*,ii,0] = input[*,slave]
     team1[*,ii,1] = teamSize[master]
     team2[*,ii,1] = teamSize[slave]
     rcomp[*,jj] = input[*,slave] - input[*,master]
     rerr[*,jj]  = sqrt(input[*,master] + input[*,slave])
     tmp0        = total(input[*,master])
     tmp1        = total(input[*,slave])
     rtcomp[jj]  = (tmp1-tmp0)
     rterr[jj]   = sqrt(tmp1+tmp0)

     meanCounts[*,jj] = (input[*,master] + input[*,slave])/2.
     meanTots[jj]   = (tmp0+tmp1)/2.

     bt = max(teamSize[hit[0:1]], bti)
     st = min(teamSize[hit[0:1]], sti)
     bti = bti[0]
     sti = sti[0]
     if bti eq sti OR teamsize(hit[0]) eq -1 OR teamsize(hit[1]) eq -1 then $
        eqteam[jj] = 1 

     if bti ne sti then $
        print, f = '(%"Tract: %7.2f; Large Tream: %s; Small Team: %s")', $
               data[hit[bti]].TRACT, data[hit[bti]].COUNTER, data[hit[sti]].COUNTER
     
     bigTeam[jj] = ([tmp0,tmp1])[bti]
     if bti eq sti then sti++
     smallTeam[jj] = ([tmp0,tmp1])[sti]
     
     if n_elements(hit) eq 3 then begin

        tmp2 = total(input[*,trip])

        team3[*,ii,0]   = input[*,trip]
        team3[*,ii,1]   = teamSize[trip]
        rcomp[*,jj+1] = team3[*,jj] - team1[*,jj]
        rcomp[*,jj+2] = team3[*,jj] - team2[*,jj]
        rerr[*,jj+1]  = sqrt(team3[*,jj] + team1[*,jj])
        rerr[*,jj+2]  = sqrt(team3[*,jj] + team2[*,jj])
        rtcomp[jj+1]  = tmp2 - tmp0
        rtcomp[jj+2]  = tmp2 - tmp1
        rterr[jj+1]   = sqrt(tmp2+tmp0)
        rterr[jj+2]   = sqrt(tmp2+tmp1)

        meanCounts[*,jj+1] = (input[*,master]+input[*,trip])/2.
        meanCounts[*,jj+2] = (input[*,slave]+input[*,trip])/2.
        meanTots[jj+1]   = (tmp0+tmp2)/2.
        meanTots[jj+2]   = (tmp1+tmp2)/2.
        
        bt = max(teamSize[hit[[0,2]]], bti)
        st = min(teamSize[hit[[0,2]]], sti)

        bti = bti[0]
        sti = sti[0]
        if bti eq sti OR teamsize(hit[0]) eq -1 OR teamsize(hit[2]) eq -1 then $
           eqteam[jj+1] = 1 

        if bti eq sti then sti++
        bigTeam[jj+1] = ([tmp0,tmp2])[bti]
        smallTeam[jj+1] = ([tmp0,tmp2])[sti]

        bt = max(teamSize[hit[[1,2]]], bti)
        st = min(teamSize[hit[[1,2]]], sti)

        bti = bti[0]
        sti = sti[0]
        if bti eq sti OR teamsize(hit[1]) eq -1 OR teamsize(hit[2]) eq -1 then $
           eqteam[jj+2] = 1 

        if bti eq sti then sti++
        bigTeam[jj+2] = ([tmp1,tmp2])[bti]
        smallTeam[jj+2] = ([tmp1,tmp2])[sti]
        
     endif
     
     ii++
     if nhit eq 3 then jj+=nhit else jj++

  endwhile

  use = where(~eqteam)

  !p.multi = [0,2,0]
  
  plot, smallTeam, bigTeam, psym = 1, xr = [0,120], yr = [0,120], /iso, $
        xtitle = 'smaller team on tract', $
        ytitle = 'larger team on tract', /nodat
  one_one, col = 'ffa500'x
;  oploterror, smallTeam, bigTeam, sqrt(smallTeam), sqrt(bigTeam), $
;              errthick = 1, psym = 1, /nohat
  oploterror, smallTeam[use], bigTeam[use], sqrt(smallTeam[use]), sqrt(bigTeam[use]), $
              errthick = 1, psym = 1, errcol = long('0000ff'x), /nohat
  oplot, smallTeam[use], bigTeam[use], psym = 1, col = 255

  use = where(~eqteam and ptracts ne 1901.00)
  
  td = sqrt((bigTeam - smallTeam)^2/(bigTeam+smallTeam))
  d  = sqrt((bigTeam - smallTeam)[where(eqteam)]^2/(bigTeam+smallTeam)[where(eqTeam)])
  d2 = td[use]

  td = td[sort(td)]
  d = d[sort(d)]
  d2 = d2[sort(d2)]
  pctles = [0.05,0.16,0.50,0.84,0.95]

  tdp = td[ceil(pctles * n_elements(td))-1]
  dp = d[ceil(pctles * n_elements(d))-1]
  d2p = d2[ceil(pctles * n_elements(d2))-1]

  em = (dp[3]-dp[1])/2/sqrt(n_elements(d))
  e2m = (d2p[3]-d2p[1])/2/sqrt(n_elements(d2))
  
  plot, td, findgen(n_elements(td))/(n_elements(td)-1), $
        xr = [0,3], $
        xtitle = 'sqrt[(!18n!D1!X!N - !18n!X!D2!N)!E2!N!18/!X(!18n!X!D1!N+!18n!X!D2!N)]', $
        ytitle = '!18f!X!Dpairs!N', /nodat
  polyfill, dp[2] + em * [-1,-1,1,1], !Y.CRANGE[[0,1,1,0]], $
            /line_fill, orien = 45, col = 'cccccc'x, spacing = 0.05
  polyfill, d2p[2] + e2m * [-1,-1,1,1], !Y.CRANGE[[0,1,1,0]], $
            /line_fill, orien = 45, col = '00a5ff'x, spacing = 0.05
  oplot, td, findgen(n_elements(td))/(n_elements(td)-1), col = '777777'x, psym = 10
  oplot, d, findgen(n_elements(d))/(n_elements(d)-1), psym = 10
  oplot, d2, findgen(n_elements(d2))/(n_elements(d2)-1), col = 255, psym = 10
  oplot, dp[2]*[1,1], !Y.CRANGE, linesty = 2
  oplot, d2p[2]*[1,1], !Y.CRANGE, linesty = 2, col = 255
  
  print, (d2p-dp)[2]/sqrt(em^2 + e2m^2)
  
  stop
  
end
