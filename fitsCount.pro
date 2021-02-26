pro fitsCount, infile, outfile, $
               emptyFrac = emptyFrac

  if NOT keyword_set(emptyFrac) then $
     emptyFrac = replicate(0.,5)
  
  readcol, infile, $
           timestamp, counter, tract, $
           adult, tay, minor, $
           car, van, RV, $
           tent, makeshift, family, $
           Verify, Flag, $
           f = 'A,A,A,F,F,F,F,F,F,F,F,F,A,A', $
           delim = ',', /preserve_null, skipline = 1
  nLines = n_elements(tract)
  for ii = 0, nLines - 1 do $
     tract[ii] = strmid(tract[ii],0,4)+'.'+strmid(tract[ii], 1, /rev)
  s = sort(tract)

  tflag = bytarr(nLines)
  tflag[where(strcompress(flag) eq 'T')] = 1
  flag = tflag  
  
  zcheck = where(tract eq '0', nbad)
  if nbad gt 0 then begin
     print, f = '(%"WARNING! There are %i lines with no tract entered!")', $
            nbad
     for ii = 0, nbad - 1 do $
        print, f = '(%"See line: %i")', zcheck[ii]+1
  endif
  print, ''
  
  fcheck = where(family gt 0, ndfam)
  if fcheck gt 0 then begin
     print, f = '(%"WARNING! There are %i lines with at least 1 family!")', $
            nfam
     for ii = 0, nfam - 1 do $
        print, f = '(%"See line: %i")', zcheck[ii]+1
  endif
  savedata = {TRACT    : '0', $
              ADULT    : 0., $
              TAY      : 0., $
              MINOR    : 0., $
              CAR      : 0., $
              VAN      : 0., $
              RV       : 0., $
              TENT     : 0., $
              MAKESHIFT: 0., $
              FAMILY   : 0., $
              COUNTER  : '0', $
;              TALLYFILE: '0', $
              FLAG     : 0, $
              EAST     : 0}
  savedata = replicate(savedata, nLines)

  for ii = 0, nLines - 1 do begin
     savedata[ii].TRACT     = tract[ii]
     savedata[ii].ADULT     = adult[ii]    
     savedata[ii].TAY       = tay[ii]    
     savedata[ii].MINOR     = minor[ii]         
     savedata[ii].CAR       = car[ii]       * (1-emptyFrac[0])         
     savedata[ii].VAN       = van[ii]       * (1-emptyFrac[1])         
     savedata[ii].RV        = rv[ii]        * (1-emptyFrac[2])                   
     savedata[ii].TENT      = tent[ii]      * (1-emptyFrac[3])                   
     savedata[ii].MAKESHIFT = makeshift[ii] * (1-emptyFrac[4])            
     savedata[ii].FAMILY    = family[ii]        
     savedata[ii].COUNTER   = counter[ii]       
     savedata[ii].FLAG      = flag[ii]
     savedata[ii].EAST      = eHoLookup(tract[ii])
  endfor
  ;; cull quality flagged entries
  savedata = savedata[s] ;; sort by tract
  bad = where(flag, nbad)

  if nbad gt 0 then begin
     for ii = 0, nbad-1 do $
        print, f = '(%"Bad data detected at tract %7.2f by counter %s")', $
               savedata[bad[ii]].TRACT, savedata[bad[ii]].COUNTER
  endif
  
  mwrfits, savedata, outfile, /create
  
end
