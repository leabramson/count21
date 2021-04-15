pro mergeSplits, dataFits, outname

  data = mrdfits(dataFits, 1)
;  data = data[where(~data.FLAG)]
  nlines = n_elements(data)

  counter = data.COUNTER
  ctr = counter[sort(counter)]
  ctr = ctr[uniq(ctr)]
  nctr = n_elements(ctr)

  input = transpose([[data.ADULT],[data.TAY], [data.MINOR], $
                     [data.CAR],[data.VAN],[data.RV],[data.TENT],[data.MAKESHIFT], $
                     [data.FAMILY]])
  output = fltarr(9)
  outTract = []
  outCount = []
  for ii = 0, nctr - 1 do begin
     hit = where(counter eq ctr[ii], nhit)
     if nhit gt 1 then begin
        t = data[hit].TRACT
        prefix = strmid(t, 0, 7)
        utract = prefix[sort(prefix)]
        utract = utract[UNIQ(utract)]
        nutract = n_elements(utract)
        outTract = [outTract, uTract]
        outCount = [outCount, replicate(ctr[ii], nutract)]
        if nutract ge 1 then begin
           for jj = 0, nutract - 1 do begin
              inds = where(prefix eq utract[jj], nhit)
              if nhit gt 1 then $
                 output = [[output], [total(input[*,hit[inds]],2)]] $
              else $
                 output = [[output], [input[*,hit[inds]]]]
           endfor
        endif            
     endif else begin
        output = [[output], [input[*,hit]]]
        outTract = [outTract, data[hit].TRACT]
        outCount = [outCount, data[hit].COUNTER]
     endelse
  endfor
  output = output[*,1:*] ;; trim first line of 0s
  nout = n_elements(outTract)
  
  savedata = data[0]
  savedata = replicate(saveData, nout)

  s = sort(outTract)
  output = output[*,s]
  outTract = strcompress(outTract[s],/rem)
  outCount = outCount[s]

  outFlag = bytarr(nout)
  outTIme = strarr(nout)
  prefix = strmid(data.TRACT, 0, 7)
  for ii = 0, nout - 1 do begin
     hit = where(prefix eq outTract[ii] and counter eq outCount[ii])
     outTime[ii] = data[hit[0]].TIMESTAMP
     if total(data[hit].FLAG) then $
        outFlag[ii] = 1
  endfor
  
  for ii = 0, nout - 1 do begin
     savedata[ii].COUNTER = outCount[ii]
     savedata[ii].TRACT   = outTract[ii]
     savedata[ii].FLAG    = outFlag[ii]
     savedata[ii].ADULT   = output[0,ii]
     savedata[ii].TAY     = output[1,ii]
     savedata[ii].MINOR   = output[2,ii]
     savedata[ii].CAR     = output[3,ii]
     savedata[ii].VAN     = output[4,ii]
     savedata[ii].RV      = output[5,ii]
     savedata[ii].TENT    = output[6,ii]
     savedata[ii].MAKESHIFT = output[7,ii]
     savedata[ii].FAMILY    = output[8,ii]
     savedata[ii].TIMESTAMP = outTime[ii]
  endfor

  mwrfits, savedata, outname, /create
  
end
