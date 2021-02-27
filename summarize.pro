pro summarize, resfile

  data     = mrdfits(resfile, 1)
  eHoDat   = data[where(data.EASTFLAG, compl = Hwood, nEho)]
  hWoodDat = data[Hwood]
  nHwood   = n_elements(hWoodDat)

  ;; Print summary stats
  
  for ii = 0, 2 do begin
     case ii of
        0: begin $
           d = data
           region = 'GREATER HOLLYWOOD'
        end
        1: begin
           d = eHoDat
           region = 'EAST HOLLYWOOD'
        end
        2: begin
           d = hWoodDat
           region = 'HOLLYWOOD'
        end
     endcase

;     ;; debias
;     d = debias(d)
     
     means = mean(d.COUNTS, dim = 2)
     errs  = stddev(d.COUNTS, dim = 2)

     tot    = total(means)
     totErr = 2*sqrt(total(errs^2))

     s = size(means, /dim)
     if n_elements(s) gt 1 then begin
        totStruct     = total(means, 2)
        totStructErrs = 2*sqrt(total(errs^2, 2))
     endif else begin
        totStruct     = means
        totStructErrs = 2*errs
     endelse
     

;     print, " THE FOLLOWING ESTIMATES HAVE NOT BEEN DE-BIASED! BEWARE OF MEDIANS! "
;     print, ""
     
     tstring = "*** SUMMARY FOR "+region+" (95% conf.) ***"
     print, ''
     print, tstring
     print, replicate('-', strlen(tstring)/2)
     print, f = '(%"Total: %i+/-%i")', round(tot), round(totErr)
     for jj = 0, n_elements(d[0].TAGS) - 1 do $
        print, f = '(%" > %s: %i+/-%i (%4.1f\%)")', d[0].TAGS[jj], $
               round(totStruct[jj]), round(totStructErrs[jj]), $
               round(totStruct[jj]/tot*100)
     
  endfor

end
