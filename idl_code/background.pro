; **********************
; *                    *
; *   BACKGROUND.PRO   *
; *                    *
; **********************

  pro BACKGROUND, tab, fon,sigma,ny, POSITIVE=POSITIVE

; Background and sigma of a 1-D distribution of intensities.

; (INPUT)
; tab = 1-D array of intensities
; (OUTPUT)
; fon,sigma = background and fluctuations (1 sigma level)
; ny = number of values upon which fon,sigma are computed
; POSITIVE = retains only tab values > 0

  if keyword_set(POSITIVE) then begin
    test=where(tab gt 0.)
    if test(0) ne -1 then tab2=tab(test) else begin
	fon=0 & sigma=-1 & ny=0 & goto,fin
    endelse
  endif else tab2=tab
  if n_elements(tab2) gt 1 then begin
	fon=mean(tab2,/double) & sigma=stddev(tab2,/double)
  endif else begin        fon=tab2(0) & sigma=0 & ny=1 & goto,fin  endelse
  ny=1

encore:
  test=where(abs(tab2-fon) lt 3.*sigma)  ; it was 2.5*sigma before
  if n_elements(test) eq 1 then goto,fin
  ny=n_elements(test)
  moy=mean(tab2(test),/double)
  sigma=stddev(tab2(test),/double)
  if moy eq fon then goto,fin
  fon=moy
  goto,encore

fin:
return
end
