pro SumThr, x, M, thr, y
; x - input 1D array
; M - size of sliding window
; thr - specified threshold for average of sequence size M
; y - previous mask of bad pixels (1=good, 0=bad)

w=where(y eq 1)
if w(0) eq -1 then return

xx=x(w)
nxx=n_elements(xx)
if nxx lt M then return

xxa=fltarr(nxx-M+1)
for i=0,M-1 do xxa=xxa+xx(i:*)

p=bytarr(nxx)+1
wa=where(xxa gt thr*M)
if wa(0) ne -1 then for i=0,M-1 do p(wa+i)=0

y(w)=p
return
end