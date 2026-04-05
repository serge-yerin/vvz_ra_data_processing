PRO launch_sumthr, data, p
; clean rows and columns separately

; for ADR
; version 26 Feb 2015
print, 'in:'
help, data ; (freq, time)
x=rotate(data,4)
p=rotate(temporary(p),4)
print, 'proc:'
help, x ; (time, freq)

sx=size(x)
nt=sx[1] & nf=sx[2]
p1=p & p2=p

; columns  - VERTICAL
M=[2, 8, 16, 128, 256]
nm=n_elements(M)
thr=fltarr(nm)
for i=0,nm-1 do thr[i]=10./(1.5^(alog(M[i])/alog(2))) ; 1.2^ ...


for i=0,nt-1 do begin
  background, x[i,*], mx, sx, ny
;  erov, x[i,*], sx, mx
  y=reform(p2[i,*])
  xx=(reform(x[i,*])-mx)/sx
  for j=0, nm-1 do Sumthr, xx, M[j], thr[j], y
  p2[i,*]=y
endfor

; rows  - HORIZONTAL
M = [1, 2, 4, 8, 64] ; 4b
nm=n_elements(M)
thr=fltarr(nm)
for i=0,nm-1 do thr[i]=10./(1.5^(alog(M[i])/alog(2)))  ; 3a

for i=0,nf-1 do begin
  background, x[*,i], mx, sx, ny
;  erov, x[*,i], sx, mx
  y=reform(p1[*,i])
  xx=(reform(x[*,i])-mx)/sx
  for j=0, nm-1 do Sumthr, xx, M[j], thr[j], y
  p1[*,i]=y
endfor

p=p1*p2

p=rotate(temporary(p),4)

print,'out:' 
help, data
help, p
return
end