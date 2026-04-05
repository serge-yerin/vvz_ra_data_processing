
PRO adr_cleaning, data, mask, pat=pat, nar=nar, wid=wid
; v 2 Mar 2015
err=0

s=size(data)
ndim=s[0]
if ndim eq 2 then begin
 wofsg=s[1]
 nofs=s[2]
endif else stop, 'Data should be 2D !'

; --------------- NaN and Inf removal --------------
mask=bytarr(wofsg,nofs)+1b
w=where(finite(data) ne 1)
if w[0] gt -1 then begin
  mask[w]=0
  data[w]=median(data)
endif
; --------------------------------------------------

if(min(data)eq max(data)) then err=1

flatarr=data

flatten, flatarr

med=0. 

if (keyword_set(pat)) then begin 
launch_sumthr, flatarr, mask
endif

if (keyword_set(wid) or keyword_set(nar)) then patrol, data, flatarr, med, mask, wid=wid, nar=nar

 data = flatarr
 data=data*mask

if (err eq 1) then data[*,*]=0.

return
END
