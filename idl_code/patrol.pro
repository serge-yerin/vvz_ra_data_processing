PRO patrol, imdat, flatdat, med, p, wid=wid, nar=nar

; for ADR data
; 2 Mar 2015

s=size(imdat)
wofsg=s[1] & nofs=s[2]

l2=4.

if (keyword_set(nar)) then begin
    
    s_sk=DBLARR(wofsg)  &  m_sk=DBLARR(wofsg)
     for j=0,wofsg-1 do begin
       testarr=where(p[j,*] eq 1)
       if(testarr[0] gt (-1)) then begin
       op=imdat[j,where(p[j,*] eq 1)]
       erov, op, st, mt
       s_sk[j]=st  &  m_sk[j]=mt
       endif
     endfor
      for j=0,wofsg-1 do if(m_sk[j] eq 0) then m_sk[j]=median(m_sk)
    op=s_sk/m_sk
    erov, op, st, mt
    testarr=where(abs(op-mt) gt (st*l2))
    if(testarr[0] gt (-1)) then begin
;      imdat[where(ispr gt 0.),*]=med ;mean(imdat[where(cle gt 0.)])
      p[testarr,*]=0
    endif
endif

if (keyword_set(wid)) then begin
    op=fltarr(nofs)
    for j=0, nofs-1 do begin
    ;op=total(imdat,1)/wofsg
    w=where(p[*,j] eq 1)
    if w[0] gt (-1) then op[j]=total(imdat[w,j])/total(p[*,j])
    endfor
    
    erov, op, st, mt
    testarr=where(abs(op-mt) gt (st*l2))
    if(testarr[0] gt (-1)) then begin 
      ;imdat[*,where(ispr gt 0.)]=med ;mean(imdat[where(cle gt 0.)])
      p[*,testarr]=0
    endif
endif

return
end