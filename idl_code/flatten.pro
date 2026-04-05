pro flatten, imdat
s=size(imdat) & nofs=s[2] & wofsg=s[1]
op=DBLARR(nofs) & s_sk=DBLARR(wofsg)  &  m_sk=DBLARR(wofsg)
 for j=0,wofsg-1 do begin
   op[*]=imdat[j,*]
   erov, op, st, mt
   s_sk[j]=st  &  m_sk[j]=mt
   imdat[j,*]=(op[*]-mt)/mt
 endfor
return
end