PRO erov, op, st, mt
count=n_elements(op) & m_in=lindgen(count) & m=op & sr=mean(m)
repeat begin
  srkv=stddev(m) & count_pr=count & sr_pr=sr & m=m*(abs(m-sr) le srkv*3.)
  m_in=where(m,count) 
  if m_in[0] ne -1 then begin
  m=m[m_in] & sr=mean(m) & ster=stddev(m) 
  endif else begin
  sr=0. & ster=0. & return
  endelse
endrep until(abs(sr_pr/sr-1) lt 1e-5) or (count_pr eq count)
mt=sr & st=ster
end
