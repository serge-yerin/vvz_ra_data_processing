PRO defDMforDSPZ,dt,fmax,fmin,DM,tot_chan,df   ;
dt=fltarr(tot_chan)
for i=0,tot_chan-1 do dt[i]=(DM/2.4103d0)*((1e4/(fmin+df*(i+1))^2)-(1e4/fmax^2))
end

PRO IndSearch, CleanedFileName, DM_Const; -----!MAIN PROGRAM!-----

HStr={Headsrtuct,Sname:BYTARR(32),Stime:BYTARR(32),Sgmtt:BYTARR(32),$
       Ssysn:BYTARR(32),Ssyst:BYTARR(32),Splace:BYTARR(96),Sdesc:BYTARR(256), $
       Sdspp:ULONARR(128)}   & HeadOffset=16 & ShiftAB=1

Name= CleanedFileName
;Name= DIALOG_PICKFILE(/READ, FILTER = '*.ucd',/ MULTIPLE_FILES)
OPENR,1, name
readu,1,HStr
print,'File name:',string(HStr.Sname),'   Local
time:',string(HStr.Stime),'   UTC time:',string(HStr.Sgmtt), $
      '   Record:',string(HStr.Ssysn),'
Place:',string(HStr.Splace),'   Descriptor:',string(HStr.Sdesc)
print,HStr.Sdspp[6:32]

nofs=HStr.sdspp[10+HeadOffset]        ;
Fclk=66.000000/2.
Fmin=HStr.Sdspp[12+HeadOffset]*Fclk/8192;
Fmax=HStr.Sdspp[13+HeadOffset]*Fclk/8192;
wofsg=HStr.Sdspp[14+HeadOffset]              ;
avrs=HStr.Sdspp[15+HeadOffset]               ;
TimeRes=avrs*8192./2./(Fclk*1E6)
fstatus = FSTAT(1) & n_kadr = floor((fstatus.size - 1024)/(4L*nofs*wofsg))
datdspz=assoc(1,fltarr(wofsg,nofs),1024)
data=fltarr(wofsg,nofs)
i_kadr=n_kadr   ; number of reading picts
defDMforDSPZ,dt,fmax,fmin,DM_Const,wofsg,(fmax-fmin)/wofsg
maxshift=65536L
range=256

temp=fltarr(range,maxshift+i_kadr*nofs)
DMstep=0.004
Dmstepnumb=50
picsize=65536L
accDM=fltarr(Dmstepnumb+1,maxshift+i_kadr*nofs)
accSB=fltarr(wofsg/range,maxshift+i_kadr*nofs)

for l=0,Dmstepnumb do begin
DM = DM_Const-Dmstepnumb/2*DMstep+l*DMstep
defDMforDSPZ,dt,fmax,fmin,DM,wofsg,(fmax-fmin)/wofsg
t=SYSTIME(1)
window,0,xsize=1100,ysize=550
fqbeg=0
fqend=wofsg/range
for k=fqbeg,fqend-1 do begin   ;
rbeg=k*range               ;
   temp[*,*]=0.
   for j=0,i_kadr-1 do begin
   data=datdspz[j]
   temp[0:range-1,maxshift+j*nofs:maxshift+(j+1)*nofs-1]=data[rbeg:rbeg+range-1,*]
   ;wset,0 & tvscl,-rebin(rotate(data[rbeg:rbeg+range-1,*],4),nofs,512)
   endfor; i_kadr
for i=0,range-1 do begin
  temp[i,*]=shift(temp[i,*],-round(dt[i+rbeg]/TimeRes))
endfor ;i
accSB[k,*]=total(temp,1)

print, 'DM step',l,'  from',Dmstepnumb,'  subband ' ,k, '  processing
time' ,SYSTIME(1)-t & t=SYSTIME(1)
endfor ; k - subband loop
accDM[l,*]=total(accSB,1)
endfor; l - DM loop

;smooth - add


openw, 5, Name+'.dmt'
writeu,5, Dmstepnumb+1L,picsize
writeu,5, accDM[*,maxshift:picsize+maxshift-1]
close, 5
;STOP

END
