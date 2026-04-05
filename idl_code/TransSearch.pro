FUNCTION Ellipse, xcenter, ycenter, majoraxis, minoraxis
   points = (2 * !PI / 99.0) * FINDGEN(100)
   x = xcenter + majoraxis * COS(points )
   y = ycenter + minoraxis * SIN(points )
   RETURN, TRANSPOSE([[x],[y]])
   END

PRO TransSearch, filename, DM_const; 29 May 2015
;TransSearch, 'e:\Survey_processing\Cleaned_ V614A_F46P12A121111_223030.jds.ucd.dmt', 3.676; For Fig3 to A&A
;TransSearch, 'e:\Survey_processing\Cleaned_ PSRB0834p06A141010_032001.jds.ucd.dmt', 12.874

;NAME:
;    TransSearch
; PURPOSE:
;    Analisys of individual and repeating pulses.
; EXPLANATION::
;
; CALLING SEQUENCE:
;   TransSearch, filename, DM_const;  
;
; INPUTS
;    *.dmt - datafile
;    DM_const - central DM
;  
;
; PROGRAM DO
;   Calculationg of DM +/-25 steps (DMstep=0.004 pc*cm^-3)
;   Data normalized to the standard deviation
;   The limists of the analized data can change with the aid of the left and right mouse buttons  
;   LF-filter (parameter smpar_background changes from 32 till 4096)
;   HF-filter (parameter smpar_background changes from 1 till 8)
;   Analisys of the individual pulses yes/no
;   Time interval (+/- 50 spectra) from .ucd-file is selected with the aid of the mouse pointer
;   DM, start spectrum and bandwidth can be changed with +/-, </>, Nerr/Wide buttons 
;   Analisys of the repeating pulses yes/no
;   FFT of full, a half... 1/16 of the file (step=1/16 of the file)
;
;
;
!P.BACKGROUND=0L
!P.COLOR=16777215L
;STOP
maxshift=65536 ;43135L - time gelay for 30 pc/cc
TimeRes=64*8192./66000000.
kadr=4096L
st_limit=5.5
window,1, xsize=1094, ysize=900, title='Transient/pulsar'
DMstepnumb=0L
picsize=0L
openr,5, filename
readu,5,DMstepnumb, picsize
dt=fltarr(DMstepnumb)
accDMfile=fltarr(DMstepnumb,picsize)
accDM=fltarr(DMstepnumb,picsize)
readu,5,accDMfile
close,5
n_ev_all=0 ; number of events in the picture
ev_snr_all=0  ; above st_limit
ev_DM_all=0
ev_time_all=0
;STOP
smpar=4
smpar_pow=9
DMpos=25
smpar_background=2^smpar_pow
nsframe=50L 
NS=nsframe
Rep=0
Ind=0
nofp_pow=0
nofp=2^nofp_pow
pNum=0

setparam:

for j=0, DMstepnumb-1 do accDM[j,*]=smooth(accDMfile[j,*],smpar,/EDG)-smooth(accDMfile[j,*],smpar_background,/EDG)

wset,1 
accDMnorm=(accDM-mean(accDM))/stddev(accDM); створюємо нормований масив
maxnorm=max(accDMnorm)
minnorm=min(accDMnorm)

maxscl=maxnorm
minscl=minnorm
pmin=190
pmax=445

setscale:
stY=0
firstLine=880
secondLine=856
for np=0,15 do begin
  accDMnorm[0,4096L*np]=maxnorm; встановлюємо максимальним значенням нормованого масиву значення останньої (чи максимальне значення масиву?) точки початкового масиву для підвищення співвідношення "С/Ш" 
  accdmnorm[1,4096L*np]=minnorm; встановлюємо мінімальним значенням нормованого масиву значення першої (чи мінімальне значення масиву?) точки початкового масиву для підвищення співвідношення "С/Ш" 
  tvscl,-rebin(rotate(accDMnorm[*,4096L*np:(np+1)*4096L-1]>minscl<maxscl,4),1024,51),70,stY+np*52
  endfor; np
  
;plot,accDM[26,*]
;for np=0,15 do tvscl,-rebin(rotate(accDM[*,4096L*np:(np+1)*4096L-1],4),1024,51),70,30+np*52
;for np=0,15 do tv,(40-rebin(rotate(accDM[*,4096L*np:(np+1)*4096L-1],4),1024,51))*5,70,30+np*52  ; Доробити!!!

STOP
;TRANSSEARCH,'h:\Survey_processing\Cleaned_ PSR0834+06E111111_024207.jds.ucd.dmt', 12.88
backg=fltarr(456,50) & backg[*,*]=100 & wset,1 & tv,backg,90,850 ,xsize=456,ysize=50, /device

XYOUTS, 455, firstLine, strtrim(strmid(string(maxscl),0,strpos(maxscl,'.')+2),1), charsize=1.2 ,/device, FONT=2 , ORIENTATION=0
XYOUTS, 155, firstLine, strtrim(strmid(string(minscl),0,strpos(minscl,'.')+2),1), charsize=1.2 ,/device, FONT=2 , ORIENTATION=0
scale=fltarr(256,1) & scale[*,0]=indgen(256) & wset,1 & tvscl,-rebin(scale,256,10),190,firstLine ,xsize=256,ysize=10, /device
plot,[minnorm,maxnorm], pos=[190,firstLine,256+190,firstLine+10], /device, /nodata, /noerase, xstyle=5, ystyle=5
AXIS, XAXIS=0, XRANGE = [minnorm,maxnorm], XSTYLE = 9, XTICKLEN=-0.5
PLOTS, [pmin,pmin],[firstLine,firstLine+10], /device, thick=5, color=(2L^16)*254L
PLOTS, [pmax,pmax],[firstLine,firstLine+10], /device , thick=5, color=255

startloop:

backg=fltarr(70,50) & backg[*,*]=100 & wset,1 & tv,backg,0,850 ,xsize=70,ysize=50, /device
wset,1 & XYOUTS, 20, 870, 'ESC', charsize=1.4 ,/device, FONT=2 , ORIENTATION=0

lenSmpar=52 & stSmpar=556 & smpCh=fltarr(lenSmpar,50) & smpCh[*,*]=100 & wset,1 & tv,smpCh,stSmpar,850 ,xsize=lenSmpar,ysize=50, /device
wset,1 & XYOUTS, stSmpar, firstLine, '  < '+strmid(string(smpar),6,2)+' >  ', charsize=1.4 ,/device, FONT=2 , ORIENTATION=0
wset,1 & XYOUTS, stSmpar, secondLine, ' smpar ', charsize=1.4 ,/device, FONT=2 , ORIENTATION=0

lenSmpar_b=72 & stSmpar_b=618 & smpCh_b=fltarr(lenSmpar_b,50) & smpCh_b[*,*]=100 & wset,1 & tv,smpCh_b,stSmpar_b,850 ,xsize=lenSmpar_b,ysize=50, /device
wset,1 & XYOUTS, stSmpar_b, firstLine, '  < '+strmid(string(smpar_background),4,4)+' >  ', charsize=1.4 ,/device, FONT=2 , ORIENTATION=0
wset,1 & XYOUTS, stSmpar_b, secondLine, '  smpar_b  ', charsize=1.4 ,/device, FONT=2 , ORIENTATION=0

lenDMs=62 & stDMs=700 & DMsCh=fltarr(lenDMs,50) & DMsCh[*,*]=100 & wset,1 & tv,DMsCh,stDMs,850 ,xsize=lenDMs,ysize=50, /device
wset,1 & if ((DMpos-25)lt 0) then sign='-' else sign='+' & XYOUTS, stDMs, firstLine, '  < '+sign+strmid(string(fix(DMpos-25)),6,2)+' >  ', charsize=1.4 ,/device, FONT=2 , ORIENTATION=0
wset,1 & XYOUTS, stDMs, secondLine, 'DMs 0.04', charsize=1.4 ,/device, FONT=2 , ORIENTATION=0

backg=fltarr(70,50) & backg[*,*]=100 & wset,1 & tv,backg,772,850 ,xsize=70,ysize=50, /device
if (Ind) then !P.color=16777215L else !P.color=150L*65536+150L*256+150
wset,1 & XYOUTS, 792, firstLine-10, 'IND', charsize=1.4 ,/device, FONT=2 , ORIENTATION=0
!P.color=16777215L 

backg=fltarr(70,50) & backg[*,*]=100 & wset,1 & tv,backg,852,850 ,xsize=70,ysize=50, /device
if (Rep) then !P.color=16777215L else !P.color=150L*65536+150L*256+150
wset,1 & XYOUTS, 872, firstLine-10, 'REP', charsize=1.4 ,/device, FONT=2 , ORIENTATION=0
!P.color=16777215L 

backg=fltarr(70,50) & backg[*,*]=100 & wset,1 & tv,backg,932,850 ,xsize=70,ysize=50, /device
wset,1 & XYOUTS, 940, firstLine, '  < '+strmid(string(nofp),6,2)+' >  ', charsize=1.4 ,/device, FONT=2 , ORIENTATION=0
wset,1 & XYOUTS, 940, secondLine, '  parts  ', charsize=1.4 ,/device, FONT=2 , ORIENTATION=0

backg=fltarr(70,50) & backg[*,*]=100 & wset,1 & tv,backg,1012,850 ,xsize=70,ysize=50, /device
wset,1 & XYOUTS, 1020, firstLine, '  < '+strmid(string(pNum),6,2)+' >  ', charsize=1.4 ,/device, FONT=2 , ORIENTATION=0
wset,1 & XYOUTS, 1016, secondLine, 'N of parts', charsize=1.4 ,/device, FONT=2 , ORIENTATION=0

end_p=pNum*52+(52*16)/nofp
if (end_p gt 52*16) then del_p=end_p-52*16 else del_p=0
if (end_p gt 52*16) then begin end_p=52*16 & beg_p=pNum*52-del_p & endif ELSE BEGIN beg_p=pNum*52 & endelse & print, beg_p, end_p
pViz=fltarr(1,52*16) & pViz[*,beg_p:end_p-1]=200 & wset,1 & tv,rebin(pViz,6,52*16),50,0 ,xsize=6,ysize=52*16, /device
pNum=pNum-del_p/52

backg=fltarr(70,50) & backg[*,*]=100 & wset,1 & tv,backg,1012,850 ,xsize=70,ysize=50, /device
wset,1 & XYOUTS, 1020, firstLine, '  < '+strmid(string(pNum),6,2)+' >  ', charsize=1.4 ,/device, FONT=2 , ORIENTATION=0
wset,1 & XYOUTS, 1016, secondLine, 'N of parts', charsize=1.4 ,/device, FONT=2 , ORIENTATION=0

wset,1 & cursor,x,y,/DOWN,/device & wait, 0.01 & print,x,y   ; Для першого вікна (для зображення зі смугами) запускаємо команду "cursor", яка дозволяє клацанням курсором миші обирати точку на зображення та показує її координати.

if ((x ge 90) and (x le 556) and (y gt 842)) then begin ; scale
if (x lt 190) then x=190
if (x gt 445) then x=445
  pscl=x-190
  pch=255/(maxnorm-minnorm)
  psnr=pscl/pch+minnorm
  print,psnr
  case !MOUSE.BUTTON of
    1: minscl=psnr
    4: maxscl=psnr
    else: minscl=psnr
    endcase
  case !MOUSE.BUTTON of
    1: pmin=pscl+190
    4: pmax=pscl+190
    else: pmin=pscl+190
    endcase
  goto, setscale
  endif ; scale

MP=0
if ((x gt stSmpar) and (x lt lenSmpar/2+stSmpar) and (y gt 842)) then MP=1  
if ((x gt lenSmpar/2+stSmpar) and (x lt lenSmpar+stSmpar) and (y gt 842)) then MP=2  
  
if ((x gt stSmpar_b) and (x lt lenSmpar_b/2+stSmpar_b) and (y gt 842)) then MP=3  
if ((x gt lenSmpar_b/2+stSmpar_b) and (x lt lenSmpar_b+stSmpar_b) and (y gt 842)) then MP=4  

if ((x gt stDMs) and (x lt lenDMs/2+stDMs) and (y gt 842)) then MP=5  
if ((x gt lenDMs/2+stDMs) and (x lt lenDMs+stDMs) and (y gt 842)) then MP=6  

if ((x ge 772) and (x le 842) and (y ge 850) and (y le 900)) then MP=7  
if ((x ge 852) and (x le 922) and (y ge 850) and (y le 900)) then MP=8  

if ((x ge 932) and (x le 967) and (y ge 850) and (y le 900)) then MP=9  
if ((x ge 968) and (x le 1002) and (y ge 850) and (y le 900)) then MP=10 

if ((x ge 1012) and (x le 1047) and (y ge 850) and (y le 900)) then MP=11  
if ((x ge 1048) and (x le 1082) and (y ge 850) and (y le 900)) then MP=12 

case MP of
1: begin & Smpar=Smpar-1 & if (Smpar eq 0) then Smpar=1 & print,MP & end
2: begin & Smpar=Smpar+1 & if (Smpar eq 9) then Smpar=8 & print,MP & end
3: begin & smpar_pow=smpar_pow-1 & if (Smpar_pow eq 4) then Smpar_pow=5 & smpar_background=2^smpar_pow & print,MP & end
4: begin & smpar_pow=smpar_pow+1 & if (Smpar_pow eq 13) then Smpar_pow=12 & smpar_background=2^smpar_pow & print,MP & end 
5: begin & DMpos=DMpos-1 & if (DMpos eq (-1)) then DMpos=0 & print,MP & end  
6: begin & DMpos=DMpos+1 & if (DMpos eq 52) then DMpos=51 & print,MP & end 
7: begin & if (Ind) then Ind=0 else Ind=1 & print, MP & end 
8: begin & if (Rep) then Rep=0 else Rep=1 & print, MP & end 
9: begin & nofp_pow=nofp_pow-1 & if (nofp_pow eq (-1)) then nofp_pow=0 & nofp=2^nofp_pow & print,MP & end 
10:begin & nofp_pow=nofp_pow+1 & if (nofp_pow eq (5)) then nofp_pow=4 & nofp=2^nofp_pow & print,MP & end
11:begin & pNum=pNum-1 & if (pNum eq (-1)) then pNum=0 & print,MP & end 
;12:begin & pNum=pNum+1 & if (pNum eq (nofp)) then pNum=nofp-1 & print,MP & end
12:begin & pNum=pNum+1 & if (pNum eq (16)) then pNum=15 & print,MP & end


    else: MP=0
endcase   
if ((MP ge 1) and (MP le 4)) then  goto, setparam
if ((MP ge 5) and (MP le 12)) then  goto, startloop

if ((x lt 70) or (x gt 1094) or (y lt stY)) then goto, endloop ; Встановлюємо умови, за яких (щоб не вилізти за межі файлу) потрібно повернутися назад (?)
print,x-70,y-stY & x=x-70 & y=y-stY ; "Обрізаємо" вікно, лишаючи лише інформативну частину.
x=(x)>0<1023 & y=((y)>0<(52*16)) ; Встановлюємо межі, в яких лежать значення x та y.
print,floor((y)/52.) & fr=floor((y)/52.) ; Визначаємо номер смуги, в якій лежить задана точка і називаємо це значення "fr".
;DMpos=y-(fr)*52
NS=4096L*fr+x*4 ; Переводимо "відносне" значення x в "абсолютне", тобто за координатами у вікні визначаємо положення даної точки в масиві даних. 
print,ns ; Відображаємо "абсолютне" значення
window,2,xsize=400,ysize=400,xpos=0,ypos=208, title='Pulse selection' ; Створюємо нове вікно (400*500 точок)
nsframe=50L ; Задаємо к-ть точок зліва і справа від нашої точки, яка буде відображатися
if (ns lt nsframe) then ns=nsframe 
if (ns ge (picsize-nsframe-1)) then ns=picsize-nsframe-1 & print, ns
wset,2 & SURFACE,rotate(accDMnorm[*,ns-nsframe:ns+nsframe]>minscl<maxscl,4),pos=[40,200,360,400],/device ; У новому вікні (2) за допомогою тривимірного графіка ми зображаємо обране місце в масиві (центральна точка і по 50 зліва та справа)
wset,2 & tvscl,-rebin(rotate(accDMnorm[*,ns-nsframe:ns+nsframe]>minscl<maxscl,4),2*(2*nsframe+1),2*DMstepnumb),200-2*nsframe,30 ,xsize=2*(2*nsframe+1),ysize=2*DMstepnumb,/device ; У новому вікні (2) за допомогою тривимірного графіка ми зображаємо обране місце в масиві (центральна точка і по 50 зліва та справа)
;window,3,xsize=400,ysize=300
;wset,3 & plot, accDMnorm[DMpos,*]>minscl<maxscl, psym=4
;;wset,1 & backg=fltarr(200,20) & tv,backg,785,870 ,xsize=200,ysize=20, /device
;;XYOUTS, 785, 875, 'DMstep = '+strtrim(string(DMpos-25),1), charsize=1.2 ,/device, FONT=2 , ORIENTATION=0
;window,5,xsize=500,ysize=300
;fftline=abs(fft(accDMnorm[DMpos,*]>minscl<maxscl))^2
;wset,5 & plot,(fftline[40:4000])

;fftarr=fltarr(picsize,61)
;for iSNR=0,60 do begin 
;  fftarr[*,iSNR]=abs(fft(accDMnorm[DMpos,*]>(iSNR/10.-3)<maxscl))^2
;  fftarr[*,iSNR]=(fftarr[*,iSNR]-max(fftarr[41:4040,iSNR]))/stddev(fftarr[41:4040,iSNR])
;  endfor; iSNR
;window,6,xsize=500,ysize=400,xpos=1280-516,ypos=208
;wset,6 & tvscl,-rebin(fftarr[41:4040,*],400,61*4),20,20
;wset,6 & surface,rebin(fftarr[41:840,*],400,61*4)

;fftsum=fltarr(picsize*5L)
;fftsum[0:picsize/2]=fftline[0:picsize/2]
;for iGarm=1L,picsize/2-1 do fftsum[iGarm]=fftsum[iGarm]+fftsum[iGarm*2]+fftsum[iGarm*3]+fftsum[iGarm*4]+fftsum[iGarm*5]+fftsum[iGarm*6]
;wset,6 & plot,(fftsum[40:4000])
if Rep then RepeatingAnalisys,accDMnorm,DMpos,minscl,maxscl, pViz
if Ind then ShowPulse, filename, DM_const, DMpos, NS, picsize, smpar, accDMnorm, DMstepnumb ;
wset,1 & XYOUTS, 20, 870, 'ESC', charsize=1.4 ,/device, FONT=2 , ORIENTATION=0, color=16777215L
goto, startloop ; Частина циклу 
endloop:
wset,1 & XYOUTS, 20, 870, 'ESC', charsize=1.4 ,/device, FONT=2 , ORIENTATION=0, color=150L*65536+150L*256+150 & !P.COLOR=16777215L
STOP

;  for np=0L,15L do begin
;  n_ev=0 ; number of events in the picture
;  ev_snr=0  ; above st_limit
;  ev_DM=0
;  ev_time=0
;  wset,1 & plot, indgen(521)/16.+kadr*np*TimeRes,intarr(521), /XSTYLE,/YSTYLE, POS=[0.063,0.05,0.996,0.986],charsize=1.8,/normal, yrange=[0.,30.], font=1,$
;        YTITLE = 'DM, pc cm^-3',  XTITLE = 'Time from file start, s', /nodata
;;XYOUTS, 0.0375, 0.54, ' pc cm!U -3!N', charsize=1.8 ,/normal, FONT=1 , ORIENTATION=90
;  tvscl,-rebin(rotate(accDM[1:3000,4096*np:(np+1)*4096-1],4),1024,1000),70,55 
;  plot, indgen(521)/16.+kadr*np*TimeRes,intarr(521), /XSTYLE,/YSTYLE, POS=[0.063,0.05,0.996,0.986],charsize=1.8,/normal, yrange=[0.,30.], font=1,$
;        YTITLE = 'DM, pc cm^-3',  XTITLE = 'Time from file start, s', /nodata, /noerase       
;  for j=0,DMstepnumb-1 do begin
;  ;indDM=accDM[j,kadr*np:(np+1)*kadr-1]
;  ;indDM=accDM[j,kadr*np:(np+1)*kadr-1]-smooth(accDM[j,kadr*np:(np+1)*kadr-1],smpar,/EDG)
;  indDM=accDM[j,kadr*np:(np+1)*kadr-1]
;  op=indDM
;  erov, op, st, mt
;  indDM=indDM-mt
;  ;print, max(indDM)/st,j
;  
;  for i=0L,kadr-1 do begin
;  if((indDM[i]/st_limit/st) gt 1.) then begin
;    n_ev=n_ev+1
;    ev_snr= [ev_snr,indDM[i]/st]
;    ev_DM= [ev_DM, j/100.]
;    ev_time= [ev_time,(np*timeres*kadr+i*timeres)]
;    n_ev_all=n_ev_all+1
;    ev_snr_all= [ev_snr_all,indDM[i]/st]
;    ev_DM_all= [ev_DM_all, j/100.]
;    ev_time_all= [ev_time_all,(n_file*picsize*timeres+np*timeres*kadr+i*timeres)]
;    
;    PLOTS, Ellipse(0.063+i/4096.*0.933, 0.05+j/3000.*0.936, 0.002*round(indDM[i]/4./st),0.002*round(indDM[i]/4./st)), /normal, color=(alog10(indDM[i]/4.)*600.)<255
;  endif
;  endfor ; i
;  ;for i=0,1024 do PLOTS, CIRCLE(i, j*30+100, round(indDM[i]/4/st)), /Device
;  endfor ; j
;  if (n_ev gt 0) then begin
;  ev_snr=ev_snr[1:*]
;  ev_DM=ev_DM[1:*]
;  ev_time=ev_time[1:*]
;  endif; n_ev
;  print, 'finished pict 4096x3000 number: ',  strcompress(string(np+1), /REMOVE_ALL), '  in  ', string(systime(1)-t2), 'sec'
;  Im=tvrd(/TRUE)
;  write_jpeg, AddPath + shortname[n_file] + strmid(string(np+1),6)+'.jpg', Im, /TRUE
;  ;STOP;
;endfor; np
;close,1
;print, 'finished file of 65536 spectra number: ', n_file+1 , '  in  ', string(systime(1)-t1), 'sec' 
;
;print, 'All files finished in', (systime(1)-tst),' sec'
;close, 5
;ev_snr_all=ev_snr_all[1:*]
;ev_DM_all=ev_DM_all[1:*]
;ev_time_all=ev_time_all[1:*]
;openw, 4, Addname+'.dat'
;writeu,4,n_ev_all
;writeu,4,ev_snr_all,ev_DM_all,ev_time_all
;close, 4
;print,'max(ev_snr_all[where(ev_DM_all gt 1.)]) = ',max(ev_snr_all[where(ev_DM_all gt 1.)])
;STOP
END

