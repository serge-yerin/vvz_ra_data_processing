PRO defDMforDSPZ,dt,fmax,fmin,DM,tot_chan,df        ; ����������������� DSPZ - �������� � ���-�� ��������
dt=fltarr(tot_chan)
for i=0,tot_chan-1 do dt[i]=(DM/2.4103d0)*((1e4/(fmin+df*(i+1))^2)-(1e4/fmax^2))
end

PRO ShowPulse, filename, DM_const, DMpos, NS, picsize, smpar, accDMnorm, DMstepnumb

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

wset,1 & XYOUTS, 20, 870, 'ESC', charsize=1.4 ,/device, FONT=2 , ORIENTATION=0, color=150L*65536+150L*256+150 & !P.COLOR=16777215L
window,7,xsize=500,ysize=800,xpos=1280-516,ypos=208, title='Individual'
window,8,xsize=1050,ysize=600,xpos=0,ypos=208, title='Cleaned data'
window,9,xsize=300,ysize=300,xpos=1280-836,ypos=208, title='Spectrum of pulse'
TimeRes=64*8192./66000000.
maxnofs=44000L
wofsg=4096L
NSframe=50L
DMstep=0.002
ShT=0
smfreq=8
pulse=fltarr(wofsg,2*NSframe+1)
filename=strmid(filename,0,strpos(filename,'.')+8)
print, filename, DM_const, DMpos, NS
DM=DM_const+0.04*(DMpos-25)
openr, inFile, filename, /get_lun
if (ns lt 2*nsframe) then ns=2*nsframe 
if (ns ge (picsize-2*nsframe-1)) then ns=picsize-2*nsframe-1 & print, ns
datucdz=assoc(inFile,FLTARR(wofsg,maxnofs),1024+4L*wofsg*(NS-2L*NSframe))
datucd=datucdz[0]
free_lun, inFile

startloop:
defDMforDSPZ,dt,33.,16.5,DM,wofsg,16.5/wofsg
  for j=0,wofsg-1 do begin
    stsp=round(dt[j]/TimeRes)+ShT+NSframe
    pulse[j,*]=datucd[j,stsp:stsp+2*NSframe]
    endfor; j
wset,8 & tvscl,-rotate(rebin(datucd[*,0:44000-1],512,1000),4),70,30
op=total(pulse[*,0:20],2)
erov, op, st, mt
nofch=2^(smfreq)
wset,9 & plot, 16.5+indgen(nofch)/(nofch-1.)*16.5,rebin((total(pulse[*,nsframe-10:nsframe+10],2)-mt)/st/sqrt(nofch/float(wofsg)),nofch), $
xrange=[16.5,33], /xstyle, pos=[50,20,310,280], /device, psym=10
XYOUTS, 10, 255, 'S/N', charsize=1.2 ,/device, FONT=2 , ORIENTATION=0
screen=tvrd() & tv,rotate(screen,1)
for j=0, wofsg-1 do pulse[j,*]=smooth(pulse[j,*],smpar,/EDG)
wset,7 & plot, (total(pulse,1)- mean(total(pulse[*,0:30],1)))/stddev(total(pulse[*,0:30],1)), pos=[100,420,403,550], /device
pulsemark=pulse
;! uncomment nrxt line for underline pulse with clear strip
;pulsemark[*,nsframe-10:nsframe+10]=pulsemark[*,nsframe-10:nsframe+10]-0.1*max(pulsemark)
wset,7 & tvscl,-rotate(rebin(pulsemark,256,(2*NSframe+1)*3),4),100,554

wset,7 & plot, (total(pulse[0:1024,*],1)- mean(total(pulse[0:1024,0:30],1)))/stddev(total(pulse[0:1024,0:30],1)), pos=[100,20,403,100], /device, /noerase
wset,7 & plot, (total(pulse[1024:2048,*],1)-mean(total(pulse[1024:2048,0:30],1)))/stddev(total(pulse[1024:2048,0:30],1)), pos=[100,120,403,200], /device, /noerase
wset,7 & plot, (total(pulse[2048:3072,*],1)-mean(total(pulse[2048:3072,0:30],1)))/stddev(total(pulse[2048:3072,0:30],1)), pos=[100,220,403,300], /device, /noerase
wset,7 & plot, (total(pulse[3072:4095,*],1)-mean(total(pulse[3072:4095,0:30],1)))/stddev(total(pulse[3072:4095,0:30],1)), pos=[100,320,403,400], /device, /noerase

backg=fltarr(60,40) & backg[*,*]=100 & wset,7 & tv,backg,0,760 ,xsize=40,ysize=40, /device
wset,7 & XYOUTS, 15, 775, 'ESC', charsize=1.4 ,/device, FONT=2 , ORIENTATION=0

backg=fltarr(100,93) & backg[*,*]=100 & wset,7 & tv,backg,403,554 ,xsize=100,ysize=93, /device
wset,7 & XYOUTS, 436, 589, ' - ', charsize=3.4 ,/device, FONT=1 , ORIENTATION=0

backg=fltarr(100,93) & backg[*,*]=100 & wset,7 & tv,backg,403,657 ,xsize=100,ysize=93, /device
wset,7 & XYOUTS, 427, 695, ' + ', charsize=3.4 ,/device, FONT=1 , ORIENTATION=0

backg=fltarr(100,93) & backg[*,*]=100 & wset,7 & tv,backg,0,554 ,xsize=100,ysize=93, /device
wset,7 & XYOUTS, 30, 589, ' > ', charsize=3.4 ,/device, FONT=1 , ORIENTATION=0

backg=fltarr(100,93) & backg[*,*]=100 & wset,7 & tv,backg,0,657 ,xsize=100,ysize=93, /device
wset,7 & XYOUTS, 30, 695, ' < ', charsize=3.4 ,/device, FONT=1 , ORIENTATION=0

backg=fltarr(60,40) & backg[*,*]=100 & wset,7 & tv,backg,14,504 ,xsize=60,ysize=40, /device
wset,7 & XYOUTS, 31, 517, 'Narr', charsize=1.4 ,/device, FONT=2 , ORIENTATION=0

backg=fltarr(60,40) & backg[*,*]=100 & wset,7 & tv,backg,14,454 ,xsize=60,ysize=40, /device
wset,7 & XYOUTS, 30, 467, 'Wide', charsize=1.4 ,/device, FONT=2 , ORIENTATION=0

wset,7 & XYOUTS, 406, 770, '  DM= '+strmid(string(DM),6,6)+'', charsize=1.4 ,/device, FONT=1 , ORIENTATION=0

wset,7 & XYOUTS, 25, 430, 'Band,', charsize=1.4 ,/device, FONT=2 , ORIENTATION=0

wset,7 & XYOUTS, 25, 415, 'kHz', charsize=1.4 ,/device, FONT=2 , ORIENTATION=0

wset,7 & XYOUTS, 25, 400, strmid(string(33000L*wofsg/float(nofch)/8192.),6,6)+'', charsize=1.4 ,/device, FONT=3 , ORIENTATION=0

wset,7 & XYOUTS, 406, 520, ' f, MHz ', charsize=1.7 ,/device, FONT=1 , ORIENTATION=0
wset,7 & XYOUTS, 406, 470, ' 16.50-33.00 ', charsize=1.7 ,/device, FONT=1 , ORIENTATION=0
wset,7 & XYOUTS, 406, 360, ' 28.88-33.00  ', charsize=1.7 ,/device, FONT=1 , ORIENTATION=0
wset,7 & XYOUTS, 406, 260, ' 24.75-28.88 ', charsize=1.7 ,/device, FONT=1 , ORIENTATION=0
wset,7 & XYOUTS, 406, 160, ' 20.63-24.75 ', charsize=1.7 ,/device, FONT=1 , ORIENTATION=0
wset,7 & XYOUTS, 406, 60, ' 16.50-20.63 ', charsize=1.7 ,/device, FONT=1 , ORIENTATION=0


wset,7 & cursor,x,y,/DOWN,/device & wait, 0.01 & print,x,y
MP=0
if ((x ge 0) and (x le 40) and (y ge 760) and (y le 800)) then goto, endloop
if ((x ge 403) and (x le 500) and (y ge 554) and (y le 647)) then MP=1  
if ((x ge 403) and (x le 500) and (y ge 657) and (y le 770)) then MP=2
if ((x ge 0) and (x le 100) and (y ge 554) and (y le 647)) then MP=3  
if ((x ge 0) and (x le 100) and (y ge 657) and (y le 770)) then MP=4
if ((x ge 14) and (x le 74) and (y ge 504) and (y le 544)) then MP=5; N
if ((x ge 14) and (x le 74) and (y ge 454) and (y le 494)) then MP=6; W

case MP of
1: begin & DM=DM-DMstep & print,MP & end
2: begin & DM=DM+DMstep & print,MP & end
3: begin & if (ShT gt -nsframe) then ShT=ShT-1 & print, MP & end 
4: begin & if (ShT lt nsframe) then ShT=ShT+1 & print, MP & end 
5: begin & if (smfreq lt 10) then smfreq=smfreq+1 & print, MP & end
6: begin & if (smfreq gt 2) then smfreq=smfreq-1 & print, MP & end
    else: MP=0
endcase
goto, startloop
endloop:
wset,7 & XYOUTS, 15, 775, 'ESC', charsize=1.4 ,/device, FONT=2 , ORIENTATION=0, color=150L*65536+150L*256+150
!P.COLOR=16777215L
STOP
set_plot,'ps'
device, filename ='d:\Survey_processing\PSRB0834.ps', /Portrait
;device, filename =filein+'_3.ps', /Portrait
DEVICE, SET_FONT='Times', /TT_FONT , BITS_PER_PIXEL=8
DEVICE, COLOR=256
;ShpsNorm=ShPs
;ShiftPar=-19

;window, 0 , xsize=800, ysize=900
plot,indgen(800),16.5+indgen(165)*0.1, /XSTYLE,/YSTYLE,charsize=1.0,/normal,THICK=5, POS=[0.1,0.75,0.55,0.95], /NODATA,font=1,ytitle='Frequency, MHz', XTICKLEN=-0.05, YTICKLEN=-0.02
tvscl,-rotate(rebin(pulsemark,256,(2*NSframe+1)*3),4), 0.1, 0.75, xsize=0.45,ysize=0.2, /normal
plot, indgen(800), (total(pulse,1)- mean(total(pulse[*,0:30],1)))/stddev(total(pulse[*,0:30],1)),font=1, /normal, /noerase, /nodata, pos=[0.1,0.4,0.55,0.65], ytitle='SNR', /XSTYLE, charsize=1.0, XTICKLEN=-0.04, YTICKLEN=-0.02
plot, (total(pulse,1)- mean(total(pulse[*,0:30],1)))/stddev(total(pulse[*,0:30],1)), font=1,/normal, pos=[0.1,0.4,0.55,0.65], /noerase, ytitle='SNR', XSTYLE=4, YSTYLE=4, charsize=1.0, XTICKLEN=-0.04, YTICKLEN=-0.02
plot,indgen(800),indgen(51)*0.004-0.1, /XSTYLE,/YSTYLE,charsize=1.0,/normal,THICK=5, POS=[0.1,0.1,0.55,0.30], /noerase, /NODATA,font=1,ytitle='    DM', xtitle='Time, msec', XTICKLEN=-0.05, YTICKLEN=-0.02
XYOUTS, 0.20, 0.98, strtrim('DM =')+strtrim(string(format='(D9.4)',DM))+strtrim('  pc cm!U -3!N'), charsize=1.4 ,/normal, FONT=1 , ORIENTATION=0
XYOUTS, 0.03, 0.17, '!7D!N', charsize=1.2 ,/normal, ORIENTATION=90
tvscl,-rebin(rotate(accDMnorm[*,ns+ShT-nsframe:ns+ShT+nsframe],4),2*(2*nsframe+1),2*DMstepnumb), 0.1, 0.1, xsize=0.45, ysize=0.2, /normal; 2*(2*nsframe+1),2*DMstepnumb),200-2*nsframe,30 ,xsize=2*(2*nsframe+1),ysize=2*DMstepnumb,/device 

device, /close
set_plot, 'win'

;for i=0,wofsg-1 do ShpsNorm[i,*]=ShPs[i,*]-mean(ShPs[i,bi:ei])
;PulseProfile=(smooth(shift(total(ShPsNorm[*,*],1),ShiftPar),1, /edg)-mean(total(ShPsNorm[*,bi:ei],1)))/STDDEV(total(ShPsNorm[*,bi:ei],1))*sqrt(1.)
;;maxScale=max(PulseProfile)
;maxScale=max(shift(AvrPulse, ShiftPar)/STDDEV(AvrPulse[bi:ei]))
;plot,indgen(Phinper)/float(Phinper-1), (smooth(shift(total(ShPsNorm[*,*],1),ShiftPar),SmPar, /edg)-mean(total(ShPsNorm[*,bi:ei],1)))/STDDEV(total(ShPsNorm[*,bi:ei],1))*sqrt(SmPar),$
;         yrange=[-5,maxScale*1.1],/ystyle,xrange=[0,1.],/xstyle, /normal, pos=[0.05,0.55,0.45,0.95], xtickname=[' ',' ',' ',' ',' ',' '], FONT=1, charsize=1.6, YTICKLEN=0.02,XTICKLEN=0.03, thick=3, ytitle='S/N'
;sigm=stddev(total(ShPsNorm[*,*],2))
;mask=fltarr(wofsg)
;mask=(abs(total(ShPsNorm[*,*],2)) le 6*sigm)
;for i=0,wofsg-1 do ShPsNorm[i,*]=ShPsNorm[i,*]*mask[i]
;plot,indgen(Phinper)/float(Phinper-1), yrange=[Fmin,Fmax],/ystyle,xrange=[0,1.],/xstyle, /normal, pos=[0.05,0.05,0.45,0.5],/nodata, /noerase, $
;         XTICKLEN=-0.03, YTICKLEN=-0.02, FONT=1, charsize=1.6, xtitle='Rotational phase', ytitle='Frequency, MHz' 
;;LOADCT,3
;tvscl,-rotate(rebin(shift(ShPsNorm,0,ShiftPar)<5>(-5),256,Phinper*9),4),0.05,0.05 ,xsize=0.4,ysize=0.45, /normal
;LOADCT,0
;tvscl,-rotate(smooth(shift(AvrPs,0,ShiftPar),[SmPar,5],/EDGE_TRUNCATE),4), 0.6,0.05 ,xsize=0.3,ysize=0.9,/normal
;plot,indgen(100)*0.01,indgen(100)*0.01, /XSTYLE,/YSTYLE,charsize=1.6,/normal,THICK=5, $
;/NOERASE,POS=[0.6,0.05,0.90,0.95],YRANGE=[-0.5,0.5],XRANGE=[0.0,1.0],/NODATA,font=1,xtitle='Rotational phase',ytitle='DM'
;
;XYOUTS, 0.086, 0.98, 'PSR '+PSRname+'    P  ='+strtrim(string(format='(D9.3)',P))+strtrim('  s'), charsize=1.8 ,/normal, FONT=1 , ORIENTATION=0
;XYOUTS, 0.66, 0.97, strtrim('DM =')+strtrim(string(format='(D9.4)',DM_const))+strtrim('  pc cm!U -3!N'), charsize=1.6 ,/normal, FONT=1 , ORIENTATION=0
;XYOUTS, 0.537, 0.45, '!7D!N', charsize=1.6 ,/normal, ORIENTATION=90
;
;scale=fltarr(1,256) & scale[0,*]=indgen(256) & tvscl,-rebin(scale,10,256),0.98,0.05 ,xsize=0.02,ysize=0.9, /normal
;;LOADCT,0
;plot,[0,1],[0,maxScale],pos=[0.98,0.05,1.,0.95],/XSTYLE,/YSTYLE, /NOERASE, /normal, THICK=5,/nodata,charsize=1.6,$
;XTICKLEN=0.0001, YTICKLEN=0.0001,FONT=1,xtickname=[' ',' ',' ',' ',' ',' '],xtitle='S/N'
;
;device, /close
;set_plot, 'win'
;STOP
END