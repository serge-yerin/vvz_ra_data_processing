PRO RepeatingAnalisys,accDMnorm,DMpos,minscl,maxscl,pViz

NAME:
;    RepeatingAnalisys
; PURPOSE:
;    Analisys of repeating pulses.
; EXPLANATION:
;
; CALLING SEQUENCE:
;   RepeatingAnalisys,accDMnorm,DMpos,minscl,maxscl,pViz;  
;
; INPUTS
;    accDMnorm - array 51x65536
;    DMpos - DMstep
;    minscl - minimal limit level
;    maxscl - maximal limit level
;    pViz - array 1x(16*52) with the information about selected time window
;      
;
; PROGRAM DO
;   Plots graph S/N vs Time and FFT of the sselected time interval

timeline=accDMnorm[DMpos,*]>minscl<maxscl

winLong=fltarr(65536L)
winFun=fltarr(16L)
winFun[*]=rebin(pViz[0,*],1,16)
for i=0,15L do winLong[i*4096L:(i+1L)*4096L-1L]=winFun[i]/200.
plotmin=min(where(winlong gt 0)) & plotmax=max(where(winlong gt 0))

window,3,xsize=400,ysize=300, title='S/N vs Time'
wset,3 & plot,plotmin+lindgen(plotmax-plotmin), timeline(where(winlong gt 0)), psym=4

window,5,xsize=500,ysize=300, title='FFT'
timeline=accDMnorm[DMpos,*]>minscl<maxscl
timeline=timeline*winLong
;wset,5 & plot, timeline
;STOP
fftline=abs(fft(timeline>minscl<maxscl))^2
fftline[0:40]=fftline[8000:8040]
wset,5 & plot,(fftline[0:4000]), /xstyle, xrange=[0,4000]

;STOP
END