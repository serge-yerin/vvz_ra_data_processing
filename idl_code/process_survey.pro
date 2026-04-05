PRO filelist, listdatfile, shortname, nunb_in_list       ; 
listdatfile = DIALOG_PICKFILE(/READ, FILTER = '*.jds',/ MULTIPLE_FILES)
nunb_in_list=n_elements(listdatfile)
shortname = listdatfile
for i=0,nunb_in_list-1 do begin
  PathName = Str_Sep(listdatfile[i],'\')
  shortname[i] = PathName[n_elements(PathName)-1]
 endfor
print,listdatfile, shortname, nunb_in_list, ' files in list'
end

PRO PROCESS_SURVEY ; ___________START______________ 
;NAME:
;    PROCESS_SURVEY
; PURPOSE:
;    processing two sucsessive files *.jds: cleaning, dedisp, individual (repetitive) pulse search, FFT in 30-sec - 8-min sequences  
; EXPLANATION:
;
; CALLING SEQUENCE:
;   PROCESS_SURVEY ... and
;    load 2 files: first with SAME NUMBER as file with pulses over 5.5 sigma  
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
;pathdat = 'E:\Survey_processing\' & DM = 5.752 & P = 1.29224132384D
pathdat = 'd:\Survey_processing\' & DM = 5.752 & P = 1.29224132384D
filelist, listdatfile, shortname, n_in_list

offset=[70,50]  ; offset for picture
mode=1       ; 0 - multiplication, 1 - sum, 2 - difference
nofs=1024     ;
;DM_Const = 2.82 & P = 1.292292258039D & AddTitle=' V137A_F55P11' ; change here
;DM_Const = 3.676 & P = 1.292292258039D & AddTitle=' V614A_F46P12' ; change here F46_12_ DM_3,676_V614_A_Nov2011_t_23810,1
;DM_Const = 5.752 & P = 1.292292258039D & AddTitle=' PSRB0809p74' ; change here F46_12_ DM_3,676_V614_A_Nov2011_t_23810,1
;DM_Const = 15.3185 & P = 1.292292258039D & AddTitle=' PSRB0943p10' ; change here F46_12_ DM_3,676_V614_A_Nov2011_t_23810,1
DM_Const = 12.872 & P = 1.292292258039D & AddTitle=' PSRB0834p06' ; change here F46_12_ DM_3,676_V614_A_Nov2011_t_23810,1
;DM_Const = 4.970 & P = 1.292292258039D & AddTitle=' V683D_F57P10' ; change here   F57_10_DM_4,970_V683_D_Nov2011_t_29461,6
;DM_Const = 24.450001 & P = 1.292292258039D & AddTitle='J1316+2033_21Nov2020_F5P4' ; change here

;   Format           DSPZ AB ( 2008 year.)
HStr={Headsrtuct,Sname:BYTARR(32),Stime:BYTARR(32),Sgmtt:BYTARR(32),$
       Ssysn:BYTARR(32),Ssyst:BYTARR(32),Splace:BYTARR(96),Sdesc:BYTARR(256), $
       Sdspp:ULONARR(128)}   & HeadOffset=16 & ShiftAB=1

if (mode eq 0) then foutname='Clean '+ ' Mul.ucd '
if (mode eq 1) then foutname='Clean '+ ' Ch1.ucd'
if (mode eq 2) then foutname='Clean '+ ' Ch2.ucd'
namestr=''

for i=0,n_elements(shortname)-1 do begin
dotrem=strpos(shortname[i], '.')
namestr=namestr+'-'+STRMID(shortname[i],0,dotrem)
endfor ; i

OPENR, 1, listdatfile[0]
readu,1,HStr
Close,1

CleanedFileName=pathdat+'Cleaned_'+addtitle+shortname[0]+'.ucd'
OPENW, 2,CleanedFileName
HStr.sdspp[10+HeadOffset]=ULONG(nofs)     ; put into header ucd value NOFS
writeu,2,HStr

for n_file = 0, n_in_list-1 do begin

OPENR, 1, listdatfile[n_file]
shortname[n_file]=shortname[n_file]+AddTitle
readu,1,HStr
print,'File name:',string(HStr.Sname),'   Local time:',string(HStr.Stime),'   UTC time:',string(HStr.Sgmtt), $
      '   Record:',string(HStr.Ssysn),'   Place:',string(HStr.Splace),'   Descriptor:',string(HStr.Sdesc)
print,HStr.Sdspp[6:32]
if (HStr.Sdspp[8+HeadOffset] eq 0) then print,'waveform'
if (HStr.Sdspp[8+HeadOffset] eq 1) then print,'spectra'
if (HStr.Sdspp[8+HeadOffset] eq 2) then print,'correlation'
Fmin=HStr.Sdspp[12+HeadOffset]*33.000000/8192; number of frequency channel LOW
Fmax=HStr.Sdspp[13+HeadOffset]*33.000000/8192; number of frequency channel HIGH
wofsg=HStr.Sdspp[14+HeadOffset]              ; working frequency range
avrs=HStr.Sdspp[15+HeadOffset]               ; average number
TimeRes=avrs*8192./66000000.
;STOP
FDatTim={Day:string(2),Mon:string(2),Year:string(2),Hour:string(2),Mins:string(2),Sec:string(2)}
FDatTim.Day=STRMID(shortname[n_file],0+ShiftAB,2)
FDatTim.Mon=STRMID(shortname[n_file],2+ShiftAB,2)
FDatTim.Year=STRMID(shortname[n_file],4+ShiftAB,2)
FDatTim.Hour=STRMID(shortname[n_file],7+ShiftAB,2)
FDatTim.Mins=STRMID(shortname[n_file],9+ShiftAB,2)
FDatTim.Sec=STRMID(shortname[n_file],11+ShiftAB,2)
print,FDatTim
StartTime=FDatTim.Hour*3600.+FDatTim.Mins*60+FDatTim.Sec   ; 6 ����� - ��� ����������
Xdat=indgen(nofs)
Ydat=intarr(fix(nofs*TimeRes)) & for l=0,fix(nofs*TimeRes)-1 do Ydat[l]=fmax
date_time = TIMEGEN(fix(nofs*TimeRes), UNITS = 'Seconds', $
START = JULDAY(FDatTim.Mon, FDatTim.Day, FDatTim.Year, FIX(FDatTim.Hour), FIX(FDatTim.Mins), FIX(FDatTim.Sec)))
date_label = LABEL_DATE(DATE_FORMAT = ['%I:%S'])
;STOP
displaySize=[nofs,512]

if (mode eq 0) then WTitle='Spectrograms - '+'Dspz ' + shortname[n_file]+ ' -Mul- sDM'+strmid(string(sDM),3)+ ' file '+ n_file + ' from '+ n_in_list
;if (mode eq 1) then WTitle='Spectrograms - '+'Dspz ' + shortname[n_file]+ ' -Ch1- sDM'+strmid(string(sDM),3) + ' file '+ n_file + ' from '+ n_in_list
if (mode eq 1) then WTitle=n_file
if (mode eq 2) then WTitle='Spectrograms - '+'Dspz ' + shortname[n_file]+ ' -Ch2- sDM'+strmid(string(sDM),3) + ' file '+ n_file + ' from '+ n_in_list

window, 1, xsize=displaySize[0]+offset[0]+20, ysize=displaySize[1]+offset[1]+40, TITLE=WTitle ; window of Spectrogram
wshow, 1, 0   ; hide this windows
datdspz=assoc(1,ULONARR(2,wofsg,nofs),1024)     ; structure of "spectra" data

fstatus = FSTAT(1) & nframe = floor((fstatus.size - 1024)/(4L*nofs*wofsg*2L))
TT = SYSTIME(1)
for i=0,nframe-1 do begin
 dataz=datdspz[i]

 chek0=reform(dataz[0,*,*],wofsg,nofs)
 chek1=reform(dataz[1,*,*],wofsg,nofs)
 ;wset,0 & plot,chek0[wofsg-1,*],/YNOZ
 ;wset,2 & plot,chek1[wofsg-1,*],/YNOZ

;STOP
 mant=(dataz and 4294967232)  ; mantissa DSPZ
 expn=(dataz and 31)          ; exponenta DSPZ
 data=double(mant)/double(long64(2)^expn)*4*2*1024./4294967296./avrs

 if (mode eq 0) then imdat=reform(data[0,*,*]-data[1,*,*],wofsg,nofs)
 if (mode eq 1) then imdat=reform(data[0,*,*],wofsg,nofs)
 if (mode eq 2) then imdat=reform(data[1,*,*],wofsg,nofs)

!P.COLOR=0
!P.BACKGROUND=16777215
!P.FONT=-1

 pat=1 & nar=1 & wid=1


wset,1 & PLOT, date_time, Ydat, /XSTYLE, /YSTYLE,  $
XRANGE=[date_time[0]+i*double(nofs)*TimeRes/86400.,date_time[0]+(i+1)*double(nofs)*TimeRes/86400.], YRANGE=[fmin,fmax],$
TITLE = 'DSPZ_' + shortname[n_file]+ strmid(string(i+1),5), XTITLE = 'Time, sec', $
YTITLE = 'Frequency, MHz', /NODATA, /DEVICE, XTICKLEN=-0.01, YTICKLEN=-0.01,$
POSITION = [offset[0]-1, offset[1]-1, $
displaySize[0] + offset[0], $
displaySize[1] + offset[1]], CHARSIZE = 1.5,$
XTICKFORMAT = 'LABEL_DATE', $
XTICKUNITS = 'Seconds', $
XTICKINTERVAL = 4
wset,1 & tvscl,-rebin(10*alog10(rotate(imdat,4))>(-12)<(-3),nofs,512),offset[0],offset[1]
print,'Frame', i+1, '    from ', nframe, ',    Time =',SYSTIME(1)-TT & TT=SYSTIME(1)

adr_cleaning, imdat, mask, pat=pat, nar=nar, wid=wid
wset,1 & tvscl,-rebin(rotate(imdat,4),nofs,512),offset[0],offset[1]
DspzIm=tvrd()
write_jpeg, pathdat +'DSPZ_'+ shortname[n_file]+ '_'+Addtitle + '_piece_of_'+strtrim(string(nofs),2)+'_spectra_#'+strtrim(string(i+1),2)+'_2Cleaned.jpg', DspzIm
wset,1 & tvscl,-rebin(rotate(mask,4),nofs,512),offset[0],offset[1]
DspzIm=tvrd()
write_jpeg, pathdat +'DSPZ_'+ shortname[n_file]+ '_'+Addtitle + '_piece_of_'+strtrim(string(nofs),2)+'_spectra_#'+strtrim(string(i+1),2)+'_3Bad_pixels.jpg', DspzIm

WRITEU, 2, float(imdat)
endfor; kadr
close,3
close,1
endfor; file
close,2
IndSearch, CleanedFileName, DM_Const
TransSearch, CleanedFileName+'.dmt', DM_Const

stop

END