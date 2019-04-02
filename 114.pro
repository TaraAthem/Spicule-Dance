;----------------------------------------------------------------------------------------
; EMD CODE MODIFIED BY T. MEHTA, DEC 2018
;
; THIS BIT HERE SHOULD BE THE ONLY PART YOU EDIT (EXPECT IF YOU NEED TO ADD
; SECONDARY MODES LATER IN THE CODE) BUT FOR THE MOST PART, DON'T EDIT BELOW
; THE STARRED LINE. 
; THE PROGRAM LOOPS WITH AN INPUT SF TO CONTINUE UNTIL A SIG. MODE IS FOUND FOR 
; A 95 % CI WITH THE GIVEN SF. OBVIOUSLY THE SF IS NOT NECESSAIRILY A GOOD CHOICE
; SO I WOULD STOP THE LOOP UNTIL YOU FIND A REASONABLE RANGE FOR A SF, AND THEN LOOP
; AND WATCH HOW MANY ITERATIONS ARE NEEDED UNTIL THE MODE IS SIG, AND IN THIS WAY 
; FIND AN APPROPRIATE SF. 

; THINGS TO WATCH OUT FOR:
; ->  THE MAX FOURIER PEAK MAY NOT NECESSAIRILY BE THE STATISTICALLY SIGNIFICANT 
; ONE KEEP AN EYE ON THE FOURIER OUTPUT AND CHECK BY EYE IF NEEDED.
; CODE WILL HAVE TO BE MODIFIED FOR NOT-GOES RESULTS
;
; -> QUALITY FACTOR IS DEFINED AS # OF PEAKS ABOVE THE EMDMAX/3
;
; -> SOMETIMES THE 95 AND 99 CI CROSS IF SAMPLENUM IS NOT GREAT ENOUGH
;
; GREAT! GO FORTH AND EMD RESPONSIBLY! GOOD LUCK! 
;
; -----------------------------------------------------------------------------------------


; Has an appropriate start time been found?
Answer='Y'

datafile='114_goes_2014-10-27T07:12:00_2014-10-27T07:30:00.dat'	
channel=4

;FOR GOES DATA
; 3 for 0.1-0.8nm flux
; 4 for 0.05-0.4nm flux

samplenum=1


; MUST KEEP THIS TIME AS THE ENTIRE TIME REGION
tstart='07:12:00 27-OCT-14'
tend = '07:30:00 27-OCT-14'


IF (Answer eq 'Y') THEN GOTO, LABEL3 
;***************************************************************************************

; CODE TO SPLIT UP SAMPLE INTO BITESIZE CHUNKS.

path='/home/space/phrsnz/Desktop/Academic/Programs/PhD/Comparison'  			; SPECIFIES THE DIRECTORY FOR YOUR DATA/PROGRAM ETC.
fname = path + '/dat_files/'+datafile
flarenum= STRSPLIT(datafile, '_', /EXTRACT)						; EXTRACTS THE FLARE NUMBER E.G. '004'
flaresav=STRMID(datafile,0,28)								; EXTRACTS THE FIRST HALF OF THE .DAT FILE FOR THE .SAV FILENAME
cols = 5 
data = readfile(fname, cols)
time = TRANSPOSE(data[0,*]+' '+data[1,*]) 						; DATE AND TIME
dat = TRANSPOSE(FLOAT(data[channel,*])) 						; DATASET 
tai = anytim2tai(time)	
						
ind = WHERE(tai GE anytim2tai(tstart) AND tai LE anytim2tai(tend))	

sets=1.										;Determines how many evenly 												; sized samples you want to break 												;the signal into. Default is 10. 
Print,'The following will produce 20 outputs: One time signal, and its corresponding Fourier Power spectrum, as the entire sample is split into even sized sub samples. These can be quickly manually analysed for evidence of oscillations and any samples showing indications of oscillatory bhevaiour can be put into the below section (Label3). '

set_plot, 'PS'
device, filename=path + '/Spectra/Time_plots/'+ flarenum[0]+'_'+flarenum[1]+'_TIMEPLOTS.eps', encapsulated=1,/color,bits=24 
loadct,39
!p.BACKGROUND =255
!p.color = 0
!p.charsize = 1
!p.charthick = 1

maxind=max(ind)
step=round(maxind/sets)

FOR n=0, sets-1  DO BEGIN
	
x=dat[n*step:(n+1)*step]
t=time[n*step:(n+1)*step]
dt=tai[ n_elements(x)/2]-tai[ n_elements(x)/2 -1]									
X_Deriv= DERIV(FINDGEN(N_elements(x)), x)						
x= X_Deriv/max(X_Deriv)

!p.multi = [0,1,2]

utplot, t, x, ystyle = 1								; LIGHTCURVE OF 'INPUT' DATA 

fit_fft = fft_alpha(x,dt,fap =0.05d)							; FAST FOURIER TRANSFORM WITH 95% CI
period = 1d/fit_fft.frequency
MAXFOURIERFLARE=Max(fit_fft.power,max_p)

plot,period, fit_fft.power,/xlog,/ylog,xstyle =1, ystyle =1$				; FOURIER POWER SPECTRUM
  ,yrange = minmax(fit_fft.power)*[0.1d,20d], xrange = minmax(period)*[0.5d,2d],$
  xtitle = 'Period [s]', ytitle = 'Log Power'
oplot,period, fit_fft.expectation, color = 64, thick =2
oplot,period,fit_fft.confidence_level,  color =250, linestyle= 2

ENDFOR

!p.multi = 0
device, /close
set_plot, 'X'

stop

;***************************************************************************************
LABEL3:
Print,'An appropriate time region has been selected, so the code will run on this region only'
LABEL1:								    			; LOOPING COMMAND

;Shiftfactors 
SF1=.02
SF2=0.09

;Time start and end points
tstart='07:12:00 27-OCT-14'
tend = '07:30:00 27-OCT-14'



;Change to 'Y' if you want the program to loop until a sig IMF is found for a given SF1, SF2
LOOP='N'

;***************************************************************************************

;---------------------------------------------------------------------------------------
; DEFINITIONS
;---------------------------------------------------------------------------------------

path='/home/space/phrsnz/Desktop/Academic/Programs/PhD/Comparison'  			; SPECIFIES THE DIRECTORY FOR YOUR DATA/PROGRAM ETC.
fname = path + '/dat_files/'+datafile
flarenum= STRSPLIT(datafile, '_', /EXTRACT)						; EXTRACTS THE FLARE NUMBER E.G. '004'
flaresav=STRMID(datafile,0,28)								; EXTRACTS THE FIRST HALF OF THE .DAT FILE FOR THE .SAV FILENAME
cols = 5 
data = readfile(fname, cols)
time = TRANSPOSE(data[0,*]+' '+data[1,*]) 						; DATE AND TIME
dat = TRANSPOSE(FLOAT(data[channel,*])) 						; DATASET 
tai = anytim2tai(time)	
						
ind = WHERE(tai GE anytim2tai(tstart) AND tai LE anytim2tai(tend))
x=dat[ind]
t=time[ind]
dt=tai[ n_elements(x)/2]-tai[ n_elements(x)/2 -1]					; CADENCE FOUND AT HALF WAY THROUGH DATASET INCASE CAD CHANGES DUE TO
Flare_Amp= max(x)									; THE EVENT INITIALLY (E.G. FERMI)

X_Deriv= DERIV(FINDGEN(N_elements(x)), x)						; DIFFERENTIATE TO REDUCE IMPACT OF TREND AND RESCALE DATA SO MAX=1
x= X_Deriv/max(X_Deriv)
sig_orig = x

;---------------------------------------------------------------------------------------
;COMPUTATION
;---------------------------------------------------------------------------------------
modes = emd(x - mean(x), shiftfactor=SF1, maxsiftings=1d5)				; EMD TO FIND TREND 
trend_emd = mean(x) + modes[*,-1] + modes[*,-2];+ modes[*,-3] 
set_plot, 'PS'
device, filename=path + '/Spectra/'+ flarenum[0]+'_'+flarenum[1]+'.eps', encapsulated=1,/color,bits=24, xsize=100, ysize=55

loadct, 39										
!p.BACKGROUND =255
!p.color = 0
!p.charsize = 4
!p.charthick = 1
!p.multi = [0,2,2]
utplot, t, x, ystyle = 1, ytitle = 'GOES Derivative'					; LIGHTCURVE OF 'INPUT' DATA AND TREND
outplot, t, trend_emd, color=250

Trendmax=max(trend_emd)
x=x-trend_emd
fit_fft = fft_alpha(x,dt,fap =0.05d)							; FAST FOURIER TRANSFORM WITH 95% CI
period = 1d/fit_fft.frequency
MAXFOURIERFLARE=Max(fit_fft.power,max_p)

plot,period, fit_fft.power,/xlog,/ylog,xstyle =1, ystyle =1$				; FOURIER POWER SPECTRUM
  ,yrange = minmax(fit_fft.power)*[0.1d,20d], xrange = minmax(period)*[0.5d,2d],$
  xtitle = 'Period [s]', ytitle = 'Log Power'
oplot,period, fit_fft.expectation, color = 64, thick =2
oplot,period,fit_fft.confidence_level,  color =250, linestyle= 2

modes = emd(x - mean(x), shiftfactor=SF2, maxsiftings=1d5)				; EMD TO FIND SIGNIFICANT IMF
alpha = fit_fft.pl_index								; ALPHA CORRESPONDING TO COLOURED NOISE (S=CF^ALPHA)
											; CONFIDENCE INTERVALS FOR WHITE, COLOURED, 95 AND 99% 
conf_c_95 = emd_noise_conf(n_elements(x),alpha,nsampl =samplenum, energy = fit_fft.color_energy, fap = 0.025d)
conf_c_99 = emd_noise_conf(n_elements(x),alpha,nsampl =samplenum, energy = fit_fft.color_energy, fap = 0.005d) 
conf_w_95 = emd_noise_conf(n_elements(x),0d,nsampl =samplenum, energy = fit_fft.white_energy, fap = 0.025d)
conf_w_99 = emd_noise_conf(n_elements(x),0d,nsampl =samplenum, energy = fit_fft.white_energy, fap = 0.005d)
sp = emd_energy_spectrum(modes)								; CALCULATES PERIOD USING GLOBAL WAVELET SPECTRUM
conf_up_95 = conf_c_95.up + conf_w_95.up
conf_up_99 = conf_c_99.up + conf_w_99.up
conf_period_95 = conf_c_95.period
conf_period_99 = conf_c_99.period

emd_clean = emd_reconstruct_signal(modes, conf_period_95, conf_up_95)			; RECONSTRUCTS ANY SIGNIFICANT IMFS ABOVE 95% CI
;emd_clean= modes[*,2]
maxemd=max(emd_clean)
sz=size(modes)

emax=MAKE_ARRAY(sz[2], /FLOAT)								; SZ[2] GIVES THE TOTAL NUMBER OF IMFS OUTPUTTED
emin=MAKE_ARRAY(sz[2], /FLOAT)								; EMAX, EMIN GIVE THE ERRORS ASSOCIATED WITH PERIODS
FOR N=0, sz[2]-1 DO BEGIN								; VIA A GAUSSIAN FITTING OF THE GLOBAL WAVELET SPECTRUM
yfit = GAUSSFIT(sp.per_n[N,*]*dt, sp.wave_n[N,*], coeff, NTERMS=4)			
pmax= coeff[1] + coeff[2]
pmin= coeff[1] - coeff[2]
emax[N]= pmax - sp.period[N]*dt
emin[N]= sp.period[N]*dt - pmin
ENDFOR

;---------------------------------------------------------------------------------------
; PLOTTING AND OUTPUTS
;---------------------------------------------------------------------------------------
Print,'************************************************************'
print, Format='("Alpha = ",F10.2)',+strcompress(alpha)
periods=sp.period*dt
FOR N=0, sz[2]-1 DO BEGIN
print, Format= '("EMD period (seconds) :",F10.2, " (+/-",F10.2,F10.2,")")', periods[N] , emax[N], emin[N]
ENDFOR

conf95=conf_c_95.up > conf_w_95.up
conf99=conf_c_99.up > conf_w_99.up
avg_energy=conf_c_95.mean_energy > conf_w_95.mean_energy
per_arr=conf_period_95*dt

!p.multi = [1,2,2]
plot,sp.period*dt,sp.energy, /xlog, /ylog, psym = 6, ytitle = 'Log Energy', xtitle = 'Period [s]', $
  yrange = minmax(sp.energy)*[1d/10d,20d], thick =5					; VISUALISATION OF THE EMD SPECTRUM

oplot,per_arr ,avg_energy, color =64, thick =3
oplot,per_arr,conf95, color =250, linestyle= 1
oplot,per_arr,conf99, color =250, linestyle= 5

!p.multi = [2,2,2]

utplot, t,x, ytitle = 'Detrended Signal', ystyle=1					; VISUALISATION OF DETRENDED SIGNAL AND IMFS
outplot, t, emd_clean, color =64, thick =3.5
!p.multi = 0
device, /close
set_plot, 'X'

maxfinder=EXTREMA(emd_clean, MAXIMA=maxima)
localmax=N_Elements(maxima)
TH=maxemd/3.										; THRESHOLD OF Q IS GIVEN AS 1/3 AS IT APPROXIMATES TO E^-1
mod_depth= 100.*(maxemd/Trendmax)
dur=N_elements(t)*dt
Q=0
FOR I=1, localmax DO BEGIN
	IF (emd_clean[maxima(I-1)] GT TH) THEN Q=Q+1
ENDFOR

print, 'Quality factor = ', Q
print, Format= '("Flare Maximum Amplitude = ", G10.3)',Flare_Amp
print,  Format= '("Flare modulation in % =",F10.2)', mod_depth
print, Format='("Duration (seconds) =", F10.2)', dur

IF (channel eq 4) THEN print,'Channel 0.05-0.4nm was used' ELSE print, 'Channel 0.1-0.8nm was used'

IF (MAXFOURIERFLARE gt [fit_fft.confidence_level[max_p]]) THEN print, 'There is a FFT peak of 95% sig. at', period[max_p] ELSE print,'No statistically significant peaks in the Fourier spectrum'

Print,'************************************************************'


;---------------------------------------------------------------------------------------
; CREATION OF SAV FILE
;---------------------------------------------------------------------------------------

IF (maxemd eq 0) and (LOOP eq 'Y') THEN GOTO, LABEL1					; THIS ENSURES THE PROGRAM CONTINUES RUNNING UNTIL A SIG IMF IS FOUND 
						
											; FOR A GIVEN SF-THIS MAY NOT NECESSAIRLY BE AN APPROPRIATE CHOICE OF SF.
mode1=emd_clean
trend=trend_emd
time=t
cad=dt
energy=sp.energy
save,mode1, trend,sig_orig,time,cad,periods,Q,mod_depth,dur,energy,conf95,conf99,avg_energy,per_arr,emin, emax,FILENAME=path + '/Sav_files/'+flaresav+'.sav'


end											


;---------------------------------------------------------------------------------------
; FIN
;---------------------------------------------------------------------------------------

