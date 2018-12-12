
LABEL1:
fname = '008_goes_2014-09-23T23:05:00_2014-09-24T00:20:00.dat'

;FOR GOES DATA
; 4 for 0.05-0.4nm flux
; 3 for 0.1-0.8nm flux
channel=4

;Shiftfactors 
SF1=0.25
SF2=0.21258

;Time start and end points
tstart = '23:08:00 23-SEPT-14'
tend = '23:13:45 23-SEPT-14'

;************************************************************

flarenum= STRSPLIT(fname, '_', /EXTRACT)
cols = 5 

data = readfile(fname, cols)
time_norh = TRANSPOSE(data[0,*]+' '+data[1,*]) ;Date and times
norh = TRANSPOSE(FLOAT(data[channel,*])) ;NoRH correlation values

tai_norh = anytim2tai(time_norh)

ind_norh = WHERE(tai_norh GE anytim2tai(tstart) AND tai_norh LE anytim2tai(tend))
time = time_norh[ind_norh]
norh = norh[ind_norh]

x=norh
t=time
dt=tai_norh[30]-tai_norh[29]

Flare_Amp= max(x)
Time_array=FINDGEN(N_elements(x))

;window, 1
;utplot, t, x

X_Deriv= DERIV(Time_array, x)
X_Deriv= X_Deriv/max(X_Deriv)
x=X_Deriv
x_original = x

;------FIRST EMD FOR TREND------
modes = emd(x - mean(x), shiftfactor=SF1, maxsiftings=1d5)
trend_emd = mean(x) + modes[*,-1]

;-------PLOTTING------------
CD, 'Spectra'

set_plot, 'PS'
device, filename='Spectra_'+ flarenum[0]+'_'+flarenum[1]+'.eps', encapsulated=1,$	
/color,bits=24, xsize=100, ysize=55

loadct, 39
!p.BACKGROUND =255
!p.color = 0
!p.charsize = 4
!p.charthick = 1
!p.multi = [0,2,2]
utplot, t, x, ystyle = 1, ytitle = 'GOES Derivative'
outplot, t, trend_emd, color=250
;--------------------------

Trendmax=max(trend_emd)
x=x-trend_emd
fit_fft = fft_alpha(x,dt,fap =0.05d)
period = 1d/fit_fft.frequency
MAXFOURIERFLARE=Max(fit_fft.power,max_p)

;------PLOTTING--------------
plot,period, fit_fft.power,/xlog,/ylog,xstyle =1, ystyle =1$
  ,yrange = minmax(fit_fft.power)*[0.1d,20d], xrange = minmax(period)*[0.5d,2d],$
  xtitle = 'Period [s]', ytitle = 'Log Power'
oplot,period, fit_fft.expectation, color = 64, thick =2
oplot,period,fit_fft.confidence_level,  color =250, linestyle= 2
;oplot, [period[max_p]], [fit_fft.confidence_level[max_p]],  psym = 6, thick=10
;----------------------------


;--------SECOND EMD FOR MODES ---

modes = emd(x - mean(x), shiftfactor=SF2, maxsiftings=1d5)
alpha = fit_fft.pl_index
; Confidence intervals for coloured noise
conf_c_95 = emd_noise_conf(n_elements(x),alpha,nsampl =50, energy = fit_fft.color_energy, fap = 0.025d) ;95% CI
conf_c_99 = emd_noise_conf(n_elements(x),alpha,nsampl =50, energy = fit_fft.color_energy, fap = 0.005d) ;99% CI
; Confidence intervals for the white noies
conf_w_95 = emd_noise_conf(n_elements(x),0d,nsampl =50, energy = fit_fft.white_energy, fap = 0.025d)
conf_w_99 = emd_noise_conf(n_elements(x),0d,nsampl =50, energy = fit_fft.white_energy, fap = 0.005d)
sp = emd_energy_spectrum(modes)
;Upper confidence interval for the combined noises
conf_up_95 = conf_c_95.up + conf_w_95.up
conf_up_99 = conf_c_99.up + conf_w_99.up
conf_period_95 = conf_c_95.period
conf_period_99 = conf_c_99.period
emd_clean = emd_reconstruct_signal(modes, conf_period_95, conf_up_95)
;emd_clean= modes[*,2];+modes[*,3]
maxemd=max(emd_clean)
Print,'************************************************************'
print, Format='("Alpha=",F10.2)',+strcompress(alpha)
periods=sp.period*dt
print, Format= '("Mean Periods are:",F10.2)', periods 
;----------------------------


;------PLOTTING--------------
!p.multi = [1,2,2]
plot,sp.period*dt,sp.energy, /xlog, /ylog, psym = 6, ytitle = 'Log Energy', xtitle = 'Period [s]', $
  yrange = minmax(sp.energy)*[1d/10d,20d], thick =5

oplot, conf_period_95*dt,conf_c_95.mean_energy > conf_w_95.mean_energy, color =64, thick =3
oplot,conf_period_95*dt,conf_c_95.up > conf_w_95.up, color =250, linestyle= 1
oplot,conf_period_99*dt,conf_c_99.up > conf_w_99.up, color =250, linestyle= 5

!p.multi = [2,2,2]

utplot, t,x, ytitle = 'Detrended Signal', ystyle=1
outplot, t, emd_clean, color =64, thick =2
!p.multi = 0
device, /close
set_plot, 'X'
CD, '..'
CD, 'Sav_files'
;-----------CREATION OF SAV FILE--------
maxfinder=EXTREMA(emd_clean, MAXIMA=maxima)
localmax=N_Elements(maxima)
TH=maxemd/3.						;Threshold, defined as 1/3 of flare max
flare_mod= 100.*(maxemd/Trendmax)
dur=N_elements(t)*dt					;given in seconds. 
Q=0
FOR I=1, localmax DO BEGIN
	IF (emd_clean[maxima(I-1)] GT TH) THEN Q=Q+1
ENDFOR

print, 'Q-factor = ', Q
print, Format= '("Flare Maximum Amplitude = ", G10.3)',Flare_Amp
print,  Format= '("Flare modulation in % =",F10.2)', flare_mod
print, Format='("duration (seconds) =", F10.2)', dur

IF (channel eq 4) THEN print,'Channel 0.05-0.4nm was used' ELSE print, 'Channel 0.1-0.8nm was used'

IF (MAXFOURIERFLARE gt [fit_fft.confidence_level[max_p]]) THEN print, 'There is a FFT peak of 95% sig. at', period[max_p] ELSE print,'No statistically significant peaks in the Fourier spectrum'

Print,'************************************************************'

IF (maxemd eq 0) THEN GOTO, LABEL1

emd_mode=emd_clean
emd_trend=trend_emd
x=x_original
time=t
cad=dt

;emd_mode, emd_trend,x,time,cad,periods,Q,flare_mod,dur
save,emd_mode, emd_trend,x,time,cad,periods,Q,flare_mod,dur,FILENAME=flarenum[0]+'_SAV.sav'

CD, '..'
end
