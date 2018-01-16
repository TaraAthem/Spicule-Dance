pro broken_power_law,x,p,model
  model = -(x LT p[1])*p[0]*x + $
    (x GT p[1])*(-(p[0]-p[3])*p[1] - p[3]*x) + p[2]
  return
end

pro linfun, x, a, f
  f = a[0] + a[1] * x
  return
end


flarename= 'flare106440.fits'
data = mrdfits(flarename, 1, header)
flux=data.flux
tim=data.time
cad=tim[1]-tim[0]
FLARETITLE=STRSPLIT(flarename, '.fits', /EXTRACT, /REGEX)
FLARETITLE=STRSPLIT(FLARETITLE, 'flare', /EXTRACT)

ns=126
nf=201
device, decomposed=0 & loadct,39


;*****************************************************************************************
;---------------------------------------PREFLARE------------------------------------------
;*****************************************************************************************

print,' **********************************'
print,' *********** PRE FLARE ************'
print,' **********************************'

flux=data.flux
tim=data.time
cad=tim[1]-tim[0]
flux_preflare=flux[0:ns]
tim_preflare=tim[0:ns]


N_preflare=n_elements(flux_preflare)
sh=0.007
modes=emd(flux_preflare-mean(flux_preflare),shi=sh, epsilon=1d-6, maxsiftings=1d8)
sz=size(modes)
Print, 'Size of (modes)', sz
Print, 'SZ[2] -1', sz[2]-1
Print, 'SZ[2] -2', sz[2]-2
trend_preflare=mean(flux_preflare)+modes[*,sz[2]-1]+modes[*,sz[2]-2]
;emd_mode=modes[*,sz[2]-2]


p = [2.D, -1.D, 0.1D, 0.D] ;Initial guesses of model parameters
result = LNP_TEST(tim_preflare, flux_preflare-trend_preflare, /DOUBLE, wk1=freq, wk2=pow, ofac=1)
npow = N_ELEMENTS(pow)
freqPRE = ALOG(freq)
powPRE = ALOG(pow)

Periodogram, tim_preflare,flux_preflare-trend_preflare, fap=0.05, power_sign=FFT95PRE			;FFT CONDIFENCE INTERVALS USING FALSE ALARM PROBABILITY
Periodogram, tim_preflare,flux_preflare-trend_preflare, fap=0.01, power_sign=FFT99PRE			;FOR 95 AND 99%

err = REPLICATE(1., npow)
pow_fitPRE = CURVEFIT(freqPRE[1:*], powPRE[1:*], err[1:*], p, FUNCTION_NAME='broken_power_law',fita=[1,1,1,0],/noderivative)
result=p

alpha_1=result[3]
alpha_2=result[0]
break_point=exp(-result[1])

print, 'Alpha_est_1 is', alpha_1
print, 'Alpha_est_2 is', alpha_2
print, 'Breaking point is',break_point, ' seconds'


sh=0.02
modes=emd(flux_preflare-trend_preflare,shi=sh, epsilon=1d-6, maxsiftings=1d8)
sz=size(modes)
trend_preflare=trend_preflare+modes[*,sz[2]-1]+modes[*,sz[2]-2]
emd_mode_preflare=modes[*,sz[2]-5]

Print, '*************************'
Print, 'Size of (modes)', sz
Print, 'Trend preflare', trend_preflare
Print, 'EMD MODE is', sz[2]-5


nmodep = 0.
pp = 0.
ep = 0.

modes/=sqrt(total((flux_preflare-mean(flux_preflare))^2))

for i=0,sz[2]-3 do begin
  p=(emd_period_energy(modes[*,i])).period
  e=total((modes[*,i])^2)
  pp = [pp, p]
  ep = [ep, e]
  nmodep = [nmodep,i+1]
endfor

trend_x=trend_preflare/sqrt(total((flux_preflare-mean(flux_preflare))^2));+modes[*,sz[2]-1]+modes[*,sz[2]-2];+modes[*,sz[2]-3]+modes[*,sz[2]-4]+modes[*,sz[2]-5]
p=(emd_period_energy(trend_x)).period
e=total((trend_x)^2)
pp = [pp, p]
ep = [ep, e]
nmodep = [nmodep,sz[2]-1]

period_preflare = pp[1:*]*cad
energy_preflare = ep[1:*]
nmode_preflare = nmodep[1:*]

print, 'periods',period_preflare



ind_fit_1=where(nmode_preflare ne sz[2] and period_preflare lt break_point)	;ALLOWS A COMPONANT OF NOISE TO BE NEGLECTED IF THERE ARE INSUFFIENT 
if (n_elements(ind_fit_1) lt 2) then begin			;MODES IN THIS REGION TO ALLOW A FITTING TO BE CARRIED OUT
  print,'Portion of the noise with alpha = ',$
  string(alpha_1,format='(F3.1)'), ' is negligible.'
  goto, label_1
endif


lin1 = [0d,alpha_1-1d]
print, 'lin1', lin1
weights1=replicate(1.0, N_ELEMENTS(period_preflare[ind_fit_1]))
fit=CURVEFIT( alog(period_preflare[ind_fit_1]), alog(energy_preflare[ind_fit_1]),fita=[1,0],$
weights1, lin1, FUNCTION_NAME= 'linfun',/noderivative)

confn_99_preflare1=emd_noise_conf(N_preflare, alpha_1, nsamples=100l, confidence_level = 0.01d)
confn_95_preflare1=emd_noise_conf(N_preflare, alpha_1, nsamples=100l, confidence_level = 0.05d)

meanenergy1=confn1_99.meanenergy[0]-alog(cad)*(alpha_1-1d)	;FIND MEAN VALUE ENERGIES USING ALPHA_1 AND 99% CI
norma_pre1=exp(meanenergy1)/exp(lin1[0])

print, meanenergy1, 'Mean energy1'
print, 'lin1', lin1
print, 'norma_pre1', norma_pre1

				
z_pre1 = alog(period_preflare)
y_pre1=(meanenergy1)+z_pre1*(alpha_1-1d)

;---------------------------LABEL 1----------------------------------------------------------------------------------------------
label_1:
sizeind=size(ind_fit_1)
ind_fit_2=  nmode_preflare[sizeind[1]:*]     ;from nmode_preflare, select all of the ones that aren't in ind_fit_1. That is, the ones with indices ranging from len(ind_fit_1) - starts 
if (n_elements(ind_fit_2) lt 2) then begin   ; from 0 - to the end of the array
  print,'Portion of the noise with alpha = ',string(alpha_2,format='(F3.1)'), ' is negligible.'
  goto, label_2
endif

weights2=replicate(1.0, N_ELEMENTS(period_preflare[ind_fit_2])) & lin2 = [0d,alpha_2-1d]
fit2=CURVEFIT( alog(period_preflare[ind_fit_2]), alog(energy_preflare[ind_fit_2]),fita=[1,0],$
weights2, lin2, FUNCTION_NAME= 'linfun',/noderivative)

confn_95_preflare2=emd_noise_conf(N_preflare, alpha_2, nsamples=100l, confidence_level = 0.05d)
confn_99_preflare2=emd_noise_conf(N_preflare, alpha_2, nsamples=100l, confidence_level = 0.01d)

meanenergy2=confn_99_preflare2.meanenergy[0]-alog(cad)*(alpha_2-1d)	;FIND MEAN VALUE ENERGIES USING ALPHA_1 AND 99% CI
norma_pre2=exp(meanenergy2)/exp(lin2[0])				;NORMALISING COEFFICIENT
z_pre2 = alog(period_preflare)
y_pre2=(meanenergy2)+z_pre2*(alpha_2-1d)







;*****************************************************************************************
;------------------------------------------FLARE------------------------------------------
;*****************************************************************************************

print,' **********************************'
print,' ************* FLARE **************'
print,' **********************************'

label_2:
flux=data.flux
tim=data.time
cad=tim[1]-tim[0]
flux_flare=flux[ns:nf]
tim_flare=tim[ns:nf]



N_flare=n_elements(flux_flare)

modes=emd(flux_flare-mean(flux_flare),shi=0.0033, epsilon=1d-6, maxsiftings=1d8)
sz=size(modes)
trend_flare=mean(flux_flare)+modes[*,sz[2]-1];+modes[*,sz[2]-2]

flux_flare=flux_flare-trend_flare

LNP_freq_finder = LNP_TEST(tim_flare, flux_flare-trend_flare, OFAC=1, WK1 = freq, WK2 = power)
ZZ=FFT(flux_flare,-1)
PowerZZ=(ABS(ZZ))^2d
PowerZZ= PowerZZ[1: N_elements(freq) -1]

p = [2.D, -1.D, 0.1D, 0.D] ;Initial guesses of model parameters
result = LNP_TEST(tim_flare, flux_flare-trend_flare, /DOUBLE, wk1=freq, wk2=pow, ofac=1)
npow = N_ELEMENTS(pow)
freqFLARE = ALOG(freq)
powFLARE = ALOG(pow)


Periodogram, tim_flare,flux_flare-trend_flare, fap=0.05, power_sign=FFT95FLARE			;FFT CONDIFENCE INTERVALS USING FALSE ALARM PROBABILITY
Periodogram, tim_flare,flux_flare-trend_flare, fap=0.01, power_sign=FFT99FLARE			;FOR 95 AND 99%



err = REPLICATE(1., npow)
;result = mpfitfun('broken_power_law', freq[1:*], pow[1:*], err[1:*], p, yfit=pow_fit, /quiet, /NAN)
pow_fitFLARE = CURVEFIT(freqFLARE[1:*], powFLARE[1:*], err[1:*], p, FUNCTION_NAME='broken_power_law',fita=[1,1,1,0],/noderivative)
result=p
print, 'Alpha_est_1 is', result[3]
print, 'Alpha_est_2 is', result[0]
print, 'Breaking point is',exp(-result[1]), ' seconds'

alpha_1=result[3]
alpha_2=result[0]
break_point=exp(-result[1])

modes=emd(flux_flare,shi=0.3, epsilon=1d-6, maxsiftings=1d8)
sz=size(modes)
modes=[[modes],[reform(trend_flare, n_elements(trend_flare),1)]]







nmodep = 0.
pp = 0.
ep = 0.

for i=0,sz[2]-1 do begin
  length=float(n_elements(flux_postflare));
  n_ext=float((size(extrema(modes[*,i])))[1])
  ;p=2d*length / n_ext
  p=(emd_period_energy(modes[*,i])).period
  e=total((modes[*,i])^2)
  pp = [pp, p]
  ep = [ep, e]
  nmodep = [nmodep,i+1]
endfor


period_flare = pp[1:*]*cad
energy_flare = ep[1:*]
nmode_flare = nmodep[1:*]

print, 'Periods:',period_flare

emd_mode=modes[*,2]

; window, 5 ,xsize=0.8*d[0],ysize=0.9*d[1]
; 
; !p.multi=[0,2,2]
; 
; plot,tim_flare, flux_flare+trend,thick=2,charsize=ch,/xsty,/ysty,title='!5 Original signal',$  ;ORIGINAL SIGNAL SCREEN
;   xtitle='Time',ytitle='Amplitude'
; oplot,tim_flare,trend,color=150,thick=2
; 
; plot,tim_flare,flux_flare,/ysty,/xsty,charsize=ch,title='Detrended signal',$        ; DETRENDED SIGNAL SCREEN
;   xtitle='Time',ytitle='Amplitude'
; oplot,tim_flare,emd_mode,color=250,thick=3
; 
; 
; pow_fit95=pow_fit + ALOG(FFT95)
; pow_fit99=pow_fit + ALOG(FFT99)
; 
; plot,EXP(-freq), pow,title='Power Spectrum', xtitle=$		;FOURIER POWER SPECTRUM
; 'Period [s]', ytitle='Log Fourier power [a.u.]', /xstyle, $
; /ystyle, charsize=2, /xlog
; oplot, EXP(-freq), pow_fit, linestyle=2, thick=3
; oplot, EXP(-freq), pow_fit95, color=250
; oplot, EXP(-freq), pow_fit99, color=150

; 
; plot,[0],[0],charsize=ch,/ylog,/xlog,xrange=[2d*cad,2d*N_flare*cad],$	;EMD SPECTRUM FOR ALPHA_0,1, VIA CREATION OF A BLANK SPECTRUM
; yrange=[0.5*min(energy_flare),2d*max(energy_flare)],xsty=1,ysty=1,title=$	
; 'EMD spectrum (FLARE)',YTICKFORMAT='(g6.0)',xtitle=$
; 'Period [s]',ytitle='Normalised modal energy [a.u.]'


ind_fit_3=where(nmode_flare ne sz[2] and period_flare lt break_point)	;ALLOWS A COMPONANT OF NOISE TO BE NEGLECTED IF THERE ARE INSUFFIENT 
if (n_elements(ind_fit_3) lt 2) then begin			;MODES IN THIS REGION TO ALLOW A FITTING TO BE CARRIED OUT
  print,'Portion of the noise with alpha = ',$
  string(alpha_1,format='(F3.1)'), ' is negligible.'
  goto, label_3
endif

;ANALYISING THE CONFIDENCE INTERVALS OF THE EMD MODES USING EMD_NOISE_CONF
;
; SHOWS HOW THE SIGNIFICANCE OF EMD MODES CAN BE ASSESSED WITH THE USE OF THE FUNCTION EMD_NOISE_CONF USING THE CALLING SEQUENCE
; confn=emd_noise_conf(N, alpha, nsamples = nsamples, confidence_level = confidence_level). HERE CONF IS SET TO 95 AND 99%


weights1=replicate(1.0, N_ELEMENTS(period_flare[ind_fit_1])) & lin1 = [0d,alpha_1-1d]
fit=CURVEFIT( alog(period_flare[ind_fit_1]), alog(energy_flare[ind_fit_1]),fita=[1,0],$
weights1, lin1, FUNCTION_NAME= 'linfun',/noderivative)

; confn_95_flare1=emd_noise_conf(N_flare, alpha_1, nsamples=100l, confidence_level = 0.05d)
; confn_99_flare1=emd_noise_conf(N_flare, alpha_1, nsamples=100l, confidence_level = 0.05d)


; oplot,period_flare,exp(y1)/norma1,color=200, thick=2.5

;---------------------------LABEL 3----------------------------------------------------------------------------------------------
label_3:

ind_fit_4=where(nmode_flare ne sz[2] and period_flare gt break_point)
if (n_elements(ind_fit_4) lt 2) then begin
  print,'Portion of the noise with alpha = ',string(alpha_2,format='(F3.1)'), ' is negligible.'
  goto, label_4
endif

weights2=replicate(1.0, N_ELEMENTS(period_flare[ind_fit_2])) & lin2 = [0d,alpha_2-1d]
fit2=CURVEFIT( alog(period_flare[ind_fit_2]), alog(energy_flare[ind_fit_2]),fita=[1,0],$
weights2, lin2, FUNCTION_NAME= 'linfun',/noderivative)

confn_95_flare1=emd_noise_conf(N_flare, alpha_1, nsamples=100l, confidence_level = 0.05d)
confn_99_flare1=emd_noise_conf(N_flare, alpha_1, nsamples=100l, confidence_level = 0.01d)

confn_95_flare2=emd_noise_conf(N_flare, alpha_2, nsamples=100l, confidence_level = 0.05d)
confn_99_flare2=emd_noise_conf(N_flare, alpha_2, nsamples=100l, confidence_level = 0.01d)

meanenergy1=confn_99_flare1.meanenergy[0]-alog(cad)*(alpha_1-1d)	;FIND MEAN VALUE ENERGIES USING ALPHA_1 AND 99% CI
;norma_flare1=exp(meanenergy1)/exp(lin1[0])				;NORMALISING COEFFICIENT





; oplot,(2d*N_flare*cad)/(confn1_95.k),confn1_95.up/norma1,color=20	;95% CONF. INTERVAL
; oplot,(2d*N_flare*cad)/(confn1_95.k),confn1_95.down/norma1,color=20
; 
; oplot,(2d*N_flare*cad)/(confn1_95.k),confn1_99.up/norma1,color=50	;99% CONF. INTERVAL
; oplot,(2d*N_flare*cad)/(confn1_95.k),confn1_99.down/norma1,color=50

z_flare1 = alog(period_flare)
y_flare1=(meanenergy1)+z_flare1*(alpha_1-1d)



meanenergy2=confn_99_flare2.meanenergy[0]-alog(cad)*(alpha_2-1d)	;FIND MEAN VALUE ENERGIES USING ALPHA_1 AND 99% CI
norma_flare2=exp(meanenergy2)/exp(lin2[0])				;NORMALISING COEFFICIENT


z_flare2 = alog(period_flare)
y_flare2=(meanenergy2)+z_flare2*(alpha_2-1d)
; oplot,period_flare,exp(y2)/norma2,color=200, thick=2.5

;---------------------------LABEL 4----------------------------------------------------------------------------------------------
label_4:
; oplot, period_flare, energy_flare, psym=8, color=250, symsize=1.1, thick=2	;PLOTS THE EMD MODES 
; !p.multi=[0,1,1]







;*****************************************************************************************
;---------------------------------------POSTFLARE------------------------------------------
;*****************************************************************************************

print,' **********************************'
print,' ********** POST FLARE ************'
print,' **********************************'


flux=data.flux
tim=data.time
cad=tim[1]-tim[0]
flux_postflare=flux[nf:299]
tim_postflare=tim[nf:299]



N_postflare=n_elements(flux_postflare)

modes=emd(flux_postflare-mean(flux_postflare),shi=0.003, epsilon=1d-6, maxsiftings=1d8)
sz=size(modes)
trend_postflare=mean(flux_postflare)+modes[*,sz[2]-1];+modes[*,sz[2]-2]

flux_postflare=flux_postflare-trend_postflare

LNP_freq_finder = LNP_TEST(tim_postflare, flux_postflare-trend_postflare, OFAC=1, WK1 = freq, WK2 = power)
ZZ=FFT(flux_postflare,-1)
PowerZZ=(ABS(ZZ))^2d
PowerZZ= PowerZZ[1: N_elements(freq) -1]

p = [2.D, -1.D, 0.1D, 0.D] ;Initial guesses of model parameters
result = LNP_TEST(tim_postflare, flux_postflare-trend_postflare, /DOUBLE, wk1=freq, wk2=pow, ofac=1)
npow = N_ELEMENTS(pow)
freqPOST = ALOG(freq)
powPOST = ALOG(pow)


Periodogram, tim_postflare,flux_postflare-trend_postflare, fap=0.05, power_sign=FFT95POST			;FFT CONDIFENCE INTERVALS USING FALSE ALARM PROBABILITY
Periodogram, tim_postflare,flux_postflare-trend_postflare, fap=0.01, power_sign=FFT99POST			;FOR 95 AND 99%



err = REPLICATE(1., npow)
;result = mpfitfun('broken_power_law', freq[1:*], pow[1:*], err[1:*], p, yfit=pow_fit, /quiet, /NAN)
pow_fitPOST = CURVEFIT(freqPOST[1:*], powPOST[1:*], err[1:*], p, FUNCTION_NAME='broken_power_law',fita=[1,1,1,0],/noderivative)
result=p
print, 'Alpha_est_1 is', result[3]
print, 'Alpha_est_2 is', result[0]
print, 'Breaking point is',exp(-result[1]), ' seconds'

alpha_1=result[3]
alpha_2=result[0]
break_point=exp(-result[1])

modes=emd(flux_postflare,shi=0.08, epsilon=1d-6, maxsiftings=1d8)
sz=size(modes)

modes=[[modes],[reform(trend_postflare, n_elements(trend_postflare),1)]]
sz=size(modes)

nmodep = 0.
pp = 0.
ep = 0.

for i=0,sz[2]-1 do begin
  length=float(n_elements(flux_postflare));
  n_ext=float((size(extrema(modes[*,i])))[1])
  ;p=2d*length / n_ext
  p=(emd_period_energy(modes[*,i])).period
  e=total((modes[*,i])^2)
  pp = [pp, p]
  ep = [ep, e]
  nmodep = [nmodep,i+1]
endfor

period_postflare = pp[1:*]*cad
energy_postflare = ep[1:*]
nmode_postflare = nmodep[1:*]

print, 'Periods:',period_postflare

emd_mode=modes[*,2]
emd_mode_postflare=modes[*,sz[2]-3]

; window, 6 ,xsize=0.8*d[0],ysize=0.9*d[1]
; 
; !p.multi=[0,2,2]
; 
; plot,tim_postflare, flux_postflare+trend,thick=2,charsize=ch,/xsty,/ysty,title='!5 Original signal',$  ;ORIGINAL SIGNAL SCREEN
;   xtitle='Time',ytitle='Amplitude'
; oplot,tim_postflare,trend,color=150,thick=2
; 
; plot,tim_postflare,flux_postflare,/ysty,/xsty,charsize=ch,title='Detrended signal',$        ; DETRENDED SIGNAL SCREEN
;   xtitle='Time',ytitle='Amplitude'
; oplot,tim_postflare,emd_mode,color=250,thick=3
; 
; 
; pow_fit95=pow_fit + ALOG(FFT95)
; pow_fit99=pow_fit + ALOG(FFT99)

; plot,EXP(-freq), pow,title='Power Spectrum', xtitle=$		;FOURIER POWER SPECTRUM
; 'Period [s]', ytitle='Log Fourier power [a.u.]', /xstyle, $
; /ystyle, charsize=2, /xlog
; oplot, EXP(-freq), pow_fit, linestyle=2, thick=3
; oplot, EXP(-freq), pow_fit95, color=250
; oplot, EXP(-freq), pow_fit99, color=150


; plot,[0],[0],charsize=ch,/ylog,/xlog,xrange=[2d*cad,2d*N_postflare*cad],$	;EMD SPECTRUM FOR ALPHA_0,1, VIA CREATION OF A BLANK SPECTRUM
; yrange=[0.5*min(energy_postflare),2d*max(energy_postflare)],xsty=1,ysty=1,title=$	
; 'EMD spectrum (POST)',YTICKFORMAT='(g6.0)',xtitle=$
; 'Period [s]',ytitle='Normalised modal energy [a.u.]'


ind_fit_5=where(nmode_postflare ne sz[2] and period_postflare lt break_point)	;ALLOWS A COMPONANT OF NOISE TO BE NEGLECTED IF THERE ARE INSUFFIENT 
if (n_elements(ind_fit_5) lt 2) then begin			;MODES IN THIS REGION TO ALLOW A FITTING TO BE CARRIED OUT
  print,'Portion of the noise with alpha = ',$
  string(alpha_1,format='(F3.1)'), ' is negligible.'
  goto, label_5
endif

;ANALYISING THE CONFIDENCE INTERVALS OF THE EMD MODES USING EMD_NOISE_CONF
;
; SHOWS HOW THE SIGNIFICANCE OF EMD MODES CAN BE ASSESSED WITH THE USE OF THE FUNCTION EMD_NOISE_CONF USING THE CALLING SEQUENCE
; confn=emd_noise_conf(N, alpha, nsamples = nsamples, confidence_level = confidence_level). HERE CONF IS SET TO 95 AND 99%


weights1=replicate(1.0, N_ELEMENTS(period_postflare[ind_fit_1])) & lin1 = [0d,alpha_1-1d]
fit=CURVEFIT( alog(period_postflare[ind_fit_1]), alog(energy_postflare[ind_fit_1]),fita=[1,0],$
weights1, lin1, FUNCTION_NAME= 'linfun',/noderivative)


; oplot,period_postflare,exp(y1)/norma1,color=200, thick=2.5

;---------------------------LABEL 5----------------------------------------------------------------------------------------------
label_5:

ind_fit_20=where(nmode_postflare ne sz[2] and period_postflare gt break_point)
if (n_elements(ind_fit_20) lt 2) then begin
  print,'Portion of the noise with alpha = ',string(alpha_2,format='(F3.1)'), ' is negligible.'
  goto, label_6
endif

weights2=replicate(1.0, N_ELEMENTS(period_postflare[ind_fit_20])) & lin2 = [0d,alpha_2-1d]
fit2=CURVEFIT( alog(period_postflare[ind_fit_20]), alog(energy_postflare[ind_fit_20]),fita=[1,0],$
weights2, lin2, FUNCTION_NAME= 'linfun',/noderivative)

confn_95_postflare2=emd_noise_conf(N_postflare, alpha_2, nsamples=100l, confidence_level = 0.05d)
confn_99_postflare2=emd_noise_conf(N_postflare, alpha_2, nsamples=100l, confidence_level = 0.01d)


confn_95_postflare1=emd_noise_conf(N_postflare, alpha_1, nsamples=100l, confidence_level = 0.05d)
confn_99_postflare1=emd_noise_conf(N_postflare, alpha_1, nsamples=100l, confidence_level = 0.01d)

meanenergy1=confn_99_postflare1.meanenergy[0]-alog(cad)*(alpha_1-1d)	;FIND MEAN VALUE ENERGIES USING ALPHA_1 AND 99% CI
;norma_post1=exp(meanenergy1)/exp(lin1[0])				;NORMALISING COEFFICIENT

; oplot,(2d*N_postflare*cad)/(confn1_95.k),confn1_95.up/norma1,color=20	;95% CONF. INTERVAL
; oplot,(2d*N_postflare*cad)/(confn1_95.k),confn1_95.down/norma1,color=20
; 
; oplot,(2d*N_postflare*cad)/(confn1_95.k),confn1_99.up/norma1,color=50	;99% CONF. INTERVAL
; oplot,(2d*N_postflare*cad)/(confn1_95.k),confn1_99.down/norma1,color=50

z_post1 = alog(period_postflare)
y_post1=(meanenergy1)+z_post1*(alpha_1-1d)

meanenergy2=confn_99_postflare2.meanenergy[0]-alog(cad)*(alpha_2-1d)	;FIND MEAN VALUE ENERGIES USING ALPHA_1 AND 99% CI
norma_post2=exp(meanenergy2)/exp(lin2[0])				;NORMALISING COEFFICIENT

; oplot,(2d*N_postflare*cad)/(confn2_95.k),confn2_95.up/norma2,color=180	;95% CONF. INTERVAL
; oplot,(2d*N_postflare*cad)/(confn2_95.k),confn2_95.down/norma2,color=180
; 
; oplot,(2d*N_postflare*cad)/(confn2_99.k),confn2_99.up/norma2,color=160	;99% CONF. INTERVAL
; oplot,(2d*N_postflare*cad)/(confn2_99.k),confn2_99.down/norma2,color=160

z_post2 = alog(period_postflare)
y_post2=(meanenergy2)+z_post2*(alpha_2-1d)
; oplot,period_postflare,exp(y2)/norma2,color=200, thick=2.5

;---------------------------LABEL 6----------------------------------------------------------------------------------------------
label_6:
; oplot, period_postflare, energy_postflare, psym=8, color=250, symsize=1.1, thick=2	;PLOTS THE EMD MODES 
; !p.multi=[0,1,1]






pow_fit95PRE=pow_fitPRE + ALOG(FFT95PRE)
pow_fit95FLARE=pow_fitFLARE + ALOG(FFT95FLARE)
pow_fit95POST=pow_fitPOST + ALOG(FFT95POST)

pow_fit99PRE=pow_fitPRE + ALOG(FFT99PRE)
pow_fit99FLARE=pow_fitFLARE + ALOG(FFT99FLARE)
pow_fit99POST=pow_fitPOST + ALOG(FFT99POST)


flux_preflare_detrended= flux_preflare-trend_preflare - emd_mode_preflare -5.
flux_flare_detrended= flux_flare - emd_mode -5.
flux_postflare_detrended= flux_postflare-trend_postflare - emd_mode_postflare -5.




M=2						;make m=/=1 to enable visualisation screens 
if (M eq 1) then begin
  goto, label_10
endif



;*****************************************************************************************
;-----------------------------------VISUALISATION-----------------------------------------
;*****************************************************************************************

d=get_screen_size() & ch=2.5
A = FIndGen(16) * (!PI*2/16.) & UserSym, cos(A), sin(A), /fill

window,1,xsize=0.5*d[0],ysize=0.7*d[1]

plot,tim, flux-25d,thick=1,charsize=ch,/xsty,/ysty,title=FLARETITLE, yrange=[-8,20],$
  xtitle='Elapsed time [s]',ytitle='Normalised amplitude [a.u.]' 
;oplot,tim_flare, trend-25d,color=250,thick=3
oplot,tim_preflare, trend_preflare-25d,color=250,thick=2.5
;oplot,tim_postflare, trend_postflare-25d,color=250,thick=3

oplot,tim_flare,flux_flare
;oplot,tim_flare,emd_mode,color=250,thick=3

oplot,tim_preflare,flux_preflare-trend_preflare
oplot,tim_preflare,emd_mode_preflare,color=200,thick=2.5

;oplot,tim_postflare,flux_postflare-trend_postflare
;oplot,tim_postflare,emd_mode_postflare,color=200,thick=3

oplot, tim_preflare, flux_preflare_detrended
;oplot, tim_flare, flux_flare_detrended
;oplot, tim_postflare, flux_postflare_detrended


;----------Preflare EMD and FFT---------------------------------------

window,2,xsize=0.6*d[0],ysize=0.7*d[1]
!p.multi=[0,2,3]

plot,[0],[0],charsize=ch,/ylog,/xlog,xrange=[2d*cad,2d*N_preflare*cad],$
  yrange=[1d-4,1d4],xsty=1,ysty=1,title='EMD Spectrum, Preflare',$
  xtitle='Period [s]',ytitle='EMD Energy [a.u.]';, position=[0.7,0.74,0.99,0.96],

if (n_elements(ind_fit_1) lt 2) then begin
  oplot, period_preflare[ind_fit_1], norma_pre2*energy_preflare(ind_fit_1),psym=8,color=50,symsize=1.1,thick=1
  goto, label_7
endif

oplot,(2d*N_preflare*cad)/(confn_99_preflare1.k),confn_99_preflare1.up,color=200,thick=1
oplot,(2d*N_preflare*cad)/(confn_99_preflare1.k),confn_99_preflare1.down,color=200,thick=1

oplot,(2d*N_preflare*cad)/(confn_95_preflare1.k),confn_95_preflare1.up,color=150,thick=1
oplot,(2d*N_preflare*cad)/(confn_95_preflare1.k),confn_95_preflare1.down,color=150,thick=1

oplot, period_preflare[ind_fit_1], norma_pre1*energy_preflare(ind_fit_1),psym=8,color=50,symsize=1.1,thick=1
oplot,period_preflare,exp(y_pre1),linesty=2

label_7:

oplot,(2d*N_preflare*cad)/(confn_99_preflare2.k),confn_99_preflare2.up,color=70,thick=1
oplot,(2d*N_preflare*cad)/(confn_99_preflare2.k),confn_99_preflare2.down,color=70,thick=1

oplot,(2d*N_preflare*cad)/(confn_95_preflare2.k),confn_95_preflare2.up,color=180,thick=1
oplot,(2d*N_preflare*cad)/(confn_95_preflare2.k),confn_95_preflare2.down,color=180,thick=1

;y=meanenergy_preflare+z*(alpha_pre-1d)

oplot,period_preflare,exp(y_pre2),linesty=2

oplot, period_preflare[ind_fit_2], norma_pre2*energy_preflare(ind_fit_2),psym=8,color=250,symsize=1.1,thick=1

plot,EXP(-freqPRE),powPRE, charsize=ch,  /xlog, /xstyle,xrange=[2d*cad,60], /ystyle,title='Fourier Spectrum, Preflare', xtitle='Period [s]', ytitle='Fourier Energy [a.u.]'
oplot, EXP(-freqPRE), pow_fitPRE, linestyle=1, thick=3
oplot, EXP(-freqPRE), pow_fit95PRE, color=150
oplot, EXP(-freqPRE), pow_fit99PRE, color=50


; ;----------Flare EMD and FFT---------------------------------------
; plot,[0],[0],charsize=ch,/ylog,/xlog,xrange=[2d*cad,2d*N_flare*cad],$
;   yrange=[1d-4,1d4],xsty=1,ysty=1,title='EMD Spectrum, Flare',$
;   xtitle='Period [s]',ytitle='EMD Energy [a.u.]';,position=[0.7,0.4,0.99,0.63]
; 
; ind_fit_4=where(nmode_postflare ne sz[2] and period_postflare gt break_point)
; if (n_elements(ind_fit_4) lt 2) then begin
;   goto, label_8
; endif
; 
; oplot,(2d*N_flare*cad)/(confn_99_flare1.k),confn_99_flare1.up,color=50,thick=1
; oplot,(2d*N_flare*cad)/(confn_99_flare1.k),confn_99_flare1.down,color=50,thick=1
; 
; oplot,(2d*N_flare*cad)/(confn_95_flare1.k),confn_95_flare1.up,color=150,thick=1
; oplot,(2d*N_flare*cad)/(confn_95_flare1.k),confn_95_flare1.down,color=150,thick=1
; oplot,period_preflare,exp(y_flare1),linesty=2
; 
; label_8:
; oplot,(2d*N_flare*cad)/(confn_99_flare2.k),confn_99_flare2.up,color=50,thick=1
; oplot,(2d*N_flare*cad)/(confn_99_flare2.k),confn_99_flare2.down,color=50,thick=1
; 
; oplot,(2d*N_flare*cad)/(confn_95_flare2.k),confn_95_flare2.up,color=150,thick=1
; oplot,(2d*N_flare*cad)/(confn_95_flare2.k),confn_95_flare2.down,color=150,thick=1
; 
; ; z = alog(period_flare)
; ; y=meanenergy_flare+z*(alpha_flare-1d)
; 
; oplot,period_preflare,exp(y_flare1),linesty=2
; oplot,period_preflare,exp(y_flare2),linesty=2
; 
; oplot, period_preflare[ind_fit_1], norma_flare1*energy_flare,psym=8,color=250,symsize=1.1,thick=2
; oplot, period_preflare[ind_fit_2], norma_flare2*energy_flare,psym=8,color=250,symsize=1.1,thick=2
; 
; plot,EXP(-freqFLARE), powFLARE, charsize=ch,  /xlog,xrange=[2d*cad,2d*N_postflare*cad], /xstyle, /ystyle,title='Fourier Spectrum, Flare', xtitle='Period [s]', ytitle='Fourier Energy [a.u.]'
; oplot, EXP(-freqFLARE), pow_fitFLARE, linestyle=2, thick=3
; oplot, EXP(-freqFLARE), pow_fit95FLARE, color=150
; oplot, EXP(-freqFLARE), pow_fit99FLARE, color=50
; 
; 
; ;----------Postflare EMD and FFT---------------------------------------
; 
; plot,[0],[0],charsize=ch,/ylog,/xlog,xrange=[2d*cad,2d*N_postflare*cad],$
;   yrange=[1d-4,1d4],xsty=1,ysty=1,title='EMD Spectrum, Postflare',$
;   xtitle='Period [s]',ytitle='EMD Energy [a.u.]';,position=[0.7,0.6,0.99,0.63]
; 
; ind_fit_5=where(nmode_postflare ne sz[2] and period_postflare gt break_point)
; if (n_elements(ind_fit_5) lt 2) then begin
;   goto, label_8
; endif
; 
; oplot,(2d*N_postflare*cad)/(confn_99_postflare1.k),confn_99_postflare1.up,color=50,thick=1
; oplot,(2d*N_postflare*cad)/(confn_99_postflare1.k),confn_99_postflare1.down,color=50,thick=1
; 
; oplot,(2d*N_postflare*cad)/(confn_95_postflare1.k),confn_95_postflare1.up,color=150,thick=1
; oplot,(2d*N_postflare*cad)/(confn_95_postflare1.k),confn_95_postflare1.down,color=150,thick=1
; oplot,period_preflare,exp(y_post1),linesty=2
; 
; label_9:
; 
; oplot,(2d*N_postflare*cad)/(confn_99_postflare2.k),confn_99_postflare2.up,color=50,thick=1
; oplot,(2d*N_postflare*cad)/(confn_99_postflare2.k),confn_99_postflare2.down,color=50,thick=1
; 
; oplot,(2d*N_postflare*cad)/(confn_95_postflare2.k),confn_95_postflare2.up,color=150,thick=1
; oplot,(2d*N_postflare*cad)/(confn_95_postflare2.k),confn_95_postflare2.down,color=150,thick=1
; 
; oplot,period_preflare,exp(y_post2),linesty=2
; 
; 
; oplot, period_postflare, norma_post1*energy_postflare,psym=8,color=250,symsize=1.1,thick=2
; 
; plot,EXP(-freqPOST), powPOST, charsize=ch, xrange=[2d*cad,2d*N_postflare*cad], /xlog, /xstyle, /ystyle,title='Fourier Spectrum, Postflare', xtitle='Period [s]', ytitle='Fourier Energy [a.u.]'
; oplot, EXP(-freqPOST), pow_fitPOST, linestyle=2, thick=3
; oplot, EXP(-freqPOST), pow_fit95POST, color=150
; oplot, EXP(-freqPOST), pow_fit99POST, color=50
; 
; !p.multi=[0,1,1]

label_10:







END