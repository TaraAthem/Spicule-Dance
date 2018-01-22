pro broken_power_law,x,p,model
  model = -(x LT p[1])*p[0]*x + $
    (x GT p[1])*(-(p[0]-p[3])*p[1] - p[3]*x) + p[2]
  return
end

pro linfun, x, a, f
  f = a[0] + a[1] * x
  return
end

;CD, '/Users/dmitriikolotkov/Dropbox/Tishtrya/flare45851'
 
A = FindGen(16) * (!PI*2/16.) & UserSym, cos(A), sin(A), /fill
device, decomposed=0 & loadct,39
d=get_screen_size()
flarename= 'flare45851.fits'
data = mrdfits(flarename, 1, header)

ns=60
nf=160

flux=data.flux
tim=data.time
cad=tim[1]-tim[0]


;for SF =0.001, 0.3001, 0.01 DO BEGIN   
;Print, SF, 'Should be 0.001 + n SF'
set_plot, 'PS'
;device, filename='Flare45851_' + StrTrim(SF,Format='(I2.2)')+'.eps', encapsulated=1,$		;DEFINES FILENAME AND IMAGE PROPERTIES
device, filename='45851_2.eps', encapsulated=1,$
/color,bits=24, xsize=75, ysize=60

;*****************************************************************************************
;---------------------------------------PREFLARE------------------------------------------
;*****************************************************************************************

print,' **********************************'
print,' ********** PREFLARE **************'
print,' **********************************'

flux_pre=flux[0:ns]
tim_pre=tim[0:ns]

x_pre=flux_pre

N_pre=n_elements(x_pre)

;EMD EXPANSION OF THE SYNTHETIC SIGNAL.
;
; VARIABLES: SHIFTFACTOR, MAXIMUM RELATIVE AMPLITUDE, AND MAXIMUM NUMBER OF SIFTINGS (SH, EPSILON, MAXSIFTINGS, RESPECTIVELY)
; CONTROLLING SENSITIVITY OF THE EXPANSION (FOR MORE DETAILS REFER TO THE MANUAL GIVEN IN EMD.PRO), EMD EXPANSION (MODES,
; SHOULD REMOVE A MEAN VALUE PRIOR TO THE EXPANSION), NUMBER OF DETECTED MODES (SZ).



;Real sf= 0.005

modes=emd(x_pre-mean(x_pre), shiftfactor=0.005, epsilon=1d-4, maxsiftings=1d5)
sz=size(modes)
trend_emd_pre=mean(x_pre)+modes[*,sz[2]-1]; +modes[*, sz[2]-2]   +modes[*, sz[2]-3] +modes[*, sz[2]-4]   ;FIND TREND FROM LAST MODE IN EMD SPECTRUM
x_pre=x_pre-trend_emd_pre             ;DETREND SIGNAL

p = [2.D, -1.D, 0.1D, 0.D]          ;INITIAL GUESSES OF MODEL PARAMETERS
result = LNP_TEST(tim_pre, x_pre, /DOUBLE, wk1=freq_pre, wk2=pow_pre, ofac=1) ;COMPUTES THE LOMB NORMALIZED PERIODOGRAM OF THE TIME AND
npow = N_ELEMENTS(pow_pre)            ;DETRENDED SINGAL SIGNAL TO FIND THE FREQUENCY ARRAY
freq_pre = ALOG(freq_pre)
pow_pre = ALOG(pow_pre)


Periodogram, tim_pre,x_pre, fap=0.05, power_sign=FFT95      ;FFT CONDIFENCE INTERVALS USING FALSE ALARM PROBABILITY
Periodogram, tim_pre,x_pre, fap=0.01, power_sign=FFT99      ;FOR 95 AND 99%

err = REPLICATE(1., npow)
pow_fit_pre = CURVEFIT(freq_pre[1:*], pow_pre[1:*], err[1:*], p,$   ;FITS A CURVE TO PRESCRIBED PARAMETRES
  FUNCTION_NAME='broken_power_law',fita=[1,1,1,0],/noderivative)
result=p
alpha_1=result[3]           ;SET AS CONSTANT GRADIENT IN A LOG-LOG PLOT
alpha_2=result[0]
break_point=exp(-result[1])

print, 'Alpha_est_1 is', alpha_1
print, 'Alpha_est_2 is', alpha_2
print, 'Breaking point is',break_point,' seconds'


;CALCULATE TOTAL ENERGY AND MAIN PERIOD OF EACH INTRINSIC MODE.
;
; SHOULD NORMALISE ALL EMD MODES FOR THE WHOLE SIGNAL'S TOTAL ENERGY TO BE UNITY
; VARIABLES: EMPTY ARRAYS TO BE FILLED IN WITH THE NUMBER, MAIN PERIOD, AND TOTAL ENERGY OF EACH EMPIRICAL MODE (NMODEP, PP, EP),
; LENGTH OF EACH MODE (LENGTH), NUMBER OF EXTREMA IN EACH MODE (N_EXT), MAIN PERIOD OF EACH MODE VIA THE WHOLE LENGTH AND NUMBER
; OF EXTREMA (P), TOTAL ENERGY OF EACH MODE AS A SUMMATION OF ALL DATA POINTS SQUARED (E).
; CHOOSE THE EMPIRICAL MODE (EMD_MODE) WITH THE PERIOD CLOSEST TO THAT OF THE SYNTHETIC MODE
; AND THE EMPIRICAL TREND (TREND_EMD) AS THE LAST (LONGEST-PERIOD) EMPIRICAL MODE (OR A COMBINATION OF A FEW EMPIRICAL MODES,
; IDEALLY IDENTICAL TO THE SYNTHETIC ONE).

modes=emd(x_pre-mean(x_pre), shiftfactor=0.07, epsilon=1d-4, maxsiftings=1d5)
modes=[[modes],[reform(trend_emd_pre,n_elements(trend_emd_pre),1)]]
sz=size(modes)

nmodep = 0.
pp = 0.
ep = 0.

for i=0,sz[2]-1 do begin
  length=float(n_elements(x_pre));
  n_ext=float((size(extrema(modes[*,i])))[1])
  ;p=2d*length / n_ext
  p=(emd_period_energy(modes[*,i])).period
  e=total((modes[*,i])^2)
  pp = [pp, p]
  ep = [ep, e]
  nmodep = [nmodep,i+1]
endfor

period_pre = pp[1:*]*cad
energy_pre = ep[1:*]
nmode = nmodep[1:*]

print, 'The periods of the extracted IMFs:',period_pre

emd_mode_pre=modes[*,2] ;+ modes[*,3]


pow_fit95_pre=pow_fit_pre + ALOG(FFT95)
pow_fit99_pre=pow_fit_pre + ALOG(FFT99)

y1_pre=0. & y2_pre=0.

ind_fit_1=where(nmode ne sz[2] and period_pre lt break_point) ;ALLOWS A COMPONANT OF NOISE TO BE NEGLECTED IF THERE ARE INSUFFIENT
if (n_elements(ind_fit_1) lt 2) then begin      ;MODES IN THIS REGION TO ALLOW A FITTING TO BE CARRIED OUT
  print,'PORTION OF THE NOISE WITH ALPHA = ',$
    string(alpha_1,format='(F10.3)'), ' IS NEGLIGIBLE.'
  goto, label_1
endif

weights1=replicate(1.0, N_ELEMENTS(period_pre[ind_fit_1])) & lin1 = [0d,alpha_1-1d]
fit=CURVEFIT( alog(period_pre[ind_fit_1]), alog(energy_pre[ind_fit_1]),fita=[1,0],$
  weights1, lin1, FUNCTION_NAME= 'linfun',/noderivative)

confn1_95_pre=emd_noise_conf(N_pre, alpha_1, nsamples=500l, confidence_level = 0.05d)
confn1_99_pre=emd_noise_conf(N_pre, alpha_1, nsamples=500l, confidence_level = 0.01d)

meanenergy1=confn1_99_pre.meanenergy[0]-alog(cad)*(alpha_1-1d)  ;FIND MEAN VALUE ENERGIES USING ALPHA_1 AND 99% CI
norma1_pre=exp(meanenergy1)/exp(lin1[0])        ;NORMALISING COEFFICIENT

z = alog(period_pre)
y1_pre=(meanenergy1)+z*(alpha_1-1d)

;---------------------------LABEL 1----------------------------------------------------------------------------------------------
label_1:
ind_fit_2=where(nmode ne sz[2] and period_pre gt break_point)
if (n_elements(ind_fit_2) lt 2) then begin
  print,'Portion of the noise with alpha = ',string(alpha_2,format='(F10.3)'), ' is negligible.'
  goto, label_2
endif

weights2=replicate(1.0, N_ELEMENTS(period_pre[ind_fit_2])) & lin2 = [0d,alpha_2-1d]
fit2=CURVEFIT( alog(period_pre[ind_fit_2]), alog(energy_pre[ind_fit_2]),fita=[1,0],$
  weights2, lin2, FUNCTION_NAME= 'linfun',/noderivative)

confn2_95_pre=emd_noise_conf(N_pre, alpha_2, nsamples=500l, confidence_level = 0.05d)
confn2_99_pre=emd_noise_conf(N_pre, alpha_2, nsamples=500l, confidence_level = 0.01d)

meanenergy2=confn2_99_pre.meanenergy[0]-alog(cad)*(alpha_2-1d)  ;FIND MEAN VALUE ENERGIES USING ALPHA_1 AND 99% CI
norma2_pre=exp(meanenergy2)/exp(lin2[0])        ;NORMALISING COEFFICIENT

z = alog(period_pre)
y2_pre=(meanenergy2)+z*(alpha_2-1d)
label_2:




;*****************************************************************************************
;------------------------------------------FLARE------------------------------------------
;*****************************************************************************************

print,' **********************************'
print,' ************* FLARE **************'
print,' **********************************'

flux_flare=flux[ns:nf]
tim_flare=tim[ns:nf]

x_flare=flux_flare

N_flare=n_elements(x_flare)

;EMD EXPANSION OF THE SYNTHETIC SIGNAL.
;
; VARIABLES: SHIFTFACTOR, MAXIMUM RELATIVE AMPLITUDE, AND MAXIMUM NUMBER OF SIFTINGS (SH, EPSILON, MAXSIFTINGS, RESPECTIVELY)
; CONTROLLING SENSITIVITY OF THE EXPANSION (FOR MORE DETAILS REFER TO THE MANUAL GIVEN IN EMD.PRO), EMD EXPANSION (MODES,
; SHOULD REMOVE A MEAN VALUE PRIOR TO THE EXPANSION), NUMBER OF DETECTED MODES (SZ).


;SF= 0.0075
;0.131 and 0.141
modes=emd(x_flare-mean(x_flare), shiftfactor=0.135, epsilon=1d-4, maxsiftings=1d5)
sz=size(modes)
trend_emd_flare=mean(x_flare)+modes[*,sz[2]-1]+modes[*,sz[2]-2]; +modes[*,sz[2]-3]  +modes[*,sz[2]-4]  ;FIND TREND FROM LAST MODE IN EMD SPECTRUM
x_flare=x_flare-trend_emd_flare             ;DETREND SIGNAL

p = [2.D, -1.D, 0.1D, 0.D]          ;INITIAL GUESSES OF MODEL PARAMETERS
result = LNP_TEST(tim_flare, x_flare, /DOUBLE, wk1=freq_flare, wk2=pow_flare, ofac=1) ;COMPUTES THE LOMB NORMALIZED PERIODOGRAM OF THE TIME AND
npow = N_ELEMENTS(pow_flare)            ;DETRENDED SINGAL SIGNAL TO FIND THE FREQUENCY ARRAY
freq_flare = ALOG(freq_flare)
pow_flare = ALOG(pow_flare)


Periodogram, tim_flare,x_flare, fap=0.05, power_sign=FFT95      ;FFT CONDIFENCE INTERVALS USING FALSE ALARM PROBABILITY
Periodogram, tim_flare,x_flare, fap=0.01, power_sign=FFT99      ;FOR 95 AND 99%

err = REPLICATE(1., npow)
pow_fit_flare = CURVEFIT(freq_flare[1:*], pow_flare[1:*], err[1:*], p,$   ;FITS A CURVE TO PRESCRIBED PARAMETRES
  FUNCTION_NAME='broken_power_law',fita=[1,1,1,0],/noderivative)
result=p
alpha_1=result[3]           ;SET AS CONSTANT GRADIENT IN A LOG-LOG PLOT
alpha_2=result[0]
break_point=exp(-result[1])

print, 'Alpha_est_1 is', alpha_1
print, 'Alpha_est_2 is', alpha_2
print, 'Breaking point is',break_point,' seconds'


;CALCULATE TOTAL ENERGY AND MAIN PERIOD OF EACH INTRINSIC MODE.
;
; SHOULD NORMALISE ALL EMD MODES FOR THE WHOLE SIGNAL'S TOTAL ENERGY TO BE UNITY
; VARIABLES: EMPTY ARRAYS TO BE FILLED IN WITH THE NUMBER, MAIN PERIOD, AND TOTAL ENERGY OF EACH EMPIRICAL MODE (NMODEP, PP, EP),
; LENGTH OF EACH MODE (LENGTH), NUMBER OF EXTREMA IN EACH MODE (N_EXT), MAIN PERIOD OF EACH MODE VIA THE WHOLE LENGTH AND NUMBER
; OF EXTREMA (P), TOTAL ENERGY OF EACH MODE AS A SUMMATION OF ALL DATA POINTS SQUARED (E).
; CHOOSE THE EMPIRICAL MODE (EMD_MODE) WITH THE PERIOD CLOSEST TO THAT OF THE SYNTHETIC MODE
; AND THE EMPIRICAL TREND (TREND_EMD) AS THE LAST (LONGEST-PERIOD) EMPIRICAL MODE (OR A COMBINATION OF A FEW EMPIRICAL MODES,
; IDEALLY IDENTICAL TO THE SYNTHETIC ONE).

modes=emd(x_flare-mean(x_flare), shiftfactor=0.003, epsilon=1d-4, maxsiftings=1d5)
modes=[[modes],[reform(trend_emd_flare,n_elements(trend_emd_flare),1)]]
sz=size(modes)

nmodep = 0.
pp = 0.
ep = 0.

for i=0,sz[2]-1 do begin
  length=float(n_elements(x_flare));
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
nmode = nmodep[1:*]

print, 'The periods of the extracted IMFs:',period_flare

emd_mode_flare=modes[*,3];+ modes[*,3]


pow_fit95_flare=pow_fit_flare + ALOG(FFT95)
pow_fit99_flare=pow_fit_flare + ALOG(FFT99)

y1_flare=0. & y2_flare=0.

ind_fit_1=where(nmode ne sz[2] and period_flare lt break_point) ;ALLOWS A COMPONANT OF NOISE TO BE NEGLECTED IF THERE ARE INSUFFIENT
if (n_elements(ind_fit_1) lt 2) then begin      ;MODES IN THIS REGION TO ALLOW A FITTING TO BE CARRIED OUT
  print,'PORTION OF THE NOISE WITH ALPHA = ',$
    string(alpha_1,format='(F10.3)'), ' IS NEGLIGIBLE.'
  goto, label_3
endif

weights1=replicate(1.0, N_ELEMENTS(period_flare[ind_fit_1])) & lin1 = [0d,alpha_1-1d]
fit=CURVEFIT( alog(period_flare[ind_fit_1]), alog(energy_flare[ind_fit_1]),fita=[1,0],$
  weights1, lin1, FUNCTION_NAME= 'linfun',/noderivative)

confn1_95_flare=emd_noise_conf(N_flare, alpha_1, nsamples=500l, confidence_level = 0.05d)
confn1_99_flare=emd_noise_conf(N_flare, alpha_1, nsamples=500l, confidence_level = 0.01d)

meanenergy1=confn1_99_flare.meanenergy[0]-alog(cad)*(alpha_1-1d)  ;FIND MEAN VALUE ENERGIES USING ALPHA_1 AND 99% CI
norma1_flare=exp(meanenergy1)/exp(lin1[0])        ;NORMALISING COEFFICIENT

z = alog(period_flare)
y1_flare=(meanenergy1)+z*(alpha_1-1d)

;---------------------------LABEL 3----------------------------------------------------------------------------------------------
label_3:
ind_fit_2=where(nmode ne sz[2] and period_flare gt break_point)
if (n_elements(ind_fit_2) lt 2) then begin
  print,'Portion of the noise with alpha = ',string(alpha_2,format='(F10.3)'), ' is negligible.'
  goto, label_4
endif

weights2=replicate(1.0, N_ELEMENTS(period_flare[ind_fit_2])) & lin2 = [0d,alpha_2-1d]
fit2=CURVEFIT( alog(period_flare[ind_fit_2]), alog(energy_flare[ind_fit_2]),fita=[1,0],$
  weights2, lin2, FUNCTION_NAME= 'linfun',/noderivative)

confn2_95_flare=emd_noise_conf(N_flare, alpha_2, nsamples=500l, confidence_level = 0.05d)
confn2_99_flare=emd_noise_conf(N_flare, alpha_2, nsamples=500l, confidence_level = 0.01d)

meanenergy2=confn2_99_flare.meanenergy[0]-alog(cad)*(alpha_2-1d)  ;FIND MEAN VALUE ENERGIES USING ALPHA_1 AND 99% CI
norma2_flare=exp(meanenergy2)/exp(lin2[0])        ;NORMALISING COEFFICIENT

z = alog(period_flare)
y2_flare=(meanenergy2)+z*(alpha_2-1d)

label_4:



;*****************************************************************************************
;---------------------------------------POSTFLARE------------------------------------------
;*****************************************************************************************

print,' **********************************'
print,' ********** POSTFLARE **************'
print,' **********************************'
flux_post=flux[nf:299]
tim_post=tim[nf:299]

x_post=flux_post

N_post=n_elements(x_post)

;EMD EXPANSION OF THE SYNTHETIC SIGNAL.
;
; VARIABLES: SHIFTFACTOR, MAXIMUM RELATIVE AMPLITUDE, AND MAXIMUM NUMBER OF SIFTINGS (SH, EPSILON, MAXSIFTINGS, RESPECTIVELY)
; CONTROLLING SENSITIVITY OF THE EXPANSION (FOR MORE DETAILS REFER TO THE MANUAL GIVEN IN EMD.PRO), EMD EXPANSION (MODES,
; SHOULD REMOVE A MEAN VALUE PRIOR TO THE EXPANSION), NUMBER OF DETECTED MODES (SZ).



;SF= 0.008
modes=emd(x_post-mean(x_post), shiftfactor=0.008, epsilon=1d-4, maxsiftings=1d5)
sz=size(modes)
trend_emd_post=mean(x_post)+modes[*,sz[2]-1] +modes[*, sz[2]-2]; +modes[*, sz[2]-3]    ;FIND TREND FROM LAST MODE IN EMD SPECTRUM
x_post=x_post-trend_emd_post             ;DETREND SIGNAL

p = [2.D, -1.D, 0.1D, 0.D]          ;INITIAL GUESSES OF MODEL PARAMETERS
result = LNP_TEST(tim_post, x_post, /DOUBLE, wk1=freq_post, wk2=pow_post, ofac=1) ;COMPUTES THE LOMB NORMALIZED PERIODOGRAM OF THE TIME AND
npow = N_ELEMENTS(pow_post)            ;DETRENDED SINGAL SIGNAL TO FIND THE FREQUENCY ARRAY
freq_post = ALOG(freq_post)
pow_post = ALOG(pow_post)


Periodogram, tim_post,x_post, fap=0.05, power_sign=FFT95      ;FFT CONDIFENCE INTERVALS USING FALSE ALARM PROBABILITY
Periodogram, tim_post,x_post, fap=0.01, power_sign=FFT99      ;FOR 95 AND 99%

err = REPLICATE(1., npow)
pow_fit_post = CURVEFIT(freq_post[1:*], pow_post[1:*], err[1:*], p,$   ;FITS A CURVE TO PRESCRIBED PARAMETRES
  FUNCTION_NAME='broken_power_law',fita=[1,1,1,0],/noderivative)
result=p
alpha_1=result[3]           ;SET AS CONSTANT GRADIENT IN A LOG-LOG PLOT
alpha_2=result[0]
break_point=exp(-result[1])

print, 'Alpha_est_1 is', alpha_1
print, 'Alpha_est_2 is', alpha_2
print, 'Breaking point is',break_point,' seconds'


;CALCULATE TOTAL ENERGY AND MAIN PERIOD OF EACH INTRINSIC MODE.
;
; SHOULD NORMALISE ALL EMD MODES FOR THE WHOLE SIGNAL'S TOTAL ENERGY TO BE UNITY
; VARIABLES: EMPTY ARRAYS TO BE FILLED IN WITH THE NUMBER, MAIN PERIOD, AND TOTAL ENERGY OF EACH EMPIRICAL MODE (NMODEP, PP, EP),
; LENGTH OF EACH MODE (LENGTH), NUMBER OF EXTREMA IN EACH MODE (N_EXT), MAIN PERIOD OF EACH MODE VIA THE WHOLE LENGTH AND NUMBER
; OF EXTREMA (P), TOTAL ENERGY OF EACH MODE AS A SUMMATION OF ALL DATA POINTS SQUARED (E).
; CHOOSE THE EMPIRICAL MODE (EMD_MODE) WITH THE PERIOD CLOSEST TO THAT OF THE SYNTHETIC MODE
; AND THE EMPIRICAL TREND (TREND_EMD) AS THE LAST (LONGEST-PERIOD) EMPIRICAL MODE (OR A COMBINATION OF A FEW EMPIRICAL MODES,
; IDEALLY IDENTICAL TO THE SYNTHETIC ONE).

modes=emd(x_post-mean(x_post), shiftfactor=0.010, epsilon=1d-4, maxsiftings=1d5)
modes=[[modes],[reform(trend_emd_post,n_elements(trend_emd_post),1)]]
sz=size(modes)

nmodep = 0.
pp = 0.
ep = 0.

for i=0,sz[2]-1 do begin
  length=float(n_elements(x_post));
  n_ext=float((size(extrema(modes[*,i])))[1])
  ;p=2d*length / n_ext
  p=(emd_period_energy(modes[*,i])).period
  e=total((modes[*,i])^2)
  pp = [pp, p]
  ep = [ep, e]
  nmodep = [nmodep,i+1]
endfor

period_post = pp[1:*]*cad
energy_post = ep[1:*]
nmode = nmodep[1:*]

print, 'The periods of the extracted IMFs:',period_post

emd_mode_post=modes[*,3];+ modes[*,3]


pow_fit95_post=pow_fit_post + ALOG(FFT95)
pow_fit99_post=pow_fit_post + ALOG(FFT99)
y1_post=0. & y2_post=0.

ind_fit_1=where(nmode ne sz[2] and period_post lt break_point) ;ALLOWS A COMPONANT OF NOISE TO BE NEGLECTED IF THERE ARE INSUFFIENT
if (n_elements(ind_fit_1) lt 2) then begin      ;MODES IN THIS REGION TO ALLOW A FITTING TO BE CARRIED OUT
  print,'PORTION OF THE NOISE WITH ALPHA = ',$
    string(alpha_1,format='(F10.3)'), ' IS NEGLIGIBLE.'
  goto, label_5
endif

weights1=replicate(1.0, N_ELEMENTS(period_post[ind_fit_1])) & lin1 = [0d,alpha_1-1d]
fit=CURVEFIT( alog(period_post[ind_fit_1]), alog(energy_post[ind_fit_1]),fita=[1,0],$
  weights1, lin1, FUNCTION_NAME= 'linfun',/noderivative)

confn1_95_post=emd_noise_conf(N_post, alpha_1, nsamples=500l, confidence_level = 0.05d)
confn1_99_post=emd_noise_conf(N_post, alpha_1, nsamples=500l, confidence_level = 0.01d)

meanenergy1=confn1_99_post.meanenergy[0]-alog(cad)*(alpha_1-1d)  ;FIND MEAN VALUE ENERGIES USING ALPHA_1 AND 99% CI
norma1_post=exp(meanenergy1)/exp(lin1[0])        ;NORMALISING COEFFICIENT

z = alog(period_post)
y1_post=(meanenergy1)+z*(alpha_1-1d)

;---------------------------LABEL 5----------------------------------------------------------------------------------------------
label_5:
ind_fit_2=where(nmode ne sz[2] and period_post gt break_point)
if (n_elements(ind_fit_2) lt 2) then begin
  print,'Portion of the noise with alpha = ',string(alpha_2,format='(F10.3)'), ' is negligible.'
  goto, label_6
endif

weights2=replicate(1.0, N_ELEMENTS(period_post[ind_fit_2])) & lin2 = [0d,alpha_2-1d]
fit2=CURVEFIT( alog(period_post[ind_fit_2]), alog(energy_post[ind_fit_2]),fita=[1,0],$
  weights2, lin2, FUNCTION_NAME= 'linfun',/noderivative)

confn2_95_post=emd_noise_conf(N_post, alpha_2, nsamples=500l, confidence_level = 0.05d)
confn2_99_post=emd_noise_conf(N_post, alpha_2, nsamples=500l, confidence_level = 0.01d)

meanenergy2=confn2_99_post.meanenergy[0]-alog(cad)*(alpha_2-1d)  ;FIND MEAN VALUE ENERGIES USING ALPHA_1 AND 99% CI
norma2_post=exp(meanenergy2)/exp(lin2[0])        ;NORMALISING COEFFICIENT

z = alog(period_post)
y2_post=(meanenergy2)+z*(alpha_2-1d)
label_6:


;*****************************************************************************************
;-----------------------------------VISUALISATION-----------------------------------------
;*****************************************************************************************

;window, 7 ,xsize=0.8*d[0],ysize=0.9*d[1]
!p.multi=[0,2,3]
d=get_screen_size() & ch=3.5
A = FindGen(16) * (!PI*2/16.) & UserSym, cos(A), sin(A), /fill
device, decomposed=0 & loadct,39

plot,EXP(-freq_pre), pow_pre,title='Fourier Spectrum - Preflare', xtitle=$   ;FOURIER POWER SPECTRUM
  'Period [s]', ytitle='Log Fourier power [a.u.]', /xstyle, $
  /ystyle, charsize=ch, /xlog, xrange=[2,100], thick=2,yrange=[-3.5,2.5]
oplot, EXP(-freq_pre), pow_fit_pre, linestyle=2, thick=3
oplot, EXP(-freq_pre), pow_fit95_pre, color=150, thick=2
oplot, EXP(-freq_pre), pow_fit99_pre, color=250, thick=2

MAXFOURIERPRE=Max(pow_pre,I_pre)
Print, 'The max preflare Fourier value corresponds to a period of', exp(-freq_pre[I_pre])

plot,[0],[0],charsize=ch,/xlog,xrange=[2d*cad,2d*N_pre*cad],$ ;EMD SPECTRUM FOR ALPHA_0,1, VIA CREATION OF A BLANK SPECTRUM
  yrange=[min(alog(energy_pre))-2d*abs(min(alog(energy_pre))),max(alog(energy_pre))+0.1*abs(max(alog(energy_pre)))],xsty=1,ysty=1,title=$
  'EMD Spectrum - Preflare',xtitle='Period [s]',ytitle='Log EMD modal energy [a.u.]'


if (n_elements(y1_pre) ne 1 and n_elements(y2_pre) ne 1) then begin

  oplot,(2d*N_pre*cad)/(confn2_99_pre.k),alog(confn2_99_pre.up/norma2_pre)>alog(confn1_99_pre.up/norma1_pre),color=250  ;99% CONF. INTERVAL
  oplot,(2d*N_pre*cad)/(confn2_95_pre.k),alog(confn2_95_pre.up/norma2_pre)>alog(confn1_95_pre.up/norma1_pre),color=150  ;99% CONF. INTERVAL
  oplot,period_pre,alog(exp(y2_pre)/norma2_pre)>alog(exp(y1_pre)/norma1_pre), thick=2, linesty=2

endif else begin

  if (n_elements(y1_pre) ne 1) then begin
    oplot,(2d*N_pre*cad)/(confn1_99_pre.k),alog(confn1_99_pre.up/norma1_pre),color=250  ;99% CONF. INTERVAL
    oplot,(2d*N_pre*cad)/(confn1_95_pre.k),alog(confn1_95_pre.up/norma1_pre),color=150  ;95% CONF. INTERVAL
    oplot,period_pre,alog(exp(y1_pre)/norma1_pre), thick=2, linesty=2
  endif

  if (n_elements(y2_pre) ne 1) then begin
    oplot,(2d*N_pre*cad)/(confn2_99_pre.k),alog(confn2_99_pre.up/norma2_pre),color=250  ;99% CONF. INTERVAL
    oplot,(2d*N_pre*cad)/(confn2_95_pre.k),alog(confn2_95_pre.up/norma2_pre),color=150  ;95% CONF. INTERVAL
    oplot,period_pre,alog(exp(y2_pre)/norma2_pre), thick=2, linesty=2
  endif

endelse

oplot, period_pre, alog(energy_pre), psym=8, color=50, symsize=1.1, thick=2  ;PLOTS THE EMD MODES


plot,EXP(-freq_flare), pow_flare,title='Fourier Spectrum - Flare', xtitle=$   ;FOURIER POWER SPECTRUM
  'Period [s]', ytitle='Log Fourier power [a.u.]', /xstyle, $
  /ystyle, charsize=ch, /xlog, xrange=[2,100], thick=2,yrange=[-4.,2.1]
oplot, EXP(-freq_flare), pow_fit_flare, linestyle=2, thick=3
oplot, EXP(-freq_flare), pow_fit95_flare, color=150, thick=2
oplot, EXP(-freq_flare), pow_fit99_flare, color=250, thick=2

MAXFOURIERFLARE=Max(pow_flare,I_flare)
Print, 'The max flare Fourier value corresponds to a period of',exp(-freq_flare[I_flare])

plot,[0],[0],charsize=ch,/xlog,xrange=[2d*cad,2d*N_flare*cad],$ ;EMD SPECTRUM FOR ALPHA_0,1, VIA CREATION OF A BLANK SPECTRUM
  yrange=[min(alog(energy_flare))-2d*abs(min(alog(energy_flare))),max(alog(energy_flare))+0.1*abs(max(alog(energy_flare)))],xsty=1,ysty=1,title=$
  'EMD Spectrum - Flare',xtitle='Period [s]',ytitle='Log EMD modal energy [a.u.]'


if (n_elements(y1_flare) ne 1 and n_elements(y2_flare) ne 1) then begin

  oplot,(2d*N_flare*cad)/(confn2_99_flare.k),alog(confn2_99_flare.up/norma2_flare)>alog(confn1_99_flare.up/norma1_flare),color=250  ;99% CONF. INTERVAL
  oplot,(2d*N_flare*cad)/(confn2_95_flare.k),alog(confn2_95_flare.up/norma2_flare)>alog(confn1_95_flare.up/norma1_flare),color=150  ;99% CONF. INTERVAL
  oplot,period_flare,alog(exp(y2_flare)/norma2_flare)>alog(exp(y1_flare)/norma1_flare), thick=2, linesty=2

endif else begin

  if (n_elements(y1_flare) ne 1) then begin
    oplot,(2d*N_flare*cad)/(confn1_99_flare.k),alog(confn1_99_flare.up/norma1_flare),color=250  ;99% CONF. INTERVAL
    oplot,(2d*N_flare*cad)/(confn1_95_flare.k),alog(confn1_95_flare.up/norma1_flare),color=150  ;95% CONF. INTERVAL
    oplot,period_flare,alog(exp(y1_flare)/norma1_flare), thick=2, linesty=2
  endif

  if (n_elements(y2_flare) ne 1) then begin
    oplot,(2d*N_flare*cad)/(confn2_99_flare.k),alog(confn2_99_flare.up/norma2_flare),color=250  ;99% CONF. INTERVAL
    oplot,(2d*N_flare*cad)/(confn2_95_flare.k),alog(confn2_95_flare.up/norma2_flare),color=150  ;95% CONF. INTERVAL
    oplot,period_flare,alog(exp(y2_flare)/norma2_flare), thick=2, linesty=2
  endif

endelse

oplot, period_flare, alog(energy_flare), psym=8, color=50, symsize=1.1, thick=2  ;PLOTS THE EMD MODES

plot,EXP(-freq_post), pow_post,title='Fourier Spectrum - Postflare', xtitle=$   ;FOURIER POWER SPECTRUM
  'Period [s]', ytitle='Log Fourier power [a.u.]', /xstyle, $
  /ystyle, charsize=ch, /xlog, xrange=[2,100], yrange=[-4, 2.7], thick=2
oplot, EXP(-freq_post), pow_fit_post, linestyle=2, thick=3
oplot, EXP(-freq_post), pow_fit95_post, color=150, thick=2
oplot, EXP(-freq_post), pow_fit99_post, color=250, thick=2

MAXFOURIERPOST=Max(pow_post,I_post)
Print, 'The max postflare Fourier value corresponds to a period of', exp(-freq_post[I_post])


plot,[0],[0],charsize=ch,/xlog,xrange=[2d*cad,2d*N_post*cad],$ ;EMD SPECTRUM FOR ALPHA_0,1, VIA CREATION OF A BLANK SPECTRUM
  yrange=[min(alog(energy_post))-2d*abs(min(alog(energy_post))),max(alog(energy_post))+0.1*abs(max(alog(energy_post)))],xsty=1,ysty=1,title=$
  'EMD Spectrum - Postflare',xtitle='Period [s]',ytitle='Log EMD modal energy [a.u.]'


if (n_elements(y1_post) ne 1 and n_elements(y2_post) ne 1) then begin

  oplot,(2d*N_post*cad)/(confn2_99_post.k),alog(confn2_99_post.up/norma2_post)>alog(confn1_99_post.up/norma1_post),color=250  ;99% CONF. INTERVAL
  oplot,(2d*N_post*cad)/(confn2_95_post.k),alog(confn2_95_post.up/norma2_post)>alog(confn1_95_post.up/norma1_post),color=150  ;99% CONF. INTERVAL
  oplot,period_post,alog(exp(y2_post)/norma2_post)>alog(exp(y1_post)/norma1_post), thick=2, linesty=2

endif else begin

  if (n_elements(y1_post) ne 1) then begin
    oplot,(2d*N_post*cad)/(confn1_99_post.k),alog(confn1_99_post.up/norma1_post),color=250  ;99% CONF. INTERVAL
    oplot,(2d*N_post*cad)/(confn1_95_post.k),alog(confn1_95_post.up/norma1_post),color=150  ;95% CONF. INTERVAL
    oplot,period_post,alog(exp(y1_post)/norma1_post), thick=2, linesty=2
  endif

  if (n_elements(y2_post) ne 1) then begin
    oplot,(2d*N_post*cad)/(confn2_99_post.k),alog(confn2_99_post.up/norma2_post),color=250  ;99% CONF. INTERVAL
    oplot,(2d*N_post*cad)/(confn2_95_post.k),alog(confn2_95_post.up/norma2_post),color=150  ;95% CONF. INTERVAL
    oplot,period_post,alog(exp(y2_post)/norma2_post), thick=2, linesty=2
  endif

endelse

oplot, period_post, alog(energy_post), psym=8, color=50, symsize=1.1, thick=2  ;PLOTS THE EMD MODES

!p.multi=[0,1,1]



;window, 6 ,xsize=0.8*d[0],ysize=0.9*d[1]

plot,[0],[0],charsize=ch,/xsty,/ysty,yrange=[-10,25],xrange=[0,300],$   ;ORIGINAL SIGNAL WINDOW
  title='Profile of Flare 45851', xtitle='Elapsed Time [s]',ytitle='Amplitude'
oplot,tim_flare,trend_emd_flare -50. ,color=130,thick=4
oplot,tim_flare, x_flare+trend_emd_flare -50.,thick=3

oplot,tim_flare,x_flare
;oplot,tim_flare,emd_mode_flare,color=250,thick=4

oplot,tim_flare,x_flare-emd_mode_flare -7. 

oplot,tim_pre, x_pre+trend_emd_pre - 50.
oplot,tim_pre,trend_emd_pre -50. ,color=50,thick=4	 ;200 is yelow, 150 is lime green, 250 red, 50 lbue

oplot,tim_pre,x_pre
;oplot,tim_pre,emd_mode_pre,color=130,thick=4

oplot,tim_pre,x_pre -7. 

oplot,tim_post, x_post+trend_emd_post - 50.,thick=3
oplot,tim_post,trend_emd_post -50. ,color=50,thick=4

oplot,tim_post,x_post
;oplot,tim_post,emd_mode_post,color=150,thick=4

oplot,tim_post,x_post -7. ;-emd_mode_post -7. 


device, /close
set_plot, 'X'


;ENDFOR
end