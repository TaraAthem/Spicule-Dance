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

flux_preflare=flux[0:ns]
tim_preflare=tim[0:ns]

N_preflare=n_elements(flux_preflare)
sh=0.007
modes=emd(flux_preflare-mean(flux_preflare),shi=sh, epsilon=1d-6, maxsiftings=1d8)
sz=size(modes)
trend_preflare=mean(flux_preflare)+modes[*,sz[2]-1]+modes[*,sz[2]-2]
;emd_mode=modes[*,sz[2]-2]

LNP_freq_finder = LNP_TEST(tim_preflare, flux_preflare, OFAC=1, WK1 = freqPRE, WK2 = power)
ZZ=FFT(flux_preflare,-1)
PowerZZ=(ABS(ZZ))^2d
PowerZZ= PowerZZ[1: N_elements(freqPRE) -1]

p = [2.D, -1.D, 0.1D, 0.D] ;Initial guesses of model parameters
result = LNP_TEST(tim_preflare, flux_preflare, /DOUBLE, wk1=freqPRE, wk2=powPRE, ofac=1)
npow = N_ELEMENTS(powPRE)
freqPRE = ALOG(freqPRE)
powPRE = ALOG(powPRE)

resultPRE = ladfit(freqPRE[1:*], powPRE[1:*])
 ; print, 'full result', resultPRE
 ; print, 'a value is gradient', resultPRE[1]

alpha_pre=-resultPRE[1]

print, 'Alpha Pre-Flare is =', alpha_pre

Periodogram, tim_preflare,flux_preflare, fap=0.05, power_sign=FFT95PRE			;FFT CONDIFENCE INTERVALS USING FALSE ALARM PROBABILITY
Periodogram, tim_preflare,flux_preflare, fap=0.01, power_sign=FFT99PRE			;FOR 95 AND 99%
err = REPLICATE(1., npow)
;result = mpfitfun('broken_power_law', freqPRE[1:*], pow[1:*], err[1:*], p, yfit=pow_fit, /quiet, /NAN)

; pow_fitPRE = CURVEFIT(freqPRE[1:*], powPRE[1:*], err[1:*], p, FUNCTION_NAME='broken_power_law',fita=[1,1,1,0],/noderivative)
; ind=where(nmode ne max(nmode))
; weights=replicate(1.0, N_ELEMENTS(period[ind]))
; lin = [0d,alpha-1d]
; pow_fitPRE=CURVEFIT(freqPRE[1:*], powPRE[1:*],fita=[1,0], weigths, lin, FUNCTION_NAME= 'linfun',/noderivative)
; 
; 
; err = REPLICATE(1., npow)
; CURVE_FIT_PRE = mpfitfun('linfun', freqPRE[1:*],  powPRE[1:*], err[1:*], p, yfit=pow_fitPRE, /quiet, /NAN)

pow_fitPRE=-alpha_pre*freqPRE + resultPRE[0]

sh=0.02
modes=emd(flux_preflare-trend_preflare,shi=sh, epsilon=1d-6, maxsiftings=1d8)
sz=size(modes)
trend_preflare=trend_preflare+modes[*,sz[2]-1]+modes[*,sz[2]-2]
emd_mode_preflare=modes[*,sz[2]-5]

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

; alpha=0d
; ind=where(nmode_preflare ne max(nmode_preflare))
; weights=replicate(1.0, N_ELEMENTS(period_preflare[ind]))
lin = [0d,alpha_pre-1d]
; fit=CURVEFIT( alog(period_preflare[ind]), alog(energy_preflare[ind]),fita=[1,0],$
;   weigths, lin, FUNCTION_NAME= 'linfun',/noderivative)
; 
; alpha_preflare=1d + lin[1]
; print,'Alpha Preflare is',alpha_preflare

confn_99_preflare=emd_noise_conf(N_preflare,alpha_pre,nsamples=300l, conf=0.01)
confn_95_preflare=emd_noise_conf(N_preflare,alpha_pre,nsamples=300l, conf=0.05)

meanenergy_preflare=confn_99_preflare.meanenergy[0]-alog(cad)*(alpha_pre-1d)
norma_preflare=exp(meanenergy_preflare)/exp(lin[0])




;*****************************************************************************************
;------------------------------------------FLARE------------------------------------------
;*****************************************************************************************

flux=data.flux
tim=data.time
cad=tim[1]-tim[0]
flux_flare=flux[ns:nf]
tim_flare=tim[ns:nf]


;window, 0 & plot,tim, flux
;oplot,tim[ns:nf],flux_preflare,color=250

N_flare=n_elements(flux_flare)


LNP_freq_finder = LNP_TEST(tim_flare, flux_flare, OFAC=1, WK1 = freqFLARE, WK2 = power)
ZZ=FFT(flux_flare,-1)
PowerZZ=(ABS(ZZ))^2d
PowerZZ= PowerZZ[1: N_elements(freqFLARE) -1]

p = [2.D, -1.D, 0.1D, 0.D] ;Initial guesses of model parameters
result = LNP_TEST(tim_flare, flux_flare, /DOUBLE, wk1=freqFLARE, wk2=powFLARE, ofac=1)
npow = N_ELEMENTS(powFLARE)
freqFLARE = ALOG(freqFLARE)
powFLARE = ALOG(powFLARE)

Periodogram, tim_flare,flux_flare, fap=0.05, power_sign=FFT95FLARE			;FFT CONDIFENCE INTERVALS USING FALSE ALARM PROBABILITY
Periodogram, tim_flare,flux_flare, fap=0.01, power_sign=FFT99FLARE			;FOR 95 AND 99%
; err = REPLICATE(1., npow)
; ;result = mpfitfun('broken_power_law', freqFLARE[1:*], pow[1:*], err[1:*], p, yfit=pow_fit, /quiet, /NAN)
; pow_fitFLARE = CURVEFIT(freqFLARE[1:*], powFLARE[1:*], err[1:*], p, FUNCTION_NAME='broken_power_law',fita=[1,1,1,0],/noderivative)
; 
; err = REPLICATE(1., npow)
; CURVE_FIT_FLARE = mpfitfun('broke_pow_law2', freqFLARE[1:*],  powFLARE[1:*], err[1:*], p, yfit=pow_fitFLARE, /quiet, /NAN)

resultFLARE = ladfit(freqFLARE[1:*], powFLARE[1:*])
alpha_flare=-resultFLARE[1]

Print, 'Alpha Flare=', alpha_flare

pow_fitFLARE=-alpha_flare*freqFLARE + resultFLARE[0]

sh=0.0033
modes=emd(flux_flare-mean(flux_flare),shi=sh, epsilon=1d-6, maxsiftings=1d8)
sz=size(modes)
trend=mean(flux_flare)+modes[*,sz[2]-1];+modes[*,sz[2]-2]
;emd_mode=modes[*,sz[2]-2]

sh=0.3
modes=emd(flux_flare-trend,shi=sh, epsilon=1d-6, maxsiftings=1d8)
sz=size(modes)
trend=trend+modes[*,sz[2]-1]+modes[*,sz[2]-2]
emd_mode=modes[*,sz[2]-5]

nmodep = 0.
pp = 0.
ep = 0.

modes/=sqrt(total((flux_flare-mean(flux_flare))^2))

for i=0,sz[2]-3 do begin
  p=(emd_period_energy(modes[*,i])).period
  e=total((modes[*,i])^2)
  pp = [pp, p]
  ep = [ep, e]
  nmodep = [nmodep,i+1]
endfor

trend_x=trend/sqrt(total((flux_flare-mean(flux_flare))^2));+modes[*,sz[2]-1]+modes[*,sz[2]-2];+modes[*,sz[2]-3]+modes[*,sz[2]-4]+modes[*,sz[2]-5]
p=(emd_period_energy(trend_x)).period
e=total((trend_x)^2)
pp = [pp, p]
ep = [ep, e]
nmodep = [nmodep,sz[2]-1]

period_flare = pp[1:*]*cad
energy_flare = ep[1:*]
nmode_flare = nmodep[1:*]

print, 'periods',period_flare

; 
; alpha=0d
; ind=where(nmode ne max(nmode))
; weights=replicate(1.0, N_ELEMENTS(period[ind]))
 lin = [0d,alpha_flare-1d]
; ;fit=CURVEFIT( alog(period[ind]), alog(energy[ind]),fita=[1,0],$
;   weigths, lin, FUNCTION_NAME= 'linfun',/noderivative)
; 
; alpha=1d + lin[1]
; ;print, 'Alpha flare is', alpha

confn_99=emd_noise_conf(N_flare,alpha_flare,nsamples=300l, conf=0.01)
confn_95=emd_noise_conf(N_flare,alpha_flare,nsamples=300l, conf=0.05)

meanenergy_flare=confn_99.meanenergy[0]-alog(cad)*(alpha_flare-1d)
norma_flare=exp(meanenergy_flare)/exp(lin[0])




;*****************************************************************************************
;---------------------------------------POSTFLARE------------------------------------------
;*****************************************************************************************


flux_postflare=flux[nf:299]
tim_postflare=tim[nf:299]

N_postflare=n_elements(flux_postflare)


LNP_freq_finder = LNP_TEST(tim_postflare, flux_postflare, OFAC=1, WK1 = freqPOST, WK2 = power)
ZZ=FFT(flux_postflare,-1)
PowerZZ=(ABS(ZZ))^2d
PowerZZ= PowerZZ[1: N_elements(freqPOST) -1]

p = [2.D, -1.D, 0.1D, 0.D] ;Initial guesses of model parameters
result = LNP_TEST(tim_postflare, flux_postflare, /DOUBLE, wk1=freqPOST, wk2=powPOST, ofac=1)
npow = N_ELEMENTS(powPOST)
freqPOST = ALOG(freqPOST)
powPOST = ALOG(powPOST)


Periodogram, tim_postflare,flux_postflare, fap=0.05, power_sign=FFT95POST			;FFT CONDIFENCE INTERVALS USING FALSE ALARM PROBABILITY
Periodogram, tim_postflare,flux_postflare, fap=0.01, power_sign=FFT99POST			;FOR 95 AND 99%

; err = REPLICATE(1., npow)
; ; ;result = mpfitfun('broken_power_law', freqPOST[1:*], powPOST[1:*], err[1:*], p, yfit=pow_fit, /quiet, /NAN)
;  pow_fitPOST = CURVEFIT(freqPOST[1:*], powPOST[1:*], err[1:*], p, FUNCTION_NAME='broken_power_law',fita=[1,1,1,0],/noderivative)
; err = REPLICATE(1., npow)
; CURVE_FIT_POST = mpfitfun('broke_pow_law2', freqPOST[1:*],  powPOST[1:*], err[1:*], p, yfit=pow_fitPOST, /quiet, /NAN)

resultPOST = ladfit(freqPOST[1:*], powPOST[1:*])
alpha_post=-resultPOST[1]

Print, 'Alpha Post Flare =', alpha_post

pow_fitPOST=-alpha_post*freqPOST + resultPOST[0]

sh=0.003
modes=emd(flux_postflare-mean(flux_postflare),shi=sh, epsilon=1d-6, maxsiftings=1d8)
sz=size(modes)
trend_postflare=mean(flux_postflare)+modes[*,sz[2]-1]+modes[*,sz[2]-2]
;emd_mode=modes[*,sz[2]-2]

sh=0.08
modes=emd(flux_postflare-trend_postflare,shi=sh, epsilon=1d-6, maxsiftings=1d8)
sz=size(modes)
trend_postflare=trend_postflare+modes[*,sz[2]-1];+modes[*,sz[2]-2]
emd_mode_postflare=modes[*,sz[2]-3]

nmodep = 0.
pp = 0.
ep = 0.

modes/=sqrt(total((flux_postflare-mean(flux_postflare))^2))

for i=0,sz[2]-2 do begin
  p=(emd_period_energy(modes[*,i])).period
  e=total((modes[*,i])^2)
  pp = [pp, p]
  ep = [ep, e]
  nmodep = [nmodep,i+1]
endfor

trend_x=trend_postflare/sqrt(total((flux_postflare-mean(flux_postflare))^2));+modes[*,sz[2]-1]+modes[*,sz[2]-2];+modes[*,sz[2]-3]+modes[*,sz[2]-4]+modes[*,sz[2]-5]
p=(emd_period_energy(trend_x)).period
e=total((trend_x)^2)
pp = [pp, p]
ep = [ep, e]
nmodep = [nmodep,sz[2]]

period_postflare = pp[1:*]*cad
energy_postflare = ep[1:*]
nmode_postflare = nmodep[1:*]

print, 'periods',period_postflare

; alpha=0d
; ind=where(nmode_postflare ne max(nmode_postflare))
; weights=replicate(1.0, N_ELEMENTS(period_postflare[ind]))
lin = [0d,alpha_post-1d]
; fit=CURVEFIT( alog(period_postflare[ind]), alog(energy_postflare[ind]),fita=[1,0],$
;   weigths, lin, FUNCTION_NAME= 'linfun',/noderivative)
; 
; alpha_post=1d + lin[1]
; print, 'Alpha postflare is',alpha_post

confn_99_postflare=emd_noise_conf(N_postflare,alpha_post,nsamples=300l, conf=0.01)
confn_95_postflare=emd_noise_conf(N_preflare,alpha_post,nsamples=300l, conf=0.05)

meanenergy_postflare=confn_99_postflare.meanenergy[0]-alog(cad)*(alpha_post-1d)
norma_postflare=exp(meanenergy_postflare)/exp(lin[0])


pow_fit95PRE=pow_fitPRE + ALOG(FFT95PRE)
pow_fit95FLARE=pow_fitFLARE + ALOG(FFT95FLARE)
pow_fit95POST=pow_fitPOST + ALOG(FFT95POST)

pow_fit99PRE=pow_fitPRE + ALOG(FFT99PRE)
pow_fit99FLARE=pow_fitFLARE + ALOG(FFT99FLARE)
pow_fit99POST=pow_fitPOST + ALOG(FFT99POST)


flux_preflare_detrended= flux_preflare-trend_preflare - emd_mode_preflare -5.
flux_flare_detrended= flux_flare-trend - emd_mode -5.
flux_postflare_detrended= flux_postflare-trend_postflare - emd_mode_postflare -5.


M=2						;make m=/=1 to enable visualisation screens 
if (M eq 1) then begin
  goto, label_5
endif



;*****************************************************************************************
;-----------------------------------VISUALISATION-----------------------------------------
;*****************************************************************************************

d=get_screen_size() & ch=2.5
A = FIndGen(16) * (!PI*2/16.) & UserSym, cos(A), sin(A), /fill

window,1,xsize=0.4*d[0],ysize=0.7*d[1]



plot,tim, flux-25d,thick=1,charsize=ch,/xsty,/ysty,title=FLARETITLE, yrange=[-8,20],$
  xtitle='Elapsed time [s]',ytitle='Normalised amplitude [a.u.]' 
oplot,tim_flare, trend-25d,color=250,thick=3
oplot,tim_preflare, trend_preflare-25d,color=250,thick=3
oplot,tim_postflare, trend_postflare-25d,color=250,thick=3

oplot,tim_flare,flux_flare-trend
oplot,tim_flare,emd_mode,color=250,thick=3

oplot,tim_preflare,flux_preflare-trend_preflare
oplot,tim_preflare,emd_mode_preflare,color=200,thick=3

oplot,tim_postflare,flux_postflare-trend_postflare
oplot,tim_postflare,emd_mode_postflare,color=200,thick=3

oplot, tim_preflare, flux_preflare_detrended
oplot, tim_flare, flux_flare_detrended
oplot, tim_postflare, flux_postflare_detrended

; 
; ;----------Preflare EMD and FFT---------------------------------------
; 
; window,2,xsize=0.6*d[0],ysize=0.7*d[1]
; !p.multi=[0,2,3]
; 
; 
; plot,[0],[0],charsize=ch,/ylog,/xlog,xrange=[2d*cad,2d*N_preflare*cad],$
;   yrange=[1d-4,1d4],xsty=1,ysty=1,title='EMD Spectrum, Preflare',$
;   xtitle='Period [s]',ytitle='EMD Energy [a.u.]';, position=[0.7,0.74,0.99,0.96],
; 
; oplot,(2d*N_preflare*cad)/(confn_99_preflare.k),confn_99_preflare.up,color=50,thick=1
; oplot,(2d*N_preflare*cad)/(confn_99_preflare.k),confn_99_preflare.down,color=50,thick=1
; 
; oplot,(2d*N_preflare*cad)/(confn_95_preflare.k),confn_95_preflare.up,color=150,thick=1
; oplot,(2d*N_preflare*cad)/(confn_95_preflare.k),confn_95_preflare.down,color=150,thick=1
; 
; z = alog(period_preflare)
; y=meanenergy_preflare+z*(alpha_pre-1d)
; oplot,period_preflare,exp(y),linesty=2
; oplot, period_preflare, norma_preflare*energy_preflare,psym=8,color=250,symsize=1.1,thick=2
; 
; 
; plot,EXP(-freqPRE),powPRE, charsize=ch,  /xlog, /xstyle,xrange=[2d*cad,2d*N_postflare*cad], /ystyle,title='Fourier Spectrum, Preflare', xtitle='Period [s]', ytitle='Fourier Energy [a.u.]'
; oplot, EXP(-freqPRE), pow_fitPRE, linestyle=2, thick=3
; oplot, EXP(-freqPRE), pow_fit95PRE, color=150
; oplot, EXP(-freqPRE), pow_fit99PRE, color=50
; 
; 
; ;----------Flare EMD and FFT---------------------------------------
; plot,[0],[0],charsize=ch,/ylog,/xlog,xrange=[2d*cad,2d*N_preflare*cad],$
;   yrange=[1d-4,1d4],xsty=1,ysty=1,title='EMD Spectrum, Flare',$
;   xtitle='Period [s]',ytitle='EMD Energy [a.u.]';,position=[0.7,0.4,0.99,0.63]
; 
; oplot,(2d*N_flare*cad)/(confn_99.k),confn_99.up,color=50,thick=1
; oplot,(2d*N_flare*cad)/(confn_99.k),confn_99.down,color=50,thick=1
; 
; oplot,(2d*N_flare*cad)/(confn_95.k),confn_95.up,color=150,thick=1
; oplot,(2d*N_flare*cad)/(confn_95.k),confn_95.down,color=150,thick=1
; 
; z = alog(period_flare)
; y=meanenergy_flare+z*(alpha_flare-1d)
; oplot,period_flare,exp(y),linesty=2
; oplot, period_flare, norma_flare*energy_flare,psym=8,color=250,symsize=1.1,thick=2
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
; oplot,(2d*N_postflare*cad)/(confn_99_postflare.k),confn_99_postflare.up,color=50,thick=1
; oplot,(2d*N_postflare*cad)/(confn_99_postflare.k),confn_99_postflare.down,color=50,thick=1
; 
; oplot,(2d*N_postflare*cad)/(confn_95_postflare.k),confn_95_postflare.up,color=150,thick=1
; oplot,(2d*N_postflare*cad)/(confn_95_postflare.k),confn_95_postflare.down,color=150,thick=1
; 
; z = alog(period_postflare)
; y=meanenergy_postflare+z*(alpha_post-1d)
; oplot,period_postflare,exp(y),linesty=2
; oplot, period_postflare, norma_postflare*energy_postflare,psym=8,color=250,symsize=1.1,thick=2
; 
; plot,EXP(-freqPOST), powPOST, charsize=ch, xrange=[2d*cad,2d*N_postflare*cad], /xlog, /xstyle, /ystyle,title='Fourier Spectrum, Postflare', xtitle='Period [s]', ytitle='Fourier Energy [a.u.]'
; oplot, EXP(-freqPOST), pow_fitPOST, linestyle=2, thick=3
; oplot, EXP(-freqPOST), pow_fit95POST, color=150
; oplot, EXP(-freqPOST), pow_fit99POST, color=50
; 
; !p.multi=[0,1,1]

label_5:

END