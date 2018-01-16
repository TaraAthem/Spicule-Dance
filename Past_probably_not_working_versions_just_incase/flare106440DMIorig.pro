pro broken_power_law,x,p,model
  model = -(x LT p[1])*p[0]*x + $
    (x GT p[1])*(-(p[0]-p[3])*p[1] - p[3]*x) + p[2]
  return
end

pro linfun, x, a, f
  f = a[0] + a[1] * x
  return
end

A = FindGen(16) * (!PI*2/16.) & UserSym, cos(A), sin(A), /fill
device, decomposed=0 & loadct,39
d=get_screen_size()
flarename= 'flare106440.fits'
data = mrdfits(flarename, 1, header)
ns=126
nf=201

flux=data.flux
tim=data.time
cad=tim[1]-tim[0]
flux_flare=flux[ns:nf]
tim_flare=tim[ns:nf]



N=n_elements(flux_flare)

modes=emd(flux_flare-mean(flux_flare),shi=0.0033, epsilon=1d-6, maxsiftings=1d8)
sz=size(modes)
trend=mean(flux_flare)+modes[*,sz[2]-1];+modes[*,sz[2]-2]

flux_flare=flux_flare-trend

LNP_freq_finder = LNP_TEST(tim_flare, flux_flare, OFAC=1, WK1 = freq, WK2 = power)
ZZ=FFT(flux_flare,-1)
PowerZZ=(ABS(ZZ))^2d
PowerZZ= PowerZZ[1: N_elements(freq) -1]

p = [2.D, -1.D, 0.1D, 0.D] ;Initial guesses of model parameters
result = LNP_TEST(tim_flare, flux_flare, /DOUBLE, wk1=freq, wk2=pow, ofac=1)
npow = N_ELEMENTS(pow)
freq = ALOG(freq)
pow = ALOG(pow)
err = REPLICATE(1., npow)
;result = mpfitfun('broken_power_law', freq[1:*], pow[1:*], err[1:*], p, yfit=pow_fit, /quiet, /NAN)
pow_fit = CURVEFIT(freq[1:*], pow[1:*], err[1:*], p, FUNCTION_NAME='broken_power_law',fita=[1,1,1,0],/noderivative)
result=p
print, 'Alpha_est_1 is', result[3]
print, 'Alpha_est_2 is', result[0]
print, 'Breaking point is',exp(-result[1]), ' seconds'

alpha_1=result[3]
alpha_2=result[0]
break_point=exp(-result[1])

modes=emd(flux_flare,shi=0.3, epsilon=1d-6, maxsiftings=1d8)
sz=size(modes)

modes=[[modes],[reform(trend, n_elements(trend),1)]]
sz=size(modes)

nmodep = 0.
pp = 0.
ep = 0.

for i=0,sz[2]-1 do begin
  length=float(n_elements(flux_flare));
  n_ext=float((size(extrema(modes[*,i])))[1])
  ;p=2d*length / n_ext
  p=(emd_period_energy(modes[*,i])).period
  e=total((modes[*,i])^2)
  pp = [pp, p]
  ep = [ep, e]
  nmodep = [nmodep,i+1]
endfor

period = pp[1:*]*cad
energy = ep[1:*]
nmode = nmodep[1:*]

print, 'Periods:',period

emd_mode=modes[*,2]

window, 4 ,xsize=0.8*d[0],ysize=0.9*d[1]

!p.multi=[0,2,2]

plot,tim_flare, flux_flare+trend,thick=2,charsize=ch,/xsty,/ysty,title='!5 Original signal',$  ;ORIGINAL SIGNAL SCREEN
  xtitle='Time',ytitle='Amplitude'
oplot,tim_flare,trend,color=150,thick=2

plot,tim_flare,flux_flare,/ysty,/xsty,charsize=ch,title='Detrended signal',$        ; DETRENDED SIGNAL SCREEN
  xtitle='Time',ytitle='Amplitude'
oplot,tim_flare,emd_mode,color=250,thick=3

plot,freq, pow,title='Power Spectrum', xtitle='Log Freq [Hz]', ytitle='Log Fourier power [a.u.]'
oplot, freq, pow_fit, linestyle=2, thick=3

plot,[0],[0],charsize=ch,/ylog,/xlog,xrange=[2d*cad,2d*N*cad],$
  yrange=[0.5*min(energy),2d*max(energy)],xsty=1,ysty=1,title='EMD spectrum for Alpha_1,2',YTICKFORMAT='(g6.0)',$
  xtitle='Mean period',ytitle='Normalised modal energy'

ind_fit_1=where(nmode ne sz[2] and period lt break_point)
if (n_elements(ind_fit_1) lt 2) then begin
  print,'Portion of the noise with alpha = ',string(alpha_1,format='(F3.1)'), ' is negligible.'
  goto, label_1
endif
weights1=replicate(1.0, N_ELEMENTS(period[ind_fit_1])) & lin1 = [0d,alpha_1-1d]
fit=CURVEFIT( alog(period[ind_fit_1]), alog(energy[ind_fit_1]),fita=[1,0],$
  weights1, lin1, FUNCTION_NAME= 'linfun',/noderivative)
confn1=emd_noise_conf(N, alpha_1, nsamples=500l, confidence_level = 0.05d)
meanenergy1=confn1.meanenergy[0]-alog(cad)*(alpha_1-1d)
norma1=exp(meanenergy1)/exp(lin1[0])
oplot,(2d*N*cad)/(confn1.k),confn1.up/norma1,color=50,thick=2
oplot,(2d*N*cad)/(confn1.k),confn1.down/norma1,color=50,thick=2
z = alog(period)
y1=(meanenergy1)+z*(alpha_1-1d)
oplot,period,exp(y1)/norma1,color=200

label_1:
ind_fit_2=where(nmode ne sz[2] and period gt break_point)
if (n_elements(ind_fit_2) lt 2) then begin
  print,'Portion of the noise with alpha = ',string(alpha_2,format='(F3.1)'), ' is negligible.'
  goto, label_2
endif
weights2=replicate(1.0, N_ELEMENTS(period[ind_fit_2])) & lin2 = [0d,alpha_2-1d]
fit2=CURVEFIT( alog(period[ind_fit_2]), alog(energy[ind_fit_2]),fita=[1,0],$
  weights2, lin2, FUNCTION_NAME= 'linfun',/noderivative)
confn2=emd_noise_conf(N, alpha_2, nsamples=500l, confidence_level = 0.05d)
meanenergy2=confn2.meanenergy[0]-alog(cad)*(alpha_2-1d)
norma2=exp(meanenergy2)/exp(lin2[0])
oplot,(2d*N*cad)/(confn2.k),confn2.up/norma2,color=40,thick=2
oplot,(2d*N*cad)/(confn2.k),confn2.down/norma2,color=40,thick=2
z = alog(period)
y2=(meanenergy2)+z*(alpha_2-1d)
oplot,period,exp(y2)/norma2,color=150

label_2:
oplot, period, energy, psym=8, color=250, symsize=1.1, thick=2

!p.multi=[0,1,1]

END