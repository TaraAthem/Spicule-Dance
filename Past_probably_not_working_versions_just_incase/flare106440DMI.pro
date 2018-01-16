pro linfun, x, a, f
  f = a[0] + a[1] * x
  return
end

flarename= 'flare106440.fits'
data = mrdfits(flarename, 1, header)
ns=126
nf=201

flux=data.flux
tim=data.time
cad=tim[1]-tim[0]
flux_flare=flux[ns:nf]
tim_flare=tim[ns:nf]

device, decomposed=0 & loadct,39
;window, 0 & plot,tim, flux
;oplot,tim[ns:nf],flux_preflare,color=250

N=n_elements(flux_flare)
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

period = pp[1:*]*cad
energy = ep[1:*]
nmode = nmodep[1:*]

print, 'periods',period

alpha=0d
ind=where(nmode ne max(nmode))
weights=replicate(1.0, N_ELEMENTS(period[ind]))
lin = [0d,alpha-1d]
fit=CURVEFIT( alog(period[ind]), alog(energy[ind]),fita=[1,0],$
  weigths, lin, FUNCTION_NAME= 'linfun',/noderivative)

alpha=1d + lin[1]
print, alpha

confn_99=emd_noise_conf(N,alpha,nsamples=500l, conf=0.01)
confn_95=emd_noise_conf(N,alpha,nsamples=500l, conf=0.05)

meanenergy=confn_99.meanenergy[0]-alog(cad)*(alpha-1d)
norma=exp(meanenergy)/exp(lin[0])




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

ns=0
nf=125

flux_preflare=flux[ns:nf]
tim_preflare=tim[ns:nf]

N_preflare=n_elements(flux_preflare)
sh=0.007
modes=emd(flux_preflare-mean(flux_preflare),shi=sh, epsilon=1d-6, maxsiftings=1d8)
sz=size(modes)
trend_preflare=mean(flux_preflare)+modes[*,sz[2]-1]+modes[*,sz[2]-2]
;emd_mode=modes[*,sz[2]-2]

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

alpha=0d
ind=where(nmode_preflare ne max(nmode_preflare))
weights=replicate(1.0, N_ELEMENTS(period_preflare[ind]))
lin = [0d,alpha-1d]
fit=CURVEFIT( alog(period_preflare[ind]), alog(energy_preflare[ind]),fita=[1,0],$
  weigths, lin, FUNCTION_NAME= 'linfun',/noderivative)

alpha_preflare=1d + lin[1]
print, alpha_preflare

confn_99_preflare=emd_noise_conf(N_preflare,alpha_preflare,nsamples=500l, conf=0.01)
confn_95_preflare=emd_noise_conf(N_preflare,alpha_preflare,nsamples=500l, conf=0.05)

meanenergy_preflare=confn_99_preflare.meanenergy[0]-alog(cad)*(alpha-1d)
norma_preflare=exp(meanenergy_preflare)/exp(lin[0])

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

ns=202
nf=299

flux_postflare=flux[ns:nf]
tim_postflare=tim[ns:nf]

N_postflare=n_elements(flux_postflare)
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

alpha=0d
ind=where(nmode_postflare ne max(nmode_postflare))
weights=replicate(1.0, N_ELEMENTS(period_postflare[ind]))
lin = [0d,alpha-1d]
fit=CURVEFIT( alog(period_postflare[ind]), alog(energy_postflare[ind]),fita=[1,0],$
  weigths, lin, FUNCTION_NAME= 'linfun',/noderivative)

alpha_postflare=1d + lin[1]
print, alpha_postflare

confn_99_postflare=emd_noise_conf(N_postflare,alpha_postflare,nsamples=500l, conf=0.01)
confn_95_postflare=emd_noise_conf(N_preflare,alpha_postflare,nsamples=500l, conf=0.05)

meanenergy_postflare=confn_99_postflare.meanenergy[0]-alog(cad)*(alpha-1d)
norma_postflare=exp(meanenergy_postflare)/exp(lin[0])

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

d=get_screen_size() & ch=2.5
A = FIndGen(16) * (!PI*2/16.) & UserSym, cos(A), sin(A), /fill
window,1,xsize=0.6*d[0],ysize=0.7*d[1]

!p.multi=[0,4,1]

plot,tim, flux-25d,thick=1,charsize=ch,/xsty,/ysty,title='!5 Flare 106440', yrange=[-5,20],$
  xtitle='Elapsed time [s]',ytitle='Normalised amplitude [a.u.]', position=[0.08,0.12,0.6,0.94]
oplot,tim_flare, trend-25d,color=250,thick=3
oplot,tim_preflare, trend_preflare-25d,color=250,thick=3
oplot,tim_postflare, trend_postflare-25d,color=250,thick=3

oplot,tim_flare,flux_flare-trend
oplot,tim_flare,emd_mode,color=250,thick=3

oplot,tim_preflare,flux_preflare-trend_preflare
oplot,tim_preflare,emd_mode_preflare,color=200,thick=3

oplot,tim_postflare,flux_postflare-trend_postflare
oplot,tim_postflare,emd_mode_postflare,color=200,thick=3

plot,[0],[0],charsize=ch,/ylog,/xlog,xrange=[2d*cad,2d*N_preflare*cad],$
  yrange=[1d-4,1d4],xsty=1,ysty=1,title='Spectrum, Preflare', position=[0.7,0.74,0.99,0.96],$
  xtitle='Period [s]',ytitle='Energy [a.u.]'

oplot,(2d*N_preflare*cad)/(confn_99_preflare.k),confn_99_preflare.up,color=50,thick=1
oplot,(2d*N_preflare*cad)/(confn_99_preflare.k),confn_99_preflare.down,color=50,thick=1

oplot,(2d*N_preflare*cad)/(confn_95_preflare.k),confn_95_preflare.up,color=150,thick=1
oplot,(2d*N_preflare*cad)/(confn_95_preflare.k),confn_95_preflare.down,color=150,thick=1

z = alog(period_preflare)
y=meanenergy_preflare+z*(alpha_preflare-1d)
oplot,period_preflare,exp(y),linesty=2

oplot, period_preflare, norma_preflare*energy_preflare,psym=8,color=200,symsize=1.1,thick=2

plot,[0],[0],charsize=ch,/ylog,/xlog,xrange=[2d*cad,2d*N*cad],$
  yrange=[1d-4,1d4],xsty=1,ysty=1,title='Spectrum, Flare',position=[0.7,0.4,0.99,0.63],$
  xtitle='Period [s]',ytitle='Energy [a.u.]'

oplot,(2d*N*cad)/(confn_99.k),confn_99.up,color=50,thick=1
oplot,(2d*N*cad)/(confn_99.k),confn_99.down,color=50,thick=1

oplot,(2d*N*cad)/(confn_95.k),confn_95.up,color=150,thick=1
oplot,(2d*N*cad)/(confn_95.k),confn_95.down,color=150,thick=1

z = alog(period)
y=meanenergy+z*(alpha-1d)
oplot,period,exp(y),linesty=2

oplot, period, norma*energy,psym=8,color=250,symsize=1.1,thick=2

plot,[0],[0],charsize=ch,/ylog,/xlog,xrange=[2d*cad,2d*N_postflare*cad],$
  yrange=[1d-4,1d4],xsty=1,ysty=1,title='Spectrum, Postflare', position=[0.7,0.06,0.99,0.3],$
  xtitle='Period [s]',ytitle='Energy [a.u.]'

oplot,(2d*N_postflare*cad)/(confn_99_postflare.k),confn_99_postflare.up,color=50,thick=1
oplot,(2d*N_postflare*cad)/(confn_99_postflare.k),confn_99_postflare.down,color=50,thick=1

oplot,(2d*N_postflare*cad)/(confn_95_postflare.k),confn_95_postflare.up,color=150,thick=1
oplot,(2d*N_postflare*cad)/(confn_95_postflare.k),confn_95_postflare.down,color=150,thick=1

z = alog(period_postflare)
y=meanenergy_postflare+z*(alpha_postflare-1d)
oplot,period_postflare,exp(y),linesty=2

oplot, period_postflare, norma_postflare*energy_postflare,psym=8,color=200,symsize=1.1,thick=2

!p.multi=[0,1,1]

END