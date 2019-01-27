
;----------------------------------------------------
;
; GOES - Geostationary Operational Environmental Satellite
;
; 008.pro is a good model for the GOES EMD (1 sig mode)
; 052.pro is a good model for the GOES EMD (2 sig modes)
;
;
; NoRH - Nobeyama RadioHeliograph
;
; 022.pro is a good model for the NoRH EMD
;
;
; VERNOV
; 010_2 is good. Note the exclusion of the [ind]. 
; THERE IS PROBABLY A WAY TO FIX THE IND ISSUE EASILY, BUT DUE TO TIME CONSTRAINTS AND MUCH CRYING,
; IVE JUST DELETED THE DATA TO ONLY READ IN IN THE REGION OF CHOICE. 
;
;
;
; Fermi - 
; 056.pro is a good model for the FERMI EMD program
; 104_2.pro is multimodal
;
;
;
; EVE
; 129.pro works. But I dont think EMD works at ALL well with EVE data (both failed to compile 642957693 times) :( 
;----------------------------------------------------




; USE BELOW CODE FOR MULTIMODAL OUTPUTS

;-----multimodal
mode1=  modes[*,%%%%'CHANGEMECHANGEME']
mode2= modes[*,%%%%%&'X']
;--------------
trend=trend_emd
time=t
cad=dt
energy=sp.energy
save,mode1, mode2, trend,sig_orig,time,cad,periods,Q,mod_depth,dur,energy,conf95,conf99,avg_energy,per_arr,emin, emax,FILENAME=path + '/Sav_files/'+flaresav+'.sav'

Pugh paper: https://arxiv.org/abs/1709.09472

