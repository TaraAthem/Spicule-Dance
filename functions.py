"""
Created on Fri 18 Jan 2019

@author: TMehta
"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mtpl
import scipy.fftpack
from math import exp
import glob
import os
from matplotlib.axes import Axes                           
from matplotlib.lines import Line2D    
from cycler import cycler  
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
import pylab
import datetime
import csv
import re
import math
import pandas as pd
import time, os
import sys


#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# LOWER TURNING POINT FUNCTION FROM N AND L (GONG)
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# THIS FUNCTION CALLS dvshiftgong AND GIVES AN OUTPUT OF A LOWER TURNING
# POINT FOR MODES DETECTED WITH GONG USING AN INPUT FROM N AND L BY 
# PASSING N AND L THROUGH dvshiftgong TO PRODUCE AN OUTPUTTED FREQUENCY
# AND THEN TAKES THAT FREQUENCY WITH THE l VALUE TO OUTPUT A LOWER
# TURNING POINT AS ltpval. 
# IF THE WEIGHTED MEAN FREQUENCY IS nan BECAUSE THE INPUT DATA SET WAS 
# INCOMPLETE THEN THIS VALUE IS SKIPPED
#
# OUTPUT:
# [0] = ltpval, (SCALAR), CORRESPONDS TO r/R (<= 1)
#
# [1] = n, (SCALAR), RADIAL DEGREE
#
# [2] = l, (SCALAR), HARMONIC DEGREE
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

def ltpfromnl(nl):
		n,l=nl
		freq = dvshiftgong(n,l)[0]					# [1] gives the weighted mean freq for a given l, n												
		if np.isnan(freq) == True:
			print('Weighted mean frequency was NaN so output for l=',l,'and n=',n,'has been skipped. ')
			ltpval=np.nan
			return ltpval,n,l
		else:
			ltpval=ltp(freq, l)
			#print('the ltp for l=',l, 'and n=', n, 'is', ltpval)
			return ltpval, n, l 


#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# LOWER TURNING POINT FUNCTION
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# THIS PROGRAM TAKES IN A VALUE FOR FREQUENCY, NU, AND ANGULAR FREQUENCY
# L AND GIVES A LOWER TURNING POINT AS A VALUE OF r/R. 
#
# FUTURE WORK; TURN THIS PROGRAM INTO A FUNCTION SO IT CAN BE CALLED INTO
# OTHER PYTHON PROGRAMS. 
#
# NOTE, THE SUB_INTER HERE TAKES ON A DIFFERENT FORM TO THAT IN THE IDL
# PROGRAM, WHERE SUB_Y GIVES THE INDEX DIRECTLY, WHEREAS IN PYTHON, THE 
# INDEX HAS TO BE EXPLICITELY FOUND VIA MAXDAT.INDEX(...)
# * http://owww.phys.au.dk/~jcd/solar_models/
#
# OUTPUT:
# ltpval, (SCALAR), CORRESPONDS TO r/R (<= 1)
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

def ltp(nu, l):
	path="/home/space/phrsnz/Desktop/Academic/Programs/PhD/Frequency_shift/misc_data/"
	R=6.955e10
	omega=1e-6*nu*2.*np.pi
	const=omega**2./(l*(l+1))
	dat=[]
	maxdat=[]
	with open(path+"ltp_txtfile.txt", "r") as f:							# READ IN THE DATA FROM THE TXT FILE FROM				
		for x in f:															# CHRISTENSEN - DALSGAARD'S SOLAR MODEL * 
			dat.append(x.split())											
	dat0 =  [float(item[0]) for item in dat]								# DAT0 GIVES r/R
	dat1 =  [float(item[1]) for item in dat]								# DAT1 GIVES c (cm/sec)
	ydat =  [float(a**2.)/float((R*b)**2.) for a,b in zip(dat1, dat0)]		# THIS IS NEEDED TO FIND THE INDEX OF WHERE														
	for i in ydat :															# THE MAX VALUE OF YDAT (WHICH IS LESS THAN 
		if i < const:														# THE CONSTANT) IS LOCATED
			maxdat.append(i)
	sub_y=max(i for i in ydat if i < const) 								# GIVES THE ACTUAL VALUE OF THE MAX (YDAT < CONST)
	sub_inter = maxdat.index(max(maxdat)) + (const - sub_y)/(ydat[maxdat.index(max(maxdat))+1]- sub_y)	
																			# FINDS THE VALUE TO INTERPOLATE ON, WHERE THIS MAX OCCURS
	fp=np.arange(0,len(dat0),1)												# CREATES AN ARRAY SO THE INDEX OF THE MAX CAN BE USED
	result= np.interp(sub_inter, fp, dat0)									# WITH DAT0 TO FIND THE LTP
	print('The lower turning point of a wave with nu =', nu, ' and l=', l, ' is', result, 'in units of r/R')
	return result



#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# FREQUENCY SHIFT CODE (GONG)
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# THIS FUNCTION TAKES IN A VALUE FOR N, L, AND GIVES FOUR OUTPUTS. 
# READS IN DATA FROM 152 DATASETS SEPARATED BY 36 DAYS, EACH GIVING A FIT
# OF FREQUENCY FROM A LORENTZIAN ON A FOURIER SPECTRUM TO GIVE NU, AND
# AN ASSOCIATED WEIGHTING VALUE, AMONST OTHER PARAMETERS. NOT ALL VALUES
# OF N AND L HAVE BEEN FOUND, SO SOME DATA SETS MAY BE INCOMPLETE. IF 
# A DATASET IS INCOMPLETE, THE COMPUTATION WILL AUTOMATICALLY TERMINATE. 
#
# OUTPUT:
# [0] = avg_freq, (SCALAR),   A WEIGHTED AVERAGE FREQUENCY READ IN FROM 
#                             INPUT DATA SETS.
# [1] = dv_shift, (1D ARRAY), A TIME SERIES OF THE FREQUENCIES NORMALISED
#							  ABOUT ZERO, BY SUBTRACTING THE avg_freq
# [2] = weight, (1D ARRAY),   A TIME SERIES OF THE WEIGHTINGS READ IN FROM
#							  THE INPUT DATA SETS
# [3] = mean_shift (SCALAR),  A MEASURE OF THE MAX AMPLITUDE OF dv_shift
# 							  FROM MAX TO MIN OVER THE ENTIRE DATASET. 
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

path="/home/space/phrsnz/Desktop/Academic/Programs/PhD/Frequency_shift/"
os.chdir(path+"GONG_v1z")
fname=sorted(glob.glob('mrv1*.txt'))										# SORTS THE FILES TO BE READ IN
dnum=152																	# NUMBER OF INPUT FILES

def dvshiftgong(n,l):
	
	def avgfreq_one(n,l,a):													# SUB FUNCTION THAT CALCULATES NU FOR EACH N,L
		modeldata= pd.read_csv(fname[a],sep="     |   |  ", skiprows=24, header=None, error_bad_lines=False, engine='python')
		nvals= modeldata.loc[:,0]											# READ_CVS READS DATA FASTER
		lvals= modeldata.loc[:,1]											# CODE FINDS INTERSECTION OF ALL THE INDEXES CORRESPONDING
		nuvals= modeldata.loc[:,2]											# TO L=l AND N=n TO GIVE A SPECIFIC N,L FROM WHICH NU AND 
		wvals= modeldata.loc[:,3]											# WEIGHT ARE READ. 
		indexl= [i for i, j in enumerate(lvals) if j == l]
		indexn= [i for i, j in enumerate(nvals) if j == n]
		index=list(set(indexl).intersection(indexn))
		if not index:
			return np.nan, np.nan, np.nan, np.nan
		else:																# FUNCTION AUTOMATICALLY QUITS IF THE DATASET IS INCOMPLETE
			return n,l,nuvals[index[0]], wvals[index[0]]
			
	def avg_freq_mult(n,l):													# THIS FUNCTION APPENDS EACH ITERATION OF AVGFREQ_ONE
		result=[]															# UNLESS A NAN IS ENCOUNTERED IN WHICH CASE IT RETURNS
		for a in range(1,152):												# A BOOLEAN NONE 
			y=avgfreq_one(n,l,a)
			if np.isnan(y[2]) == False:
				a+= a
				result.append([y[0],y[1],y[2],y[3]])
			else:
				result= None
				break
		return result
		
	result= avg_freq_mult(n,l)
	if result is not None:													# ONLY IF THE DATASET IS COMPLETE WILL THE COMPUTATION BE PERFORMED
		freq=[float(item[2])*float(item[3]) for item in result]				# DIFFERENT DATAPOINTS HAVE DIFFERENT WEIGHTINGS WHICH NEED
		weight=[float(item[3]) for item in result]							# TO BE CONSIDERED WHEN FINDING AN AVG VALUE
		avg_freq= sum(freq)/sum(weight)										# WEIGHTED ARITHMETIC MEAN (en.wikipedia.org/wiki/Weighted_arithmetic_mean)
		dv_shift= [float(item[2])-avg_freq for item in result]				# NORMALISE ABOUT ZERO.
		mean_shift= np.max(dv_shift)-np.min(dv_shift)						# 'MEAN SHIFT' IS A MISNOMER; ACTUALLY DEFINED AS THE MAXIMUM AMPLITUDE OF SIGNAL
		return avg_freq, dv_shift, weight, mean_shift		
	else:
		 print("Dataset for l=",l,"n=",n,"was incomplete so the weighted mean frequency cannot be calculated.")
		 return np.nan, np.nan, np.nan, np.nan	




#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# INERTIA RATIO FROM N, L (GONG)
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# THIS FUNCTION TAKES IN A VALUE OF N AND L (INPUTTED AS (n,l) E.G. 
# meanshiftfromnl((11,55)) FOR N=11, L=55 AND GIVES FIVE OUTPUTS:
# 
# OUTPUTS:
#
# [0] = MEANSHIFT,(SCALAR), DIRECTLY TAKEN FROM THE OUTPUT OF THE 
#							dvshiftgong CODE (SEE ABOVE).
#
# [1] = FREQ, (SCALAR), 	WEIGHTED MEAN FREQ. TAKEN FROM dvshiftgong
#
# [2] = n, (SCALAR),		RADIAL DEGREE
#
# [3] = l, (SCALAR), 		HARMONIC DEGREE  
#
# [4] = Q, (SCALAR), 		RATIO OF MASS OF MODE DIVIDED BY THE MASS OF
#							AN L = 0 MODE INTERPOLATED AT THE SAME FREQ.
# * http://owww.phys.au.dk/~jcd/solar_models/
#
# TO RUN THIS CODE OVER A RANGE OF VALUES USE THE BELOW:
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# nlarray=[]
# for n in range(nmin,nmax+1,1):											# CREATE ARRAY OF [(N1,L1),(N1,L2),...,(NM,LM)] #
#	for l in range(lmin,lmax+1,1):											# WHICH CAN BE READ INTO meanshiftfromnl FUNCTION		
#		nlarray.append((n,l))
#
# pool=mp.Pool() 															# USE MULTIPLE PROCESSORS TO SPEED UP COMPUTATION
# output=pool.map(meanshiftfromnl,nlarray)
# pool.close()
# pool.join()
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

def meanshiftfromnl(nl):
	path="/home/space/phrsnz/Desktop/Academic/Programs/PhD/Frequency_shift/misc_data/"
	dat=[]
	with open(path+"masses.txt", "r") as f:									# READ IN THE DATA FROM THE TXT FILE FROM				
		for x in f:															# CHRISTENSEN - DALSGAARD'S SOLAR MODEL * 
			dat.append(x.split())											
			ntxt =  [int(item[0]) for item in dat]							
			ltxt =  [int(item[1]) for item in dat]
			freqtxt = [float(item[2]) for item in dat]
			masstxt = [float(item[3]) for item in dat]

	n,l=nl
	x=dvshiftgong(n,l)
	meanshift=np.nan														# MEAN SHIFT IS GIVEN BY MAX TO MIN HEIGHT OF SHIFT
	freq=np.nan																# AND IS NOT AN ACTUAL 'MEAN' AS THE NAME SUGGESTS
	Q=np.nan
	
	if np.isnan(x[0]) == False:												# CHECKS IF THE DATASET IS EMPTY => MEAN= O EXACTLY. 
		meanshift= x[3]
		freq=x[0]
		indexl= [i for i, j in enumerate(ltxt) if j == l]
		indexn= [i for i, j in enumerate(ntxt) if j == n]
		#print('index=', list(set(indexl).intersection(indexn))[0]) 		# FINDS THE INDEXES IN masses.txt FOR GIVEN L AND N
		mass_orig=masstxt[list(set(indexl).intersection(indexn))[0]]		# FINDS THE UNIQUE INDEX FOR THE COMBINATION OF L,N
		mass_mod= np.interp(freq, freqtxt[0:35], masstxt[0:35])				# INTERPOLATES THE MODE MASS FOR THE GIVEN FREQUENCY 
		Q = mass_orig/mass_mod												# FOR AN L=0 MODE (FIRST 35 ENTRIES OF THE TXT FILE)
		print('Dataset for l=',l,'n=',n,'gives Q= %.2g, with a weighted mean frequency of %.2f \u03bcHz.' %(Q, x[0]))
	return meanshift,freq, l, n, Q



#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# FREQUENCY SHIFT CODE (BiSON)
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# THIS FUNCTION TAKES IN A VALUE FOR N, L. IT AUTOMATICALLY
# SCALES THROUGH FROM NMIN TO NMAX IN STEPS OF 1 (AND SIMILARLY FOR L)
# AND GIVES AN OUTPUT OF A WEIGHTED MEAN FREQUENCY IN MICROHZ AS A PRINTED
# RESULT.
# NOTE- OUTPUTS ARE EXLUDED TO HAVE A PEAK FREQ BETWEEN USER DEFINED
# FREQMIN AND FREQMAX (DEFAULT IS 1000,4000 microHz) AND ARE NOT USED IF 
# THE DATA SET IS INCOMPLETE (i.e HAS MISSING VALUES) 
#
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

		
def dvshiftbison(nmin,nmax,lmin, lmax,plot):
	ndiff=nmax-nmin
	path="/home/space/phrsnz/Desktop/Academic/Programs/PhD/Frequency_shift/"
	os.chdir(path+"BiSON_Data")
	fname=glob.glob('Main_Fits_waverage_fill_365d_*.dat')
	fname2='10_7_NOAA_FLUX.csv'
	fname3= 'sunspot_area.csv'

	solarmax1989= datetime.datetime(1989, 11, 1) 		# Dates of solar maxima
	solarmax2001= datetime.datetime(2001, 11, 1)	
	solarmax2014= datetime.datetime(2014, 4, 1)
	func= lambda arg:int(arg[:-4].split("_")[5])		# Sorts files 0,1,2...dnum
	fname=sorted(fname, key=func)						# by splitting filename 
														# and selecting integer
	dnum=127											# Total number of files
	cad=91.25 											#days
	x=np.zeros(shape=(lmax,dnum,ndiff))					# Create zero arrays
	w=np.zeros(shape=(lmax,dnum,ndiff))					# To fill with data
	flux=np.zeros(shape=(20000))
	area=np.zeros(shape=(20000))
	y=np.zeros(shape=(70,dnum))							# 70 is the upper limit 
	dv_avg=np.zeros(shape=(dnum))						# of modes to plot

	x_avg_w=np.zeros(shape=(lmax,ndiff))
	x_avg_nw=np.zeros(shape=(lmax,ndiff))
	w_sum=np.zeros(shape=(lmax,ndiff))
	
	freqmin=1000
	freqmax=4000

#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# COMPUTATION
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

	for a in range(0,dnum):								# Loops to read
		modeldata= np.loadtxt(fname[a], dtype=float)	# through all values
		#print(fname[a])
		totalcolumn=np.shape(modeldata)[0]				# of l and n over
		for b in range(0,totalcolumn):					# all files. Reads 
			l = modeldata[b][0]							# peak frequency into
			n = modeldata[b][1]							# x array, and weights
			for c in range(lmin,lmax):						# into w array. 
				lnum = c								# x_avg is partially
				for d in range(0,ndiff):				# computed in this loop
					nnum=nmin + d
					if l == lnum and n == nnum : 
						x[c][a][d]=modeldata[b][2]
						w[c][a][d]=modeldata[b][3]
						x_avg_w[c][d]=x_avg_w[c][d] + x[c][a][d]*w[c][a][d]
						w_sum[c][d]=w_sum[c][d] + w[c][a][d]

#index1=np.argwhere(x != 0)[:,0]						# Uncomment to read
#index2=np.argwhere(x != 0)[:,2]						# indexes of null data
	
	for c in range(lmin,lmax):
		for d in range(0,ndiff):
			x_avg_w[c][d]=x_avg_w[c][d]/w_sum[c][d]		# Weighted average  
			x_avg_nw[c][d]=np.mean(x[c,:,d])			# Not- weighted
			x[c,:,d] -= x_avg_w[c][d]					# Normalised about (weighted) mean
			if (min(x[c,:,d])>-530 and (freqmin< x_avg_w[c][d] <= freqmax) ) :
				print(u'Weighted average frequency for the l=',lmin+c,'n=',d+nmin,'mode is = %.2f \u03bcHz.' %(x_avg_w[c][d]))
			else:
				print ('No valid output for l=',lmin+c,'and n=',d+nmin,'as the dataset was incomplete')
			#print('Average frequency for the l=',lmin+c,'n=',d,'mode is =',x_avg_nw[c][d],'micro Hz.')
			nu=x_avg_w[c][d]


	n=0
	with open(fname2) as csvfile:						# 10.7cm flux cvs file
		readCSV = csv.reader(csvfile, delimiter=',')
		next(readCSV, None)								# Skips reading the header
		for row in readCSV:								# Create a flux array
			flux[n]=row[1]
			n+=1
		
	n=0
	with open(fname3) as csvfile:						# Sunspot area cvs file
		readCSV = csv.reader(csvfile, delimiter=' ')
		next(readCSV, None)								# Skips reading the header
		for row in readCSV:								# Create an array of areas
			area[n]=row[-1]
			n+=1

	flux=flux[0:91*dnum]								# To match cadence of t2 array
	flux=flux[(flux>=0)]/max(flux)						# Array is normalised
	flux=flux[0::91]

	area=area[1320:1699]								# Select elements starting at 1/1/85
	area=area[0::3]/max(area)							# and ending 6/1/16. Then selecting 
														# every 3rd element to match cadence 

#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# VISUALISATION
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

	startdate = datetime.datetime(1985, 1, 1)			# Cadence of data is actually 91.25 days
	enddate = datetime.datetime(2016, 6, 1)				# not 91, but needed an interger input. Revise?
	t=[startdate + datetime.timedelta(days=i) for i in range(0,127*91,91)]

	e=0
	if plot == 'Y':  
		plt.figure(100)

		ax= plt.gca()
		ax.grid()
		ax.set_xlabel('Year')
		ax.set_ylabel(r'$ \delta \nu\ (\mu Hz)$')
		ax.set_title(r'Change in $\delta \nu $ over time')
		ax.set_ylim(-2.75, 1.2)
		plt.xlim(left=startdate, right= enddate)
		plt.plot([solarmax1989,solarmax1989],[-5,5], 'ko--')
		plt.plot([solarmax2001,solarmax2001],[-5,5], 'ko--')
		plt.plot([solarmax2014,solarmax2014],[-5,5], 'ko--')              
		
		for c in range(lmin,lmax):
			for d in range(0,ndiff):						
				if (min(x[c,:,d])>-530) :					# Exclude data with zero values
					plt.plot(t, x[c,:,d])					
					y[e,:]=x[c,:,d]							# Valid data in a new y array
					e+=1									 
												
		print(e,'modes have been plotted and averaged')
		#y=y[0:e-1,:]										# Only select the non-zero arrays		
		#dv_avg= np.mean(y, axis=0)							# Averaged to give the mean change in nu		
		#plt.plot( t, dv_avg, linewidth=3.0)
		plt.plot(t, flux-1.5, color='red')					# -1.2 shifts the data
		plt.plot(t, area - 2.5, color='blue')				# -2.4 shifts data
		e=0
		plt.show()
	
	
	# NON PLOTTING OPTION- STILL GIVES Y[0,:] BUT DOESN'T PLOT OUTPUT. 
	else:
		for c in range(lmin,lmax):
			for d in range(0,ndiff):							
				if (min(x[c,:,d])>-530 and (freqmin< x_avg_w[c][d] <= freqmax) ) :					
					y[e,:]=x[c,:,d]							# Valid data in a new y array
					e+=1	
	print(e,'of the', ndiff*(lmax-lmin),'(%.1f %%) input modes have full datasets.' %(100*float(e/(ndiff*(lmax-lmin)))) )
	return y[e-1,:]
		
