import csv
import numpy as np
import scipy
#import fish
import pandas as pd
import bces
#import seaborn
import time
#import fancy_plots
#import nmmn.plots
#import plots
import misc
#import pylab
import xplot
import stats
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib.ticker import FormatStrFormatter
from matplotlib import rc
from matplotlib.ticker import LogLocator

plt.rcParams["font.family"] = "Times New Roman"

filename='New_ETG_LTG3.csv'
df=pd.read_csv(filename)


df=pd.read_csv(filename)

#df=df[df.ETG == 'no']
#df=df[df.galaxy != 'Cygnus A'] #trim outlier, Elliptical?
#df=df[df.galaxy != 'NGC 4395'] #trim outlier, bulgeless
#df=df[df.Galaxy != 'n5128']
#df=df[df.Galaxy != 'n1316']
#df=df[df.Galaxy != 'ic4296']
#df=df[df.Galaxy != 'm89']
#df=df[df.galaxy != 'NGC 4699']

df=df[df.in_out != 'o']
df=df[df.disk == 'E']  # E type
df=df[df.Gal_type == 'ETG']
#df=df[df.Bar == 'nb']


#names = ["names","y","perr","nerr","y_error","x","x_error"]
#data= np.genfromtxt("M_BH-M_sph.csv", dtype=None, delimiter=',', names=names)

#median = np.median(df['Log_Mass_gal_New'])
median = 10.699
print( 'Median = ')
print(median)

ydata=df['Log_Mass_BH']
erry=df['Err_Log_Mbh_dex']

xdata=df['Log_Sph_mass_new']-median
errx=df['Err_Log_mass_sph_quad']

#xdata=df['Log_Gal_mass_new']-median
#errx=df['Err_log_mass_gal_dex']

labels=df['Galaxy']
cov = np.zeros(len(xdata))

a,b,aerr,berr,covab=bces.bces(xdata,errx,ydata,erry,cov)

delta0=0.0
delta0 = delta0 + (ydata - (a[0]*xdata + b[0]))**2
delta0 = np.sqrt(np.sum(delta0)/len(xdata))

delta1=0.0
delta1 = delta1 + (ydata - (a[1]*xdata + b[1]))**2
delta1 = np.sqrt(np.sum(delta1)/len(xdata))

delta2=0.0
delta2 = delta2 + (ydata - (a[2]*xdata + b[2]))**2
delta2 = np.sqrt(np.sum(delta2)/len(xdata))

delta3=0.0
delta3 = delta3 + (ydata - (a[3]*xdata + b[3]))**2
delta3 = np.sqrt(np.sum(delta3)/len(xdata))

print('\nWITHOUT BOOTSTRAPPING')
print('E-Type Galaxies')
print ('\nNumber Galaxies = ')
print (len(xdata))
print ('\nFitting form: Y=AX+B.')
print ('Output = A, A error, B, B error, cov(AB), RMS Scatter\n')
print ('OLS(Y|X)')
print (a[0], aerr[0], b[0], berr[0], covab[0], delta0)
print ('')
print ('OLS(X|Y)')
print (a[1], aerr[1], b[1], berr[1], covab[1], delta1)
print ('')
print ('Bisector')
print (a[2], aerr[2], b[2], berr[2], covab[2], delta2)
print ('')
print ('Orthogonal')
print (a[3], aerr[3], b[3], berr[3], covab[3], delta3)
print ('')

i=0
r,rho,rchisq,scat,iscat,r2=stats.fitstats(xdata,ydata,errx,erry,a[i],b[i])
#r,rho,rchisq,scat,iscat,r2=stats.fitstats(xdata,ydata,errx,erry,1.3026649175678038,7.1939262508806125)

print ('BCES (Y|X)\n')

print ('(Pearson r, p-value) =')
print( r)
print ('(Spearman rho, p-value) =')
print (rho)
print ('Reduced chi-squared = ')
print (rchisq)
print ('Orthogonal Raw scatter = ')
print (scat)
print ('(Orthogonal Median intrinsic scatter, MAD) = ')
print (iscat)
print ('Coefficient of determination R^2 = ')
print (r2)

sd=stats.scatterfit(xdata,ydata,a[i],b[i])
#sd=stats.scatterfit(xdata,ydata,1.3026649175678038,7.1939262508806125)

print ('Traditional scatter = ')
print (sd)

isd=stats.myintscat(xdata,ydata,errx,erry,a[i],b[i])
#isd=stats.myintscat(xdata,ydata,errx,erry,1.3026649175678038,7.1939262508806125)

print ('Traditional intrinsic scatter = ')
print (isd)

i=1
r,rho,rchisq,scat,iscat,r2=stats.fitstats(xdata,ydata,errx,erry,a[i],b[i])
#r,rho,rchisq,scat,iscat,r2=stats.fitstats(xdata,ydata,errx,erry,1.3026649175678038,7.1939262508806125)

print ('\nBCES (X|Y)\n')

print ('(Pearson r, p-value) =')
print( r)
print ('(Spearman rho, p-value) =')
print (rho)
print ('Reduced chi-squared = ')
print (rchisq)
print ('Orthogonal Raw scatter = ')
print (scat)
print ('(Orthogonal Median intrinsic scatter, MAD) = ')
print (iscat)
print ('Coefficient of determination R^2 = ')
print (r2)

sd=stats.scatterfit(xdata,ydata,a[i],b[i])
#sd=stats.scatterfit(xdata,ydata,1.3026649175678038,7.1939262508806125)

print ('Traditional scatter = ')
print (sd)

isd=stats.myintscat(xdata,ydata,errx,erry,a[i],b[i])
#isd=stats.myintscat(xdata,ydata,errx,erry,1.3026649175678038,7.1939262508806125)

print ('Traditional intrinsic scatter = ')
print (isd)

i=2
r,rho,rchisq,scat,iscat,r2=stats.fitstats(xdata,ydata,errx,erry,a[i],b[i])

print ('\nBCES Bisector\n')

print ('(Pearson r, p-value) =')
print (r)
print ('(Spearman rho, p-value) =')
print (rho)
print ('Reduced chi-squared = ')
print (rchisq)
print ('Raw scatter = ')
print (scat)
print ('(Median intrinsic scatter, MAD) = ')
print (iscat)
print ('Coefficient of determination R^2 = ')
print (r2)

sd=stats.scatterfit(xdata,ydata,a[i],b[i])
#sd=stats.scatterfit(xdata,ydata,1.800237,7.235332)

print ('Traditional scatter = ')
print (sd)

isd=stats.myintscat(xdata,ydata,errx,erry,a[i],b[i])
#isd=stats.myintscat(xdata,ydata,errx,erry,1.800237,7.235332)

print ('Traditional intrinsic scatter = ')
print (isd)

i=2

lcb,ucb,xcb=stats.confband(xdata,ydata,a[i],b[i],conf=0.682689492,x=np.linspace(-15,15,100))
plt.fill_between(xcb+median, lcb, ucb, alpha=0.3, facecolor='red')
lpb,upb,xpb=stats.predband(xdata,ydata,a[i],b[i],conf=0.682689492,x=np.linspace(-15,15,100))
plt.fill_between(xpb+median, lpb, upb, alpha=0.10, facecolor='red')


ybces=a[i]*(xdata)+b[i]
####################################################################
####################################################################


df2=pd.read_csv(filename)
#df2=df2[df2.Galaxy != 'NGC 7052']
#df2=df2[df2.Core2 == 'yes'] #trim core-Sersic

df2=df2[df2.in_out != 'o']
df2=df2[df2.Gal_type == 'ETG']
df2=df2[df2.disk == 'S'] #ES

median2 = 10.699
print( 'Median = ')
print(median2)

ydata2=df2['Log_Mass_BH']
erry2=df2['Err_Log_Mbh_dex']

xdata2=df2['Log_Sph_mass_new']-median
errx2=df2['Err_Log_mass_sph_quad']

#xdata2=df2['Log_Gal_mass_new']-median
#errx2=df2['Err_log_mass_gal_dex']

labels2=df2['Galaxy']
cov2 = np.zeros(len(xdata2))

## number of bootstrapping trials
#nboot=2**15 #1000000 takes 305.585825 s USE 2**20

##time
## Performs the BCES fit (parallel = .bcesp)
#a,b,aerr,berr,covab=bces.bcesp(xdata2,errx2,ydata2,erry2,cov2,nboot)

a,b,aerr,berr,covab=bces.bces(xdata2,errx2,ydata2,erry2,cov2)

delta0=0.0
delta0 = delta0 + (ydata2 - (a[0]*xdata2 + b[0]))**2
delta0 = np.sqrt(np.sum(delta0)/len(xdata2))

delta1=0.0
delta1 = delta1 + (ydata2 - (a[1]*xdata2 + b[1]))**2
delta1 = np.sqrt(np.sum(delta1)/len(xdata2))

delta2=0.0
delta2 = delta2 + (ydata2 - (a[2]*xdata2 + b[2]))**2
delta2 = np.sqrt(np.sum(delta2)/len(xdata2))

delta3=0.0
delta3 = delta3 + (ydata2 - (a[3]*xdata2 + b[3]))**2
delta3 = np.sqrt(np.sum(delta3)/len(xdata2))

print('\nWITHOUT BOOTSTRAPPING')
print('ES-Type Galaxies')
print ('\nNumber Galaxies = ')
print (len(xdata2))
print ('\nFitting form: Y=AX+B.')
print ('Output = A, A error, B, B error, cov(AB), RMS Scatter\n')
print ('OLS(Y|X)')
print (a[0], aerr[0], b[0], berr[0], covab[0], delta0)
print ('')
print ('OLS(X|Y)')
print (a[1], aerr[1], b[1], berr[1], covab[1], delta1)
print ('')
print ('Bisector')
print (a[2], aerr[2], b[2], berr[2], covab[2], delta2)
print ('')
print ('Orthogonal')
print (a[3], aerr[3], b[3], berr[3], covab[3], delta3)
print ('')

i=0
r,rho,rchisq,scat,iscat,r2=stats.fitstats(xdata2,ydata2,errx2,erry2,a[i],b[i])
#r,rho,rchisq,scat,iscat,r2=stats.fitstats(xdata,ydata,errx,erry,1.3026649175678038,7.1939262508806125)

print ('BCES (Y|X)\n')

print ('(Pearson r, p-value) =')
print( r)
print ('(Spearman rho, p-value) =')
print (rho)
print ('Reduced chi-squared = ')
print (rchisq)
print ('Orthogonal Raw scatter = ')
print (scat)
print ('(Orthogonal Median intrinsic scatter, MAD) = ')
print (iscat)
print ('Coefficient of determination R^2 = ')
print (r2)

sd=stats.scatterfit(xdata2,ydata2,a[i],b[i])
#sd=stats.scatterfit(xdata,ydata,1.3026649175678038,7.1939262508806125)

print ('Traditional scatter = ')
print (sd)

isd=stats.myintscat(xdata2,ydata2,errx2,erry2,a[i],b[i])
#isd=stats.myintscat(xdata,ydata,errx,erry,1.3026649175678038,7.1939262508806125)

print ('Traditional intrinsic scatter = ')
print (isd)

i=1
r,rho,rchisq,scat,iscat,r2=stats.fitstats(xdata2,ydata2,errx2,erry2,a[i],b[i])
#r,rho,rchisq,scat,iscat,r2=stats.fitstats(xdata,ydata,errx,erry,1.3026649175678038,7.1939262508806125)

print ('\nBCES (X|Y)\n')

print ('(Pearson r, p-value) =')
print( r)
print ('(Spearman rho, p-value) =')
print (rho)
print ('Reduced chi-squared = ')
print (rchisq)
print ('Orthogonal Raw scatter = ')
print (scat)
print ('(Orthogonal Median intrinsic scatter, MAD) = ')
print (iscat)
print ('Coefficient of determination R^2 = ')
print (r2)

sd=stats.scatterfit(xdata2,ydata2,a[i],b[i])
#sd=stats.scatterfit(xdata,ydata,1.3026649175678038,7.1939262508806125)

print ('Traditional scatter = ')
print (sd)

isd=stats.myintscat(xdata2,ydata2,errx2,erry2,a[i],b[i])
#isd=stats.myintscat(xdata,ydata,errx,erry,1.3026649175678038,7.1939262508806125)

print ('Traditional intrinsic scatter = ')
print (isd)

i=2
r,rho,rchisq,scat,iscat,r2=stats.fitstats(xdata2,ydata2,errx2,erry2,a[i],b[i])

print ('\nBCES Bisector\n')

print ('(Pearson r, p-value) =')
print (r)
print ('(Spearman rho, p-value) =')
print (rho)
print ('Reduced chi-squared = ')
print (rchisq)
print ('Raw scatter = ')
print (scat)
print ('(Median intrinsic scatter, MAD) = ')
print (iscat)
print ('Coefficient of determination R^2 = ')
print (r2)

sd=stats.scatterfit(xdata2,ydata2,a[i],b[i])
#sd=stats.scatterfit(xdata,ydata,1.800237,7.235332)

print ('Traditional scatter = ')
print (sd)

isd=stats.myintscat(xdata2,ydata2,errx2,erry2,a[i],b[i])
#isd=stats.myintscat(xdata,ydata,errx,erry,1.800237,7.235332)

print ('Traditional intrinsic scatter = ')
print (isd)



#Selector for which line to plot
i=2

lcb2,ucb2,xcb2=stats.confband(xdata2,ydata2,a[i],b[i],conf=0.682689492,x=np.linspace(-15,15,100))
plt.fill_between(xcb+median2, lcb2, ucb2, alpha=0.3, facecolor='green')
lpb2,upb2,xpb2=stats.predband(xdata2,ydata2,a[i],b[i],conf=0.682689492,x=np.linspace(-15,15,100))
plt.fill_between(xpb+median2, lpb2, upb2, alpha=0.10, facecolor='green')



ybces2=a[i]*(xdata2)+b[i]


##############################################################################################
##############################################################################################

df3=pd.read_csv(filename)
#df2=df2[df2.Galaxy != 'NGC 7052']
#df2=df2[df2.Core2 == 'yes'] #trim core-Sersic

df3=df3[df3.in_out != 'o']
df3=df3[df3.Gal_type == 'LTG'] #ES

median3 = 10.699
print( 'Median = ')
print(median3)

ydata3=df3['Log_Mass_BH']
erry3=df3['Err_Log_Mbh_dex']

xdata3=df3['Log_Sph_mass_new']-median
errx3=df3['Err_Log_mass_sph_quad']

#xdata3=df3['Log_Gal_mass_new']-median
#errx3=df3['Err_log_mass_gal_dex']

labels3=df3['Galaxy']
cov3 = np.zeros(len(xdata3))

## number of bootstrapping trials
#nboot=2**15 #1000000 takes 305.585825 s USE 2**20

##time
## Performs the BCES fit (parallel = .bcesp)
#a,b,aerr,berr,covab=bces.bcesp(xdata2,errx2,ydata2,erry2,cov2,nboot)

a,b,aerr,berr,covab=bces.bces(xdata3,errx3,ydata3,erry3,cov3)

delta0=0.0
delta0 = delta0 + (ydata3 - (a[0]*xdata3 + b[0]))**2
delta0 = np.sqrt(np.sum(delta0)/len(xdata3))

delta1=0.0
delta1 = delta1 + (ydata3 - (a[1]*xdata3 + b[1]))**2
delta1 = np.sqrt(np.sum(delta1)/len(xdata3))

delta2=0.0
delta2 = delta2 + (ydata3 - (a[2]*xdata3 + b[2]))**2
delta2 = np.sqrt(np.sum(delta2)/len(xdata3))

delta3=0.0
delta3 = delta3 + (ydata3 - (a[3]*xdata3 + b[3]))**2
delta3 = np.sqrt(np.sum(delta3)/len(xdata3))

print('\nWITHOUT BOOTSTRAPPING')
print('S0-Type Galaxies')
print ('\nNumber Galaxies = ')
print (len(xdata3))
print ('\nFitting form: Y=AX+B.')
print ('Output = A, A error, B, B error, cov(AB), RMS Scatter\n')
print ('OLS(Y|X)')
print (a[0], aerr[0], b[0], berr[0], covab[0], delta0)
print ('')
print ('OLS(X|Y)')
print (a[1], aerr[1], b[1], berr[1], covab[1], delta1)
print ('')
print ('Bisector')
print (a[2], aerr[2], b[2], berr[2], covab[2], delta2)
print ('')
print ('Orthogonal')
print (a[3], aerr[3], b[3], berr[3], covab[3], delta3)
print ('')

i=0
r,rho,rchisq,scat,iscat,r2=stats.fitstats(xdata3,ydata3,errx3,erry3,a[i],b[i])
#r,rho,rchisq,scat,iscat,r2=stats.fitstats(xdata,ydata,errx,erry,1.3026649175678038,7.1939262508806125)

print ('BCES (Y|X)\n')

print ('(Pearson r, p-value) =')
print( r)
print ('(Spearman rho, p-value) =')
print (rho)
print ('Reduced chi-squared = ')
print (rchisq)
print ('Orthogonal Raw scatter = ')
print (scat)
print ('(Orthogonal Median intrinsic scatter, MAD) = ')
print (iscat)
print ('Coefficient of determination R^2 = ')
print (r2)

sd=stats.scatterfit(xdata3,ydata3,a[i],b[i])
#sd=stats.scatterfit(xdata,ydata,1.3026649175678038,7.1939262508806125)

print ('Traditional scatter = ')
print (sd)

isd=stats.myintscat(xdata3,ydata3,errx3,erry3,a[i],b[i])
#isd=stats.myintscat(xdata,ydata,errx,erry,1.3026649175678038,7.1939262508806125)

print ('Traditional intrinsic scatter = ')
print (isd)

i=1
r,rho,rchisq,scat,iscat,r2=stats.fitstats(xdata3,ydata3,errx3,erry3,a[i],b[i])
#r,rho,rchisq,scat,iscat,r2=stats.fitstats(xdata,ydata,errx,erry,1.3026649175678038,7.1939262508806125)

print ('\nBCES (X|Y)\n')

print ('(Pearson r, p-value) =')
print( r)
print ('(Spearman rho, p-value) =')
print (rho)
print ('Reduced chi-squared = ')
print (rchisq)
print ('Orthogonal Raw scatter = ')
print (scat)
print ('(Orthogonal Median intrinsic scatter, MAD) = ')
print (iscat)
print ('Coefficient of determination R^2 = ')
print (r2)

sd=stats.scatterfit(xdata3,ydata3,a[i],b[i])
#sd=stats.scatterfit(xdata,ydata,1.3026649175678038,7.1939262508806125)

print ('Traditional scatter = ')
print (sd)

isd=stats.myintscat(xdata3,ydata3,errx3,erry3,a[i],b[i])
#isd=stats.myintscat(xdata,ydata,errx,erry,1.3026649175678038,7.1939262508806125)

print ('Traditional intrinsic scatter = ')
print (isd)

i=2
r,rho,rchisq,scat,iscat,r2=stats.fitstats(xdata3,ydata3,errx3,erry3,a[i],b[i])

print ('\nBCES Bisector\n')

print ('(Pearson r, p-value) =')
print (r)
print ('(Spearman rho, p-value) =')
print (rho)
print ('Reduced chi-squared = ')
print (rchisq)
print ('Raw scatter = ')
print (scat)
print ('(Median intrinsic scatter, MAD) = ')
print (iscat)
print ('Coefficient of determination R^2 = ')
print (r2)

sd=stats.scatterfit(xdata3,ydata3,a[i],b[i])
#sd=stats.scatterfit(xdata,ydata,1.800237,7.235332)

print ('Traditional scatter = ')
print (sd)

isd=stats.myintscat(xdata3,ydata3,errx3,erry3,a[i],b[i])
#isd=stats.myintscat(xdata,ydata,errx,erry,1.800237,7.235332)

print ('Traditional intrinsic scatter = ')
print (isd)



#Selector for which line to plot
i=2

lcb3,ucb3,xcb3=stats.confband(xdata3,ydata3,a[i],b[i],conf=0.682689492,x=np.linspace(-15,15,100))
plt.fill_between(xcb+median3, lcb3, ucb3, alpha=0.3, facecolor='blue')
lpb3,upb3,xpb3=stats.predband(xdata3,ydata3,a[i],b[i],conf=0.682689492,x=np.linspace(-15,15,100))
plt.fill_between(xpb+median3, lpb3, upb3, alpha=0.10, facecolor='blue')

#lpb,upb,xpb=stats.predband(xdata,ydata,a[i],b[i],conf=0.8664,x=np.linspace(-15,15,100))
#plt.fill_between(xpb+median, lpb, upb, alpha=0.3, facecolor='red')
#lpb,upb,xpb=stats.predband(xdata,ydata,a[i],b[i],conf=0.954499736,x=np.linspace(-15,15,100))
#plt.fill_between(xpb+median, lpb, upb, alpha=0.3, facecolor='yellow')

ybces3=a[i]*(xdata3)+b[i]


#############################################################################################
df4=pd.read_csv(filename)

df4=df4[df4.Gal_type == 'ETG']
df4=df4[df4.in_out == 'o']

#df4=df4[df4.Core2 == 'no']

ydata4=df4['Log_Mass_BH']
erry4=df4['Err_Log_Mbh_dex']

xdata4=df4['Log_Sph_mass_new']-median
errx4=df4['Err_Log_mass_sph_quad']

#xdata4=df4['Log_Gal_mass_new']-median
#errx4=df4['Err_log_mass_gal_dex']

labels4=df4['Galaxy']

ybces4=a[i]*(xdata4)+b[i]
#############################################################################################

#-----------#-----------#-----------#-----------#-----------#-----------#-----------
Log_Msph_ETG1= np.linspace(8.50,12.50,100)
Log_Msph_ETG=pd.Series(Log_Msph_ETG1)
Log_Mbh_ETG= 1.27*(Log_Msph_ETG-np.log10(5e10))+8.41

#############################################################################################

#plt.plot(xdata+median,ybces,color='black',linestyle='-',zorder=15,linewidth=2,label='ETG (E+ES+S0)')

plt.plot(xdata+median,ybces,color='red',linestyle='-',zorder=15,linewidth=2,label='E-Type')
plt.plot(np.linspace(min(xdata2+median2),max(xdata2+median2)),np.linspace(min(ybces2),max(ybces2)),color='green',linestyle='-',zorder=15,linewidth=2,label='(ES+S0)-Type')

plt.plot(np.linspace(min(xdata3+median3),max(xdata3+median3)),np.linspace(min(ybces3),max(ybces3)),color='blue',linestyle='-',zorder=15,linewidth=2,label='LTG (S)')


plt.xlabel(r'$\log({M_{\rm *,sph}/{\rm M}_{\odot}})$', fontsize=22)
#plt.xlabel(r'$\log({M_{\rm *,gal}/{\rm M}_{\odot}})$', fontsize=22)
plt.ylabel(r'$\log({M_{\rm BH}/{\rm M}_{\odot}})$', fontsize=22)

#plt.legend(loc='best', fontsize=13)
plt.tick_params(axis='x', labelsize=24,direction='in')
plt.tick_params(axis='y', labelsize=24,direction='in')

ax=plt.gca()
#plt.tick_params(labeltop=True, labelright=True)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')

#ax.set_ylim([4.7,11.30])
#ax.set_xlim([8.4,12.7])
#ax.set_ylim([3.60,10.99])
#ax.set_xlim([8.8,12.7])
#ax.set_xlim([7.6,12.7])

#new range
#ax.set_xlim([9.14,12.26])
#ax.set_ylim([5.72,10.78])

ax.set_xlim([8.84,12.69]) #Msph
#ax.set_xlim([8.44,12.7])#Msph
#ax.set_xlim([9.24,12.29]) #Mgal 
ax.set_ylim([5.42,10.82]) 


fig=plt.gcf()
fig.subplots_adjust(bottom=0.15)
fig.set_size_inches(11,8)

plt.errorbar(xdata+median,ydata,xerr=errx,yerr=erry,ls='None',linewidth=0.5,ecolor='r',zorder=20,mew=0.5,label=None)
#plt.errorbar(xdata+median,ydata,xerr=errx,yerr=erry,ls='None',linewidth=0.2,ecolor='r',zorder=2,mew=0.2,label=None)
plt.plot(xdata+median,ydata,'r^',markersize=5,zorder=25,label='E-Type')   

plt.errorbar(xdata2+median2,ydata2,xerr=errx2,yerr=erry2,ls='None',linewidth=0.5,ecolor='green',zorder=20,mew=0.5,label=None)
plt.plot(xdata2+median2,ydata2,'gs',markersize=5,zorder=25,label='(ES+S0)-Type')

#plt.errorbar(xdata+median,ydata,xerr=errx,yerr=erry,ls='None',linewidth=0.2,ecolor='r',zorder=2,mew=0.2,label=None)
#plt.plot(xdata+median,ydata,'ro',markersize=5,zorder=25,label='E-Type')   

#plt.errorbar(xdata2+median2,ydata2,xerr=errx2,yerr=erry2,ls='None',linewidth=0.2,ecolor='b',zorder=2,mew=0.2,label=None)
#plt.plot(xdata2+median2,ydata2,'b^',markersize=5,zorder=25,label='(ES+S0)-Type')


plt.errorbar(xdata3+median3,ydata3,xerr=errx3,yerr=erry3,ls='None',linewidth=0.5,ecolor='b',zorder=20,mew=0.5,label=None)
plt.plot(xdata3+median3,ydata3,'bo',markersize=5,zorder=25,label='S-Type')

plt.errorbar(xdata4+median,ydata4,xerr=errx4,yerr=erry4,ls='None',linewidth=0.5,ecolor='k',zorder=20,mew=0.5,label=None)
plt.plot(xdata4+median,ydata4,'k*',markersize=12,zorder=25, label= 'Excluded')


#plt.plot(Log_Msph_ETG,Log_Mbh_ETG,color='black',linestyle=':',zorder=7,linewidth=1,label='ETG (E+ES+S0)')

plt.legend(loc='lower right', fontsize=16)


#for label, x, y in zip(labels, xdata+median, ydata):
#   plt.annotate(
#        label,
#       xy=(x, y), xytext=(-10, 125),
#        textcoords='offset points', ha='right', va='bottom',
#        bbox=dict(boxstyle='round,pad=0.4', fc='yellow', alpha=0.4),
#       arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'),fontsize=7)
#for label, x, y in zip(labels2, xdata2+median2, ydata2):
#   plt.annotate(
#        label,
#       xy=(x, y), xytext=(-10, 125),
#        textcoords='offset points', ha='right', va='bottom',
#        bbox=dict(boxstyle='round,pad=0.4', fc='yellow', alpha=0.4),
#       arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'),fontsize=7)  
for label, x, y in zip(labels4, xdata4+median, ydata4):
   plt.annotate(label,xy=(x, y), xytext=(-5, 0), zorder=30,
        textcoords='offset points', ha='right', va='bottom', 
        bbox=dict(boxstyle='round,pad=0.3', fc='yellow', alpha=0.3),
         fontsize=7)

plt.show()
