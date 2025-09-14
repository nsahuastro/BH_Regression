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

##############################################################################

########################### ETG ###########################

filename='New_ETG_LTG3.csv'

df=pd.read_csv(filename)

#df=df[df.ETG == 'no']

df=df[df.Galaxy != 'n5128']
df=df[df.Galaxy != 'n1316']
df=df[df.Galaxy != 'n4486B']
df=df[df.Galaxy != 'n404']
df=df[df.Galaxy != 'n4342']
df=df[df.Galaxy != 'n2787']
df=df[df.Galaxy != 'n1277']




df=df[df.Gal_type == 'ETG'] 

#median = np.median(df['Log_Mass_sph_New'])
median = 10.699
print( 'Median = ')
print(median)


ydata=df['Log_Mass_BH']

#xdata=df['Log_Sph_mass_new']-median
#errx=df['Err_Log_mass_sph_quad']

xdata=df['Log_Gal_mass_new']-median
errx=df['Err_log_mass_gal_dex']



erry=df['Err_Log_Mbh_dex']

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




#####################################################


################measurement error################
#k=[0,0,0,0]
#k1=[0,0,0,0]

#k[0]= erry**2 + (a[0]*errx)**2
#k1[0]= np.sqrt(sum(k[0])/len(xdata)) 
 
#k[1]= erry**2 + (a[1]*errx)**2
#k1[1]= np.sqrt(sum(k[1])/len(xdata)) 
 
#k[2]= erry**2 + (a[2]*errx)**2
#k1[2]= np.sqrt(sum(k[2])/len(xdata)) 
 
#k[3]= erry**2 + (a[3]*errx)**2
#k1[3]= np.sqrt(sum(k[3])/len(xdata)) 
 

################intrinsic scatter################
#k2=[0,0,0,0]

#k2[0]= np.sqrt((delta0)**2 - k1[0]**2) 
#k2[1]= np.sqrt((delta1)**2 - k1[1]**2) 
#k2[2]= np.sqrt((delta2)**2 - k1[2]**2) 
#k2[3]= np.sqrt((delta3)**2 - k1[3]**2) 
 
###############################################



print('\nWITHOUT BOOTSTRAPPING')
print('ETGs')
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

#isd=stats.myintscat(xdata,ydata,errx,erry,a[i],b[i])
#print ('Traditional intrinsic scatter = ')
#print (isd)
isd2=stats.fitexyiscat(xdata,ydata,errx,erry,a[i],b[i])
print ('Intrinsic scatter = ')
print (isd2)


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


print ('Traditional scatter = ')
print (sd)

#isd=stats.myintscat(xdata,ydata,errx,erry,a[i],b[i])
#print ('Traditional intrinsic scatter = ')
#print (isd)

isd2=stats.fitexyiscat(xdata,ydata,errx,erry,a[i],b[i])
print ('Intrinsic scatter = ')
print (isd2)


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

print ('Traditional scatter = ')
print (sd)

#isd=stats.myintscat(xdata,ydata,errx,erry,a[i],b[i])

#print ('Traditional intrinsic scatter = ')
#print (isd)


isd2=stats.fitexyiscat(xdata,ydata,errx,erry,a[i],b[i])
print ('Intrinsic scatter = ')
print (isd2)


i=2

lcb,ucb,xcb=stats.confband(xdata,ydata,a[i],b[i],conf=0.682689492,x=np.linspace(-15,15,100))
plt.fill_between(xcb+median, lcb, ucb, alpha=0.3, facecolor='red')
lpb,upb,xpb=stats.predband(xdata,ydata,a[i],b[i],conf=0.682689492,x=np.linspace(-15,15,100))
plt.fill_between(xpb+median, lpb, upb, alpha=0.3, facecolor='red')


#x_asymmetric_error=[nerrx,perrx]
#y_asymmetric_error=[nerry,perry]

ybces=a[i]*(xdata)+b[i]



#######################################################################################

################### LTG ###################


df2=pd.read_csv(filename)

df2=df2[df2.Galaxy != 'LEDA 87300']

df2=df2[df2.Gal_type == 'LTG'] 


#median2 = np.median(df2['Log_Mass_sph_New'])

median2 = 10.699
print( 'Median = ')
print(median2)


ydata2=df2['Log_Mass_BH']

#xdata2=df2['Log_Sph_mass_new']-median2
#errx2=df2['Err_Log_mass_sph_quad']

xdata2=df2['Log_Gal_mass_new']-median2
errx2=df2['Err_log_mass_gal_dex']


erry2=df2['Err_Log_Mbh_dex']

labels2=df2['Galaxy']
cov2 = np.zeros(len(xdata2))



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





################measurement error################
#k=[0,0,0,0]
#k1=[0,0,0,0]

#k[0]= erry2**2 + (a[0]*errx2)**2
#k1[0]= np.sqrt(sum(k[0])/len(xdata2)) 
 
#k[1]= erry2**2 + (a[1]*errx2)**2
#k1[1]= np.sqrt(sum(k[1])/len(xdata2)) 
 
#k[2]= erry2**2 + (a[2]*errx2)**2
#k1[2]= np.sqrt(sum(k[2])/len(xdata2)) 
 
#k[3]= erry2**2 + (a[3]*errx2)**2
#k1[3]= np.sqrt(sum(k[3])/len(xdata2)) 
 

################intrinsic scatter################
#k2=[0,0,0,0]



#k2[0]= np.sqrt((delta0)**2 - k1[0]**2) 
#k2[1]= np.sqrt((delta1)**2 - k1[1]**2) 
#k2[2]= np.sqrt((delta2)**2 - k1[2]**2)
#k2[3]= np.sqrt((delta3)**2 - k1[3]**2) 
 
###############################################




print('\nWITHOUT BOOTSTRAPPING')
print('LTGs')
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
print ('Traditional scatter = ')
print (sd)

#isd=stats.myintscat(xdata2,ydata2,errx2,erry2,a[i],b[i])
#print ('Traditional intrinsic scatter = ')
#print (isd)

isd2=stats.fitexyiscat(xdata2,ydata2,errx2,erry2,a[i],b[i])
print ('Intrinsic scatter = ')
print (isd2)



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
print ('Traditional scatter = ')
print (sd)

#isd=stats.myintscat(xdata2,ydata2,errx2,erry2,a[i],b[i])
#print ('Traditional intrinsic scatter = ')
#print (isd)

isd2=stats.fitexyiscat(xdata2,ydata2,errx2,erry2,a[i],b[i])
print ('Intrinsic scatter = ')
print (isd2)


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
print ('Traditional scatter = ')
print (sd)

#isd=stats.myintscat(xdata2,ydata2,errx2,erry2,a[i],b[i])
#print ('Traditional intrinsic scatter = ')
#print (isd)

isd2=stats.fitexyiscat(xdata2,ydata2,errx2,erry2,a[i],b[i])
print ('Intrinsic scatter = ')
print (isd2)



#Selector for which line to plot
i=2

lcb2,ucb2,xcb2=stats.confband(xdata2,ydata2,a[i],b[i],conf=0.682689492,x=np.linspace(-15,15,100))
plt.fill_between(xcb2+median2, lcb2, ucb2, alpha=0.3, facecolor='blue')
lpb2,upb2,xpb2=stats.predband(xdata2,ydata2,a[i],b[i],conf=0.682689492,x=np.linspace(-15,15,100))
plt.fill_between(xpb2+median2, lpb2, upb2, alpha=0.3, facecolor='blue')
#lpb,upb,xpb=stats.predband(xdata,ydata,a[i],b[i],conf=0.8664,x=np.linspace(-15,15,100))
#pylab.fill_between(xpb+median, lpb, upb, alpha=0.3, facecolor='red')
#lpb,upb,xpb=stats.predband(xdata,ydata,a[i],b[i],conf=0.954499736,x=np.linspace(-15,15,100))
#pylab.fill_between(xpb+median, lpb, upb, alpha=0.3, facecolor='yellow')

#x_asymmetric_error2=[nerrx2,perrx2]
#y_asymmetric_error2=[nerry2,perry2]



ybces2=a[i]*(xdata2)+b[i]
#plt.plot(xdata+15,ybces0,label='BCES (Y|X)',color='yellow')
#ymarc=8.44-0.076*(xdata+15)
#yjoel=8.21-0.062*(xdata+15)
#plt.plot(np.linspace(min(xdata+15),max(xdata+15),100),np.linspace(max(ymarc),min(ymarc),100),label='Seigar et al. (2008)',color='magenta',zorder=5,linewidth=3,linestyle=':')
#plt.plot(np.linspace(min(xdata+15),max(xdata+15),100),np.linspace(max(yjoel),min(yjoel),100),label='Berrier et al. (2013)',color='cyan',zorder=10,linewidth=3,linestyle='--')

######################################################################################################


#### %%%%%%  Excluded From Regression %%%%%% #####



df3=pd.read_csv(filename)

df3=df3[df3.Gal_type == 'ETG']
df3=df3[df3.in_out == 'o']

#df3=df3[df3.Core2 == 'no']

ydata3=df3['Log_Mass_BH']

erry3=df3['Err_Log_Mbh_dex']


#xdata3=df3['Log_Sph_mass_new']-median
#errx3=df3['Err_Log_mass_sph_quad']

xdata3=df3['Log_Gal_mass_new']-median
errx3=df3['Err_log_mass_gal_dex']


labels3=df3['Galaxy']

ybces3=a[i]*(xdata3)+b[i]

#######################################################################################################
######################## %%%%%%%% Peculiar Galaxies Sersic %%%%% ################################

df4=pd.read_csv(filename)


df4=df4[df4.Gal_type == 'ETG']
df4=df4[df4.Pec == 'pec']

#df4=df4[df4.Core == 'no']

ydata4=df4['Log_Mass_BH']

#xdata4=df4['Log_Sph_mass_new']-median
#errx4=df4['Err_Log_mass_sph_quad']

xdata4=df4['Log_Gal_mass_new']-median
errx4=df4['Err_log_mass_gal_dex']


erry4=df4['Err_Log_Mbh_dex']

labels4=df4['Galaxy']

ybces4=a[i]*(xdata4)+b[i]

######################################################################################################

#df5=pd.read_csv(filename)

#df5=df5[df5.bulge == 'c']
#df5=df5[df5.Galaxy == 'NGC 224']
#df5=df5[df5.Galaxy == 'NGC 1398']
#df5=df5[df5.Galaxy == 'NGC 2974']
#df5=df5[df5.Galaxy == 'NGC 3031']
#df5=df5[df5.Galaxy == 'NGC 4151']

#ydata5=df5['Log_Mass_BH']

#xdata5=df5['Log_Sph_mass_new']-median2

#errx5=df5['Err_Log_mass_sph_quad']

#errx5=df5['Err_log_mass_gal_dex']


#erry5=df5['Err_Log_Mbh_dex']

#labels5=df5['Galaxy']

#print(xdata5)
#ybces5=a[i]*(xdata5)+b[i]



###################################################################################

#filename2='s16.csv'
#filename2='Sh16_3.csv'
#df5=pd.read_csv(filename2)

#xdata5= df5['x'] 
#xdata5= df5['x_cor']
#ydata5= df5['y_cor'] 


#################################################################################

#filename2='Sh16_3.csv'
#df6=pd.read_csv(filename2)

#xdata5= df5['x'] 
#xdata6= df6['x_cor']
#ydata6= df6['y'] 

######################################################################################################


plt.plot(xdata+median,ybces,color='red',zorder=15,linewidth=3,label='ETG')
plt.plot(np.linspace(min(xdata2+median2),max(xdata2+median2)),np.linspace(min(ybces2),max(ybces2)),color='blue',ls='--', zorder=15,linewidth=3,label='LTG')

#plt.plot(xdata5,ydata5,color='black',zorder=15,linewidth=3,label='Shankar16_cor')
#plt.plot(xdata6,ydata6,color='c',zorder=15,linewidth=3,label='Shankar16')

plt.xlabel(r'$\log({M_{\rm *,gal}/{\rm M}_{\odot}})$', fontsize=24)
plt.ylabel(r'$\log({M_{\rm BH}/{\rm M}_{\odot}})$', fontsize=24)

#plt.xlabel(r'$\log(\, \rm Spheroid \, Stellar \, Mass \,)$', fontsize=24)
#plt.ylabel(r'$\log(\, \rm Black \, Hole \, Mass \,)$', fontsize=24)

plt.tick_params(axis='x', labelsize=24, direction='in')
plt.tick_params(axis='y', labelsize=24,direction='in')

ax=plt.gca()

#plt.tick_params(labeltop=True, labelright=True)

ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')

ax.set_ylim([3.60,10.99])
ax.set_xlim([8.8,12.7]) #Mgal
#ax.set_xlim([7.6,12.7])

fig=plt.gcf()
fig.subplots_adjust(bottom=0.15)
fig.set_size_inches(11,8)

#plt.errorbar(xdata+median,ydata,xerr=errx,yerr=erry,ls='None',linewidth=0.5,ecolor='r',zorder=20,mew=0.5,label=None)
#plt.plot(xdata+median,ydata,'r^',markersize=7,zorder=25,label='ETG')

#plt.errorbar(xdata2+median2,ydata2,xerr=errx2,yerr=erry2,ls='None',linewidth=0.5,ecolor='b',zorder=20,mew=0.5,label=None)
#plt.plot(xdata2+median2,ydata2,'bs',markersize=6,zorder=25,label='LTG')

#plt.errorbar(xdata4+median,ydata4,xerr=errx4,yerr=erry4,ls='None',linewidth=0.5,ecolor='c',zorder=20,mew=0.5,label=None)
#plt.plot(xdata4+median,ydata4,'co',markersize=8,zorder=25, label= 'Pec.ETG')

##plt.errorbar(xdata5+median2,ydata5,xerr=errx5,yerr=erry5,ls='None',linewidth=0.5,ecolor='m',zorder=20,mew=0.5,label=None)
##plt.plot(xdata5+median2,ydata5,'m*',markersize=12,zorder=25, label='Pec.Core-Se'+u'\u0301''rsic')


#plt.errorbar(xdata3+median,ydata3,xerr=errx3,yerr=erry3,ls='None',linewidth=0.5,ecolor='k',zorder=20,mew=0.5,label=None)
#plt.plot(xdata3+median,ydata3,'k*',markersize=12,zorder=25, label= 'Excluded_ETG')

#plt.errorbar(xdata5+median2,ydata5,xerr=errx5,yerr=erry5,ls='None',linewidth=0.5,ecolor='m',zorder=20,mew=0.5,label=None)
#plt.plot(xdata5+median2,ydata5,'m*',markersize=12,zorder=35, label= 'Classical_bulge')


plt.legend(loc='lower right', fontsize=15)


#for label, x, y in zip(labels3, xdata3+median, ydata3):
#   plt.annotate(label,xy=(x, y), xytext=(-5, 0), zorder=30,
#        textcoords='offset points', ha='right', va='bottom', 
#        bbox=dict(boxstyle='round,pad=0.3', fc='yellow', alpha=0.3),
#         fontsize=7)




plt.show()
