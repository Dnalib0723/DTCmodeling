#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 15:02:13 2022

@author: cherri
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
#from scipy.integrate import odeint
from scipy.interpolate import interp1d
import math
#model wild type


# import the steady-state profile simulated under constant LIN-42 levels

X_Lin42=np.array([9.06666666666667,18.1333333333333,27.2,36.2666666666667, \
                  45.3333333333333,54.4,63.4666666666667,72.5333333333333,\
                81.6,90.6666666666667,	99.7333333333333,108.8,\
                117.866666666667,126.933333333333,136,145.066666666667,	\
                154.133333333333,163.2,172.266666666667,181.333333333333,\
                190.4,199.466666666667,	208.533333333333,217.6,226.666666666667,\
                235.733333333333,244.8,253.866666666667,	262.933333333333,272])
      
 
# steady state UNC-5 level simulated with a HIGH UNC-5 level in Circuit I  
Y_unc5_as_g=np.array([187.724883784531,	187.673021377725,187.527610618038,\
      187.225766457742,186.676582601996,185.742438835436,\
      184.218217302627,181.817882228353,178.180990681039,\
      172.903909243840,165.590310347614,155.914309436989,\
      143.694283395829,128.974670205889,111.919450919407,\
      88.7202510705562,36.8424443385704,33.1320524347998,\
      30.4280599248321,28.4232939320258,26.9161499310741,\
      25.7682855974081,24.8828743745947,24.1913447593757,\
      23.6446642217679,23.2074375361411,22.8538613839326,\
      22.5649323990428,22.3265089524354,22.1279569965509])
    
 
    

# steady state UNC-5 level simulated with a LOW UNC-5 level in Circuit I  
Y_unc5_de_r=np.array([187.724883784531,187.673021377725,187.527610618038,\
     187.225766457741,186.676582601996,185.742438835436,\
    184.218217302627,181.817882228353,178.180990681039,\
    172.903909243840,165.590310347614,155.914309436989,\
    143.694283395829,128.974670205887,50.1370116367671,\
    42.0854986221808,36.8424443385696,33.1320524347992,\
    30.4280599248318,28.4232939320255,26.9161499310738,\
    25.7682855974078,24.8828743745945,24.1913447593755,\
    23.6446642217677,23.2074375361408,22.8538613839323,\
    22.5649323990426,22.3265089524351,22.1279569965506])
 
      
 
# steady state UNC-5 level simulated with a HIGH UNC-5 level in Circuit II   
Y_unc5_noAC_as_g=np.array([180.007096397800,179.821950014411,179.312772856926,\
                           178.292311622865,176.523736016318,173.684708887614,\
                          169.332818364685,162.874476596037,153.463293736831,\
                          138.871546345664,82.0966629674532,67.8566977844166,\
                          57.4691571952955,49.3301930435774,42.9266174170897,\
                          37.9469383991530,34.1235447703947,31.2133141668375,\
                          29.0058060218587,27.3292379631204,26.0496412057633,\
                          25.0656919212357,24.3021538732462,23.7036797676562,\
                          23.2296701741846,22.8502964996226,22.5435425943005,\
                         22.2930514634669,22.0865752674345,21.9148643391358])
    
 
# steady state UNC-5 level simulated with a LOW UNC-5 level in Circuit II
Y_unc5_noAC_de_r=np.array([180.007096397799,179.821950014410,179.312772856925,\
                           178.292311622864,176.523736016317,173.684708887612,\
                           169.332818364683,162.874476596034,153.463293736830,\
                           138.871546345662,82.0966629674514,67.8566977844149,\
                           57.4691571952946,49.3301930435766,42.9266174170892,\
                           37.9469383991527,34.1235447703944,31.2133141668372,\
                           29.0058060218584,27.3292379631201,26.0496412057630,\
                           25.0656919212354,24.3021538732460,23.7036797676560,\
                           23.2296701741843,22.8502964996223,22.5435425943004,\
                          22.2930514634668,22.0865752674342,21.9148643391357])
    

    

Lin42=np.array([136])



period=3.5
omega=2.*math.pi/period  # period is about 3.5 hours
b1_lin42=5  #basal transcription
b2_lin42=50 #max transcription

ndata=1000
t=np.linspace(0,30,ndata)
delta_t=t[1]-t[0]

ndataa=200
tt=np.linspace(0,4,ndataa)
delta_tt=tt[1]-tt[0]

#burst=40*0.075 # UNC-5 protein burst size
burst=1.0 # UNC-5 protein burst size

LIN42_burst=40*0.075 # LIN-42 protein burst size


lin42_scale=1.0 #tuning for burst frequency through manipulating the protein degradation
unc5_scale=1.0

titlelabel='WT'


beta_lin42=1.0*14.36*lin42_scale #translation rate


gamma_lin42=7.92*lin42_scale # protein degradation

gamma_unc5=1.98*unc5_scale

Yscale_unc5_as_g=Y_unc5_as_g*gamma_unc5/burst
f_as_g = interp1d(X_Lin42, Yscale_unc5_as_g,fill_value='extrapolate')
 
Yscale_unc5_de_r=Y_unc5_de_r*gamma_unc5/burst
f_de_r = interp1d(X_Lin42, Yscale_unc5_de_r,fill_value='extrapolate')


Yscale_unc5_noAC_as_g=Y_unc5_noAC_as_g*gamma_unc5/burst
f_noAC_as_g = interp1d(X_Lin42, Yscale_unc5_noAC_as_g,fill_value='extrapolate')
 
Yscale_unc5_noAC_de_r=Y_unc5_noAC_de_r*gamma_unc5/burst
f_noAC_de_r = interp1d(X_Lin42, Yscale_unc5_noAC_de_r,fill_value='extrapolate')



lin42p_max=LIN42_burst*beta_lin42*(b1_lin42+b2_lin42)/gamma_lin42 #max(X_Lin42) #
lin42p_min=LIN42_burst*beta_lin42*(b1_lin42)/gamma_lin42 #min(X_Lin42) #
print('lin42p,max and min: ', lin42p_max,lin42p_min)
unc5_max=burst*f_de_r(lin42p_min)/gamma_unc5
unc5_min=burst*f_de_r(lin42p_max)/gamma_unc5
print('unc5, max and min: ',unc5_max,unc5_min)
kappa=1/2*(lin42p_max-lin42p_min)

unc5_mean=np.zeros(ndata)
unc5_std=np.zeros(ndata)
unc5_det=np.zeros(ndata)

#%%
unc5_noDAF_end=np.zeros(len(X_Lin42))
unc5_noAC_end=np.zeros(len(X_Lin42))
unc5_noDAFnoAC_end=np.zeros(len(X_Lin42))
unc5_end=np.zeros(len(X_Lin42))


ntraj=100
lin42p=np.zeros((ndata,ntraj))

lin42p_det=np.zeros(ndata)
lin42p_mean=np.zeros(ndata)
lin42p_std=np.zeros(ndata)

unc5H=np.zeros((ndata,ntraj))

h=np.zeros((ndata,ntraj))

h_det=np.zeros(ndata)
h_mean=np.zeros(ndata)
h_std=np.zeros(ndata)

start_time=11+1.5*period
end_time=30


sdata=1
scale=np.linspace(0.4,0.4,sdata)
Lin42_mean=np.zeros((ndata,sdata))
Hill_mean=np.zeros((ndata,sdata))
Hill_std=np.zeros((ndata,sdata))
Hill_det=np.zeros((ndata,sdata))
Hill_Lin42_mean=np.zeros((ndata,sdata))

Unc5p_mean=np.zeros((ndata,sdata))
Unc5p_std=np.zeros((ndata,sdata))
Unc5p_CV=np.zeros((ndata,sdata))

unc5_det_all=np.zeros((ndata,sdata))

index=0
for q in scale:

  def lin42mH(t):
   return b1_lin42+b2_lin42*q    

#  print (q)

  for i in range(1,ndata):
# simulation of noisy LIN-42 input
     lin42p[0]=Lin42[index]
     lin42p_det[0]=Lin42[index]
     lin42_prod=beta_lin42*lin42mH(t[i])*delta_t
     lin42_deg=[gamma_lin42*lin42p[i-1,j]*delta_t for j in range(ntraj)]

     noise=[np.random.normal(size=(ntraj))[j]*LIN42_burst*np.sqrt(lin42_prod)+\
             np.random.normal(size=(ntraj))[j]*np.sqrt(abs(lin42_deg[j])) 
               for j in range(ntraj)]

     lin42p[i]=[lin42p[i-1,j]+LIN42_burst*lin42_prod-lin42_deg[j]+noise[j] \
               for j in range(ntraj)]
     lin42p_mean[i]=np.mean(lin42p[i])
     lin42p_std[i]=np.std(lin42p[i])
     lin42p_det[i]=lin42p_det[i-1]\
         +LIN42_burst*lin42_prod\
             -gamma_lin42*lin42p_det[i-1]*delta_t

     for j in range(ntraj):
         if lin42p[i,j] > lin42p[i-1,j]:
 # switch for circuit I or II

          h[i,j]=f_noAC_as_g(lin42p[i-1,j])  
 #         h[i,j]=f_as_g(lin42p[i-1,j]) 
 #
         else:
           h[i,j]=f_noAC_de_r(lin42p[i-1,j])
 #         h[i,j]=f_de_r(lin42p[i-1,j])
         
     h_det[i]=f_noAC_de_r(lin42p_det[i-1]) 
 #    h_det[i]=f_de_r(lin42p_det[i-1])   

 # simulation of UNC-5 output
    
     unc5_prod=[h[i,j]*delta_t for j in range(ntraj)] #leaving out burst
     unc5_deg=[gamma_unc5*unc5H[i-1,j]*delta_t for j in range(ntraj)] 
     noise=[np.random.normal(size=(ntraj))[j]*burst*np.sqrt(unc5_prod[j])+\
            np.random.normal(size=(ntraj))[j]*np.sqrt(abs(unc5_deg[j])) \
               for j in range(ntraj)]
     unc5H[i]=[unc5H[i-1,j]+unc5_prod[j]*burst-unc5_deg[j]+noise[j] for j in range(ntraj)]
     unc5_mean[i]=np.mean(unc5H[i])
     unc5_std[i]=np.std(unc5H[i])
     unc5_det[i]=unc5_det[i-1]\
         +burst*h_det[i]*delta_t\
             -gamma_unc5*unc5_det[i-1]*delta_t
 
     Lin42_mean[i,index]=np.mean(lin42p[i])     


     Unc5p_mean[i,index]=np.mean(unc5H[i])
     Unc5p_std[i,index]=np.std(unc5H[i])
     Unc5p_CV[i,index]=np.std(unc5H[i])/np.mean(unc5H[i])
     unc5_det_all[i,index]=unc5_det[i]
  
  
  index=index+1           
  print (q)

#plot stochastic UNC-5 output
f6=plt.plot(t, Unc5p_mean)
plt.xlabel('Time', fontsize=18)
plt.ylabel('Stochasic UNC-5', fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.show()

#plot stochastic LIN-42 of a constant mean
f7=plt.plot(t,  lin42p[:,0])
plt.xlabel('Time', fontsize=18)
plt.ylabel('LIN-42', fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.xlim([0,20]);
plt.ylim([0,600]);
plt.show()

#plot deterministic UNC-5 output
f8=plt.plot(t,  unc5_det_all)
plt.xlabel('Time', fontsize=18)
plt.ylabel('Deterministic UNC-5', fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.show()



#f4=plt.plot(lin42p_det,unc5_det) #[625:800
#plt.show()
   
    
