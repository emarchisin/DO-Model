#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 16:06:23 2024

@author: emmamarchisin
"""

## modeling DO

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import datetime
import sys

#Import data
MEdata_all=pd.read_csv("/Users/emmamarchisin/Research/Data/NTL LTER 1995-2022.csv")
#MEdata_allclean=MEdata_all.dropna() #drop Nan rows
ME_1meter=MEdata_all[MEdata_all['depth']==1]
ME_1meter=ME_1meter.dropna(subset=['wtemp','shortwaverad'])
ME_1meter.reset_index(inplace=True)
Date=pd.to_datetime(ME_1meter['sampledate']) #biweekly
DO_actual=ME_1meter['o2']
Date_daily=pd.date_range(start=Date.min(),end=Date.max(),freq='D') #created daily values
#Relevant columns
T=ME_1meter['wtemp'] #(degC)
I=ME_1meter['shortwaverad']

#Interpolate I and T
ForInterpdata=pd.DataFrame({'Date':Date,'T':T,'I':I})
T_daily=np.interp(Date_daily.astype(np.int64),ForInterpdata['Date'].astype(np.int64),ForInterpdata['T'])
I_daily=np.interp(Date_daily.astype(np.int64),ForInterpdata['Date'].astype(np.int64),ForInterpdata['I'])
daily_inputs=pd.DataFrame({
    'Date':Date_daily,
    'T daily':T_daily,
    'I daily':I_daily
    })                

Tw_abs=T_daily+273.15 #Absolute Water Temp (K)

# Parameters
ThetaAr=1.08 #Ahrenius constant
ThetaArR=1.04 #Ahrenius constant for Respiration
T20=20 #(degC)
IP=(3e10^-5) #madeup 
nDays=365 #(days)
dt=1 #timestep (days)
z=5 #thermocline depth (m) 
cDOC=.003 #Resp DOC coeff
cPOC=.15 #Resp POC coeff
c1=.2#Respiration Coeff
c2=0.01 #respiration Coeff
cLEC1=.02
cLEC2=.7
OC= 5 # organic carbon (mgL^-1) made up
k=5 #made up (md^-1)
elev=.256 #elevation (m)
LA=39*10**7 #(m^2) Lake Area of Mendota
LV=LA*z #volume
Sal=0 #Salinity (in freshwter lake=0)

Qin=38/365*LA #estimation g/yr streaminflow
Qout=Qin*.8 #estimation stream outflow


#Setup
nt=len(Date_daily) #number of steps
ArGPP=np.empty(nt+1)
ArR=np.empty(nt+1)
GPP=np.empty(nt+1)
R=np.empty(nt+1)
Fatm=np.empty(nt+1)
Fatm2=np.empty(nt+1)
Calculated_DO=np.empty(nt+1)
DOCr=np.empty(nt+1)
POC=np.empty(nt+1)
RDOCr=np.empty(nt+1)
RPOC=np.empty(nt+1)
Rtot=np.empty(nt+1)
Rboth=np.empty(nt+1)
Calculated_DO2=np.empty(nt+1)
LEC=np.empty(nt+1)
Secchi=np.empty(nt+1)

Calculated_DO[0]=13.1*LV #inital DO condition
Calculated_DO2[0]=13.1*LV 
DOCr[0]=5*LV
POC[0]=2.0*LV

#DO Sat Calculations
PsiSF=-139.34411+((1.575701*10**5)/Tw_abs)-((6.642308*10**7)/(Tw_abs**2))+((1.243800*10**10)/(Tw_abs**3))-((8.621949*10**11)/(Tw_abs**4))#DO sat of sealevel freshwater
PsiSal=(2.718281828459045)**(-Sal*(1.7674*10**-2+(10.754/Tw_abs)-(2140.7/(Tw_abs**2)))) #fractional reduction of Sat from Salinity
PsiElev=1-(0.11988*elev)+(6.10834*(10**-3)*(elev**2))-1.60747*(10**-4)*(elev**3)#fractional reduction of freshwater sat sealevel to elevation
DOsat=PsiElev*PsiSal*(2.718281828459045)**(PsiSF)*LV


for i in range(1,nt): 
    ArGPP[i]=(ThetaAr**(T_daily[i]-T20)) #Ahrenius value
    ArR[i]=(ThetaArR**(T_daily[i]-T20))
    GPP[i]=IP*I_daily[i]*ArGPP[i]*LV #DO variables
    Fatm[i]=k*((DOsat[i]-Calculated_DO[i-1])/z)
    Fatm2[i]=k*((DOsat[i]-Calculated_DO2[i-1])/z)
    RDOCr[i]=DOCr[i-1]*cDOC*ArR[i] #Resp for DOCr scaled up for LV already
    RPOC[i]=POC[i-1]*cPOC*ArR[i] #Resp for POC scaled up for LV already
    Rtot[i]=(c1*GPP[i]+c2*OC*LV)*ArR[i] #respiration calcualted as a function of GPP
    Rboth[i]=(32/12)*(RDOCr[i]+RPOC[i]) #Combined respiration calcualted from respiration of DOC and POC
    #everything below is mass balance equation
    DOCr[i]=DOCr[i-1]-RDOCr[i]+Qin-Qout #dDOCr/dt
    POC[i]=POC[i-1]+(12/32)*GPP[i]-RPOC[i]-0.4*POC[i-1]
    Calculated_DO[i]=Calculated_DO[i-1]+GPP[i]-Rtot[i]+Fatm[i] #dDO/dt
    Calculated_DO2[i]=Calculated_DO2[i-1]+GPP[i]-Rboth[i]+Fatm2[i]
    
    LEC[i]=.125+cLEC1*(DOCr[i]/LV)+cLEC2*(POC[i]/LV)
    Secchi[i]=1.7/LEC[i]

#results
respiration=pd.DataFrame({
    'Date':Date_daily[:nt],
    'R from DO calc':Rtot[:nt],
    'Radded':Rboth[:nt],
    'RDOCr':RDOCr[:nt],
    'RPOC':RPOC[:nt],
    })

results=pd.DataFrame({
    'Date':Date_daily[:nt],
    'Calculated_DO': Calculated_DO[:nt],
    'GPP':GPP[:nt],
    'R':Rtot[:nt],
    'Fatm':Fatm[:nt],
    'ArGPP':ArGPP[:nt],
    'ArR':ArR[:nt],
    })


#print(GPP/LV)
#print(R/LV)
#print(Fatm/LV)
print(Calculated_DO/LV) #back in mg/L units
#print (DOsat/LV)


#Graphing
plt.figure(figsize=(12,6))
plt.plot(Date_daily[237:1332], (Calculated_DO[237:1332]/LV), label='Calculated DO', color='blue',marker='o')
plt.plot(Date_daily[237:1332],(DOsat[237:1332]/LV),label='Calculated DO', color='green',marker='o')
plt.plot(Date[15:61], (DO_actual[15:61]), label='Calculated DO', color='red',marker='o')
plt.plot(Date_daily[237:1332], (Calculated_DO2[237:1332]/LV), label='Calculated DO', color='purple',marker='o')
plt.xlabel('Date')
plt.ylabel('Calculated DO (mg/L')
plt.show()

plt.plot(Date_daily[237:1332], GPP[237:1332], label='GPP', color='green',marker='o')
plt.xlabel('Date')
plt.xticks(rotation=45)
plt.ylabel('GPP')
plt.show()

plt.plot(Date_daily[237:1332], Rtot[237:1332], label='R', color='red',marker='o')
plt.plot(Date_daily[237:1332], Rboth[237:1332], label='R', color='blue',marker='o')
plt.xlabel('Date')
plt.xticks(rotation=45)
plt.ylabel('Respiration')
plt.show()

plt.plot(Date_daily[237:1332], Fatm[237:1332], label='Fatm', color='purple',marker='o')
plt.xlabel('Date')
plt.xticks(rotation=45)
plt.ylabel('Fatm')
plt.show()

plt.plot(Date_daily[237:1332],(DOCr[237:1332]/LV), label='DOCr', color='purple',marker='o')
plt.xlabel('Date')
plt.xticks(rotation=45)
plt.ylabel('DOCr')
plt.show()

plt.plot(Date_daily[237:1332], (POC[237:1332]/LV), label='POC', color='purple',marker='o')
plt.xlabel('Date')
plt.xticks(rotation=45)
plt.ylabel('POC')
plt.show()


plt.scatter(T_daily[237:1332],(Calculated_DO/LV)[237:1332],linewidths=(.1))
plt.xlabel('Temp (Deg C)')
plt.ylabel('Calculated DO (mg/L')
plt.show()

plt.plot(Date_daily[237:1332], Secchi[237:1332], label='Secchi Depth', color='purple',marker='o')
plt.xlabel('Date')
plt.xticks(rotation=45)
plt.ylabel('Secchi Depth')
plt.show()

secchi_resanddate=pd.DataFrame({
    'Date':Date_daily[:nt],
    'Secchi Depth':Secchi[:nt]})

secchi_resanddate.to_csv('~/Downloads/SecchiandDate')